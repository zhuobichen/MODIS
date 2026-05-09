#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 MODIS MCD12Q1 (LC_Type5) 最近邻归并到 CMAQ 规则网格，并统计每格各 PFT 像元个数。
- 只读取 ncl_convert2nc 转出的 .nc（MCD12Q1.AYYYY001.hXXvYY....nc）
- 自动根据 h/v 计算经纬度（MODIS Sinusoidal → 地理坐标，反投影公式正确）
- 处理 byte/NaN 缺测，限定有效类别 0..11
"""

import os
import glob
import numpy as np
import xarray as xr
from scipy.spatial import cKDTree

# ============== 用户参数 ==============
indir   = "./Input"                                  # .nc 文件所在目录
grid_nc = "GRIDCRO2D_2000121_GuangDongD3"         # CMAQ 网格
YEAR    = 2000
out_nc  = f"MODIS_PFT_{YEAR}_GuangDong_auto.nc"
hlines  = [27,28]                                   # e.g. 27,28
vlines  = [6]                                        # e.g. 06
ntype   = 12
chunk_n = 2_000_000

# ============== 目标网格 ==============
ds_grid = xr.open_dataset(grid_nc)
dstlats = ds_grid["LAT"].isel(TSTEP=0, LAY=0).values  # (nlat, nlon)
dstlons = ds_grid["LON"].isel(TSTEP=0, LAY=0).values
nlat, nlon = dstlats.shape

def ll_to_unitvec(lat_deg, lon_deg):
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    cl  = np.cos(lat)
    return np.stack([cl*np.cos(lon), cl*np.sin(lon), np.sin(lat)], axis=-1)

dst_xyz = ll_to_unitvec(dstlats.ravel(), dstlons.ravel())
tree = cKDTree(dst_xyz)

# ============== Sinusoidal 反投影（正确版） ==============
def modis_tile_latlon(h, v, nrows=2400, ncols=2400, pixel_size=463.31271653):
    """
    输入：h、v（两位数，含前导0）、像元大小 463.31271653 m
    返回：该瓦片 2D 的 (lat, lon)（度）
    """
    R = 6371007.181        # 球半径 (m)
    TILE = 1200000.0       # 单瓦片边长 (m)

    # 瓦片左上角（以整球为基准）
    x_ul = -20015109.354 + h * TILE
    y_ul =  10007554.677 - v * TILE

    # 像元中心坐标
    x = x_ul + (np.arange(ncols) + 0.5) * pixel_size
    y = y_ul - (np.arange(nrows) + 0.5) * pixel_size
    xv, yv = np.meshgrid(x, y)

    # 反投影：Sinusoidal → geographic
    lat_rad = yv / R
    # 避免极区 cos(lat)=0 的除零
    coslat = np.cos(lat_rad)
    coslat[coslat == 0] = np.finfo(np.float64).eps
    lon_rad = xv / (R * coslat)

    lat = np.rad2deg(lat_rad)
    lon = np.rad2deg(lon_rad)
    return lat, lon

# ============== 读取 LC_Type5 并生成经纬度 ==============
def read_modis_tile_from_nc(ncfile):
    """
    只读 .nc（ncl_convert2nc 输出），提取 LC_Type5，自动生成经纬度
    返回：lat1d, lon1d, pft1d
    """
    base = os.path.basename(ncfile)

    # 从文件名解析 h/v
    try:
        h = int(base.split(".h")[1][:2])
        v = int(base.split("v")[1][:2])
    except Exception:
        raise ValueError(f"无法从文件名解析 h/v: {base}")

    ds = xr.open_dataset(ncfile, mask_and_scale=False, decode_times=False)
    if "LC_Type5" not in ds.variables:
        ds.close()
        raise RuntimeError(f"{base} 中缺少 LC_Type5 变量")

    arr = ds["LC_Type5"].values
    ds.close()

    # arr 可能是带 NaN 的浮点或 int8；统一成 int16，并把 NaN 当缺测(-1)
    if np.issubdtype(arr.dtype, np.floating):
        miss = ~np.isfinite(arr)
        arr = arr.astype(np.int16, copy=False)
        arr[miss] = -1
    else:
        arr = arr.astype(np.int16, copy=False)

    # byte→uint8（-1→255），随后再做有效类别筛选
    arr = np.where(arr < 0, arr + 256, arr)

    # 生成该瓦片的经纬度
    lat2d, lon2d = modis_tile_latlon(h, v)

    # 仅保留真实 PFT 类别（0..11）
    mask = (arr >= 0) & (arr <= 11)
    if not np.any(mask):
        return None, None, None

    pft1d = arr[mask].ravel().astype(np.int16)
    lat1d = lat2d[mask].ravel()
    lon1d = lon2d[mask].ravel()
    return lat1d, lon1d, pft1d

# ============== 主循环：最近邻 + 计数 ==============
pft_counts = np.zeros((ntype, nlat, nlon), dtype=np.int32)

for v in vlines:
    for h in hlines:
        patt = os.path.join(indir, f"MCD12Q1.A{YEAR}001.h{h:02d}v{v:02d}.061.*.nc")
        files = sorted(glob.glob(patt))
        if not files:
            print(f"[WARN] 未找到 .nc 瓦片: h{h:02d}v{v:02d}")
            continue

        ncfile = files[0]
        print(f"[INFO] 处理: {os.path.basename(ncfile)}")

        try:
            lat1d, lon1d, pft1d = read_modis_tile_from_nc(ncfile)
            if lat1d is None:
                print(f"[INFO] {os.path.basename(ncfile)} 无有效像元，跳过")
                continue
        except Exception as e:
            print(f"[ERROR] 打开 {ncfile} 失败: {e}")
            continue

        src_xyz = ll_to_unitvec(lat1d, lon1d)
        N = src_xyz.shape[0]
        for i0 in range(0, N, chunk_n):
            i1 = min(i0 + chunk_n, N)
            # SciPy >=1.6 可用 workers 并行
            try:
                dist, idx = tree.query(src_xyz[i0:i1], k=1, workers=-1)
            except TypeError:
                dist, idx = tree.query(src_xyz[i0:i1], k=1)

            ilat = idx // nlon
            ilon = idx % nlon
            cls  = pft1d[i0:i1]

            # === 新增：距离阈值判断（3000/sqrt(2) 米） ===
            R = 6371007.181
            max_dist_m = 3000.0 / np.sqrt(2)
            max_dist_rad = max_dist_m / R
            good = (cls >= 0) & (cls < ntype) & (dist < max_dist_rad)
            if np.any(good):
                np.add.at(pft_counts, (cls[good], ilat[good], ilon[good]), 1)

# ============== 输出 ==============
out = xr.Dataset(
    data_vars=dict(
        PFT=(("type", "y", "x"), pft_counts.astype(np.int32)),
        lat=(("y", "x"), dstlats.astype(np.float32)),
        lon=(("y", "x"), dstlons.astype(np.float32)),
    ),
    coords=dict(
        type=np.arange(ntype, dtype=np.int16),
        y=np.arange(nlat, dtype=np.int32),
        x=np.arange(nlon, dtype=np.int32),
    ),
    attrs=dict(
        title="MODIS MCD12Q1 LC_Type5 → CMAQ (nearest neighbor, sinusoidal inverse fixed)",
        note_cn="每格各 PFT 的像元计数；可按格内总数归一化得到比例。",
        method="KDTree (LL on unit sphere), MODIS sinusoidal inverse projection",
        year=YEAR,
    ),
)

comp = dict(zlib=True, complevel=4, shuffle=True)
out.to_netcdf(out_nc, encoding={"PFT": comp, "lat": comp, "lon": comp})
print(f"[OK] 输出结果文件 -> {out_nc}")
