#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MCD12Q1 (LC_Type5 分类变量) → CMAQ 网格
最近邻插值（KDTree），经纬度按 MODIS 正弦投影严格计算
修复：TILE_SIZE/原点常数取值错误导致的“全蓝/inf”问题
"""

import os
import glob
import numpy as np
import xarray as xr
from pyhdf.SD import SD, SDC
from scipy.spatial import cKDTree

# -----------------------------
# 用户参数
# -----------------------------
indir   = "./Input"                               # HDF 输入目录
grid_nc = "GRIDCRO2D_2000121_GuangDongD3"      # CMAQ 网格
YEAR    = 2000
hlines  = [27, 28]
vlines  = [6]
out_nc  = f"MODIS_LCType5_{YEAR}_NN.nc"

# -----------------------------
# MODIS 正弦投影常数（标准值）
# -----------------------------
R         = 6371007.181                  # 球半径 (m)
PIX_SIZE  = 463.3127165                  # 像元尺寸 (m)
NCOLS     = 2400
NROWS     = 2400
TILE_SIZE = PIX_SIZE * NCOLS             # 1,111,950.519667 m  ← 关键修复
X0        = -20015109.354                # 全球左边界 (m)
Y0        =  10007554.677                # 全球上边界 (m)

# 最近邻阈值：像元半对角线，稍放宽 10%
DIST_THRESH = (PIX_SIZE * np.sqrt(2) / 2.0) * 1.1
CHUNK = 300_000

# -----------------------------
# 工具：投影与经纬度
# -----------------------------
def modis_tile_latlon(h, v):
    """计算瓦片 h,v 的像元中心经纬度，与 NCL GridLat/GridLon 一致"""
    x_ul = X0 + h * TILE_SIZE          # 瓦片左上角 x
    y_ul = Y0 - v * TILE_SIZE          # 瓦片左上角 y（v 向下）
    x = x_ul + (np.arange(NCOLS) + 0.5) * PIX_SIZE
    y = y_ul - (np.arange(NROWS) + 0.5) * PIX_SIZE
    xv, yv = np.meshgrid(x, y)         # 以米为单位的正弦投影坐标
    lat = np.rad2deg(yv / R)
    lon = np.rad2deg(xv / (R * np.cos(yv / R)))
    return lat, lon

def geo2sinu(lat_deg, lon_deg):
    """地理经纬度 → 正弦投影米坐标"""
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    x = R * lon * np.cos(lat)
    y = R * lat
    return x, y

# -----------------------------
# 读取 LC_Type5
# -----------------------------
def read_modis_LCType5(hdf_path):
    hdf = SD(hdf_path, SDC.READ)
    arr = np.array(hdf.select("LC_Type5")[:], dtype=np.int16)
    hdf.end()
    # 过滤无效
    arr[(arr < 0) | (arr >= 254)] = -1
    return arr

# -----------------------------
# 读取 CMAQ 网格
# -----------------------------
ds_grid = xr.open_dataset(grid_nc)
lat_cmaq = ds_grid["LAT"].isel(TSTEP=0, LAY=0).values
lon_cmaq = ds_grid["LON"].isel(TSTEP=0, LAY=0).values
ny, nx = lat_cmaq.shape

x_cmaq, y_cmaq = geo2sinu(lat_cmaq, lon_cmaq)
pts_cmaq = np.stack([x_cmaq.ravel(), y_cmaq.ravel()], axis=-1)

lc_out  = np.full((ny, nx), -1, dtype=np.int16)
dist_nn = np.full((ny, nx), np.inf, dtype=np.float32)

# 目标网格经纬度包络（便于快速重叠判断）
lat_min_d, lat_max_d = np.nanmin(lat_cmaq), np.nanmax(lat_cmaq)
lon_min_d, lon_max_d = np.nanmin(lon_cmaq), np.nanmax(lon_cmaq)

print(f"[INFO] CMAQ bbox: lat[{lat_min_d:.4f},{lat_max_d:.4f}] "
      f"lon[{lon_min_d:.4f},{lon_max_d:.4f}]")

# -----------------------------
# 主循环：逐瓦片最近邻
# -----------------------------
total_assigned = 0

for v in vlines:
    for h in hlines:
        patt = os.path.join(
            indir, f"MCD12Q1.A{YEAR}001.h{h:02d}v{v:02d}.061.*.hdf"
        )
        files = sorted(glob.glob(patt))
        if not files:
            print(f"[WARN] 未找到瓦片 h{h:02d}v{v:02d}")
            continue

        hfile = files[0]
        print(f"[INFO] 处理: {os.path.basename(hfile)}")

        # 源数据与经纬度
        lc = read_modis_LCType5(hfile)
        lat, lon = modis_tile_latlon(h, v)

        # 快速重叠判断（按经纬度包络框）
        t_latmin, t_latmax = float(lat.min()), float(lat.max())
        t_lonmin, t_lonmax = float(lon.min()), float(lon.max())
        overlap = not (t_latmax < lat_min_d or t_latmin > lat_max_d or
                       t_lonmax < lon_min_d or t_lonmin > lon_max_d)
        print(f"       tile bbox: lat[{t_latmin:.4f},{t_latmax:.4f}] "
              f"lon[{t_lonmin:.4f},{t_lonmax:.4f}] -> overlap={overlap}")
        if not overlap:
            continue

        mask = (lc >= 0)
        if not np.any(mask):
            print("       [INFO] 瓦片无有效像元，跳过")
            continue

        xs, ys = geo2sinu(lat[mask], lon[mask])
        cls = lc[mask]
        tree = cKDTree(np.stack([xs, ys], axis=-1))

        assigned_here = 0
        N = pts_cmaq.shape[0]
        for i0 in range(0, N, CHUNK):
            i1 = min(i0 + CHUNK, N)
            d, idx = tree.query(pts_cmaq[i0:i1], k=1, workers=-1)

            ok = d <= DIST_THRESH
            if not np.any(ok):
                continue

            i_all = np.arange(i0, i1)[ok]
            iy, ix = i_all // nx, i_all % nx

            nearer = d[ok] < dist_nn[iy, ix]
            if np.any(nearer):
                iy2, ix2 = iy[nearer], ix[nearer]
                lc_out[iy2, ix2] = cls[idx[ok][nearer]]
                dist_nn[iy2, ix2] = d[ok][nearer]
                assigned_here += nearer.sum()

        total_assigned += assigned_here
        cover = assigned_here / (ny * nx) * 100.0
        print(f"       [OK] 本瓦片新赋值 {assigned_here} 个格点 "
              f"({cover:.2f}% of domain)")

# -----------------------------
# 收尾与输出
# -----------------------------
print(f"[SUM] 总共赋值 {total_assigned} / {ny*nx} 个格点 "
      f"({100.0*total_assigned/(ny*nx):.2f}%)")

out = xr.Dataset(
    data_vars=dict(
        LC_Type5=(("y", "x"), lc_out),
        lat=(("y", "x"), lat_cmaq.astype(np.float32)),
        lon=(("y", "x"), lon_cmaq.astype(np.float32)),
        dist_m=(("y", "x"), dist_nn.astype(np.float32)),
    ),
    attrs=dict(
        title="MCD12Q1 LC_Type5 → CMAQ 最近邻 (Sinusoidal, fixed constants)",
        note_cn="分类变量只能最近邻；常数使用 MODIS 标准值，避免整体偏移",
        year=YEAR,
    ),
)
comp = dict(zlib=True, complevel=4, shuffle=True)
out.to_netcdf(out_nc, encoding={"LC_Type5": comp, "lat": comp, "lon": comp, "dist_m": comp})
print(f"[OK] 输出结果 -> {out_nc}")
