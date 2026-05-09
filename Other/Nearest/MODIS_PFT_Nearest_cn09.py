#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MCD12Q1 (LC_Type5 分类变量) → CMAQ 网格
最近邻插值（KDTree），经纬度按 MODIS 正弦投影严格计算
增加：支持通过 skip_tiles = ["h30v06","h29v07"] 跳过指定瓦片
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
grid_nc = "GRIDCRO2D_2000121_GuangDongD2"         # CMAQ 网格文件
YEAR    = 2000
hlines  = [26, 27, 28, 29]               # 横向瓦片
vlines  = [5, 6, 7]                           # 纵向瓦片
out_nc  = f"MODIS_LCType5_{YEAR}_NN_cn09.nc"

# 👉 要跳过的瓦片列表（输入格式形如 "h30v06"）
skip_tiles = ["h26v07","h29v05"]

# -----------------------------
# MODIS 正弦投影常数（标准值）
# -----------------------------
R         = 6371007.181
PIX_SIZE  = 463.3127165
NCOLS     = 2400
NROWS     = 2400
TILE_SIZE = PIX_SIZE * NCOLS
X0        = -20015109.354
Y0        =  10007554.677
DIST_THRESH = (PIX_SIZE * np.sqrt(2) / 2.0) * 1.1
CHUNK = 300_000

# -----------------------------
# 工具函数
# -----------------------------
def modis_tile_latlon(h, v):
    """根据瓦片号(h,v)计算经纬度"""
    x_ul = X0 + h * TILE_SIZE
    y_ul = Y0 - v * TILE_SIZE
    x = x_ul + (np.arange(NCOLS) + 0.5) * PIX_SIZE
    y = y_ul - (np.arange(NROWS) + 0.5) * PIX_SIZE
    xv, yv = np.meshgrid(x, y)
    lat = np.rad2deg(yv / R)
    lon = np.rad2deg(xv / (R * np.cos(yv / R)))
    return lat, lon

def geo2sinu(lat_deg, lon_deg):
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    x = R * lon * np.cos(lat)
    y = R * lat
    return x, y

def read_modis_LCType5(hdf_path):
    """读取 MODIS LC_Type5"""
    hdf = SD(hdf_path, SDC.READ)
    arr = np.array(hdf.select("LC_Type5")[:], dtype=np.int16)
    hdf.end()
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

lat_min_d, lat_max_d = np.nanmin(lat_cmaq), np.nanmax(lat_cmaq)
lon_min_d, lon_max_d = np.nanmin(lon_cmaq), np.nanmax(lon_cmaq)

print(f"[INFO] CMAQ bbox: lat[{lat_min_d:.4f},{lat_max_d:.4f}], lon[{lon_min_d:.4f},{lon_max_d:.4f}]")

# -----------------------------
# 主循环：逐瓦片处理
# -----------------------------
total_assigned = 0

for v in vlines:
    for h in hlines:
        rname = f"h{h:02d}v{v:02d}"

        # ✅ 跳过黑名单瓦片
        if rname in skip_tiles:
            print(f"[SKIP] 跳过瓦片 {rname}")
            continue

        patt = os.path.join(indir, f"MCD12Q1.A{YEAR}001.{rname}.061.*.hdf")
        files = sorted(glob.glob(patt))
        if not files:
            print(f"[WARN] 未找到 {rname}")
            continue

        hfile = files[0]
        print(f"[INFO] 处理: {os.path.basename(hfile)}")

        lc = read_modis_LCType5(hfile)
        lat, lon = modis_tile_latlon(h, v)

        t_latmin, t_latmax = float(lat.min()), float(lat.max())
        t_lonmin, t_lonmax = float(lon.min()), float(lon.max())
        overlap = not (t_latmax < lat_min_d or t_latmin > lat_max_d or
                       t_lonmax < lon_min_d or t_lonmin > lon_max_d)
        if not overlap:
            print(f"       [SKIP] 与CMAQ网格不重叠")
            continue

        mask = (lc >= 0)
        if not np.any(mask):
            print("       [INFO] 瓦片无有效像元")
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
        print(f"       [OK] 新赋值 {assigned_here} 点 ({cover:.2f}% of domain)")

# -----------------------------
# 输出结果
# -----------------------------
print(f"[SUM] 共赋值 {total_assigned}/{ny*nx} 个格点 ({100.0*total_assigned/(ny*nx):.2f}%)")

out = xr.Dataset(
    data_vars=dict(
        LC_Type5=(("y", "x"), lc_out),
        lat=(("y", "x"), lat_cmaq.astype(np.float32)),
        lon=(("y", "x"), lon_cmaq.astype(np.float32)),
        dist_m=(("y", "x"), dist_nn.astype(np.float32)),
    ),
    attrs=dict(
        title="MCD12Q1 LC_Type5 → CMAQ 最近邻 (支持跳过指定瓦片)",
        year=YEAR,
    ),
)
comp = dict(zlib=True, complevel=4, shuffle=True)
out.to_netcdf(out_nc, encoding={"LC_Type5": comp, "lat": comp, "lon": comp, "dist_m": comp})
print(f"[OK] 输出结果 -> {out_nc}")
