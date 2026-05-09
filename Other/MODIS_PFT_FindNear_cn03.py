#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MCD12Q1 (LC_Type5 离散分类) → CMAQ 规则网格
—— 每个 CMAQ 网格中心为正方形 3 km 邻域，统计各类别百分比（0..11），保存 NetCDF 与 CSV。

优化思路：
1) 先把所有 MODIS 源点合并到一个 KDTree（正弦投影米坐标）；
2) 对每个 CMAQ 点：用 query_ball_point(半径=√2*1.5km) 快速圈选候选，然后再做正方形判定；
3) numpy.bincount(minlength=12) 直出 12 类计数 → 百分比。
"""

import os
import glob
import math
import numpy as np
import xarray as xr
import pandas as pd
from pyhdf.SD import SD, SDC
from scipy.spatial import cKDTree
from multiprocessing import Pool, cpu_count

# ===========================
# 用户参数
# ===========================
indir    = "./Input"                       # HDF 输入目录
grid_nc  = "GRIDCRO2D_2000121_GuangDongD3" # CMAQ 网格（.nc 或无后缀均可）
YEAR     = 2000
hlines   = [27, 28]                        # h 索引
vlines   = [6]                             # v 索引
exclude_tiles = set()                      # 例如 {"h30v06"}
SIDE_M   = 3000.0                          # 正方形边长（米）
HALF_M   = SIDE_M / 2.0
RADIUS   = HALF_M * math.sqrt(2.0)         # 圆半径（先圈选，再做方形过滤）
BATCH    = 7000                          # 每批处理的 CMAQ 点数
PARALLEL = False                           # 是否并行（默认否）
N_PROCS  = max(1, min(cpu_count(), 8))     # 并行进程数（启用并行时使用）
out_nc   = f"PFT_frac_{YEAR}_3km_square.nc"
out_csv  = f"PFT_frac_{YEAR}_3km_square.csv"

# ===========================
# MODIS 正弦投影常数
# ===========================
R         = 6371007.181        # 球半径 (m)
PIX_SIZE  = 463.3127165        # 像元尺寸 (m)
NCOLS     = 2400
NROWS     = 2400
TILE_SIZE = PIX_SIZE * NCOLS   # 1,111,950.519667 m
X0        = -20015109.354      # 全球左边界 (m)
Y0        =  10007554.677      # 全球上边界 (m)

# ===========================
# 官方 LC_Type5 名称（全称列名）
# ===========================
CLASS_NAMES = [
    "Water",
    "Evergreen Needleleaf trees",
    "Evergreen Broadleaf trees",
    "Deciduous Needleleaf trees",
    "Deciduous Broadleaf trees",
    "Shrub",
    "Grass",
    "Cereal crops",
    "Broad-leaf crops",
    "Urban and built-up",
    "Snow and ice",
    "Barren or sparse vegetation",
]

# ===========================
# 工具函数
# ===========================
def modis_tile_xy(h: int, v: int):
    """返回瓦片 h,v 的像元中心正弦投影坐标 (xv, yv)，单位: 米。"""
    x_ul = X0 + h * TILE_SIZE            # 瓦片左上角 x
    y_ul = Y0 - v * TILE_SIZE            # 瓦片左上角 y（v 向下）
    x = x_ul + (np.arange(NCOLS) + 0.5) * PIX_SIZE
    y = y_ul - (np.arange(NROWS) + 0.5) * PIX_SIZE
    xv, yv = np.meshgrid(x, y)
    return xv, yv

def read_modis_lc_type5(hdf_path: str) -> np.ndarray:
    """读取 HDF 里的 LC_Type5（uint8），转 int16，并把无效置为 -1。"""
    hdf = SD(hdf_path, SDC.READ)
    arr = np.array(hdf.select("LC_Type5")[:], dtype=np.int16)
    hdf.end()
    arr[(arr < 0) | (arr >= 254)] = -1
    return arr

def geo2sinu(lat_deg: np.ndarray, lon_deg: np.ndarray):
    """地理经纬度 → 正弦投影米坐标（与 MODIS Sine 一致的球面公式）。"""
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    x = R * lon * np.cos(lat)
    y = R * lat
    return x, y

# ===========================
# 读取 CMAQ 网格
# ===========================
ds_grid  = xr.open_dataset(grid_nc)
lat_cmaq = ds_grid["LAT"].isel(TSTEP=0, LAY=0).values
lon_cmaq = ds_grid["LON"].isel(TSTEP=0, LAY=0).values
ny, nx   = lat_cmaq.shape
x_cmaq, y_cmaq = geo2sinu(lat_cmaq, lon_cmaq)
pts_cmaq = np.stack([x_cmaq.ravel(), y_cmaq.ravel()], axis=0).T  # (N,2)
N_pts    = pts_cmaq.shape[0]
print(f"[INFO] CMAQ grid: {ny} x {nx} = {N_pts} points")

# ===========================
# 合并瓦片 → 源点库 + KDTree
# ===========================
xs_all, ys_all, cls_all = [], [], []
xmin_d, xmax_d = x_cmaq.min(), x_cmaq.max()
ymin_d, ymax_d = y_cmaq.min(), y_cmaq.max()

def tile_bbox_xy(h, v):
    """给出瓦片边界（米），便于快速判断是否与域重叠。"""
    x0 = X0 + h * TILE_SIZE
    y0 = Y0 - v * TILE_SIZE
    x1 = x0 + TILE_SIZE
    y1 = y0 - TILE_SIZE
    # 注意：y 向下递减
    return min(x0, x1), max(x0, x1), min(y1, y0), max(y1, y0)

for v in vlines:
    for h in hlines:
        hv_tag = f"h{h:02d}v{v:02d}"
        if hv_tag in exclude_tiles: 
            print(f"[SKIP] exclude {hv_tag}")
            continue

        patt  = os.path.join(indir, f"MCD12Q1.A{YEAR}001.{hv_tag}.061.*.hdf")
        files = sorted(glob.glob(patt))
        if not files:
            print(f"[WARN] missing {hv_tag}")
            continue

        hfile = files[0]
        print(f"[INFO] reading {os.path.basename(hfile)}")

        # 与 CMAQ 域做快速 bbox 过滤（用 x/y）
        txmin, txmax, tymin, tymax = tile_bbox_xy(h, v)
        overlap = not (txmax < xmin_d or txmin > xmax_d or tymax < ymin_d or tymin > ymax_d)
        if not overlap:
            print("       no overlap with domain, skip")
            continue

        lc  = read_modis_lc_type5(hfile)
        xv, yv = modis_tile_xy(h, v)
        mask = (lc >= 0)
        if not np.any(mask):
            print("       no valid pixels")
            continue

        xs_all.append(xv[mask].ravel())
        ys_all.append(yv[mask].ravel())
        cls_all.append(lc[mask].ravel().astype(np.int32))

if not xs_all:
    raise RuntimeError("没有任何源像元可用，请检查输入瓦片/路径/年份/索引。")

xs = np.concatenate(xs_all)
ys = np.concatenate(ys_all)
cl = np.concatenate(cls_all)                # (M,)
src_pts = np.stack([xs, ys], axis=1)        # (M,2)
print(f"[INFO] source points: {src_pts.shape[0]:,}")

tree = cKDTree(src_pts)

# ===========================
# 逐点统计（可并行）
# ===========================
def process_chunk(i0_i1):
    """处理 CMAQ 点的一个批次，返回 (rows, cols, counts12) 列表。"""
    i0, i1 = i0_i1
    local = []
    P = pts_cmaq[i0:i1]    # (K,2)
    found = tree.query_ball_point(P, r=RADIUS)  # 列表，长度 K

    for k, idxs in enumerate(found):
        # 批次内的平面索引 → 行列
        gidx = i0 + k
        iy, ix = divmod(gidx, nx)

        if not idxs:  # 候选为空
            local.append((iy, ix, np.zeros(12, dtype=np.int32)))
            continue

        # 方形过滤：|dx|≤HALF_M & |dy|≤HALF_M
        dx = np.abs(xs[idxs] - P[k, 0])
        dy = np.abs(ys[idxs] - P[k, 1])
        ok = (dx <= HALF_M) & (dy <= HALF_M)
        if not np.any(ok):
            local.append((iy, ix, np.zeros(12, dtype=np.int32)))
            continue

        cls_sel = cl[idxs][ok]   # 只会是 0..11
        # 计数（12 类）
        counts = np.bincount(cls_sel, minlength=12)[:12].astype(np.int32)
        local.append((iy, ix, counts))

    return local

# 分派任务
tasks = []
for i0 in range(0, N_pts, BATCH):
    i1 = min(i0 + BATCH, N_pts)
    tasks.append((i0, i1))

if PARALLEL and N_PROCS > 1:
    print(f"[INFO] parallel map with {N_PROCS} procs ...")
    with Pool(processes=N_PROCS) as pool:
        parts = pool.map(process_chunk, tasks)
    results = [item for part in parts for item in part]
else:
    print("[INFO] single-process map ...")
    results = []
    for t in tasks:
        results.extend(process_chunk(t))

# 聚合到 (12, ny, nx) 百分比数组
counts_3d = np.zeros((12, ny, nx), dtype=np.int32)
for (iy, ix, cnt) in results:
    counts_3d[:, iy, ix] = cnt

totals = counts_3d.sum(axis=0, keepdims=True)  # (1, ny, nx)
# 避免除零：总数=0 的格点 → 百分比全 0
with np.errstate(divide='ignore', invalid='ignore'):
    fracs = np.where(totals > 0, counts_3d * 100.0 / totals, 0.0).astype(np.float32)  # (12, ny, nx)

# ===========================
# 写 NetCDF（12类百分比分变量形式）
# ===========================
# 计算索引与坐标
ICELL = np.tile(np.arange(1, nx + 1), ny)
JCELL = np.repeat(np.arange(1, ny + 1), nx)
CELLID = np.arange(1, ny * nx + 1)

# 构建包含所有百分比的变量字典
data_vars = {
    "lat": (("y", "x"), lat_cmaq.astype(np.float32)),
    "lon": (("y", "x"), lon_cmaq.astype(np.float32)),
}

for ci, cname in enumerate(CLASS_NAMES):
    data_vars[cname] = (("y", "x"), fracs[ci].astype(np.float32))

# 构建 Dataset
ds_out = xr.Dataset(
    data_vars=data_vars,
    coords={
        "type": np.arange(12, dtype=np.int16),
        "y": np.arange(ny, dtype=np.int32),
        "x": np.arange(nx, dtype=np.int32),
    },
    attrs=dict(
        title="MCD12Q1 LC_Type5 → CMAQ 3 km Square Neighborhood Fractions",
        description="包含12个分类百分比变量（单位%），保持三维结构 (y, x)。",
        year=YEAR,
        note="Missing neighborhood → all-zero fractions",
    ),
)
comp = dict(zlib=True, complevel=4, shuffle=True)
enc = {k: comp for k in data_vars}
ds_out.to_netcdf(out_nc, encoding=enc)
print(f"[OK] Wrote NetCDF -> {out_nc}")

# ===========================
# 写 CSV（平铺格式）
# ===========================
# 创建 DataFrame
LAT = lat_cmaq.ravel()
LON = lon_cmaq.ravel()

df = pd.DataFrame({
    "CELLID": CELLID,
    "JCELL": JCELL,
    "ICELL": ICELL,
    "LAT": LAT,
    "LON": LON,
})
for ci, cname in enumerate(CLASS_NAMES):
    df[cname] = fracs[ci].ravel()

df.to_csv(out_csv, index=False, encoding="utf-8-sig")
print(f"[OK] Wrote CSV -> {out_csv}")
