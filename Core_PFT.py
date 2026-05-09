"""
Core_PFT.py - MODIS PFT 处理核心逻辑
======================================
从 MODIS MCD12Q1 LC_Type5 提取 PFT 分布到 CMAQ 规则网格。

支持两种模式:
  1. 正方形邻域百分比统计（推荐）: 每个 CMAQ 格点周边正方形区域内各类别百分比
  2. 最近邻单类提取: 每个格点取最近 MODIS 像元类别，赋值 100%

所有路径和参数通过函数签名传入，无硬编码。

提供函数:
  # 模式 1: 正方形邻域
  - run_square_neighborhood_pipeline() -> dict
      完整管道: 读取瓦片 → KDTree → 邻域统计 → 输出 NC + CSV

  # 模式 2: 最近邻
  - run_nearest_neighbor_pipeline() -> dict
      完整管道: 读取瓦片 → KDTree → 最近邻 → 输出 NC

  # 模式 3: 单类提取
  - extract_single_class_pft() -> pd.DataFrame
      从最近邻 NC 结果中提取单类 CSV（每格点一类=100%）
"""

import os
import glob
import math
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
from multiprocessing import Pool, cpu_count
from typing import List, Optional, Tuple, Dict

from modis_geo_utils import (
    modis_tile_xy, geo2sinu,
    tile_bbox_xy, check_overlap,
    CLASS_NAMES,
)
from Core_MODIS_IO import (
    read_modis_lc_type5, read_modis_tile_from_nc,
    read_cmaq_grid, ll_to_unitvec,
)


# ===========================
# 分辨率配置
# ===========================

RESOLUTION_CONFIG = {
    "3km":  {"side_m": 3000.0,  "grid_suffix": "D3"},
    "9km":  {"side_m": 9000.0,  "grid_suffix": "D2"},
    "27km": {"side_m": 27000.0, "grid_suffix": "D1"},
}


def _get_resolution_params(resolution: str) -> Tuple[float, float, float]:
    """根据分辨率名称返回 SIDE_M, HALF_M, RADIUS。"""
    if resolution in RESOLUTION_CONFIG:
        side_m = RESOLUTION_CONFIG[resolution]["side_m"]
    else:
        raise ValueError(f"Unknown resolution: {resolution}. "
                         f"Available: {list(RESOLUTION_CONFIG.keys())}")
    half_m = side_m / 2.0
    radius = half_m * math.sqrt(2.0)
    return side_m, half_m, radius


# ===========================
# 模式 1: 正方形邻域百分比统计
# ===========================

def build_source_points(
    input_dir: str,
    year: int,
    hlines: List[int],
    vlines: List[int],
    exclude_tiles: Optional[set] = None,
    skip_tiles: Optional[List[str]] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """读取所有 MODIS 瓦片，合并源点（正弦投影米坐标 + 分类）。

    参数:
        input_dir: HDF 输入目录
        year: 年份 (e.g. 2000)
        hlines: 横向瓦片编号列表
        vlines: 纵向瓦片编号列表
        exclude_tiles: 要排除的瓦片集合 (e.g. {"h30v06"})
        skip_tiles: 要跳过的瓦片列表 (等同 exclude_tiles)

    返回:
        xs: (M,) 正弦投影 x 坐标 (米)
        ys: (M,) 正弦投影 y 坐标 (米)
        cl: (M,) LC_Type5 分类值 (0..11)
    """
    if exclude_tiles is None:
        exclude_tiles = set(skip_tiles or [])

    xs_all, ys_all, cls_all = [], [], []

    for v in vlines:
        for h in hlines:
            hv_tag = f"h{h:02d}v{v:02d}"
            if hv_tag in exclude_tiles:
                print(f"[SKIP] exclude {hv_tag}")
                continue

            patt = os.path.join(input_dir, f"MCD12Q1.A{year}001.{hv_tag}.061.*.hdf")
            files = sorted(glob.glob(patt))
            if not files:
                print(f"[WARN] missing {hv_tag}")
                continue

            hfile = files[0]
            print(f"[INFO] reading {os.path.basename(hfile)}")

            lc = read_modis_lc_type5(hfile)
            xv, yv = modis_tile_xy(h, v)
            mask = (lc >= 0)
            if not np.any(mask):
                print("       no valid pixels")
                continue

            xs_all.append(xv[mask].ravel())
            ys_all.append(yv[mask].ravel())
            cls_all.append(lc[mask].ravel().astype(np.int32))

    if not xs_all:
        raise RuntimeError("No valid source pixels found. Check input tiles/path/year.")

    xs = np.concatenate(xs_all)
    ys = np.concatenate(ys_all)
    cl = np.concatenate(cls_all)
    print(f"[INFO] source points: {xs.shape[0]:,}")
    return xs, ys, cl


def _process_chunk_square(args):
    """处理 CMAQ 点的一个批次（正方形邻域统计）。用于并行 map。"""
    i0, i1, pts_cmaq, xs, ys, cl, tree, radius, half_m, nx = args
    local = []

    P = pts_cmaq[i0:i1]
    found = tree.query_ball_point(P, r=radius)

    for k, idxs in enumerate(found):
        gidx = i0 + k
        iy, ix = divmod(gidx, nx)

        if not idxs:
            local.append((iy, ix, np.zeros(12, dtype=np.int32)))
            continue

        # 方形过滤
        dx = np.abs(xs[idxs] - P[k, 0])
        dy = np.abs(ys[idxs] - P[k, 1])
        ok = (dx <= half_m) & (dy <= half_m)
        if not np.any(ok):
            local.append((iy, ix, np.zeros(12, dtype=np.int32)))
            continue

        cls_sel = cl[idxs][ok]
        counts = np.bincount(cls_sel, minlength=12)[:12].astype(np.int32)
        local.append((iy, ix, counts))

    return local


def run_square_neighborhood_pipeline(
    year: int,
    resolution: str,
    input_dir: str,
    output_dir: str,
    grid_nc: str,
    hlines: List[int],
    vlines: List[int],
    skip_tiles: Optional[List[str]] = None,
    batch: int = 6000,
    parallel: bool = True,
    n_procs: Optional[int] = None,
) -> Dict:
    """正方形邻域 PFT 百分比统计完整管道。

    参数:
        year: 年份
        resolution: 分辨率 ("3km"/"9km"/"27km")
        input_dir: MODIS HDF 输入目录
        output_dir: 输出目录
        grid_nc: CMAQ 网格 NetCDF 文件路径
        hlines: 横向瓦片编号列表
        vlines: 纵向瓦片编号列表
        skip_tiles: 要跳过的瓦片列表
        batch: 每批处理的 CMAQ 点数
        parallel: 是否并行
        n_procs: 并行进程数 (None=auto)

    返回:
        dict: {"nc": netcdf_path, "csv": csv_path, "shape": (12, ny, nx)}
    """
    side_m, half_m, radius = _get_resolution_params(resolution)

    if n_procs is None:
        n_procs = max(1, min(cpu_count(), 8))

    # 读取 CMAQ 网格
    lat_cmaq, lon_cmaq, ny, nx = read_cmaq_grid(grid_nc)
    x_cmaq, y_cmaq = geo2sinu(lat_cmaq, lon_cmaq)
    pts_cmaq = np.stack([x_cmaq.ravel(), y_cmaq.ravel()], axis=0).T
    N_pts = pts_cmaq.shape[0]
    print(f"[INFO] CMAQ grid: {ny} x {nx} = {N_pts} points, "
          f"resolution={resolution}, side={side_m}m")

    # 构建源点库 + KDTree
    xs, ys, cl = build_source_points(
        input_dir, year, hlines, vlines,
        exclude_tiles=set(skip_tiles or []),
    )
    tree = cKDTree(np.stack([xs, ys], axis=1))
    print(f"[INFO] KDTree built with {xs.shape[0]:,} source points")

    # 分派任务
    tasks = []
    for i0 in range(0, N_pts, batch):
        i1 = min(i0 + batch, N_pts)
        tasks.append((i0, i1, pts_cmaq, xs, ys, cl, tree, radius, half_m, nx))

    # 执行
    if parallel and n_procs > 1:
        print(f"[INFO] parallel map with {n_procs} procs ...")
        with Pool(processes=n_procs) as pool:
            parts = pool.map(_process_chunk_square, tasks)
        results = [item for part in parts for item in part]
    else:
        print("[INFO] single-process map ...")
        results = []
        for t in tasks:
            results.extend(_process_chunk_square(t))

    # 聚合到 (12, ny, nx) 百分比数组
    counts_3d = np.zeros((12, ny, nx), dtype=np.int32)
    for (iy, ix, cnt) in results:
        counts_3d[:, iy, ix] = cnt

    totals = counts_3d.sum(axis=0, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        fracs = np.where(totals > 0, counts_3d * 100.0 / totals, 0.0).astype(np.float32)

    # 输出 NetCDF
    os.makedirs(output_dir, exist_ok=True)
    out_nc = os.path.join(output_dir, f"PFT_frac_{year}_{resolution}_square.nc")
    out_csv = os.path.join(output_dir, f"PFT_frac_{year}_{resolution}_square.csv")

    data_vars = {
        "lat": (("y", "x"), lat_cmaq.astype(np.float32)),
        "lon": (("y", "x"), lon_cmaq.astype(np.float32)),
    }
    for ci, cname in enumerate(CLASS_NAMES):
        data_vars[cname] = (("y", "x"), fracs[ci].astype(np.float32))

    ds_out = xr.Dataset(
        data_vars=data_vars,
        coords={
            "type": np.arange(12, dtype=np.int16),
            "y": np.arange(ny, dtype=np.int32),
            "x": np.arange(nx, dtype=np.int32),
        },
        attrs=dict(
            title=f"MCD12Q1 LC_Type5 → CMAQ {resolution} Square Neighborhood Fractions",
            year=year, resolution=resolution,
            note="Missing neighborhood → all-zero fractions",
        ),
    )
    comp = dict(zlib=True, complevel=4, shuffle=True)
    enc = {k: comp for k in data_vars}
    ds_out.to_netcdf(out_nc, encoding=enc)
    print(f"[OK] Wrote NetCDF -> {out_nc}")

    # 输出 CSV
    CELLID = np.arange(1, ny * nx + 1)
    ICELL = np.tile(np.arange(1, nx + 1), ny)
    JCELL = np.repeat(np.arange(1, ny + 1), nx)

    df = pd.DataFrame({
        "CELLID": CELLID,
        "JCELL": JCELL,
        "ICELL": ICELL,
        "LAT": lat_cmaq.ravel(),
        "LON": lon_cmaq.ravel(),
    })
    for ci, cname in enumerate(CLASS_NAMES):
        df[cname] = fracs[ci].ravel()

    df.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"[OK] Wrote CSV -> {out_csv}")

    return {
        "nc": out_nc,
        "csv": out_csv,
        "shape": (12, ny, nx),
        "resolution": resolution,
    }


# ===========================
# 模式 2: 最近邻单类提取
# ===========================

def run_nearest_neighbor_pipeline(
    year: int,
    input_dir: str,
    output_dir: str,
    grid_nc: str,
    hlines: List[int],
    vlines: List[int],
    skip_tiles: Optional[List[str]] = None,
    use_sinusoidal: bool = True,
    chunk_n: int = 2_000_000,
) -> Dict:
    """最近邻 PFT 提取完整管道（单位球面或正弦投影 KDTree）。

    参数:
        year: 年份
        input_dir: MODIS HDF 或 NC 输入目录
        output_dir: 输出目录
        grid_nc: CMAQ 网格 NetCDF
        hlines, vlines: 瓦片范围
        skip_tiles: 跳过的瓦片
        use_sinusoidal: True=正弦投影KDTree, False=单位球面KDTree
        chunk_n: 批处理大小

    返回:
        dict: {"nc": path, "shape": (ny, nx)}
    """
    lat_cmaq, lon_cmaq, ny, nx = read_cmaq_grid(grid_nc)
    N_pts = ny * nx

    if use_sinusoidal:
        x_cmaq, y_cmaq = geo2sinu(lat_cmaq, lon_cmaq)
        pts_cmaq = np.stack([x_cmaq.ravel(), y_cmaq.ravel()], axis=-1)
        dist_thresh = (463.3127165 * np.sqrt(2) / 2.0) * 1.1
    else:
        pts_cmaq = ll_to_unitvec(lat_cmaq.ravel(), lon_cmaq.ravel())
        dist_thresh = 3000.0 / 6371007.181

    lc_out = np.full((ny, nx), -1, dtype=np.int16)
    dist_nn = np.full((ny, nx), np.inf, dtype=np.float32)

    total_assigned = 0
    exclude_tiles = set(skip_tiles or [])

    for v in vlines:
        for h in hlines:
            rname = f"h{h:02d}v{v:02d}"
            if rname in exclude_tiles:
                print(f"[SKIP] {rname}")
                continue

            patt = os.path.join(input_dir, f"MCD12Q1.A{year}001.{rname}.061.*.hdf")
            files = sorted(glob.glob(patt))
            if not files:
                print(f"[WARN] missing {rname}")
                continue

            hfile = files[0]
            print(f"[INFO] processing {os.path.basename(hfile)}")

            if use_sinusoidal:
                lc = read_modis_lc_type5(hfile)
                xv, yv = modis_tile_xy(h, v)
                mask = (lc >= 0)
                if not np.any(mask):
                    continue
                src_pts = np.stack([xv[mask].ravel(), yv[mask].ravel()], axis=-1)
                cls_data = lc[mask].ravel()
            else:
                lat1d, lon1d, pft1d = read_modis_tile_from_nc(
                    os.path.join(input_dir, os.path.basename(hfile))
                )
                if lat1d is None:
                    continue
                src_pts = ll_to_unitvec(lat1d, lon1d)
                cls_data = pft1d

            tree = cKDTree(src_pts)
            assigned_here = 0

            for i0 in range(0, N_pts, chunk_n):
                i1 = min(i0 + chunk_n, N_pts)
                try:
                    d, idx = tree.query(pts_cmaq[i0:i1], k=1, workers=-1)
                except TypeError:
                    d, idx = tree.query(pts_cmaq[i0:i1], k=1)

                ok = d <= dist_thresh
                if not np.any(ok):
                    continue

                i_all = np.arange(i0, i1)[ok]
                iy, ix = i_all // nx, i_all % nx
                nearer = d[ok] < dist_nn[iy, ix]
                if np.any(nearer):
                    iy2, ix2 = iy[nearer], ix[nearer]
                    lc_out[iy2, ix2] = cls_data[idx[ok][nearer]]
                    dist_nn[iy2, ix2] = d[ok][nearer]
                    assigned_here += nearer.sum()

            total_assigned += assigned_here
            cover = assigned_here / N_pts * 100.0
            print(f"       [OK] assigned {assigned_here} pts ({cover:.2f}%)")

    print(f"[SUM] total assigned {total_assigned}/{N_pts} "
          f"({100.0 * total_assigned / N_pts:.2f}%)")

    os.makedirs(output_dir, exist_ok=True)
    method = "sinusoidal" if use_sinusoidal else "unit_sphere"
    out_nc = os.path.join(output_dir, f"MODIS_LCType5_{year}_NN_{method}.nc")

    out = xr.Dataset(
        data_vars=dict(
            LC_Type5=(("y", "x"), lc_out),
            lat=(("y", "x"), lat_cmaq.astype(np.float32)),
            lon=(("y", "x"), lon_cmaq.astype(np.float32)),
            dist_m=(("y", "x"), dist_nn.astype(np.float32)),
        ),
        attrs=dict(
            title=f"MCD12Q1 LC_Type5 → CMAQ nearest neighbor ({method})",
            year=year,
        ),
    )
    comp = dict(zlib=True, complevel=4, shuffle=True)
    out.to_netcdf(out_nc, encoding={"LC_Type5": comp, "lat": comp, "lon": comp, "dist_m": comp})
    print(f"[OK] Wrote -> {out_nc}")

    return {"nc": out_nc, "shape": (ny, nx)}


# ===========================
# 模式 3: 单类提取 (NC → CSV)
# ===========================

# LC_Type5 → 全称列名映射 (用于单类提取)
PFT_EXTRACT_COLS = [
    "Evergreen_Needleleaf_Trees",
    "Evergreen_Broadleaf_Trees",
    "Deciduous_Needleleaf_Trees",
    "Deciduous_Broadleaf_Trees",
    "Shrub",
    "Grass",
    "Cereal_Crops",
    "Broadleaf_Crops",
    "Urban_and_Builtup",
    "Snow_and_Ice",
    "Barren_or_Sparse_Vegetation",
    "Water",
]

CODE2COL = {
    0: "Water",
    1: "Evergreen_Needleleaf_Trees",
    2: "Evergreen_Broadleaf_Trees",
    3: "Deciduous_Needleleaf_Trees",
    4: "Deciduous_Broadleaf_Trees",
    5: "Shrub",
    6: "Grass",
    7: "Cereal_Crops",
    8: "Broadleaf_Crops",
    9: "Urban_and_Builtup",
    10: "Snow_and_Ice",
    11: "Barren_or_Sparse_Vegetation",
}


def extract_single_class_pft(
    input_nc: str,
    output_csv: str,
) -> pd.DataFrame:
    """从最近邻 NC 结果中提取单类 CSV（每格点仅一类=100%，其余=0）。

    参数:
        input_nc: 最近邻 NC 文件路径（含 LC_Type5 变量）
        output_csv: 输出 CSV 路径

    返回:
        df: 提取后的 DataFrame
    """
    print(f"[INFO] Reading: {input_nc}")
    ds = xr.open_dataset(input_nc)
    lc = ds["LC_Type5"].values  # (ny, nx)
    ny, nx = lc.shape
    print(f"[INFO] Shape: {ny} x {nx} -> total {ny * nx} grids")

    rows = []
    cellid = 1

    for j in range(ny):
        for i in range(nx):
            val = int(lc[j, i]) if np.isfinite(lc[j, i]) else -1
            row = {c: 0 for c in PFT_EXTRACT_COLS}
            row["CELLID"] = cellid
            row["ICELL"] = i + 1
            row["JCELL"] = j + 1

            if val in CODE2COL:
                row[CODE2COL[val]] = 100

            rows.append(row)
            cellid += 1

    df = pd.DataFrame(rows, columns=["CELLID", "ICELL", "JCELL"] + PFT_EXTRACT_COLS).astype(int)
    assert len(df) == ny * nx, f"Row count ({len(df)}) != grid count ({ny * nx})"

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    df.to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"[OK] Wrote: {output_csv}  ({len(df)} rows)")
    return df
