#!/usr/bin/env python3
"""
Run_PFT_Nearest.py - PFT 最近邻单类提取
=========================================
输入: Data/Input/HDF/ (MCD12Q1 .hdf 瓦片)
输出: Data/PFT/{year}/ (NetCDF)

使用方式:
    python Run_PFT_Nearest.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Core_PFT import run_nearest_neighbor_pipeline, extract_single_class_pft

# ====== 配置 ======
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(BASE_DIR, "Data", "Input", "HDF")
OUTPUT_DIR = os.path.join(BASE_DIR, "Data", "PFT", "2000")
GRID_NC = os.path.join(BASE_DIR, "GRIDCRO2D_2000121_GuangDongD1")

YEAR = 2000
HLINES = [23, 24, 25, 26, 27, 28, 29, 30]
VLINES = [3, 4, 5, 6, 7, 8]

SKIP_TILES = [
    "h30v03", "h29v03", "h28v03", "h27v03",
    "h30v04", "h29v04",
    "h30v05", "h23v05", "h24v05",
    "h23v08", "h24v08", "h25v08", "h26v08",
    "h23v07", "h24v07", "h25v07", "h26v07",
    "h24v06", "h25v06", "h30v06",
]

# ====== 执行 ======
if __name__ == "__main__":
    # 步骤 1: 最近邻提取 (正弦投影 KDTree)
    result = run_nearest_neighbor_pipeline(
        year=YEAR,
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        grid_nc=GRID_NC,
        hlines=HLINES,
        vlines=VLINES,
        skip_tiles=SKIP_TILES,
        use_sinusoidal=True,
    )
    print(f"Nearest neighbor NC: {result['nc']}")

    # 步骤 2: 单类提取 (NC → CSV)
    output_csv = os.path.join(OUTPUT_DIR, f"pft_cn27_GuangDong_{YEAR}.csv")
    extract_single_class_pft(
        input_nc=result['nc'],
        output_csv=output_csv,
    )
    print(f"Single class CSV: {output_csv}")
