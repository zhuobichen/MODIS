#!/usr/bin/env python3
"""
Run_PFT_Percentage.py - PFT 正方形邻域百分比提取
===================================================
输入: Data/Input/HDF/ (MCD12Q1 .hdf 瓦片)
输出: Data/PFT/{year}/ (NetCDF + CSV)

支持多分辨率批量处理。

使用方式:
    python Run_PFT_Percentage.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Core_PFT import run_square_neighborhood_pipeline

# ====== 配置 ======
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(BASE_DIR, "Data", "Input", "HDF")
OUTPUT_BASE = os.path.join(BASE_DIR, "Data", "PFT")

YEAR = 2000

# CMAQ 网格文件 (按分辨率)
GRID_FILES = {
    "27km": os.path.join(BASE_DIR, "GRIDCRO2D_2000121_GuangDongD1"),
    "9km":  os.path.join(BASE_DIR, "GRIDCRO2D_2000121_GuangDongD2"),
    "3km":  os.path.join(BASE_DIR, "GRIDCRO2D_2000121_GuangDongD3"),
}

# MODIS 瓦片范围 (中国广东)
HLINES = [23, 24, 25, 26, 27, 28, 29, 30]
VLINES = [3, 4, 5, 6, 7, 8]

# 跳过的瓦片 (边界外/数据质量问题)
SKIP_TILES = [
    "h30v03", "h29v03", "h28v03", "h27v03",
    "h30v04", "h29v04",
    "h30v05", "h23v05", "h24v05",
    "h23v08", "h24v08", "h25v08", "h26v08",
    "h23v07", "h24v07", "h25v07", "h26v07",
    "h24v06", "h25v06", "h30v06",
]

# 要处理的分辨率列表
RESOLUTIONS = ["27km", "9km", "3km"]

# ====== 执行 ======
if __name__ == "__main__":
    for res in RESOLUTIONS:
        print(f"\n{'='*60}")
        print(f"Processing resolution: {res}")
        print(f"{'='*60}")

        grid_nc = GRID_FILES[res]
        output_dir = os.path.join(OUTPUT_BASE, str(YEAR))

        result = run_square_neighborhood_pipeline(
            year=YEAR,
            resolution=res,
            input_dir=INPUT_DIR,
            output_dir=output_dir,
            grid_nc=grid_nc,
            hlines=HLINES,
            vlines=VLINES,
            skip_tiles=SKIP_TILES,
            parallel=True,
        )
        print(f"  -> NC: {result['nc']}")
        print(f"  -> CSV: {result['csv']}")

    print("\nAll resolutions processed.")
