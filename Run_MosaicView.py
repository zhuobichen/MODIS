#!/usr/bin/env python3
"""
Run_MosaicView.py - MODIS 瓦片拼接可视化
==========================================
输入: Data/Input/HDF/ (MCD12Q1 / MCD15A2H .hdf 瓦片)
输出: Picture/MosaicView/ (PNG 拼接图)

使用方式:
    python Run_MosaicView.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Core_MosaicView import run_mosaic_view_pipeline

# ====== 配置 ======
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(BASE_DIR, "Data", "Input", "HDF")
OUTPUT_DIR = os.path.join(BASE_DIR, "Picture", "MosaicView")

YEAR = 2000

# PFT 拼接参数
HLINES_PFT = list(range(1, 36))  # 全球范围
VLINES_PFT = list(range(1, 11))

# LAI 拼接参数 (广东区域)
HLINES_LAI = [23, 24, 25, 26, 27, 28, 29, 30]
VLINES_LAI = [3, 4, 5, 6, 7, 8]

# ====== 执行 ======
if __name__ == "__main__":
    # 1. PFT 瓦片拼接
    print("=" * 60)
    print("Mosaicking PFT tiles...")
    print("=" * 60)
    run_mosaic_view_pipeline(
        year=YEAR,
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        hlines=HLINES_PFT,
        vlines=VLINES_PFT,
        data_type="PFT",
    )

    # 2. LAI 瓦片拼接 (示例: DOY=361)
    print("\n" + "=" * 60)
    print("Mosaicking LAI tiles (DOY=361)...")
    print("=" * 60)
    run_mosaic_view_pipeline(
        year=YEAR,
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        hlines=HLINES_LAI,
        vlines=VLINES_LAI,
        data_type="LAI",
        doy=361,
    )

    print("\nMosaic view completed.")
