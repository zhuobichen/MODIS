#!/usr/bin/env python3
"""
Run_LAI.py - LAI 时序拼接
==========================
输入: Data/Input/LAI/{resolution}_{project}_{year}/ (NetCDF 切片)
输出: Data/LAI/{year}/ (CSV)

支持多分辨率批量处理。

使用方式:
    python Run_LAI.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Core_LAI import run_lai_pipeline

# ====== 配置 ======
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
YEAR = "2000"
PROJECT = "GuangDong"

# 各分辨率参数
RESOLUTION_CONFIGS = [
    {"resolution": "cn27", "rows": 175, "cols": 124, "use_1based": False},
    {"resolution": "cn09", "rows": 157, "cols": 199, "use_1based": True},
    {"resolution": "cn03", "rows": 175, "cols": 124, "use_1based": True},
]

# ====== 执行 ======
if __name__ == "__main__":
    for cfg in RESOLUTION_CONFIGS:
        res = cfg["resolution"]
        input_dir = os.path.join(
            BASE_DIR, "Data", "Input", "LAI",
            f"{res}_{PROJECT}_{YEAR}",
        )
        output_dir = os.path.join(BASE_DIR, "Data", "LAI", YEAR)

        print(f"\n{'='*60}")
        print(f"Processing LAI: {res} ({cfg['rows']}x{cfg['cols']})")
        print(f"{'='*60}")

        outfile = run_lai_pipeline(
            year=YEAR,
            resolution=res,
            project=PROJECT,
            input_dir=input_dir,
            output_dir=output_dir,
            rows=cfg["rows"],
            cols=cfg["cols"],
            use_1based_coords=cfg["use_1based"],
        )
        print(f"  -> {outfile}")

    print("\nAll LAI resolutions processed.")
