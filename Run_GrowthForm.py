#!/usr/bin/env python3
"""
Run_GrowthForm.py - PFT → Growth Form 转换
=============================================
输入: Data/PFT/{year}/ (PFT 百分比 CSV)
输出: Data/GrowthForm/{year}/ (Growth Form CSV)

使用方式:
    python Run_GrowthForm.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Core_GrowthForm import run_growth_form_pipeline

# ====== 配置 ======
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
YEAR = 2000
INPUT_DIR = os.path.join(BASE_DIR, "Data", "PFT", str(YEAR))
OUTPUT_DIR = os.path.join(BASE_DIR, "Data", "GrowthForm", str(YEAR))

# 输入 → 输出映射
FILE_MAPPING = {
    f"PFT_frac_{YEAR}_27km_square.csv": "grid_growth_form_cn27.csv",
    f"PFT_frac_{YEAR}_9km_square.csv":  "grid_growth_form_cn09.csv",
    f"PFT_frac_{YEAR}_3km_square.csv":  "grid_growth_form_cn03.csv",
}

# ====== 执行 ======
if __name__ == "__main__":
    results = run_growth_form_pipeline(
        file_mapping=FILE_MAPPING,
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
    )
    print("\nResults:")
    for inp, out in results.items():
        print(f"  {os.path.basename(inp)} -> {os.path.basename(out)}")
