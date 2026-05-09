#!/usr/bin/env python3
"""
Run_EFP.py - MEGAN 排放因子处理
=================================
输入: Data/Input/EFP/ (排放因子表、物种组成、生态类型)
      Data/GrowthForm/{year}/ (Growth Form CSV)
输出: Data/EFP/Output/ (OutputGridEF CSV)
      Data/EFP/Database/ (SQLite 数据库)

使用方式:
    python Run_EFP.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Core_EFP import run_efp_pipeline

# ====== 配置 ======
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# 场景配置
SCEN_NAME = "GD_cn27"

# 输入目录 (EFP 输入 CSV)
INPUT_DIR = os.path.join(BASE_DIR, "Data", "Input", "EFP")

# 输出目录
OUTPUT_DIR = os.path.join(BASE_DIR, "Data", "EFP", "Output")
DATABASE_DIR = os.path.join(BASE_DIR, "Data", "EFP", "Database")

# Growth Form 和 Ecotype 输入 (如果不在 Data/Input/EFP/ 下)
# grid_growth_form 也可从 Data/GrowthForm/ 复制/软链到 Data/Input/EFP/

# ====== 执行 ======
if __name__ == "__main__":
    output_csv = run_efp_pipeline(
        scen_name=SCEN_NAME,
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        database_dir=DATABASE_DIR,
    )
    print(f"\nEFP completed. Output: {output_csv}")
