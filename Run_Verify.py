#!/usr/bin/env python3
"""
Run_Verify.py - MODIS → MEGAN 全链路质控校验
===============================================
输入: Data/PFT/、Data/GrowthForm/、Data/LAI/、Data/EFP/Output/
输出: Picture/QC/MODIS_Pipeline_QC_Report.txt

使用方式:
    python Run_Verify.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Core_Verification import run_verification_pipeline

# ====== 配置 ======
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
YEAR = 2000

# 输入文件 (根据实际输出调整年份和文件名)
PFT_CSV = os.path.join(BASE_DIR, "Data", "PFT", str(YEAR),
                       f"PFT_frac_{YEAR}_27km_square.csv")
GF_CSV = os.path.join(BASE_DIR, "Data", "GrowthForm", str(YEAR),
                      "grid_growth_form_cn27.csv")
LAI_CSV = os.path.join(BASE_DIR, "Data", "LAI", str(YEAR),
                       "LAI3_cn27.csv")
EF_CSV = os.path.join(BASE_DIR, "Data", "EFP", "Output",
                      "OutputGridEF.GD_cn27.csv")

OUTPUT_REPORT = os.path.join(BASE_DIR, "Picture", "QC",
                             "MODIS_Pipeline_QC_Report.txt")

# ====== 执行 ======
if __name__ == "__main__":
    report_path = run_verification_pipeline(
        pft_csv=PFT_CSV,
        growth_form_csv=GF_CSV,
        lai_csv=LAI_CSV,
        ef_csv=EF_CSV,
        output_report=OUTPUT_REPORT,
    )
    print(f"\nQC report: {report_path}")
