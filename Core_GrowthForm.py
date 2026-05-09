"""
Core_GrowthForm.py - PFT → Growth Form 转换
=============================================
将 MODIS PFT 12 类百分比聚合为 MEGAN Growth Form 4 类:
  - TreeFrac  = Sum(常绿针叶 + 常绿阔叶 + 落叶针叶 + 落叶阔叶)
  - CropFrac  = Sum(谷类作物 + 阔叶作物)
  - ShrubFrac = 灌木
  - HerbFrac  = 草本

所有路径通过函数签名传入，无硬编码。

提供函数:
  - pft_to_growth_form(df) -> pd.DataFrame
      单文件转换: PFT百分比 DataFrame → Growth Form DataFrame
  - run_growth_form_pipeline() -> dict
      管道函数: 批量处理多个分辨率的 PFT → Growth Form
"""

import os
import pandas as pd
from typing import List, Dict, Optional


# Growth Form 聚合列映射
TREE_COLS = [
    "Evergreen Needleleaf trees",
    "Evergreen Broadleaf trees",
    "Deciduous Needleleaf trees",
    "Deciduous Broadleaf trees",
]
CROP_COLS = [
    "Cereal crops",
    "Broad-leaf crops",
]


def pft_to_growth_form(df: pd.DataFrame) -> pd.DataFrame:
    """将 PFT 12 类百分比 DataFrame 转换为 Growth Form 4 类 DataFrame。

    参数:
        df: 含 CELLID 列 + 12 个 CLASS_NAMES 列的 DataFrame

    返回:
        DataFrame: gridID, TreeFrac, CropFrac, ShrubFrac, HerbFrac
    """
    result = pd.DataFrame()
    result["gridID"] = df["CELLID"]
    result["TreeFrac"] = df[TREE_COLS].sum(axis=1)
    result["CropFrac"] = df[CROP_COLS].sum(axis=1)
    result["ShrubFrac"] = df["Shrub"]
    result["HerbFrac"] = df["Grass"]
    return result


def run_growth_form_pipeline(
    file_mapping: Dict[str, str],
    input_dir: str = "",
    output_dir: str = "",
) -> Dict[str, str]:
    """批量 PFT → Growth Form 转换管道。

    参数:
        file_mapping: {"输入文件相对路径": "输出文件相对路径"} 映射
        input_dir: 输入根目录 (可选，加在输入路径前)
        output_dir: 输出根目录 (可选，加在输出路径前)

    返回:
        dict: {"输入路径": "输出路径"} 结果映射

    示例:
        file_mapping = {
            "PFT_frac_2000_27km_square.csv": "grid_growth_form_cn27.csv",
            "PFT_frac_2000_9km_square.csv": "grid_growth_form_cn09.csv",
        }
    """
    # 检查输入文件
    missing_files = []
    for input_rel in file_mapping:
        input_path = os.path.join(input_dir, input_rel) if input_dir else input_rel
        if not os.path.exists(input_path):
            missing_files.append(input_path)

    if missing_files:
        print("Missing input files:")
        for f in missing_files:
            print(f"  {f}")
        raise FileNotFoundError(f"{len(missing_files)} input file(s) not found")

    results = {}
    for input_rel, output_rel in file_mapping.items():
        input_path = os.path.join(input_dir, input_rel) if input_dir else input_rel
        output_path = os.path.join(output_dir, output_rel) if output_dir else output_rel

        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

        print(f"\nProcessing: {input_path}")
        df_pft = pd.read_csv(input_path)
        df_gf = pft_to_growth_form(df_pft)
        df_gf.to_csv(output_path, index=False)
        print(f"  -> {output_path}  ({len(df_gf)} rows)")
        results[input_path] = output_path

    print("\nAll growth form conversions completed.")
    return results
