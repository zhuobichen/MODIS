"""
Core_Verification.py - MODIS → MEGAN 全链路质控
=================================================
对 MODIS → MEGAN 管道各阶段输出进行质量校验。

检查项:
  1. PFT vs Growth Form 一致性
  2. Growth Form 内部合理性
  3. LAI 缺失统计
  4. EF 列极值检测
  5. EF 零值网格解释

所有路径通过函数签名传入，无硬编码。

提供函数:
  - run_verification_pipeline() -> str
      完整管道: 加载数据 → 所有检查 → 输出报告
"""

from __future__ import annotations
import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional

# 默认阈值
DEFAULT_CONFIG = {
    "pft_growth_diff_tol": 1.0,
    "growth_total_tol": 2.0,
    "zero_threshold_ef": 0.5,
    "high_std_multiplier": 5.0,
    "lai_zero_ratio_grid_flag": 0.9,
}

# PFT 列分类模式
TREE_PATTERNS = ["Evergreen", "Deciduous"]
SHRUB_PATTERN = "Shrub"
GRASS_PATTERN = "Grass"
CROP_PATTERNS = ["Cereal", "Broad-leaf", "Crop"]
EXCLUDE_PFT_PATTERNS = ["Water", "Urban", "Barren", "Snow", "Ice"]


def _identify_pft_columns(df: pd.DataFrame) -> dict:
    """分类 PFT 列: tree, shrub, grass, crops, exclude。"""
    cols = df.columns
    tree_cols = [c for c in cols if any(p in c for p in TREE_PATTERNS) and "Crop" not in c]
    shrub_cols = [c for c in cols if SHRUB_PATTERN in c]
    grass_cols = [c for c in cols if GRASS_PATTERN in c]
    crop_cols = [c for c in cols if any(p in c for p in CROP_PATTERNS)]
    exclude_cols = [c for c in cols if any(p in c for p in EXCLUDE_PFT_PATTERNS)]

    def _unique(lst):
        uniq = []
        for x in lst:
            if x not in uniq:
                uniq.append(x)
        return uniq

    return {
        "tree": _unique(tree_cols),
        "shrub": shrub_cols,
        "grass": grass_cols,
        "crop": crop_cols,
        "exclude": exclude_cols,
    }


def _compute_pft_growth_consistency(pft_df: pd.DataFrame, gf_df: pd.DataFrame) -> pd.DataFrame:
    """计算 PFT 植被总和与 Growth Form 总和差异。"""
    mapping = _identify_pft_columns(pft_df)
    pft_df["PFT_VegSum"] = (
        pft_df[mapping["tree"]].sum(axis=1)
        + pft_df[mapping["shrub"]].sum(axis=1)
        + pft_df[mapping["grass"]].sum(axis=1)
        + pft_df[mapping["crop"]].sum(axis=1)
    )
    required_gf_cols = ["TreeFrac", "CropFrac", "ShrubFrac", "HerbFrac"]
    for c in required_gf_cols:
        if c not in gf_df.columns:
            raise KeyError(f"Missing growth form column: {c}")
    gf_df["GF_VegSum"] = gf_df[required_gf_cols].sum(axis=1)
    merged = gf_df.merge(
        pft_df[["CELLID", "PFT_VegSum"]],
        left_on="gridID", right_on="CELLID", how="left"
    )
    merged["VegDiff"] = merged["GF_VegSum"] - merged["PFT_VegSum"]
    return merged


def _lai_missing_stats(lai_df: pd.DataFrame, flag_ratio: float = 0.9) -> dict:
    """LAI 缺失统计。"""
    lai_cols = [c for c in lai_df.columns if c.startswith("LAI")]
    lai_vals = lai_df[lai_cols].values
    zero_mask = (lai_vals == 0)
    total_entries = lai_vals.size
    total_zeros = int(zero_mask.sum())
    overall_zero_ratio = total_zeros / total_entries if total_entries > 0 else np.nan
    grid_all_zero = zero_mask.all(axis=1)
    grids_all_zero_ratio = float(grid_all_zero.mean())
    grid_zero_ratio = zero_mask.sum(axis=1) / lai_vals.shape[1]
    high_zero_grids_ratio = float((grid_zero_ratio > flag_ratio).mean())
    return {
        "total_entries": total_entries,
        "total_zeros": total_zeros,
        "overall_zero_ratio": overall_zero_ratio,
        "grids_all_zero_ratio": grids_all_zero_ratio,
        "high_zero_grids_ratio": high_zero_grids_ratio,
    }


def _ef_extreme_stats(ef_df: pd.DataFrame, zero_threshold: float = 0.5,
                      high_std_multiplier: float = 5.0) -> pd.DataFrame:
    """EF 列极值统计。"""
    ef_cols = [c for c in ef_df.columns if c.startswith("EF")]
    stats_rows = []
    for col in ef_cols:
        arr = ef_df[col].values
        zero_ratio = float((arr == 0).mean())
        mean = float(arr.mean())
        std = float(arr.std(ddof=0))
        p01 = float(np.nanpercentile(arr, 1))
        p99 = float(np.nanpercentile(arr, 99))
        maxv = float(arr.max())
        minv = float(arr.min())
        high_flag = (std > 0 and maxv > mean + high_std_multiplier * std)
        zero_flag = (zero_ratio > zero_threshold)
        stats_rows.append({
            "EF": col, "min": minv, "p01": p01, "mean": mean,
            "std": std, "p99": p99, "max": maxv,
            "zero_ratio": zero_ratio,
            "flag_zero_excess": zero_flag,
            "flag_high_outlier": high_flag,
        })
    return pd.DataFrame(stats_rows)


def _interpret_ef_grid_zero(ef_df: pd.DataFrame, gf_df: pd.DataFrame) -> dict:
    """解释 EF 零值网格。"""
    ef_cols = [c for c in ef_df.columns if c.startswith("EF")]
    ef_values = ef_df[ef_cols].values
    grid_all_zero_mask = (ef_values == 0).all(axis=1)
    gf_sum = gf_df[["TreeFrac", "CropFrac", "ShrubFrac", "HerbFrac"]].sum(axis=1)
    acceptable_zero = (gf_sum < 1) & grid_all_zero_mask
    suspect_zero = (gf_sum >= 1) & grid_all_zero_mask
    return {
        "total_grids": len(gf_df),
        "grids_all_zero_EF": int(grid_all_zero_mask.sum()),
        "acceptable_zero_grids": int(acceptable_zero.sum()),
        "suspect_zero_grids": int(suspect_zero.sum()),
    }


def run_verification_pipeline(
    pft_csv: str,
    growth_form_csv: str,
    lai_csv: str,
    ef_csv: str,
    output_report: str,
    config: Optional[Dict] = None,
) -> str:
    """全链路 QC 校验管道。

    参数:
        pft_csv: PFT 百分比 CSV 路径
        growth_form_csv: Growth Form CSV 路径
        lai_csv: LAI CSV 路径
        ef_csv: EF 输出 CSV 路径
        output_report: 报告输出 TXT 路径
        config: 阈值配置字典 (可选, 覆盖默认值)

    返回:
        str: 报告文件路径
    """
    cfg = {**DEFAULT_CONFIG, **(config or {})}

    # 加载数据
    pft_df = pd.read_csv(pft_csv)
    gf_df = pd.read_csv(growth_form_csv)
    lai_df = pd.read_csv(lai_csv)
    ef_df = pd.read_csv(ef_csv)

    lines = []
    lines.append("MODIS -> MEGAN Pipeline QC Report")
    lines.append("=" * 50)

    # 1. PFT vs Growth Form 一致性
    lines.append("\n1. PFT vs Growth Form Consistency:")
    consistency_df = _compute_pft_growth_consistency(pft_df, gf_df)
    diff_stats = {
        "diff_mean": consistency_df["VegDiff"].mean(),
        "diff_std": consistency_df["VegDiff"].std(ddof=0),
        "diff_min": consistency_df["VegDiff"].min(),
        "diff_max": consistency_df["VegDiff"].max(),
        "abs_diff_gt_tol": int(
            (consistency_df["VegDiff"].abs() > cfg["pft_growth_diff_tol"]).sum()
        ),
    }
    for k, v in diff_stats.items():
        lines.append(f"   {k}: {v}")
    gf_total = consistency_df["GF_VegSum"]
    gf_total_flag = int((gf_total > (100 + cfg["growth_total_tol"])).sum())
    lines.append(f"   GrowthForm total > 100+tol count: {gf_total_flag}")

    # 2. LAI 缺失统计
    lines.append("\n2. LAI Missing (Zero) Statistics:")
    lai_stats = _lai_missing_stats(lai_df, cfg["lai_zero_ratio_grid_flag"])
    for k, v in lai_stats.items():
        lines.append(f"   {k}: {v}")

    # 3. EF 列极值
    lines.append("\n3. EF Column Extremes:")
    ef_stats_df = _ef_extreme_stats(
        ef_df, cfg["zero_threshold_ef"], cfg["high_std_multiplier"]
    )
    for _, row in ef_stats_df.iterrows():
        lines.append(
            f"   {row['EF']}: min={row['min']:.4g} p01={row['p01']:.4g} "
            f"mean={row['mean']:.4g} std={row['std']:.4g} "
            f"p99={row['p99']:.4g} max={row['max']:.4g} "
            f"zero%={row['zero_ratio']*100:.1f}% "
            f"zeroFlag={row['flag_zero_excess']} highFlag={row['flag_high_outlier']}"
        )

    # 4. EF 零值网格解释
    lines.append("\n4. EF Grid Zero Interpretation:")
    ef_zero_interpret = _interpret_ef_grid_zero(ef_df, gf_df)
    for k, v in ef_zero_interpret.items():
        lines.append(f"   {k}: {v}")

    # 解读指南
    lines.append("\n" + "=" * 50)
    lines.append("Interpretation Guide:")
    lines.append(" - abs_diff_gt_tol > 0: rounding or missing PFT categories.")
    lines.append(" - GrowthForm total > 100+tol: double counting or non-vegetation not excluded.")
    lines.append(" - High LAI zero ratio: may signal missing data rather than phenology.")
    lines.append(" - EF columns zero% > threshold: check ecotype/speciation coverage.")
    lines.append(" - Max > mean + k*std: potential outlier or unit mismatch.")
    lines.append(" - Suspect zero EF grids (vegetation present): ecotype/growth form ID mismatch.")

    os.makedirs(os.path.dirname(output_report) or ".", exist_ok=True)
    Path(output_report).write_text("\n".join(lines), encoding="utf-8")
    print(f"QC report written to {output_report}")
    return output_report
