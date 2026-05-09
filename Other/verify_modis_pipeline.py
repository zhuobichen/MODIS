#!/usr/bin/env python3
"""
Quality control checks for MODIS → MEGAN pipeline outputs.

Checks implemented:
1. PFT vs Growth Form consistency:
   - Sum vegetation PFT fractions (tree components + shrub + grass + crops)
   - Compare with Growth Form (TreeFrac + CropFrac + ShrubFrac + HerbFrac)
   - Report difference stats and flag grids with |diff| > tolerance.
2. Growth Form internal sanity:
   - Each fraction in [0, 100]; total ≤ 100 + tolerance.
3. LAI missing ratio:
   - Treat value == 0 as potential missing (can refine later).
   - Compute % zeros per grid, % grids all zero, overall zero percentage.
4. EF columns extremes:
   - Basic stats (min, p01, mean, p99, max, std, zero% per EF)
   - Flag columns with excessive zeros (> zero_threshold) or high outliers (max > mean + high_std_multiplier*std).
5. EF per-grid nonzero interpretation:
   - Count grids with all EF == 0.
   - For grids all-zero EF but non-zero GrowthForm vegetation -> suspect mismatch.
   - For grids vegetation sum ≈ 0 (e.g., water/urban) and EF all zero -> acceptable.

Adjust CONFIG constants below as necessary.
"""
from __future__ import annotations
import pandas as pd
import numpy as np
from pathlib import Path

# ------------------ CONFIG ------------------
# Paths (adjust if different scenario/resolution)
PFT_FRAC_CSV = Path("PFT_frac_2000_27km_square.csv")  # neighbor square percentages
GROWTH_FORM_CSV = Path("grid_growth_form_cn27.csv")
LAI_CSV = Path("LAI3_cn27.csv")
EF_OUTPUT_CSV = Path("MEGEFP32/outputs/OutputGridEF.GD_cn27.csv")  # example location
REPORT_OUT = Path("MODIS_Pipeline_QC_Report.txt")

# Numerical tolerances & thresholds
PFT_GROWTH_DIFF_TOL = 1.0          # percent points tolerance for vegetation sum difference
GROWTH_TOTAL_TOL = 2.0             # total growth form may exceed vegetation sum by ≤2 due to rounding
ZERO_THRESHOLD_EF = 0.5            # >50% zeros across a column triggers flag
HIGH_STD_MULTIPLIER = 5.0          # max > mean + k*std triggers high value flag
LAI_ALL_ZERO_THRESHOLD = 0.0       # treat exactly all zero as missing
LAI_ZERO_RATIO_GRID_FLAG = 0.9     # if a grid has >90% zeros, flag

# Column name patterns (flexible matching)
TREE_PATTERNS = ["Evergreen", "Deciduous"]  # needleleaf/broadleaf variants
SHRUB_PATTERN = "Shrub"
GRASS_PATTERN = "Grass"
CROP_PATTERNS = ["Cereal", "Broad-leaf", "Crop"]
EXCLUDE_PFT_PATTERNS = ["Water", "Urban", "Barren", "Snow", "Ice"]

# -------------------------------------------------

def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def identify_pft_columns(df: pd.DataFrame) -> dict:
    """Classify columns into tree, shrub, grass, crops, excluded."""
    cols = df.columns
    tree_cols = [c for c in cols if any(p in c for p in TREE_PATTERNS) and "Crop" not in c]
    shrub_cols = [c for c in cols if SHRUB_PATTERN in c]
    grass_cols = [c for c in cols if GRASS_PATTERN in c]
    crop_cols = [c for c in cols if any(p in c for p in CROP_PATTERNS)]
    exclude_cols = [c for c in cols if any(p in c for p in EXCLUDE_PFT_PATTERNS)]
    # Remove overlaps
    def unique(lst):
        uniq = []
        for x in lst:
            if x not in uniq:
                uniq.append(x)
        return uniq
    tree_cols = unique(tree_cols)
    return {
        "tree": tree_cols,
        "shrub": shrub_cols,
        "grass": grass_cols,
        "crop": crop_cols,
        "exclude": exclude_cols,
    }


def compute_pft_growth_consistency(pft_df: pd.DataFrame, gf_df: pd.DataFrame) -> pd.DataFrame:
    mapping = identify_pft_columns(pft_df)
    # vegetation sum from PFT detailed
    pft_df["PFT_VegSum"] = (
        pft_df[mapping["tree"]].sum(axis=1)
        + pft_df[mapping["shrub"]].sum(axis=1)
        + pft_df[mapping["grass"]].sum(axis=1)
        + pft_df[mapping["crop"]].sum(axis=1)
    )
    # Growth form sum
    required_gf_cols = ["TreeFrac", "CropFrac", "ShrubFrac", "HerbFrac"]
    for c in required_gf_cols:
        if c not in gf_df.columns:
            raise KeyError(f"Missing growth form column: {c}")
    gf_df["GF_VegSum"] = gf_df[required_gf_cols].sum(axis=1)
    merged = gf_df.merge(pft_df[["CELLID", "PFT_VegSum"]], left_on="gridID", right_on="CELLID", how="left")
    merged["VegDiff"] = merged["GF_VegSum"] - merged["PFT_VegSum"]
    return merged


def lai_missing_stats(lai_df: pd.DataFrame) -> dict:
    lai_cols = [c for c in lai_df.columns if c.startswith("LAI")]
    lai_vals = lai_df[lai_cols].values
    zero_mask = (lai_vals == 0)
    total_entries = lai_vals.size
    total_zeros = zero_mask.sum()
    overall_zero_ratio = total_zeros / total_entries if total_entries > 0 else np.nan
    grid_all_zero = zero_mask.all(axis=1)
    grids_all_zero_ratio = grid_all_zero.mean()
    grid_zero_ratio = zero_mask.sum(axis=1) / lai_vals.shape[1]
    high_zero_grids_ratio = (grid_zero_ratio > LAI_ZERO_RATIO_GRID_FLAG).mean()
    return {
        "total_entries": total_entries,
        "total_zeros": int(total_zeros),
        "overall_zero_ratio": overall_zero_ratio,
        "grids_all_zero_ratio": grids_all_zero_ratio,
        "high_zero_grids_ratio": high_zero_grids_ratio,
    }


def ef_extreme_stats(ef_df: pd.DataFrame) -> pd.DataFrame:
    ef_cols = [c for c in ef_df.columns if c.startswith("EF")]
    stats_rows = []
    for col in ef_cols:
        series = ef_df[col]
        arr = series.values
        zero_ratio = (arr == 0).mean()
        mean = arr.mean()
        std = arr.std(ddof=0)
        p01 = np.nanpercentile(arr, 1)
        p99 = np.nanpercentile(arr, 99)
        maxv = arr.max()
        minv = arr.min()
        high_flag = False
        if std > 0 and maxv > mean + HIGH_STD_MULTIPLIER * std:
            high_flag = True
        zero_flag = zero_ratio > ZERO_THRESHOLD_EF
        stats_rows.append({
            "EF": col,
            "min": minv,
            "p01": p01,
            "mean": mean,
            "std": std,
            "p99": p99,
            "max": maxv,
            "zero_ratio": zero_ratio,
            "flag_zero_excess": zero_flag,
            "flag_high_outlier": high_flag,
        })
    return pd.DataFrame(stats_rows)


def interpret_ef_grid_zero(ef_df: pd.DataFrame, gf_df: pd.DataFrame) -> dict:
    ef_cols = [c for c in ef_df.columns if c.startswith("EF")]
    ef_values = ef_df[ef_cols].values
    grid_all_zero_mask = (ef_values == 0).all(axis=1)
    # vegetation sum from growth form
    gf_sum = gf_df[["TreeFrac", "CropFrac", "ShrubFrac", "HerbFrac"]].sum(axis=1)
    # Acceptable zero if vegetation sum very small (< 1%)
    acceptable_zero = (gf_sum < 1) & grid_all_zero_mask
    suspect_zero = (gf_sum >= 1) & grid_all_zero_mask
    return {
        "total_grids": len(gf_df),
        "grids_all_zero_EF": int(grid_all_zero_mask.sum()),
        "acceptable_zero_grids": int(acceptable_zero.sum()),
        "suspect_zero_grids": int(suspect_zero.sum()),
    }


def main():
    pft_df = load_csv(PFT_FRAC_CSV)
    gf_df = load_csv(GROWTH_FORM_CSV)
    lai_df = load_csv(LAI_CSV)
    ef_df = load_csv(EF_OUTPUT_CSV)

    # 1. PFT vs Growth Form consistency
    consistency_df = compute_pft_growth_consistency(pft_df, gf_df)
    diff_stats = {
        "diff_mean": consistency_df["VegDiff"].mean(),
        "diff_std": consistency_df["VegDiff"].std(ddof=0),
        "diff_min": consistency_df["VegDiff"].min(),
        "diff_max": consistency_df["VegDiff"].max(),
        "abs_diff_gt_tol": int((consistency_df["VegDiff"].abs() > PFT_GROWTH_DIFF_TOL).sum()),
    }

    # 2. Growth Form total sanity
    gf_total = consistency_df["GF_VegSum"]
    gf_total_flag = int((gf_total > (100 + GROWTH_TOTAL_TOL)).sum())

    # 3. LAI missing stats
    lai_stats = lai_missing_stats(lai_df)

    # 4. EF extremes
    ef_stats_df = ef_extreme_stats(ef_df)

    # 5. EF per-grid zero interpretation
    ef_zero_interpret = interpret_ef_grid_zero(ef_df, gf_df)

    # Compose report
    lines = []
    lines.append("MODIS → MEGAN Pipeline QC Report")
    lines.append("=================================================")
    lines.append("1. PFT vs Growth Form Consistency:")
    for k, v in diff_stats.items():
        lines.append(f"   {k}: {v}")
    lines.append(f"   GrowthForm total > 100+tol count: {gf_total_flag}")
    lines.append("")
    lines.append("2. LAI Missing (Zero) Statistics:")
    for k, v in lai_stats.items():
        lines.append(f"   {k}: {v}")
    lines.append("")
    lines.append("3. EF Column Extremes:")
    for _, row in ef_stats_df.iterrows():
        lines.append(
            f"   {row['EF']}: min={row['min']:.4g} p01={row['p01']:.4g} mean={row['mean']:.4g} std={row['std']:.4g} p99={row['p99']:.4g} max={row['max']:.4g} zero%={row['zero_ratio']*100:.1f}% zeroFlag={row['flag_zero_excess']} highFlag={row['flag_high_outlier']}"
        )
    lines.append("")
    lines.append("4. EF Grid Zero Interpretation:")
    for k, v in ef_zero_interpret.items():
        lines.append(f"   {k}: {v}")
    lines.append("")
    lines.append("Interpretation Guide:")
    lines.append(" - abs_diff_gt_tol > 0: Investigate rounding or missing PFT categories.")
    lines.append(" - GrowthForm total > 100+tol: Possible double counting or lack of exclusion of non-vegetation.")
    lines.append(" - High LAI zero ratios may signal missing data rather than phenology if concentrated spatially.")
    lines.append(" - EF columns with zero% > threshold: Check ecotype/speciation coverage; zeros may be expected for non-emitters.")
    lines.append(" - Max > mean + k*std: Potential outlier or unit mismatch; inspect source EF entries.")
    lines.append(" - Suspect zero EF grids (vegetation present): Likely ecotype or growth form ID mismatch.")

    REPORT_OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"QC report written to {REPORT_OUT}")

if __name__ == "__main__":
    main()
