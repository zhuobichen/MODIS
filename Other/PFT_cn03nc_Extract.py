#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
读取 MODIS LC_Type5 最近邻结果并输出为 CSV
- 每格点一行（不跳过任何网格）
- 所属类别列 = 100，其余类别列 = 0
- 列名为全称
- 输出列：CELLID, ICELL, JCELL, 各类PFT
"""

import xarray as xr
import pandas as pd
import numpy as np

# ---------- 参数 ----------
infile  = "MODIS_LCType5_2000_NN_cn03.nc"
outfile = "pft_cn03_GuangDong_2000.csv"

# LC_Type5 → 全称列名
pft_cols = [
    "Evergreen_Needleleaf_Trees",
    "Evergreen_Broadleaf_Trees",
    "Deciduous_Needleleaf_Trees",
    "Deciduous_Broadleaf_Trees",
    "Shrub",
    "Grass",
    "Cereal_Crops",
    "Broadleaf_Crops",
    "Urban_and_Builtup",
    "Snow_and_Ice",
    "Barren_or_Sparse_Vegetation",
    "Water",
]

# 分类编号对应列名
code2col = {
    0:  "Water",
    1:  "Evergreen_Needleleaf_Trees",
    2:  "Evergreen_Broadleaf_Trees",
    3:  "Deciduous_Needleleaf_Trees",
    4:  "Deciduous_Broadleaf_Trees",
    5:  "Shrub",
    6:  "Grass",
    7:  "Cereal_Crops",
    8:  "Broadleaf_Crops",
    9:  "Urban_and_Builtup",
    10: "Snow_and_Ice",
    11: "Barren_or_Sparse_Vegetation",
}

print(f"[INFO] 读取: {infile}")
ds = xr.open_dataset(infile)
lc = ds["LC_Type5"].values  # (ny, nx)
ny, nx = lc.shape
print(f"[INFO] 维度: {ny} x {nx} -> 总格点 {ny*nx}")

# 构建输出
rows = []
cellid = 1

for j in range(ny):          # 行（y）
    for i in range(nx):      # 列（x）
        val = int(lc[j, i]) if np.isfinite(lc[j, i]) else -1

        # 初始化行：所有类别=0
        row = {c: 0 for c in pft_cols}
        row["CELLID"] = cellid
        row["ICELL"]  = i + 1  # 列号（1-based）
        row["JCELL"]  = j + 1  # 行号（1-based）

        # 合法类别(0-11)：对应列置100
        if val in code2col:
            row[code2col[val]] = 100

        rows.append(row)
        cellid += 1

df = pd.DataFrame(rows, columns=["CELLID","ICELL","JCELL"] + pft_cols).astype(int)

# 自检：行数必须等于 nx*ny
assert len(df) == ny * nx, f"导出行数({len(df)}) != 网格数({ny*nx})"

# 输出
df.to_csv(outfile, index=False, encoding="utf-8-sig")
print(f"[OK] 已写出: {outfile}  共 {len(df)} 行")
