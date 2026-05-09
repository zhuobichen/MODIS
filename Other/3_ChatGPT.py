# -*- coding: utf-8 -*-
"""
任务 3：12 个土地利用变量之间的关系分析与可视化

说明：
- 读取 NetCDF 文件 PFT_frac_2000_3km_square.nc
- 整理 12 个土地利用百分比变量，展平成 (像元数, 12)
- 计算 Pearson 相关系数矩阵并绘制热力图
- 选取代表性的变量对绘制散点图
- 做 PCA 分析，绘制 PC1-PC2 散点图
- 最后根据相关性和 PCA 结果自动生成 5–10 句中文总结并打印
"""

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# -----------------------------
# 1. 数据读取与整理
# -----------------------------
# 请根据实际路径修改文件名
nc_file = "PFT_frac_2000_3km_square.nc"

# 打开 NetCDF 数据集
ds = xr.open_dataset(nc_file)

# 12 个土地利用百分比变量名
# 注意：这里的字符串要和 NetCDF 里的变量名完全一致
pft_vars = [
    "Water",
    "Evergreen Needleleaf trees",
    "Evergreen Broadleaf trees",
    "Deciduous Needleleaf trees",
    "Deciduous Broadleaf trees",
    "Shrub",
    "Grass",
    "Cereal crops",
    "Broad-leaf crops",
    "Urban and built-up",
    "Snow and ice",
    "Barren or sparse vegetation",
]

# 读取经纬度（此处主要是满足“读取 lat, lon”要求，
# 后续分析主要用 12 个土地利用变量）
lat = ds["lat"]
lon = ds["lon"]

# 将 12 个变量组合成一个 xarray.DataArray，维度为 (variable, y, x)
pft_da = ds[pft_vars].to_array("variable")  # shape: (12, y, x)

# 将 (y, x) 展平为 pixel 维度，得到 (pixel, variable)
# pixel = y * x
pft_flat = (
    pft_da
    .stack(pixel=("y", "x"))          # 维度 (variable, pixel)
    .transpose("pixel", "variable")   # 变为 (pixel, variable)
)

# 转为 numpy 数组，单位应为百分比 (0–100)
data = pft_flat.values  # shape: (n_pixel, 12)

# 同样展平 lat, lon，便于需要时使用（例如按像元对应空间位置）
lat_flat = lat.values.ravel()
lon_flat = lon.values.ravel()

# 处理缺失值：忽略含 NaN 的像元
mask_valid = ~np.isnan(data).any(axis=1)

# 处理全 0 行：这里选择“剔除全 0 像元”，因为它们不提供土地利用信息。
# 在论文或报告中可以写明这一点。
mask_nonzero = ~(np.all(data == 0, axis=1))

mask = mask_valid & mask_nonzero
data_valid = data[mask]
lat_valid = lat_flat[mask]
lon_valid = lon_flat[mask]

# 将数据放入 DataFrame，列名为变量名，便于相关性计算和后续处理
df = pd.DataFrame(data_valid, columns=pft_vars)

print(f"有效像元数量: {df.shape[0]}")
print(f"每个像元的变量数量: {df.shape[1]}")

# -----------------------------
# 2. 相关性分析与热力图
# -----------------------------
# 计算 Pearson 相关系数矩阵 (12x12)
corr_matrix = df.corr(method="pearson")

# 绘制相关系数热力图
fig, ax = plt.subplots(figsize=(10, 8))

# 使用对称色标 [-1, 1]，中性 0
im = ax.imshow(corr_matrix.values, vmin=-1, vmax=1, cmap="coolwarm")

# 轴刻度与标签
num_vars = len(pft_vars)
ax.set_xticks(np.arange(num_vars))
ax.set_yticks(np.arange(num_vars))
ax.set_xticklabels(pft_vars, rotation=45, ha="right", fontsize=8)
ax.set_yticklabels(pft_vars, fontsize=8)

# 在每个格子上标注相关系数值（保留两位小数）
for i in range(num_vars):
    for j in range(num_vars):
        value = corr_matrix.values[i, j]
        ax.text(
            j, i, f"{value:.2f}",
            ha="center", va="center", fontsize=6, color="black"
        )

# 添加颜色条
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Pearson correlation coefficient", fontsize=9)

ax.set_title("PFT Fractions Correlation Matrix (2000)", fontsize=12)
fig.tight_layout()
fig.savefig("PFT_corr_heatmap_2000.png", dpi=300)
plt.close(fig)

print("已保存相关系数热力图: PFT_corr_heatmap_2000.png")

# -----------------------------
# 3. 代表性变量对的散点图
# -----------------------------
# 计算上三角中的所有变量对及其相关系数
pairs = []
for i in range(num_vars):
    for j in range(i + 1, num_vars):
        r = corr_matrix.values[i, j]
        pairs.append(((i, j), r))

# 选取相关性最高的一对（正相关最大）
pairs_sorted_by_r = sorted(pairs, key=lambda x: x[1], reverse=True)
top_pos_pairs = [p for p in pairs_sorted_by_r if p[1] > 0]
best_pos_pair = top_pos_pairs[0] if top_pos_pairs else pairs_sorted_by_r[0]

# 选取相关性最低的一对（最负相关）
pairs_sorted_by_r_asc = sorted(pairs, key=lambda x: x[1])
top_neg_pairs = [p for p in pairs_sorted_by_r_asc if p[1] < 0]
best_neg_pair = top_neg_pairs[0] if top_neg_pairs else pairs_sorted_by_r_asc[0]

# 再选一对“接近 0 相关”的组合
pairs_sorted_by_abs = sorted(pairs, key=lambda x: abs(x[1]))
near_zero_pair = pairs_sorted_by_abs[0]

# 整理要绘制的三对变量
pairs_to_plot = [
    ("strong_positive", best_pos_pair),
    ("strong_negative", best_neg_pair),
    ("near_zero", near_zero_pair),
]

# 为了避免点太密，对像元随机下采样
np.random.seed(0)
sample_size = min(5000, len(df))
sample_idx = np.random.choice(len(df), size=sample_size, replace=False)
df_sample = df.iloc[sample_idx]

# 绘制三个散点图放在同一张图中
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for ax, (kind, ((i, j), r)) in zip(axes, pairs_to_plot):
    var_x = pft_vars[i]
    var_y = pft_vars[j]

    ax.scatter(
        df_sample[var_x],
        df_sample[var_y],
        s=5,
        alpha=0.3,
        edgecolors="none",
    )
    ax.set_xlabel(f"{var_x} (%)", fontsize=8)
    ax.set_ylabel(f"{var_y} (%)", fontsize=8)
    if kind == "strong_positive":
        subtitle = "强正相关示例"
    elif kind == "strong_negative":
        subtitle = "强负相关示例"
    else:
        subtitle = "接近零相关示例"
    ax.set_title(f"{subtitle}\n{var_x} vs {var_y} (r={r:.2f})", fontsize=9)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

fig.tight_layout()
fig.savefig("PFT_scatter_examples_2000.png", dpi=300)
plt.close(fig)

print("已保存散点图示例: PFT_scatter_examples_2000.png")

# -----------------------------
# 4. PCA 分析（基于标准化的 12 个变量）
# -----------------------------
# 标准化：去均值、除以标准差
scaler = StandardScaler()
data_scaled = scaler.fit_transform(df.values)

# 做 PCA，保留前 3 个主成分
pca = PCA(n_components=3)
scores = pca.fit_transform(data_scaled)  # shape: (n_samples, 3)

explained_ratio = pca.explained_variance_ratio_  # 长度为 3
components = pca.components_  # 形状 (3, 12)，每行对应一个主成分的载荷

# 绘制 PC1 vs PC2 的散点图
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(scores[:, 0], scores[:, 1], s=5, alpha=0.3, edgecolors="none")
ax.set_xlabel("PC1", fontsize=10)
ax.set_ylabel("PC2", fontsize=10)
ax.set_title(
    "PCA of PFT Fractions (PC1 vs PC2)\n"
    f"Explained var: PC1 {explained_ratio[0]*100:.1f}%, "
    f"PC2 {explained_ratio[1]*100:.1f}%",
    fontsize=10,
)
ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
fig.tight_layout()
fig.savefig("PFT_PCA_PC1_PC2_scatter_2000.png", dpi=300)
plt.close(fig)

print("已保存 PCA 散点图: PFT_PCA_PC1_PC2_scatter_2000.png")

# -----------------------------
# 5. 根据相关性 & PCA 自动生成中文总结
# -----------------------------
summary_lines = []

# 5.1 相关性方面的总结
# 取前三个较强正相关的变量对
top_pos_pairs = [p for p in pairs_sorted_by_r if p[1] > 0][:3]
if top_pos_pairs:
    for ((i, j), r) in top_pos_pairs:
        v1 = pft_vars[i]
        v2 = pft_vars[j]
        summary_lines.append(
            f"{v1} 与 {v2} 之间存在较强的正相关关系 (r≈{r:.2f})，说明这些土地利用类型往往在同一像元中同时出现。"
        )

# 取前三个最明显的负相关变量对
top_neg_pairs = [p for p in pairs_sorted_by_r_asc if p[1] < 0][:3]
if top_neg_pairs:
    for ((i, j), r) in top_neg_pairs:
        v1 = pft_vars[i]
        v2 = pft_vars[j]
        summary_lines.append(
            f"{v1} 与 {v2} 呈显著负相关 (r≈{r:.2f})，表明某一种类型的增加往往伴随另一种类型的减少。"
        )
else:
    summary_lines.append(
        "各土地利用类型之间几乎没有明显的强负相关关系，大多数变量对的相关性以正相关或弱相关为主。"
    )

# 接近零相关的示例
(i0, j0), r0 = near_zero_pair
summary_lines.append(
    f"{pft_vars[i0]} 与 {pft_vars[j0]} 的相关系数接近 0 (r≈{r0:.2f})，说明二者在空间上的变化相对独立。"
)

# 5.2 PCA 方面的总结
# 总方差解释率
summary_lines.append(
    f"前两个主成分 PC1 和 PC2 累计解释了约 {np.sum(explained_ratio[:2]) * 100:.1f}% 的总方差，能够概括大部分土地利用组成信息。"
)

# 分析前两个主成分的主要贡献变量
for pc_index in range(2):
    comp = components[pc_index]  # 该主成分的载荷
    # 按绝对值降序取前三个最重要的变量
    idx_sorted = np.argsort(np.abs(comp))[::-1][:3]
    var_contrib = []
    for idx in idx_sorted:
        sign = "正" if comp[idx] >= 0 else "负"
        var_contrib.append(f"{pft_vars[idx]}({sign}载荷)")
    contrib_str = "、".join(var_contrib)
    summary_lines.append(
        f"主成分 PC{pc_index + 1} 主要由 {contrib_str} 等类型共同贡献，代表了这些土地利用类型组合变化的主导梯度。"
    )

# 如果还没到 5 句，再加一两句总体性的描述（较为中性、不会违背数据）
if len(summary_lines) < 5:
    summary_lines.append(
        "从整体上看，森林类、草地/灌丛类与农田/裸地等类型在相关性和 PCA 空间中表现出不同的组合模式，反映了区域土地利用结构的差异。"
    )

# 控制在 5–10 句之间（如果过多就裁剪）
summary_lines = summary_lines[:10]

print("\n==== 相关性与 PCA 简要总结 ====\n")
for line in summary_lines:
    print("- " + line)
