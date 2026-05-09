# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # 只为激活 3D
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap, ScalarMappable

# 1. 打开 NetCDF 文件
ds = xr.open_dataset("PFT_frac_2000_3km_square.nc")

# 2. 定义要“叠放”的土地利用类型（按需修改顺序和内容）
#    示例：水体、常绿针叶林、Grass 三层
selected_vars = [
    "Water",
    "Evergreen Needleleaf trees",
    "Grass",
]
n_layers = len(selected_vars)

# 3. 取经纬度和数据
lon = ds["lon"].values   # shape (y, x)
lat = ds["lat"].values   # shape (y, x)

# 保证是 numpy 数组
X = lon
Y = lat

# 统一色标范围
vmin, vmax = 0, 100
norm = Normalize(vmin=vmin, vmax=vmax)
cmap = get_cmap("viridis")

# 4. 创建 3D 图像
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")

# 每一层在 z 方向的高度偏移
z_step = 10.0  # 层间距，可改大/改小
z_offsets = np.arange(n_layers) * z_step

# 5. 循环画每一层的“悬浮切片”
for k, (var_name, z0) in enumerate(zip(selected_vars, z_offsets)):
    data = ds[var_name].values  # shape (y, x)

    # 颜色映射
    colors = cmap(norm(data))

    # 当前层的 z 面：常数平面
    Z = np.full_like(X, z0, dtype=float)

    # 用 plot_surface，把二维场贴在 z=z0 的平面上
    surf = ax.plot_surface(
        X, Y, Z,
        rstride=1, cstride=1,
        facecolors=colors,
        linewidth=0,
        antialiased=False,
        shade=False,
        alpha=0.95
    )

    # 在每一层上方写上变量名
    ax.text(
        np.nanmean(X), np.nanmax(Y) + 0.5, z0 + 0.5,
        var_name,
        fontsize=8, ha="center", va="bottom"
    )

# 6. 设置视角和坐标轴
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_zlabel("Layer index")

# 调整视角（可以多试几个角度）
ax.view_init(elev=25, azim=-60)

# Z 轴范围稍微放宽一点
ax.set_zlim(-2, z_offsets[-1] + 5)

# 7. 加一个公共颜色条
mappable = ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array([])
cbar = fig.colorbar(mappable, ax=ax, fraction=0.03, pad=0.08)
cbar.set_label("Percentage (%)")

plt.title("Stacked PFT fractions (example layers)", fontsize=12)
plt.tight_layout()
plt.savefig("PFT_stacked_layers_3D.png", dpi=300, bbox_inches="tight")
plt.show()

ds.close()
