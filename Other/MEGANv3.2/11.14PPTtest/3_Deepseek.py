import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

# 打开NetCDF文件
ds = xr.open_dataset('PFT_frac_2000_3km_square.nc')

# 定义12个土地利用类型变量名列表
pft_vars = [
    'Water',
    'Evergreen Needleleaf trees', 
    'Evergreen Broadleaf trees',
    'Deciduous Needleleaf trees',
    'Deciduous Broadleaf trees',
    'Shrub',
    'Grass',
    'Cereal crops',
    'Broad-leaf crops',
    'Urban and built-up',
    'Snow and ice',
    'Barren or sparse vegetation'
]

# 创建图形和网格布局 - 使用4x4网格实现斜向排版
fig = plt.figure(figsize=(16, 12))
gs = GridSpec(4, 4, figure=fig, wspace=0.3, hspace=0.3)  # 调整间距实现斜向分布效果

# 斜向排版的子图位置序列：从左上到右下对角线方向
diagonal_positions = [(0, 0), (0, 1), (0, 2),  # 第一行
                      (1, 0), (1, 1), (1, 2), (1, 3),  # 第二行  
                      (2, 1), (2, 2), (2, 3),  # 第三行
                      (3, 2), (3, 3)]  # 第四行

# 统一色标范围
vmin, vmax = 0, 100

# 循环绘制12个土地利用类型
axes = []
for i, (var_name, pos) in enumerate(zip(pft_vars, diagonal_positions)):
    # 在斜向位置创建子图
    ax = fig.add_subplot(gs[pos[0], pos[1]])
    axes.append(ax)
    
    # 获取数据并绘图
    data = ds[var_name]
    im = ax.pcolormesh(ds.lon, ds.lat, data, vmin=vmin, vmax=vmax, 
                       shading='auto', cmap='viridis')
    
    # 设置小标题（简化变量名显示）
    short_name = var_name.replace(' ', '\n')  # 换行显示长名称
    ax.set_title(short_name, fontsize=10, pad=5)
    
    # 隐藏坐标轴刻度
    ax.set_xticks([])
    ax.set_yticks([])

# 关闭剩余的空白子图（4x4网格中未使用的2个位置）
unused_positions = [(0, 3), (3, 0), (3, 1)]  # 未使用的位置
for pos in unused_positions:
    ax = fig.add_subplot(gs[pos[0], pos[1]])
    ax.axis('off')

# 添加公共颜色条
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # 右侧颜色条位置
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label('Percentage (%)', fontsize=12)

# 添加总标题
plt.suptitle('PFT fractions in 2000 (12 classes)', fontsize=16, y=0.95)

# 调整布局并保存
plt.tight_layout()
plt.subplots_adjust(top=0.92, right=0.9)  # 为总标题和颜色条留出空间
plt.savefig('PFT_12_classes_2000.png', dpi=300, bbox_inches='tight')
plt.show()

# 关闭数据集
ds.close()