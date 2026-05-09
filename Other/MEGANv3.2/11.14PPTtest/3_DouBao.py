import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

# 读取NetCDF文件
ds = xr.open_dataset("PFT_frac_2000_3km_square.nc")

# 提取所需变量（纬度、经度和草地百分比）
lat = ds["lat"]
lon = ds["lon"]
grass = ds["Grass"]  # 处理带空格的变量名时使用[]访问

# 创建图形和轴，使用PlateCarree投影
# 若未安装cartopy，请先运行：pip install cartopy
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": ccrs.PlateCarree()})

# 绘制草地百分比空间分布
mesh = ax.pcolormesh(lon, lat, grass, 
                     vmin=0, vmax=100,  # 颜色范围固定为0-100
                     cmap="YlGn",       # 绿色系 colormap 适合表示植被
                     shading="auto")

# 添加海岸线和国家边界（可选）
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=':')

# 添加颜色条，设置单位
cbar = plt.colorbar(mesh, ax=ax, orientation="vertical", pad=0.05)
cbar.set_label("Percentage (%)")

# 设置标题和坐标轴标签
ax.set_title("Grass fraction (%) in 2000", fontsize=14)
ax.set_xlabel("Longitude", fontsize=12)
ax.set_ylabel("Latitude", fontsize=12)

# 调整坐标轴刻度
ax.set_xticks(np.arange(int(lon.min()), int(lon.max()) + 1, 2))
ax.set_yticks(np.arange(int(lat.min()), int(lat.max()) + 1, 2))

# 设置合适的显示范围（根据数据自动调整）
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])

# 调整布局
plt.tight_layout()

# 保存图像
plt.savefig("Grass_fraction_2000.png", dpi=300, bbox_inches="tight")

# 显示图像（可选，如需在脚本中查看）
# plt.show()

# 关闭数据集
ds.close()