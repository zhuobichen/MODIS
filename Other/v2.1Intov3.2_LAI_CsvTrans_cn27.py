import os
import glob
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# 参数设置
row = 175
col = 124
year = "2000"
grid = "cn27"
project = "GuangDong"
# 预期的日序（从001开始，每8天一次，共46个）
expected_days = [f"{1 + 8*i:03d}" for i in range(46)]
ndays = len(expected_days)

# 输入输出路径（修改输出文件名格式）
indir = f"./LAI/{grid}_{project}_{year}/"
outfile = f"./LAI3_{grid}.csv"  # 新文件名格式：LAI3_cn27.csv

# 预先检查日期序列完整性
print("检查日期序列完整性...")
if expected_days[-1] != "361":
    raise ValueError("日期序列不完整，最后一个日期应为361")
for i in range(1, ndays):
    prev = int(expected_days[i-1])
    curr = int(expected_days[i])
    if curr - prev != 8:
        raise ValueError(f"日期间隔错误: {expected_days[i-1]}到{expected_days[i]}不是8天")
print("日期序列检查通过")

# 检查输入文件是否存在
input_files = []
for day in expected_days:
    pattern = f"{indir}MODIS_LAI_{year}{day}_{project}_{grid}.nc"
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"缺少文件: {pattern}")
    input_files.append(files[0])

# 读取经纬度信息（二维网格，从第一个文件获取）
print("读取经纬度信息...")
first_file = Dataset(input_files[0], "r")
lat = first_file.variables["lat"][:]  # 二维数组 (row, col)
lon = first_file.variables["lon"][:]  # 二维数组 (row, col)
first_file.close()

# 确保经纬度维度正确（二维）
if lat.shape != (row, col) or lon.shape != (row, col):
    raise ValueError(f"经纬度维度不匹配: 预期({row},{col})，实际lat{lat.shape}, lon{lon.shape}")

# 初始化数据数组（缺失值初始化为0）
print("开始处理LAI数据...")
total_pixels = row * col
lai_data = np.full((total_pixels, ndays), 0.0)  # 用0作为填充值

# 读取每个日期的LAI数据
for d, (day, file_path) in enumerate(zip(expected_days, input_files)):
    print(f"处理第{d+1}/{ndays}个文件: {os.path.basename(file_path)}")
    with Dataset(file_path, "r") as nc:
        lai = nc.variables["LAI"][:]  # 假设LAI变量为二维 (row, col)
        
        # 检查数据维度
        if lai.shape != (row, col):
            raise ValueError(f"LAI维度不匹配: 预期({row},{col})，实际{lai.shape}")
        
        # 展平为1维并存储（按行优先）
        lai_flat = lai.flatten(order='C')
        # 替换填充值为0
        fill_value = getattr(nc.variables["LAI"], "_FillValue", -999.0)
        lai_flat = np.where(lai_flat == fill_value, 0.0, lai_flat)  # 缺失值用0填充
        lai_data[:, d] = lai_flat

# 生成CELL_ID、X、Y坐标
cell_ids = np.arange(1, total_pixels + 1)  # CELL_ID从1开始
y_coords, x_coords = np.mgrid[0:row, 0:col]  # y是行索引，x是列索引
x_flat = x_coords.flatten(order='C')
y_flat = y_coords.flatten(order='C')

# 展平经纬度
lat_flat = lat.flatten(order='C')
lon_flat = lon.flatten(order='C')

# 创建DataFrame
print("生成输出文件...")
df = pd.DataFrame({
    "CELL_ID": cell_ids,
    "X": x_flat,
    "Y": y_flat,
    "LAT": lat_flat,
    "LONG": lon_flat
})

# 添加LAI列（LAI01到LAI46）
for i in range(ndays):
    df[f"LAI{i+1:02d}"] = lai_data[:, i]

# 后处理：检查缺失值（统计0的数量）
missing_count = df.iloc[:, 5:].eq(0.0).sum().sum()
if missing_count > 0:
    print(f"警告: 数据中存在{missing_count}个缺失值（已用0标记）")

# 保存为CSV
df.to_csv(outfile, index=False)
print(f"处理完成，输出文件: {outfile}")