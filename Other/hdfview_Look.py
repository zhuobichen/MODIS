from pyhdf.SD import SD, SDC

# 文件路径
file_path = 'MCD12Q1.A2000001.h08v07.061.2022147203319.hdf'

# 打开 HDF 文件
hdf = SD(file_path, SDC.READ)

# 查看所有数据集
print("📂 文件中的数据集：")
datasets = hdf.datasets()
for name, info in datasets.items():
    print(f"{name} - 维度: {info[0]}, 类型: {info[3]}")

# 查看全局属性
print("\n🧾 全局属性：")
for attr in hdf.attributes():
    print(f"{attr}: {hdf.attributes()[attr]}")
