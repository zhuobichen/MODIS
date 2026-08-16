# MODIS → MEGAN 植被参数处理

> 基于 MODIS 卫星数据，构建 MEGAN v3.2 生物源排放模型所需的植被功能类型（PFT）、叶面积指数（LAI）与排放因子（EFP）输入。

## 流程总览

MODIS MCD12Q1 土地覆盖 + MODIS LAI → 网格化 PFT 占比 → Growth Form 映射 → LAI 时序 → MEGAN EFP 排放因子（EF1..EF19）与光依赖因子（LDF3..LDF6）生成。

## 核心脚本

| 脚本 | 功能 |
|------|------|
| `Core_MODIS_IO.py` | MODIS 数据读写 |
| `Core_PFT.py` / `Core_GrowthForm.py` | 植被功能类型与生长型处理 |
| `Core_LAI.py` | 叶面积指数时序 |
| `Core_EFP.py` | 排放因子生成 |
| `Core_MosaicView.py` | 瓦片拼接与可视化 |

## 数据准备

本仓库**不包含**输入数据（卫星 HDF/NetCDF 体积大，未入库）。运行前需自行准备以下数据，目录结构见 `MODIS_Workflow_Detailed_Documentation.md` 第 2 节：

| 数据 | 来源 | 说明 |
|------|------|------|
| MODIS MCD12Q1 土地覆盖（LC_Type5） | NASA LP DAAC（`lpdaac.usgs.gov`） | 原始 HDF 瓦片 |
| MODIS LAI（8-day，MOD15A2H 等） | NASA LP DAAC | 已重采样为区域 NetCDF |
| GRIDCRO2D 网格文件 | WRF/CMAQ 预处理输出 | 提供目标域经纬度格点 |

> 依赖 `pyhdf` 读取 MODIS HDF，官方 wheel 仅覆盖 Windows / 旧版 Python；Linux/macOS 需用 conda 安装或源码编译（见 `requirements.txt` 注释）。

## 文档

- `MODIS_Workflow_Detailed_Documentation.md` — 全流程详细文档

## 技术栈

Python（numpy / pandas / pyhdf / netCDF4），输出多分辨率（3km / 9km / 27km）格点数据。
