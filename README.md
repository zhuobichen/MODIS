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

## 文档

- `MODIS_Workflow_Detailed_Documentation.md` — 全流程详细文档

## 技术栈

Python（numpy / pandas / pyhdf / netCDF4），输出多分辨率（3km / 9km / 27km）格点数据。
