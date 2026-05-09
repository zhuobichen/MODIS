# MODIS → MEGAN v3.2 植被功能类型 (PFT) 与 LAI 处理及 EFP 排放因子生成全流程详细文档

> 文档生成时间：2025-11-20  
> 目录：`/DeepLearning/mnt/shixiansheng/MODIS`  
> 适用场景：以 MODIS MCD12Q1 土地覆盖 (LC_Type5) 与 MODIS LAI (8-day) 数据为源，构建多分辨率 (3km / 9km / 27km) 的 PFT 分布与 LAI 时序，转换为 MEGAN v3.2 预处理所需的 `grid_growth_form` / `LAI3_*` 等输入，再调用 MEGAN EFP (Emission Factor Processor) 生成格点排放因子 (EF1..EF19) 与光依赖因子 (LDF3..LDF6)。

---
## 1. 总体流程总览

```
原始 MODIS HDF (MCD12Q1 LC_Type5, LAI 8-day)             
          │
          ├─ 瓦片拼接与可视化 (hdfview_merged_PFT.py)
          │
          ├─ 最近邻或邻域统计 → 规则网格 PFT 百分比 (MODIS_PFT_FindNear_*.py / KDTree)
          │          │
          │          └─ 输出：PFT_frac_2000_{3|9|27}km_square(.csv/.nc)
          │
          ├─ 规则网格单分类最近邻提取 (PFT_cn27nc_Extract.py 等) → 每格点一个主类 (100/0)
          │
          ├─ PFT 分布 → Growth Form 映射 (v2.1Intov3.2_PFT_CsvTrans_Muti.py)
          │          └─ 输出：grid_growth_form_cn03/cn09/cn27.csv
          │
          ├─ LAI 8-day NetCDF 批量读取与拼接 (v2.1Intov3.2_LAI_CsvTrans_cn*.py)
          │          └─ 输出：LAI3_cn03 / LAI3_cn09 / LAI3_cn27.csv (46期 LAI01..LAI46)
          │
          ├─ 与生态类型 (grid_ecotype.*.csv) 合并 (外部生成)
          │
          ├─ MEGAN EFP 数据库构建 (MEGEFP32/MEGAN_EFP.py 调 run_M3EFP.py)
          │          ├─ SQLite 表：EF, SpeciationCrop/Tree/Shrub/Herb, GridEcotype, GridGrowthForm
          │          ├─ 中间表：IntermediateTreeEcoEF / Shrub / Herb / Crop
          │          └─ 最终表：OutputGridEF → OutputGridEF.<SCEN>.csv
          │
          └─ 下游：CMAQ / MEGAN 运行时使用 EF / LDF + LAI + GrowthForm + Ecotype
```

命名约定：
- `cn03` / `cn09` / `cn27`：对应区域中国 (GuangDong) 3km / 9km / 27km 分辨率网格。  
- `PFT_frac_2000_*km_square.csv`：每个格点内 12 个土地覆盖功能类型的百分比 (0–100)，来自正方形邻域统计。  
- `grid_growth_form_cn*.csv`：将 PFT 聚合为 MEGAN Growth Form 四类 (Tree / Crop / Shrub / Herb)。  
- `LAI3_cn*.csv`：2000 年 46 个 8-day LAI 切片 (LAI01..LAI46) + 空间坐标。  
- `OutputGridEF.<SCEN>.csv`：排放因子与光依赖因子汇总。

---
## 2. 数据源与坐标体系

### 2.1 MODIS 数据产品
| 产品 | 文件示例 | 关键变量 | 时间分辨率 | 投影 |
|------|----------|---------|-----------|------|
| MCD12Q1 v6 | `MCD12Q1.A2000001.h08v07.061.*.hdf` | `LC_Type5` (PFT/IGBP扩展分类) | 年度 (当年 DOY=001) | MODIS Sinusoidal |
| MODIS LAI (MOD15A2H等重采样结果，已转为区域 NetCDF) | `MODIS_LAI_2000001_GuangDong_cn27.nc` 等 | `LAI` | 8-day (46期) | 已对齐到目标规则网格 |

### 2.2 投影与网格
- MCD12Q1 使用 **正弦投影 (Sinusoidal)**：需要瓦片 h/v → 空间坐标转换：
  - 常数：`R=6371007.181`, `PIX_SIZE=463.3127165`, `NCOLS=NROWS=2400`, 全球左上角 (`X0=-20015109.354`, `Y0=10007554.677`).
  - 通过：`x_ul = X0 + h * TILE_SIZE`, `y_ul = Y0 - v * TILE_SIZE` 计算瓦片左上地理位置。
- CMAQ / 规则网格文件：`GRIDCRO2D_2000121_GuangDongD{1|2|3}`，用于获取目标域的经纬度格点集合 (LAT/LON)。
- 处理流程中有两种空间关系：
  1. **最近邻分类提取** (`PFT_cn*_Extract.py`)：每格点只有一个主类，赋值 100%。
  2. **正方形邻域比例统计** (`MODIS_PFT_FindNear_cn27.py`)：每格点一个 27km 正方形邻域内多类别比例 (0–100)。

---
## 3. 脚本功能详解

### 3.1 瓦片拼接与可视化：`hdfview_merged_PFT.py`
- 作用：批量读取指定年份 (DOY=001) 的所有可用 MCD12Q1 瓦片 HDF → 裁剪到设定 h/v 范围 → 按经纬度重排拼接为一个二维矩阵。
- 关键步骤：
  - 读取 `LC_Type5` 数据，将 `_FillValue` / >250 的值转为 NaN。
  - 从文件元属性或瓦片索引推导空间边界 (lat/lon)。
  - `pixel_per_degree = 2400 / 10` 确定拼接目标尺寸。
  - 独立绘制并叠加瓦片框线，标记是否在中国区域 (lon:73–135, lat:18–54)。
- 输出：`MCD12Q1_PFT_Tiles_A2000001.png`，用于人工核查数据覆盖情况与缺失瓦片。
- 使用场景：数据完整性检验、瓦片范围确认、后续 KDTree 最近邻处理前可视化核验。

### 3.2 最近邻 + 正方形邻域统计：`MODIS_PFT_FindNear_cn27.py`
- 目标：为目标域每个 CMAQ 格点生成 12 分类 (Water, Evergreen Needleleaf trees, ...) 百分比。
- 核心逻辑：
  1. 汇总瓦片 → 构建源点集合 (`src_pts`)：筛选有效分类 (0–11)。
  2. 空间索引：`cKDTree(src_pts)` → 为每格点中心用半径 = `HALF_M * sqrt(2)` 做圆形初筛。
  3. 方形过滤：对初筛结果再用 `|dx|≤HALF_M && |dy|≤HALF_M` 剔除圆中但不在正方形邻域的点。
  4. 分类计数：`np.bincount(cls_sel, minlength=12)` 得到 12 类计数；再除以总数转百分比 (0–100)。
  5. 并行：可选 `PARALLEL=True`，多进程分批 (BATCH=6000)。
- 输出：
  - NetCDF：`PFT_frac_2000_27km_square.nc`（变量名即分类全称 + `lat`/`lon`）
  - CSV：`PFT_frac_2000_27km_square.csv`（平坦表：CELLID, ICELL, JCELL, LAT, LON, 每类百分比）
- 优势：相比单类最近邻，更能反映格点周边土地覆盖混合程度；用于后续 GrowthForm 总和。
- 注意：空邻域 → 百分比全 0。需检查是否由于跳过瓦片列表 `skip_tiles` 过多造成。

### 3.3 单分类最近邻提取：`PFT_cn27nc_Extract.py`/同类 `pft_cn*_GuangDong_2000.csv`
- 功能：对每个目标格点仅保留一个主类 (LC_Type5 最近邻)，其对应列置 100，其余分类列置 0。
- 输出结构：`CELLID, ICELL, JCELL, Evergreen_Needleleaf_Trees, ..., Water`
- 使用场景：需要单一主土地类型的简化模型、或用于某些只接受单类输入的旧版预处理链。
- 与邻域比例文件相比：信息更简化、不含混合分布；可能丢失细分植被结构。

### 3.4 Growth Form 转换：`v2.1Intov3.2_PFT_CsvTrans_Muti.py`
- 输入：`PFT_frac_*km_square.csv`（27/9/3 km 三种分辨率）。
- 映射规则：
  - TreeFrac = Sum(Evergreen Needleleaf + Evergreen Broadleaf + Deciduous Needleleaf + Deciduous Broadleaf)
  - CropFrac = Sum(Cereal Crops + Broad-leaf Crops)
  - ShrubFrac = Shrub
  - HerbFrac = Grass
- 输出：`grid_growth_form_cn27.csv` / `grid_growth_form_cn09.csv` / `grid_growth_form_cn03.csv`：列 `gridID, TreeFrac, CropFrac, ShrubFrac, HerbFrac`。
- 注意：不做归一化重算；若源百分比已为 0–100 总和约 100，将保持一致。建议后续验证：`Tree+Crop+Shrub+Herb` ≤ 100，否则需检查分类覆盖是否缺失水体/城市等。

### 3.5 LAI 时序拼接：`v2.1Intov3.2_LAI_CsvTrans_cn03/09/27.py`
- 输入：按 8-day 间隔切分的 NetCDF：`MODIS_LAI_<year><DOY>_GuangDong_<grid>.nc` 共 46 期 (001 → 361)。
- 检查：
  - 期序列：`expected_days = [001, 009, 017, ..., 361]` 间隔固定 8 天。
  - 缺文件即报错退出；保证时间完整性。
- 处理：
  - 读取每期二维 LAI → 展平 → 缺失值 `_FillValue` 替换为 0（注意：0 表示无叶面积或缺失需下游区分）。
  - 构造 `CELL_ID (1-based), X, Y, LAT, LONG, LAI01..LAI46`。
  - 不同脚本对 X/Y 编号有细微差异：cn03/cn09 使用 1-based 行列；cn27 脚本中仍保留 0-based → 建议统一。
- 输出：`LAI3_cn03.csv`, `LAI3_cn09.csv`, `LAI3_cn27.csv`。
- 注意：如果后续 MEGAN 预处理期望 `ROW/COL` 或 `X/Y` 为 1-based，需要对 cn27 输出做兼容转换。

### 3.6 MEGAN EFP 排放因子处理：`MEGEFP32/MEGAN_EFP.py`, `src/run_M3EFP.py`, `src/M3GEFP.py`

#### 3.6.1 入口脚本：`MEGAN_EFP.py`
- 参数：
  - `scen_name = "GD_cn27"` → 决定输出文件与数据库命名。
  - 数据库路径：`./database/M3GEFP_database.<scen_name>.db`
  - 输入 CSV 目录：`./inputs/EFP/` 包含：
    - 植物物种排放因子表：`EFv210806.csv`
    - 物种组成 (各生长形态生态类型物种构成)：`SpeciationCrop/Herb/Tree/Shrub*.csv`
    - 空间生态类型：`grid_ecotype.<scen>.csv`
    - 生长形态分布：`grid_growth_form.<scen>.csv`（由前述 Growth Form 处理生成）。
  - 类别数：`EFClasses=19`, `LDFClasses0=3`, `LDFClasses1=6`（若与默认不同需同步调整中间查询 build 逻辑）。

#### 3.6.2 驱动：`run_M3EFP.py`
- 功能：创建 SQLite 数据库 + 加载所有输入表。
- 步骤：
  1. 校验物种 EF 完整性：若某个 `VegID` 在物种组成表出现但不在 EF 表，则提示是否补零。
  2. 将各 CSV 载入 SQLite：表名固定 (EF, SpeciationCrop, SpeciationShrub, SpeciationHerb, SpeciationTree, GridEcotype, GridGrowthForm)。
  3. 调用 `M3GEFP.run_M3GEFP_DB()` 执行中间表与最终表生成。

#### 3.6.3 数据库模块：`M3GEFP.py`
- 中间表：`Intermediate<Tree|Shrub|Herb|Crop>EcoEF`：对每生态类型聚合加权的物种 EF 与 LDF：
  ```sql
  Sum([EF{i}] * [TreeSpecfrac]) AS TreeEcoEF{i}
  Sum([LDF{j}] * [TreeSpecfrac]) AS TreeEcoLDF{j}
  ```
- 最终表构造：`OutputGridEF`：按格点循环聚合：
  ```
  EF{i} = Sum( EcotypeFrac * (
        CropFrac*CropEcoEF{i} + TreeFrac*TreeEcoEF{i} + HerbFrac*HerbEcoEF{i} + ShrubFrac*ShrubEcoEF{i} ) )
  LDF{j} = 同理聚合
  ```
- 输出：`OutputGridEF.<SCEN>.csv` 包含：`gridID, EF1..EF19, LDF3..LDF6`。

### 3.7 其他脚本与目录
| 文件/目录 | 说明 | 备注 |
|-----------|------|------|
| `3_ChatGPT.py`, `3_Deepseek_QingXie.py` | 说明/测试版本脚本 | 多为探索性生成，不在主链路 |
| `MOIDS_PFT.py` | 可能为旧版 PFT 处理入口 | 未阅读，建议标注 deprecated |
| `Nearest/` *MODIS_PFT_Nearest_cn*.py* | 可能是简化最近邻实现 | 与 KDTree + 正方形逻辑类似，验证性能差异 |
| `MEGEFP32/outputs/` | 多域/分辨率 EF 输出 | 如 `OutputGridEF.GD_cn09.csv` |
| `MEGEFP32/jupyter/MEGEFP32.ipynb` | 交互式测试与验证 | 可转为自动化单元测试 |
| `LAI/` | 原始 8-day LAI NetCDF 文件目录 | 输入由其它脚本生成或下载 |
| `Input/` | 原始 MODIS HDF 瓦片目录 | MCD12Q1 文件放置路径 |

---
## 4. 输入与输出字段详述

### 4.1 PFT 百分比文件 `PFT_frac_2000_27km_square.csv`
| 字段 | 含义 | 范围 |
|------|------|------|
| CELLID | 全域唯一格点 ID (1-based) | 1..N |
| JCELL / ICELL | 行 / 列号 (1-based) | |
| LAT / LON | 格点中心经纬度 | 实数 |
| 12 分类列 | 各类在正方形邻域内的百分比 | 0–100 (总和≈100 或 ≤100 若有缺失) |

### 4.2 Growth Form 文件 `grid_growth_form_cn27.csv`
| 字段 | 含义 |
|------|------|
| gridID | 继承自 PFT_frac 的 CELLID 或外部映射 | 唯一格点 |
| TreeFrac | 树木功能类型合计百分比 | 0–100 |
| CropFrac | 农作物功能类型合计百分比 | 0–100 |
| ShrubFrac | 灌木百分比 | 0–100 |
| HerbFrac | 草本百分比 | 0–100 |

### 4.3 LAI 文件 `LAI3_cn27.csv`
| 字段 | 含义 |
|------|------|
| CELL_ID | 格点 ID (1-based) |
| X / Y | 列 / 行索引 (需确认是否统一 1-based) |
| LAT / LONG | 格点经纬度 |
| LAI01..LAI46 | 每 8 天 LAI 值，缺失填 0 | 建议后续 0 区分 (原始缺失 vs 冬季叶落)

### 4.4 EFP 输出 `OutputGridEF.GD_cn27.csv`
| 字段 | 含义 |
|------|------|
| gridID | 空间格点 ID，对齐 `grid_growth_form` 与 `grid_ecotype` |
| EF1..EF19 | 各排放物或功能类型的整合排放因子 | 数值（单位依数据定义） |
| LDF3..LDF6 | 光依赖因子 (某些挥发性有机物相关) | 数值 |

---
## 5. 运行推荐顺序

1. 数据准备：放置所有 MCD12Q1 HDF 文件到 `Input/`，准备 LAI NetCDF 切片到 `LAI/`。
2. 瓦片完整性与覆盖检查：运行 `hdfview_merged_PFT.py` → 生成拼接图。
3. 生成 PFT 百分比：执行 `MODIS_PFT_FindNear_cn27.py`（按需要复制并改分辨率参数与输出文件名）。
4. 可选：生成最近邻主类提取：`PFT_cn27nc_Extract.py`（若模型需单类）。
5. 构建 Growth Form：运行 `v2.1Intov3.2_PFT_CsvTrans_Muti.py` → 获取 `grid_growth_form_cn*.csv`。
6. 拼接 LAI：依分辨率运行对应 `v2.1Intov3.2_LAI_CsvTrans_cn*.py` → 得到 `LAI3_cn*.csv`。
7. 准备生态类型 `grid_ecotype.<SCEN>.csv`（外部生成，确保与格点 ID 对齐）。
8. 运行排放因子：`python MEGEFP32/MEGAN_EFP.py` → 生成 `OutputGridEF.<SCEN>.csv`。
9. 校验：检查 EF 列是否存在极端 0 或异常高值；验证 GrowthForm 总和与 PFT 分布一致性；LAI 缺失比例。
10. 下游：MEGAN/CMAQ 输入文件构建（与 `run_prepmegan4cmaq.csh` 主脚本联动）。

---
## 6. 质量控制与常见问题

| 问题 | 可能原因 | 排查与解决 |
|------|----------|-----------|
| 某些格点 PFT 全 0 | 邻域外无有效像元 / 跳过瓦片过多 | 检查 `skip_tiles` 列表与域边界重叠；放宽筛选半径或检查瓦片缺失 |
| Growth Form 总和 > 100 | 输入 PFT 百分比含重复或未剔除水体/城市 | 在转换前先过滤非植被类别或归一化 `植被类总和` |
| LAI 全 0 或大量 0 | 填充值被替换但原始文件缺失 / 空间错配 | 校验 `_FillValue`，确认 NetCDF 维度与脚本 `row/col` 参数一致 |
| EF 缺少物种 | EF 表没有对应 `VegID` | 运行时脚本会提示并可补 0；建议补全真实因子避免后续低估 |
| EF/LDF 全部为 0 | 输入 `grid_ecotype` 或 `grid_growth_form` 未对齐 | 检查 `gridID` 一致性与行数；打印中间表行计数 |
| 运行极慢 (PFT 统计) | KDTree 查询批次过大 / 并行关闭 | 调整 `BATCH`，启用 `PARALLEL=True` 并提高 `N_PROCS` |
| 27km/9km/3km 文件混淆 | 命名不统一 | 严格使用模板：`PFT_frac_2000_<res>km_square.csv` + `grid_growth_form_cn<res>.csv` |
| LAI 不与 EF 格点对齐 | 经纬度栅格源不同 | 对比 `LAT/LON` 差异；必要时重采样或投影统一 |

---
## 7. 改进与重构建议

| 方向 | 建议 | 价值 |
|------|------|------|
| 投影与瓦片工具抽象 | 抽出 `modis_geo.py` 封装坐标、瓦片 bbox、KDTree 构建 | 降低脚本重复度，统一常数 |
| 分辨率参数化 | 通过 argparse 传递 row/col/grid/resolution 而不是复制三个 LAI 脚本 | 易扩展至其它区域年份 |
| Growth Form 正常化 | 增加总和校验与 `植被总和` 与 `非植被` 拆分 | 防止输入异常传播到 EF 计算 |
| LAI 缺失标记 | 使用 `NaN` 而非 0，后续再在 MEGAN 阶段插值或填补 | 避免误将缺失视为叶面积为 0 |
| 中间表输出缓存 | 保存 `Intermediate*EcoEF` 表到 CSV 可复用 | 便于调试 EF 与物种组合权重 |
| 单元测试 | 针对：EF 完整性 / GrowthForm 计算 / 邻域统计边界条件 | 降低维护风险 |
| 统一 ID 与坐标 | 强制使用 `ROW/COL` + `CELLID`，并记录投影坐标 X/Y | 减少空间错位风险 |
| 并行优化 | 使用 `joblib.Parallel` 或 `vectorized filtering` 替换多进程开销 | 提升大域处理性能 |

---
## 8. 依赖与环境

| 组件 | 用途 |
|------|------|
| `pyhdf` | 读取 MODIS HDF (LC_Type5) |
| `xarray` / `netCDF4` | NetCDF 读写与多维数据处理 |
| `numpy` / `pandas` | 数值与表处理 |
| `scipy.spatial.cKDTree` | 快速邻域候选检索 |
| `sqlite3` | MEGAN EFP 数据库构造 |
| `matplotlib` | 拼接图与可视化验证 |

额外建议：在 `MEGEFP32/requirements.txt` 明确版本范围，避免 Python 3.6/3.11 在类型与行为上的差异影响（例如 pandas append 已弃用）。

---
## 9. 关键脚本间依赖关系图

```
[PFT_frac_2000_*km_square.csv] ← MODIS_PFT_FindNear_cn*.py
        │
        ├──► grid_growth_form_cn*.csv (v2.1Intov3.2_PFT_CsvTrans_Muti.py)
        │
        ├──► (可选) pft_cn*_GuangDong_2000.csv (PFT_cn*_Extract.py)
        │
[LAI3_cn*.csv] ← v2.1Intov3.2_LAI_CsvTrans_cn*.py
[grid_ecotype.<SCEN>.csv] ← 外部生态类型预处理
        │
        ▼
MEGAN_EFP.py → run_M3EFP.py → M3GEFP.py → OutputGridEF.<SCEN>.csv
        │
        ▼
下游 MEGAN/CMAQ 模型运行 (排放计算) 结合：EF/LDF + LAI + GrowthForm + Ecotype
```

---
## 10. 快速核验脚本（建议新增）
可新增一个 `verify_modis_pipeline.py`：
- 检查：
  - PFT 百分比是否在 [0,100] 且总和 <= 100 + 容差。
  - Growth Form 四类与源 PFT 总和关系。
  - LAI 列数是否 46，缺失比例统计。
  - EF 输出是否每格点均有非零值（或可解释的零）。
- 统计摘要写出 `MODIS_Pipeline_QC_Report.txt`。

---
## 11. 常用命令示例

```bash
# 生成 27km PFT 百分比
python MODIS_PFT_FindNear_cn27.py

# 转换为 Growth Form（多分辨率）
python v2.1Intov3.2_PFT_CsvTrans_Muti.py

# 拼接 LAI → LAI3_cn27.csv
python v2.1Intov3.2_LAI_CsvTrans_cn27.py

# 运行 MEGAN EFP 场景 GD_cn27
python MEGEFP32/MEGAN_EFP.py
```

---
## 12. 结语
本目录完成了从遥感原始分类与 LAI 观测到 MEGAN v3.2 所需空间植被结构与排放因子输入的全链路。其关键价值在于：
1. **空间一致性保障**：通过 KDTree 最近邻与正方形邻域统计，映射分类到多尺度 CMAQ 网格。  
2. **植被功能抽象**：PFT → Growth Form → Ecotype 与物种构成的层次化整合。  
3. **排放参数生成自动化**：基于 SQLite 中间表构建与批量聚合计算 EF/LDF。  
4. **可扩展性**：分辨率参数化后可推广至不同区域与年份，增强区域排放研究的统一性。  

为未来维护与性能提升，建议优先实施：坐标与分辨率参数统一、异常值自动质控、以及公共函数模块化抽离。

---
*如需按单脚本继续逐行解释或补充 QC 工具，请告知。*
