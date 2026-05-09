# MODIS 项目 Core+Run 分层架构规范

> 目标：将 MODIS → MEGAN 遥感处理全链路重构为 Core+Run 分层架构  
> 版本: 2026-05 (v1 - 初始设计)  
> 参考基准: `/DeepLearning/mnt/shixiansheng/DataFusion_China_CleanAir/docs/DataFusion项目架构规范_Core+Run分层设计.md`  
> 基准项目: `/DeepLearning/mnt/shixiansheng/MODIS`

---

## 一、项目概述

MODIS 项目完成从原始 MODIS 卫星遥感数据（MCD12Q1 土地覆盖分类 + MODIS LAI 叶面积指数）到 MEGAN v3.2 排放模型所需输入的全链路处理：

```
原始 MODIS HDF (MCD12Q1 LC_Type5 + LAI 8-day)
    │
    ├─ 瓦片拼接与可视化 → Picture/MosaicView/
    │
    ├─ PFT 邻域统计 → Data/PFT/
    │       └─ Growth Form 转换 → Data/GrowthForm/
    │
    ├─ LAI 时序拼接 → Data/LAI/
    │
    └─ MEGAN EFP 排放因子计算 → Data/EFP/Output/
            └─ 下游: CMAQ / MEGAN 排放模型
```

---

## 二、顶层目录结构

```
MODIS/
├── Core_*.py              # 核心模块（纯逻辑，无硬编码路径，无 __main__）
├── Run_*.py               # 入口脚本（薄层，仅配置 + 调用 Core）
├── Data/                   # 所有输入/输出数据
│   ├── Input/              #   原始输入数据
│   │   ├── HDF/            #     MODIS HDF 文件 (MCD12Q1, MCD15A2H)
│   │   ├── LAI/            #     LAI NetCDF 切片 (按分辨率/年份)
│   │   └── EFP/            #     EFP 输入 CSV (EF表、物种组成、生态类型)
│   ├── PFT/                #   PFT 处理输出
│   ├── GrowthForm/         #   Growth Form 输出
│   ├── LAI/                #   LAI CSV 输出
│   └── EFP/                #   EFP 输出
│       ├── Database/       #     SQLite 数据库
│       └── Output/         #     OutputGridEF CSV
├── Picture/                # 所有可视化输出
│   ├── MosaicView/         #   瓦片拼接图
│   └── QC/                 #   质控报告图表
├── megan_efp/              # MEGAN EFP 算法子模块（独立可复用）
│   ├── src/                #   核心算法 (run_M3EFP.py, M3GEFP.py)
│   ├── jupyter/            #   交互式测试
│   └── requirements.txt
├── modis_geo_utils.py      # MODIS 地理工具（正弦投影、瓦片坐标、KDTree）
├── docs/                   # 流程文档
├── Other/                  # 废弃旧代码归档
├── .gitignore
└── MODIS_Workflow_Detailed_Documentation.md  # 全流程详解文档
```

---

## 三、Core+Run 分层架构

### 3.1 核心原则

```
┌──────────────────────────────────────────┐
│  Run_*.py                                │
│  - 硬编码路径配置                          │
│  - 硬编码分辨率/年份/瓦片范围              │
│  - 仅调用 Core_ 模块的函数                 │
│  - 不超过 50 行有效逻辑                    │
├──────────────────────────────────────────┤
│  Core_*.py                               │
│  - 纯函数，所有路径/参数通过函数签名传入     │
│  - 可被多个 Run_* 或 其他 Core_* 复用       │
│  - 包含上层管道函数（接受 years/resolution） │
│  - 无 if __name__ == "__main__"            │
├──────────────────────────────────────────┤
│  megan_efp/ , modis_geo_utils.py         │
│  - 底层工具/算法库，无项目特定路径          │
│  - 通过 import 或 sys.path 加载            │
└──────────────────────────────────────────┘
```

### 3.2 Run 入口脚本模板

```python
#!/usr/bin/env python3
"""
Run_XXX.py - 功能描述
======================
输入: Data/Input/XXX
输出: Data/XXX 或 Picture/XXX

使用方式:
    python Run_XXX.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Core_XXX import run_pipeline

# ====== 配置 ======
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(BASE_DIR, "Data", "Input", "HDF")
OUTPUT_DIR = os.path.join(BASE_DIR, "Data", "PFT")
GRID_NC = os.path.join(BASE_DIR, "GRIDCRO2D_2000121_GuangDongD1")
YEAR = 2000
RESOLUTION = "27km"
# ... 其他硬编码配置

# ====== 执行 ======
if __name__ == "__main__":
    run_pipeline(
        year=YEAR,
        resolution=RESOLUTION,
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        grid_nc=GRID_NC,
    )
```

### 3.3 Core 模块模板

```python
"""
Core_XXX.py - 模块描述
=======================
提供函数:
  - func_a(): ...
  - run_pipeline(): 上层管道函数
"""

def func_a(data: np.ndarray, resolution: str) -> np.ndarray:
    """纯逻辑函数。"""
    ...

def run_pipeline(
    year: int,
    resolution: str,
    input_dir: str,
    output_dir: str,
    grid_nc: str,
    ...
) -> dict:
    """管道函数：串联多个 func_* 完成完整流程。"""
    ...
    return summary
```

---

## 四、模块清单与职责

### 4.1 Core 模块

| Core 模块 | 职责 | 抽离自 |
|-----------|------|--------|
| `Core_MODIS_IO.py` | HDF 文件读取（LC_Type5 / LAI / QC）、MODIS 正弦投影坐标计算（经纬度与米坐标互转）、瓦片 bbox 与目标域重叠判断 | `hdfview_merged_PFT.py`, `MODIS_PFT_FindNear_cn27.py`, `hdfview_data_merged.py` |
| `Core_PFT.py` | KDTree 构建与查询、正方形邻域 PFT 百分比统计、最近邻单类提取、百分比数组转 CSV/NetCDF 输出 | `MODIS_PFT_FindNear_cn*.py`, `MOIDS_PFT.py`, `Nearest/MODIS_PFT_Nearest_cn*.py` |
| `Core_GrowthForm.py` | PFT 12 类百分比 → Growth Form 4 类（Tree/Crop/Shrub/Herb）聚合转换 | `v2.1Intov3.2_PFT_CsvTrans_Muti.py` |
| `Core_LAI.py` | LAI NetCDF 时序批量读取、日期序列校验、展平拼接为 CSV | `v2.1Intov3.2_LAI_CsvTrans_cn*.py` |
| `Core_MosaicView.py` | MODIS 瓦片拼接（PFT / LAI）、拼接图绘制与瓦片边界标注 | `hdfview_merged_PFT.py`, `hdfview_data_merged.py` |
| `Core_EFP.py` | MEGAN EFP 入口封装：构造场景配置 → 调用 `megan_efp/` 算法库生成排放因子 | `MEGEFP32/MEGAN_EFP.py` |
| `Core_Verification.py` | 全链路 QC 校验：PFT/GrowthForm 一致性、LAI 缺失统计、EF 极值检测、零值网格解释 | `verify_modis_pipeline.py` |

### 4.2 Run 入口一览

| Run 脚本 | 功能 | 调用 Core | 输入 | 输出 |
|----------|------|-----------|------|------|
| `Run_MosaicView.py` | MODIS 瓦片拼接可视化 | `Core_MosaicView` | `Data/Input/HDF/` | `Picture/MosaicView/` |
| `Run_PFT_Percentage.py` | PFT 百分比提取（正方形邻域） | `Core_MODIS_IO`, `Core_PFT` | `Data/Input/HDF/`, 网格 NC | `Data/PFT/` |
| `Run_PFT_Nearest.py` | PFT 最近邻单类提取 | `Core_MODIS_IO`, `Core_PFT` | `Data/Input/HDF/`, 网格 NC | `Data/PFT/` |
| `Run_GrowthForm.py` | PFT → Growth Form 转换 | `Core_GrowthForm` | `Data/PFT/` | `Data/GrowthForm/` |
| `Run_LAI.py` | LAI 时序拼接 | `Core_LAI` | `Data/Input/LAI/` | `Data/LAI/` |
| `Run_EFP.py` | MEGAN 排放因子计算 | `Core_EFP` | `Data/Input/EFP/`, `Data/GrowthForm/` | `Data/EFP/Output/` |
| `Run_Verify.py` | 全链路质控校验 | `Core_Verification` | `Data/PFT/`, `Data/GrowthForm/`, `Data/LAI/`, `Data/EFP/Output/` | `Picture/QC/` + TXT 报告 |

---

## 五、Data 目录规范

```
Data/
├── Input/                    # 原始输入数据（只读，不写入）
│   ├── HDF/                  #   MODIS HDF 瓦片
│   │   ├── MCD12Q1.A2000001.h23v03.061.*.hdf
│   │   └── MCD15A2H.A2000361.h23v03.*.hdf
│   ├── LAI/                  #   LAI NetCDF 切片（按分辨率/项目/年份）
│   │   └── {grid}_{project}_{year}/
│   │       └── MODIS_LAI_{year}{doy}_{project}_{grid}.nc
│   └── EFP/                  #   EFP 输入表
│       ├── EFv210806.csv
│       ├── SpeciationCrop210806.csv
│       ├── SpeciationHerb210806.csv
│       ├── SpeciationShrub210806.csv
│       ├── SpeciationTree210725.csv
│       └── grid_ecotype.{scen}.csv
├── PFT/                      # PFT 处理输出
│   └── {year}/
│       ├── PFT_frac_{year}_{resolution}_square.nc
│       └── PFT_frac_{year}_{resolution}_square.csv
├── GrowthForm/               # Growth Form 输出
│   └── {year}/
│       └── grid_growth_form_{resolution}.csv
├── LAI/                      # LAI CSV 输出
│   └── {year}/
│       └── LAI3_{resolution}.csv
└── EFP/                      # EFP 输出
    ├── Database/             #   SQLite 数据库
    │   └── M3GEFP_database.{scen}.db
    └── Output/               #   排放因子 CSV
        └── OutputGridEF.{scen}.csv
```

**规则**:
- 子目录按 **数据来源/功能** 命名
- 年份用子目录（如 `PFT/{year}/`）
- `Data/Input/` 为只读原始数据，`Data/{功能}/` 为处理后输出
- EFP 的 `inputs/EFP/` 归入 `Data/Input/EFP/`

---

## 六、Picture 目录规范

```
Picture/
├── MosaicView/               # 瓦片拼接可视化
│   ├── MCD12Q1_PFT_Tiles_A2000001.png
│   └── MCD15A2H_LAI_Mosaic_A2000361.png
└── QC/                       # 质控报告
    └── MODIS_Pipeline_QC_Report.txt
```

**规则**:
- 按 **图表类型/分析目的** 建子目录
- **仅输出 PNG**，不输出 PDF
- 根目录不再散落图片文件（当前 `*.png` 全部归入 `Picture/`）

---

## 七、耦合规则

### 7.1 依赖层次（由上到下）

```
Run_*.py
  └── Core_*.py (同级 import)
        ├── modis_geo_utils.py    (MODIS 地理工具)
        ├── megan_efp/            (MEGAN EFP 算法库)
        └── 标准库 (numpy, pandas, xarray, scipy, pyhdf, netCDF4, sqlite3, matplotlib)
```

### 7.2 禁止的依赖

- `Core_*.py` 不能 `import Run_*.py`
- `Core_*.py` 不能硬编码路径（必须通过参数传入）
- `Run_*.py` 之间不能互相 import
- `Core_*.py` 不能互相循环引用
- 不能 import `Other/` 下的文件
- `Run_*.py` 不能直接 import `megan_efp/`（应通过 `Core_EFP.py`）

### 7.3 允许的依赖

- `Core_A.py` 可以 `import Core_B.py`（如 `Core_PFT` → `Core_MODIS_IO`）
- `Run_*.py` 可以 `sys.path.insert(0, ...)` 后 `import` 任意 Core 模块
- `Core_EFP.py` 可以 `from megan_efp.src import run_M3EFP`

### 7.4 数据流方向

```
外部数据 → Data/Input/ → Core 加载 → 计算 → Data/{功能}/ → 下游模型
                                       ↓
                                   Core 绘图 → Picture/
```

---

## 八、共享工具库

### 8.1 `modis_geo_utils.py` — MODIS 地理工具

```python
"""
modis_geo_utils.py - MODIS 正弦投影与瓦片坐标工具
==================================================
提供函数:
  - modis_tile_xy(h, v) → (xv, yv)          # 瓦片正弦投影米坐标
  - modis_tile_latlon(h, v) → (lat, lon)    # 瓦片经纬度
  - geo2sinu(lat, lon) → (x, y)             # 经纬度 → 正弦投影米坐标
  - sinu2geo(x, y) → (lat, lon)             # 正弦投影米坐标 → 经纬度
  - tile_bbox_xy(h, v) → (xmin, xmax, ymin, ymax)  # 瓦片边界框
  - tile_bbox_latlon(h, v) → (lat_min, lat_max, lon_min, lon_max)
  - check_overlap(bbox1, bbox2) → bool       # 两个 bbox 是否重叠

常量:
  - R = 6371007.181          # MODIS 球半径 (m)
  - PIX_SIZE = 463.3127165   # 像元尺寸 (m)
  - NCOLS = NROWS = 2400     # 瓦片行列数
  - TILE_SIZE = 1111950.52   # 单瓦片边长 (m)
  - X0 = -20015109.354       # 全球左边界 (m)
  - Y0 = 10007554.677        # 全球上边界 (m)
"""
```

**规则**:
- 纯粹的数学/地理工具，不依赖任何项目特定路径
- 从 `MODIS_PFT_FindNear_cn27.py`、`MOIDS_PFT.py`、`Nearest/` 中重复出现的投影代码抽离
- 多处重复定义的正弦投影常数统一到此文件

### 8.2 `megan_efp/` — MEGAN EFP 算法子模块

```
megan_efp/
├── src/
│   ├── run_M3EFP.py         # 数据库创建 + EFP 驱动
│   ├── M3GEFP.py            # SQL 中间表 + 最终表生成
│   └── validation.py        # EF 完整性校验
├── jupyter/
│   └── MEGEFP32.ipynb       # 交互式测试
└── requirements.txt
```

**规则**:
- 独立的算法子模块，不包含项目特定配置
- 场景名称、数据库路径、输入文件路径全部由调用方（`Core_EFP.py`）传入
- 当前 `MEGEFP32/` 目录下的 `inputs/`、`outputs/`、`database/` 分离到 `Data/` 目录中
- `MEGEFP32/env/` 为旧的 virtualenv，不移入，在 `Other/` 归档或删除

---

## 九、当前问题分析

### 9.1 结构问题汇总

| # | 问题 | 具体表现 | 影响 |
|---|------|----------|------|
| 1 | 脚本按分辨率复制 | `MODIS_PFT_FindNear_cn27.py` / `cn09` / `cn03` 三个文件，仅 `SIDE_M`/`HALF_M` 和输出文件名不同 | 修改一处逻辑需同步 3 个文件 |
| 2 | LAI 脚本三份重复 | `v2.1Intov3.2_LAI_CsvTrans_cn03/09/27.py` 仅 `row`/`col`/`grid` 不同 | 同上 |
| 3 | 配置与逻辑混合 | 所有脚本的 `hlines`/`vlines`/`indir`/`outfile` 硬编码在脚本顶部 | 切换场景需修改源码 |
| 4 | PFT 提取多套实现 | `MOIDS_PFT.py`（旧版）、`MODIS_PFT_FindNear_cn*.py`（正方形邻域）、`Nearest/` 目录（最近邻变体）共存 | 不知道哪个是权威版本 |
| 5 | 可视化脚本重复 | `hdfview_data_merged.py` 与 `hdfview_data_look.py` 完全一致（393行） | 维护两份相同代码 |
| 6 | 数据散落根目录 | `PFT_frac_*.csv`、`LAI3_*.csv`、`*.png` 散落在项目根目录 | 难以区分输入/输出/临时文件 |
| 7 | MEGEFP32 深层嵌套 | `MEGEFP32/inputs/EFP/`、`MEGEFP32/outputs/`、`MEGEFP32/database/` 混在模块内 | 数据与代码混杂 |
| 8 | 探索性脚本未归档 | `3_ChatGPT.py`、`3_Deepseek_QingXie.py`、`MEGANv3.2/11.14PPTtest/` | 干扰主线理解 |
| 9 | 命名不统一 | 有的用 `CELLID` 有的用 `gridID`，有的 `X`/`Y` 0-based 有的 1-based | 下游对接易出错 |
| 10 | 无 `__init__.py` 或包结构 | 纯脚本堆叠，无模块化导入 | import 依赖 `sys.path` 操作 |

### 9.2 当前文件 → 目标位置对照表

| 当前文件 | → 目标位置 | 说明 |
|----------|-----------|------|
| `MODIS_PFT_FindNear_cn27.py` | → `Core_PFT.py`（逻辑）+ `Run_PFT_Percentage.py`（配置） | 三合一，resolution 参数化 |
| `MODIS_PFT_FindNear_cn09.py` | → 同上 | 合并 |
| `MODIS_PFT_FindNear_cn03.py` | → 同上 | 合并 |
| `MOIDS_PFT.py` | → `Other/` 归档 | 旧版实现 |
| `Nearest/MODIS_PFT_Nearest_cn*.py` | → `Other/` 归档 | 与正方形邻域功能重叠 |
| `PFT_cn03nc_Extract.py` | → `Core_PFT.py`（`extract_single_class` 函数） | 三合一 |
| `PFT_cn09nc_Extract.py` | → 同上 | 合并 |
| `PFT_cn27nc_Extract.py` | → 同上 | 合并 |
| `v2.1Intov3.2_PFT_CsvTrans_Muti.py` | → `Core_GrowthForm.py` + `Run_GrowthForm.py` | 分离配置与逻辑 |
| `v2.1Intov3.2_LAI_CsvTrans_cn27.py` | → `Core_LAI.py`（逻辑）+ `Run_LAI.py`（配置） | 三合一 |
| `v2.1Intov3.2_LAI_CsvTrans_cn09.py` | → 同上 | 合并 |
| `v2.1Intov3.2_LAI_CsvTrans_cn03.py` | → 同上 | 合并 |
| `hdfview_merged_PFT.py` | → `Core_MosaicView.py` + `Run_MosaicView.py` | 分离配置与逻辑 |
| `hdfview_data_merged.py` | → 同上（合并入 `Core_MosaicView.py`） | 保留一份，删 `_look` 副本 |
| `hdfview_data_look.py` | → 删除（与 `_merged` 重复） | 393 行完全相同 |
| `hdfview_merged_LAI.py` | → `Other/` 归档 或 合并入 `Core_MosaicView.py` | 视需要决定 |
| `hdfview_Look.py` | → `Core_MODIS_IO.py`（`inspect_hdf` 工具函数） | 19 行小工具 |
| `MEGEFP32/MEGAN_EFP.py` | → `Core_EFP.py` + `Run_EFP.py` | 分离场景配置与调用 |
| `MEGEFP32/src/run_M3EFP.py` | → `megan_efp/src/run_M3EFP.py` | 算法子模块 |
| `MEGEFP32/src/M3GEFP.py` | → `megan_efp/src/M3GEFP.py` | 算法子模块 |
| `MEGEFP32/src/validation.py` | → `megan_efp/src/validation.py` | 算法子模块 |
| `MEGEFP32/jupyter/` | → `megan_efp/jupyter/` | 算法子模块 |
| `MEGEFP32/inputs/EFP/` | → `Data/Input/EFP/` | 数据与代码分离 |
| `MEGEFP32/outputs/` | → `Data/EFP/Output/` | 数据与代码分离 |
| `MEGEFP32/database/` | → `Data/EFP/Database/` | 数据与代码分离 |
| `MEGEFP32/env/` | → `Other/` 归档或删除 | 旧 virtualenv |
| `MEGEFP32/11.14PPTtest/` | → `Other/` 归档 | 探索性测试 |
| `verify_modis_pipeline.py` | → `Core_Verification.py` + `Run_Verify.py` | 分离 QC 逻辑与路径配置 |
| `3_ChatGPT.py` | → `Other/` 归档 | 探索性生成 |
| `3_Deepseek_QingXie.py` | → `Other/` 归档 | 探索性生成 |
| `MEGANv3.2/` | → `Other/` 归档 | MEGAN v3.2 旧版实验 |
| `*.png` (根目录散落) | → `Picture/MosaicView/` 或 `Picture/QC/` | 图片统一管理 |
| `PFT_frac_*.csv` / `PFT_frac_*.nc` | → `Data/PFT/{year}/` | 数据归类 |
| `LAI3_*.csv` (MEGANv3.2/) | → `Data/LAI/{year}/` | 数据归类 |
| `grid_growth_form_*.csv` | → `Data/GrowthForm/{year}/` | 数据归类 |
| `lai_cn*_GuangDong_2000.csv` | → `Data/LAI/2000/` | 数据归类 |
| `pft_cn*_GuangDong_2000.csv` | → `Data/PFT/2000/` | 数据归类 |
| `GRIDCRO2D_*` | → 保留在根目录（CMAQ 网格文件） | 外部输入，不归入 Data/ |

---

## 十、重构步骤建议

### 第一阶段：目录骨架搭建

1. 创建顶层目录：`Data/Input/{HDF,LAI,EFP}/`、`Data/{PFT,GrowthForm,LAI,EFP/{Database,Output}}/`、`Picture/{MosaicView,QC}/`、`Other/`、`docs/`
2. 创建 `modis_geo_utils.py`：从现有脚本抽离所有正弦投影/瓦片坐标函数
3. 创建 `megan_efp/`：从 `MEGEFP32/src/` 迁移 `run_M3EFP.py`、`M3GEFP.py`、`validation.py`
4. 更新 `.gitignore`：排除 `Data/`（大文件）、`Picture/`（生成物）、`__pycache__/`、数据库文件

### 第二阶段：Core 模块提取

按依赖顺序提取，每提取一个模块立即编写配套 Run 入口验证：

1. `Core_MODIS_IO.py` — HDF 读取、坐标计算、瓦片 bbox（最底层，无项目依赖）
2. `Core_PFT.py` — PFT 百分比统计 + 最近邻提取（依赖 Core_MODIS_IO）
3. `Core_GrowthForm.py` — Growth Form 转换（依赖 Core_PFT 输出格式）
4. `Core_LAI.py` — LAI 时序拼接（独立，依赖标准库）
5. `Core_MosaicView.py` — 瓦片拼接可视化（依赖 Core_MODIS_IO）
6. `Core_EFP.py` — EFP 入口封装（依赖 megan_efp/）
7. `Core_Verification.py` — 全链路 QC（依赖所有 Core 输出格式）

### 第三阶段：Run 入口编写

对应每个 Core 模块创建 Run 入口：

1. `Run_MosaicView.py` — 配置瓦片范围/年份，调用 `Core_MosaicView`
2. `Run_PFT_Percentage.py` — 配置分辨率列表，批量调用 `Core_PFT.run_percentage`
3. `Run_PFT_Nearest.py` — 配置分辨率，调用 `Core_PFT.run_nearest`
4. `Run_GrowthForm.py` — 配置 PFT 输入→GrowthForm 输出映射
5. `Run_LAI.py` — 配置分辨率列表，批量调用 `Core_LAI`
6. `Run_EFP.py` — 配置场景名/分辨率，调用 `Core_EFP`
7. `Run_Verify.py` — 配置所有中间产物路径，调用 `Core_Verification`

### 第四阶段：数据迁移与归档

1. 移动 `PFT_frac_*.csv` / `*.nc` → `Data/PFT/2000/`
2. 移动 `grid_growth_form_*.csv` → `Data/GrowthForm/2000/`
3. 移动 `LAI3_*.csv` → `Data/LAI/2000/`
4. 移动 `MEGEFP32/inputs/EFP/*` → `Data/Input/EFP/`
5. 移动 `MEGEFP32/outputs/*` → `Data/EFP/Output/`
6. 移动 `MEGEFP32/database/*` → `Data/EFP/Database/`
7. 移动根目录 `*.png` → `Picture/MosaicView/`
8. 归档废弃脚本 → `Other/`
9. 验证所有 Run 入口可正常执行

---

## 十一、关键规范清单

| # | 规范 | 说明 |
|---|------|------|
| 1 | Core/Run 分离 | Core 纯函数、Run 纯配置 |
| 2 | 路径参数化 | Core 不出现任何硬编码路径 |
| 3 | 分辨率参数化 | 不复刻脚本，通过参数 `resolution` 区分 3km/9km/27km |
| 4 | 平面顶层 | 所有 Core_*/Run_* 在项目根目录 |
| 5 | Data/ 按功能 | 子目录以数据来源/用途命名 |
| 6 | Picture/ 按图类型 | 子目录以图表类型命名 |
| 7 | 仅 PNG 输出 | 不生成 PDF |
| 8 | 代码与数据分离 | megan_efp/ 中不含 inputs/outputs/database |
| 9 | Other/ 归档 | 废弃代码不移入 Core/Run |
| 10 | 单一地理工具源 | 正弦投影常数只定义在 `modis_geo_utils.py` 一处 |
| 11 | 单模块职责 | 每个 Core 不超过 500 行，超过则考虑拆分 |
| 12 | 禁止用户交互 | Core 中不出现 `input()` 等交互函数（EFP 校验逻辑改为参数控制） |
