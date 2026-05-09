"""
Core_EFP.py - MEGAN 排放因子处理入口
======================================
封装 MEGAN EFP (Emission Factor Processor) 调用，
提供参数化场景配置 → 调用 megan_efp/ 算法库。

所有路径通过函数签名传入，无硬编码。
依赖: megan_efp/ 算法子模块

提供函数:
  - run_efp_pipeline() -> str
      完整管道: 校验 → 建库 → 中间表 → 输出 CSV
"""

import os
import sys
from typing import Optional


def run_efp_pipeline(
    scen_name: str,
    input_dir: str,
    output_dir: str,
    database_dir: str,
    ef_file: str = "EFv210806.csv",
    ecotype_crop: str = "SpeciationCrop210806.csv",
    ecotype_herb: str = "SpeciationHerb210806.csv",
    ecotype_shrub: str = "SpeciationShrub210806.csv",
    ecotype_tree: str = "SpeciationTree210725.csv",
    grid_ecotype: Optional[str] = None,
    grid_growth_form: Optional[str] = None,
    ef_classes: int = 19,
    ldf_classes0: int = 3,
    ldf_classes1: int = 6,
    megan_efp_dir: Optional[str] = None,
) -> str:
    """MEGAN EFP 排放因子计算管道。

    参数:
        scen_name: 场景名称 (e.g. "GD_cn27")
        input_dir: EFP 输入 CSV 目录
        output_dir: 输出 CSV 目录
        database_dir: SQLite 数据库目录
        ef_file: 排放因子表文件名
        ecotype_crop/herb/shrub/tree: 物种组成文件名
        grid_ecotype: 生态类型网格文件名 (默认: grid_ecotype.{scen}.csv)
        grid_growth_form: Growth Form 网格文件名 (默认: grid_growth_form.{scen}.csv)
        ef_classes: EF 类别数 (默认 19)
        ldf_classes0: LDF 起始索引 (默认 3)
        ldf_classes1: LDF 结束索引 (默认 6)
        megan_efp_dir: megan_efp 目录路径 (默认: 自动检测)

    返回:
        str: 输出 CSV 文件路径
    """
    # 确定 megan_efp 路径
    if megan_efp_dir is None:
        megan_efp_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "megan_efp",
        )

    src_dir = os.path.join(megan_efp_dir, "src")
    if not os.path.isdir(src_dir):
        raise FileNotFoundError(f"megan_efp/src not found: {src_dir}")

    # 添加 megan_efp/src 到 Python 路径
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)

    import run_M3EFP as efp

    # 默认文件名
    if grid_ecotype is None:
        grid_ecotype = f"grid_ecotype.{scen_name}.csv"
    if grid_growth_form is None:
        grid_growth_form = f"grid_growth_form.{scen_name}.csv"

    # 构建路径
    database_path = os.path.join(database_dir, f"M3GEFP_database.{scen_name}.db")

    print(f"\n=== MEGAN EFP ===")
    print(f"Scenario: {scen_name}")
    print(f"Input dir: {input_dir}")
    print(f"Database: {database_path}")
    print(f"Output dir: {output_dir}")

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(database_dir, exist_ok=True)

    # 调用 MEGAN EFP 驱动
    efp.m3efp_driver(
        scen_name,
        input_dir if input_dir.endswith('/') else input_dir + '/',
        ef_file,
        ecotype_crop,
        ecotype_herb,
        ecotype_tree,
        ecotype_shrub,
        grid_ecotype,
        grid_growth_form,
        database_path,
        1,
        ef_classes,
        ldf_classes0,
        ldf_classes1,
        output_dir if output_dir.endswith('/') else output_dir + '/',
    )

    output_csv = os.path.join(output_dir, f"OutputGridEF.{scen_name}.csv")
    print(f"\nEFP output: {output_csv}")
    return output_csv
