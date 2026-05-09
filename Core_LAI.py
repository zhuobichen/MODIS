"""
Core_LAI.py - MODIS LAI 时序拼接
==================================
从 MODIS LAI 8-day NetCDF 切片批量读取并拼接为单 CSV 文件。

所有路径和参数通过函数签名传入，无硬编码。

提供函数:
  - run_lai_pipeline() -> str
      完整管道: 日期校验 → 批量读取 → 展平拼接 → 输出 CSV
"""

import os
import glob
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from typing import List, Optional


def _build_expected_days() -> List[str]:
    """生成 MODIS 8-day LAI 期望日期序列 (001, 009, ..., 361) 共 46 期。"""
    return [f"{1 + 8 * i:03d}" for i in range(46)]


def _validate_day_sequence(expected_days: List[str]) -> None:
    """校验日期序列完整性：每期间隔 8 天，末期为 361。"""
    if expected_days[-1] != "361":
        raise ValueError("Date sequence incomplete: last day should be 361")
    ndays = len(expected_days)
    for i in range(1, ndays):
        prev = int(expected_days[i - 1])
        curr = int(expected_days[i])
        if curr - prev != 8:
            raise ValueError(
                f"Date interval error: {expected_days[i-1]} to {expected_days[i]} "
                f"is not 8 days"
            )


def run_lai_pipeline(
    year: str,
    resolution: str,
    project: str,
    input_dir: str,
    output_dir: str,
    rows: int,
    cols: int,
    use_1based_coords: bool = True,
) -> str:
    """LAI NetCDF 时序拼接管道。

    参数:
        year: 年份 (e.g. "2000")
        resolution: 分辨率标识 (e.g. "cn27")
        project: 项目/区域名 (e.g. "GuangDong")
        input_dir: LAI NetCDF 切片目录
        output_dir: 输出 CSV 目录
        rows: 网格行数
        cols: 网格列数
        use_1based_coords: X/Y 是否为 1-based（默认 True）

    返回:
        str: 输出 CSV 文件路径

    输入文件命名格式:
        {input_dir}/MODIS_LAI_{year}{doy}_{project}_{resolution}.nc

    输出文件命名格式:
        {output_dir}/LAI3_{resolution}.csv
    """
    expected_days = _build_expected_days()
    ndays = len(expected_days)

    print("Validating date sequence...")
    _validate_day_sequence(expected_days)
    print("Date sequence valid.")

    # 检查输入文件
    input_files = []
    for day in expected_days:
        pattern = os.path.join(
            input_dir,
            f"MODIS_LAI_{year}{day}_{project}_{resolution}.nc",
        )
        files = glob.glob(pattern)
        if not files:
            raise FileNotFoundError(f"Missing file: {pattern}")
        input_files.append(files[0])

    # 读取经纬度（从第一个文件）
    print("Reading lat/lon info...")
    first_file = Dataset(input_files[0], "r")
    lat = first_file.variables["lat"][:]
    lon = first_file.variables["lon"][:]
    first_file.close()

    if lat.shape != (rows, cols) or lon.shape != (rows, cols):
        raise ValueError(
            f"Lat/lon dimension mismatch: expected ({rows},{cols}), "
            f"got lat{lat.shape}, lon{lon.shape}"
        )

    # 处理 LAI
    total_pixels = rows * cols
    lai_data = np.full((total_pixels, ndays), 0.0)
    print("Processing LAI data...")

    for d, (day, file_path) in enumerate(zip(expected_days, input_files)):
        print(f"  [{d+1}/{ndays}] {os.path.basename(file_path)}")
        with Dataset(file_path, "r") as nc:
            lai = nc.variables["LAI"][:]

            if lai.shape != (rows, cols):
                raise ValueError(
                    f"LAI dimension mismatch: expected ({rows},{cols}), "
                    f"got {lai.shape}"
                )

            lai_flat = lai.flatten(order='C')
            fill_value = getattr(nc.variables["LAI"], "_FillValue", -999.0)
            lai_flat = np.where(lai_flat == fill_value, 0.0, lai_flat)
            lai_data[:, d] = lai_flat

    # 生成坐标
    cell_ids = np.arange(1, total_pixels + 1)
    if use_1based_coords:
        y_coords, x_coords = np.mgrid[1:rows+1, 1:cols+1]
    else:
        y_coords, x_coords = np.mgrid[0:rows, 0:cols]
    x_flat = x_coords.flatten(order='C')
    y_flat = y_coords.flatten(order='C')

    lat_flat = lat.flatten(order='C')
    lon_flat = lon.flatten(order='C')

    # 构建 DataFrame
    print("Building output DataFrame...")
    df = pd.DataFrame({
        "CELL_ID": cell_ids,
        "X": x_flat,
        "Y": y_flat,
        "LAT": lat_flat,
        "LONG": lon_flat,
    })

    for i in range(ndays):
        df[f"LAI{i+1:02d}"] = lai_data[:, i]

    # 检查缺失
    missing_count = df.iloc[:, 5:].eq(0.0).sum().sum()
    if missing_count > 0:
        print(f"Warning: {missing_count} missing values (filled with 0)")

    # 保存
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, f"LAI3_{resolution}.csv")
    df.to_csv(outfile, index=False)
    print(f"Done. Output: {outfile}  ({len(df)} rows, {ndays} LAI columns)")
    return outfile
