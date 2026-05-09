"""
Core_MODIS_IO.py - MODIS 数据读取与网格工具
=============================================
提供 MODIS HDF/NC 文件读取及 CMAQ 目标网格加载函数。
所有路径通过函数签名传入，无硬编码。

提供函数:
  # HDF 读取
  - read_modis_lc_type5(hdf_path) -> np.ndarray
      读取 MCD12Q1 HDF 中 LC_Type5 数据（int16, 无效=-1）
  - read_modis_tile_from_nc(nc_path) -> (lat1d, lon1d, pft1d)
      读取 NCL 转换后的 MCD12Q1 NC 文件
  - read_single_mcd12q1_pft(file_path) -> dict
      读取单个 MCD12Q1 HDF 文件（含元数据/地理范围）
  - read_single_mcd15a2h(file_path) -> dict
      读取单个 MCD15A2H HDF 文件（含 LAI + QC + 地理范围）

  # CMAQ 网格
  - read_cmaq_grid(grid_nc) -> (lat, lon, ny, nx)
      读取 CMAQ 规则网格经纬度

  # 工具
  - inspect_hdf(file_path) -> None
      打印 HDF 文件的数据集与全局属性
  - ll_to_unitvec(lat, lon) -> np.ndarray
      经纬度 → 单位球面坐标 (N,3)
"""

import os
import numpy as np
import xarray as xr

from modis_geo_utils import (
    modis_tile_latlon, parse_tile_from_filename,
    R, PIX_SIZE, NCOLS, NROWS,
)


# ===========================
# MCD12Q1 LC_Type5 (PFT)
# ===========================

def read_modis_lc_type5(hdf_path: str) -> np.ndarray:
    """读取 MODIS MCD12Q1 HDF 里的 LC_Type5（uint8 → int16，无效置为 -1）。

    参数:
        hdf_path: MCD12Q1 .hdf 文件路径

    返回:
        arr: (2400, 2400) int16 数组，有效值为 0..11，无效为 -1
    """
    from pyhdf.SD import SD, SDC
    hdf = SD(hdf_path, SDC.READ)
    arr = np.array(hdf.select("LC_Type5")[:], dtype=np.int16)
    hdf.end()
    arr[(arr < 0) | (arr >= 254)] = -1
    return arr


def read_modis_tile_from_nc(nc_path: str):
    """读取 NCL 转换后的 MCD12Q1 NC 文件，提取 LC_Type5 并自动生成经纬度。

    参数:
        nc_path: .nc 文件路径（文件名含 h/v 信息）

    返回:
        lat1d, lon1d, pft1d: 一维有效像元经纬度和 PFT 分类
        无有效像元时返回 (None, None, None)
    """
    h, v = parse_tile_from_filename(nc_path)

    ds = xr.open_dataset(nc_path, mask_and_scale=False, decode_times=False)
    if "LC_Type5" not in ds.variables:
        ds.close()
        raise RuntimeError(f"{os.path.basename(nc_path)} 中缺少 LC_Type5 变量")

    arr = ds["LC_Type5"].values
    ds.close()

    # NaN/float → int16
    if np.issubdtype(arr.dtype, np.floating):
        miss = ~np.isfinite(arr)
        arr = arr.astype(np.int16, copy=False)
        arr[miss] = -1
    else:
        arr = arr.astype(np.int16, copy=False)

    # byte → uint8
    arr = np.where(arr < 0, arr + 256, arr)

    # 生成该瓦片的经纬度
    lat2d, lon2d = modis_tile_latlon(h, v)

    # 仅保留真实 PFT 类别（0..11）
    mask = (arr >= 0) & (arr <= 11)
    if not np.any(mask):
        return None, None, None

    pft1d = arr[mask].ravel().astype(np.int16)
    lat1d = lat2d[mask].ravel()
    lon1d = lon2d[mask].ravel()
    return lat1d, lon1d, pft1d


def read_single_mcd12q1_pft(file_path: str) -> dict:
    """读取单个 MCD12Q1 文件的 LC_Type5 数据及地理范围。

    参数:
        file_path: MCD12Q1 .hdf 文件路径

    返回:
        dict: {
            'tile': 瓦片标识 (e.g. "h23v03"),
            'data': (2400,2400) float32, NaN 为无效,
            'lat_range': (lat_min, lat_max),
            'lon_range': (lon_min, lon_max),
        }
    """
    from pyhdf.SD import SD, SDC

    try:
        hdf = SD(file_path, SDC.READ)
    except Exception as e:
        print(f"Cannot open {file_path}: {e}")
        return None

    try:
        ds = hdf.select('LC_Type5')
        data = ds[:].astype(np.float32)
        fill_value = ds.attributes().get('_FillValue', 255)
        data[data == fill_value] = np.nan
        data[data > 250] = np.nan
    except Exception as e:
        print(f"Read LC_Type5 failed: {e}")
        hdf.end()
        return None

    # 获取地理范围
    try:
        meta = hdf.attr('StructMetadata.0').get().decode('utf-8', errors='ignore')
        lat_min = float(meta.split('SOUTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lat_max = float(meta.split('NORTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lon_min = float(meta.split('WESTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lon_max = float(meta.split('EASTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
    except Exception:
        h, v = parse_tile_from_filename(file_path)
        lon_min = -180 + h * 10
        lon_max = lon_min + 10
        lat_max = 90 - v * 10
        lat_min = lat_max - 10

    hdf.end()
    tile = os.path.basename(file_path).split('.')[2]
    return {
        'tile': tile,
        'data': data,
        'lat_range': (lat_min, lat_max),
        'lon_range': (lon_min, lon_max),
    }


# ===========================
# MCD15A2H LAI
# ===========================

def read_single_mcd15a2h(file_path: str) -> dict:
    """读取单个 MCD15A2H 文件的 LAI 数据、QC 数据及地理范围。

    参数:
        file_path: MCD15A2H .hdf 文件路径

    返回:
        dict: {
            'lai_data': LAI 数组,
            'qc_data': QC 数组,
            'lai_attrs': LAI 属性字典,
            'qc_attrs': QC 属性字典,
            'global_attrs': 全局属性,
            'struct_metadata': 结构元数据字符串,
            'lat_range': (lat_min, lat_max),
            'lon_range': (lon_min, lon_max),
            'tile': 瓦片标识,
            'original_shape': 数据形状,
        }
    """
    from pyhdf.SD import SD, SDC

    try:
        hdf = SD(file_path, SDC.READ)
    except Exception:
        print(f"Cannot open {file_path}")
        return None

    try:
        # LAI
        lai = hdf.select('Lai_500m')
        lai_data = lai[:].astype(np.float32)
        lai_attrs = lai.attributes()
        scale_factor = lai_attrs.get('scale_factor', 0.1)
        fill_value = lai_attrs.get('_FillValue', 255)
        lai_data = lai_data * scale_factor
        lai_data[lai_data == fill_value * scale_factor] = np.nan

        # QC
        qc = hdf.select('FparLai_QC')
        qc_data = qc[:].astype(np.uint8)
        qc_attrs = qc.attributes()
        valid_mask = (qc_data & 0b11) <= 0b11
        lai_data = np.where(valid_mask, lai_data, np.nan)

        # 全局属性
        global_attrs = hdf.attributes()
        struct_meta = ""

        # 地理范围
        try:
            struct_meta_bytes = hdf.attr('StructMetadata.0').get()
            struct_meta = struct_meta_bytes.decode('utf-8', errors='ignore')
            lat_min = float(struct_meta.split('SOUTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
            lat_max = float(struct_meta.split('NORTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
            lon_min = float(struct_meta.split('WESTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
            lon_max = float(struct_meta.split('EASTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        except Exception:
            h, v = parse_tile_from_filename(file_path)
            lon_min = -180 + h * 10.0
            lon_max = lon_min + 10.0
            lat_max = 90 - v * 10.0
            lat_min = lat_max - 10.0

        hdf.end()

        h, v = parse_tile_from_filename(file_path)
        return {
            'lai_data': lai_data,
            'lai_attrs': lai_attrs,
            'qc_data': qc_data,
            'qc_attrs': qc_attrs,
            'global_attrs': global_attrs,
            'struct_metadata': struct_meta,
            'lat_range': (lat_min, lat_max),
            'lon_range': (lon_min, lon_max),
            'tile': f"h{h:02d}v{v:02d}",
            'original_shape': lai_data.shape,
        }

    except Exception as e:
        print(f"Read {file_path} data failed: {e}")
        try:
            hdf.end()
        except Exception:
            pass
        return None


# ===========================
# CMAQ 网格
# ===========================

def read_cmaq_grid(grid_nc: str):
    """读取 CMAQ 规则网格的经纬度。

    参数:
        grid_nc: GRIDCRO2D NetCDF 文件路径

    返回:
        lat: (ny, nx) 纬度数组
        lon: (ny, nx) 经度数组
        ny, nx: 行列数
    """
    ds_grid = xr.open_dataset(grid_nc)
    lat = ds_grid["LAT"].isel(TSTEP=0, LAY=0).values
    lon = ds_grid["LON"].isel(TSTEP=0, LAY=0).values
    ny, nx = lat.shape
    ds_grid.close()
    return lat, lon, ny, nx


# ===========================
# 工具
# ===========================

def inspect_hdf(file_path: str) -> None:
    """打印 HDF 文件的数据集与全局属性（快速查看工具）。

    参数:
        file_path: .hdf 文件路径
    """
    from pyhdf.SD import SD, SDC
    hdf = SD(file_path, SDC.READ)

    print("Datasets in file:")
    datasets = hdf.datasets()
    for name, info in datasets.items():
        print(f"  {name} - dims: {info[0]}, type: {info[3]}")

    print("\nGlobal attributes:")
    for attr in hdf.attributes():
        print(f"  {attr}: {hdf.attributes()[attr]}")
    hdf.end()


def ll_to_unitvec(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    """经纬度 → 单位球面上的 3D 坐标（用于 KDTree 在大圆弧上查询）。

    参数:
        lat_deg: 纬度 (度)
        lon_deg: 经度 (度)

    返回:
        xyz: (N, 3) 单位球面坐标
    """
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    cl = np.cos(lat)
    return np.stack([cl * np.cos(lon), cl * np.sin(lon), np.sin(lat)], axis=-1)
