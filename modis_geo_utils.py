"""
modis_geo_utils.py - MODIS 正弦投影与瓦片坐标工具
==================================================
提供 MODIS Sinusoidal 投影的全部坐标转换与瓦片边界计算函数。

所有函数均为纯数学/地理计算，不依赖任何项目特定路径。

常量:
    R         - MODIS 球半径 (m)
    PIX_SIZE  - 像元尺寸 (m)
    NCOLS     - 瓦片列数
    NROWS     - 瓦片行数
    TILE_SIZE - 单瓦片边长 (m)
    X0        - 全球左边界正弦投影 x (m)
    Y0        - 全球上边界正弦投影 y (m)

提供函数:
    modis_tile_xy(h, v) -> (xv, yv)
        返回瓦片 (h,v) 的像元中心正弦投影米坐标网格。

    modis_tile_latlon(h, v) -> (lat, lon)
        返回瓦片 (h,v) 的像元中心经纬度网格。

    geo2sinu(lat, lon) -> (x, y)
        经纬度 → 正弦投影米坐标。

    sinu2geo(x, y) -> (lat, lon)
        正弦投影米坐标 → 经纬度。

    tile_bbox_xy(h, v) -> (xmin, xmax, ymin, ymax)
        返回瓦片在正弦投影米坐标下的边界框。

    tile_bbox_latlon(h, v) -> (lat_min, lat_max, lon_min, lon_max)
        返回瓦片在经纬度下的边界框。

    check_overlap(bbox1, bbox2) -> bool
        判断两个 bbox 是否有重叠。
"""

import numpy as np

# ===========================
# MODIS 正弦投影常量
# ===========================
R         = 6371007.181        # 球半径 (m)
PIX_SIZE  = 463.3127165        # 像元尺寸 (m)
NCOLS     = 2400               # 瓦片列数
NROWS     = 2400               # 瓦片行数
TILE_SIZE = PIX_SIZE * NCOLS   # 1,111,950.519667 m
X0        = -20015109.354      # 全球左边界 (m)
Y0        =  10007554.677      # 全球上边界 (m)

# ===========================
# 瓦片像元中心坐标
# ===========================

def modis_tile_xy(h: int, v: int):
    """返回瓦片 h,v 的像元中心正弦投影坐标 (xv, yv)，单位: 米。

    参数:
        h: 横向瓦片编号 (0-based, e.g. 23..30)
        v: 纵向瓦片编号 (0-based, e.g. 3..8)

    返回:
        xv: (NROWS, NCOLS) x 坐标 (米)
        yv: (NROWS, NCOLS) y 坐标 (米)
    """
    x_ul = X0 + h * TILE_SIZE            # 瓦片左上角 x
    y_ul = Y0 - v * TILE_SIZE            # 瓦片左上角 y（v 向下）
    x = x_ul + (np.arange(NCOLS) + 0.5) * PIX_SIZE
    y = y_ul - (np.arange(NROWS) + 0.5) * PIX_SIZE
    xv, yv = np.meshgrid(x, y)
    return xv, yv


def modis_tile_latlon(h: int, v: int):
    """返回瓦片 (h,v) 的像元中心经纬度网格。

    参数:
        h: 横向瓦片编号
        v: 纵向瓦片编号

    返回:
        lat: (NROWS, NCOLS) 纬度 (度)
        lon: (NROWS, NCOLS) 经度 (度)
    """
    x_ul = X0 + h * TILE_SIZE
    y_ul = Y0 - v * TILE_SIZE
    x = x_ul + (np.arange(NCOLS) + 0.5) * PIX_SIZE
    y = y_ul - (np.arange(NROWS) + 0.5) * PIX_SIZE
    xv, yv = np.meshgrid(x, y)

    # 反投影：Sinusoidal → geographic
    lat_rad = yv / R
    coslat = np.cos(lat_rad)
    coslat[coslat == 0] = np.finfo(np.float64).eps
    lon_rad = xv / (R * coslat)

    lat = np.rad2deg(lat_rad)
    lon = np.rad2deg(lon_rad)
    return lat, lon


# ===========================
# 经纬度 ↔ 正弦投影互转
# ===========================

def geo2sinu(lat_deg: np.ndarray, lon_deg: np.ndarray):
    """地理经纬度 → 正弦投影米坐标（与 MODIS Sine 一致的球面公式）。

    参数:
        lat_deg: 纬度 (度)
        lon_deg: 经度 (度)

    返回:
        x: 正弦投影 x 坐标 (米)
        y: 正弦投影 y 坐标 (米)
    """
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    x = R * lon * np.cos(lat)
    y = R * lat
    return x, y


def sinu2geo(x: np.ndarray, y: np.ndarray):
    """正弦投影米坐标 → 地理经纬度。

    参数:
        x: 正弦投影 x 坐标 (米)
        y: 正弦投影 y 坐标 (米)

    返回:
        lat: 纬度 (度)
        lon: 经度 (度)
    """
    lat = np.rad2deg(y / R)
    lon = np.rad2deg(x / (R * np.cos(np.deg2rad(lat))))
    return lat, lon


# ===========================
# 瓦片边界框
# ===========================

def tile_bbox_xy(h: int, v: int):
    """返回瓦片在正弦投影米坐标下的边界框。

    参数:
        h: 横向瓦片编号
        v: 纵向瓦片编号

    返回:
        (xmin, xmax, ymin, ymax): 单位米
        注意 y 向下递减，ymin < ymax
    """
    x0 = X0 + h * TILE_SIZE
    y0 = Y0 - v * TILE_SIZE
    x1 = x0 + TILE_SIZE
    y1 = y0 - TILE_SIZE
    return (min(x0, x1), max(x0, x1), min(y1, y0), max(y1, y0))


def tile_bbox_latlon(h: int, v: int):
    """返回瓦片在经纬度下的边界框。

    参数:
        h: 横向瓦片编号
        v: 纵向瓦片编号

    返回:
        (lat_min, lat_max, lon_min, lon_max): 单位度
        基于四个角点计算，逼近实际范围。
    """
    # 取四个角点 + 中心点来计算近似范围
    corners_xy = [
        (X0 + h * TILE_SIZE,               Y0 - v * TILE_SIZE),               # 左上
        (X0 + (h + 1) * TILE_SIZE,         Y0 - v * TILE_SIZE),               # 右上
        (X0 + h * TILE_SIZE,               Y0 - (v + 1) * TILE_SIZE),         # 左下
        (X0 + (h + 1) * TILE_SIZE,         Y0 - (v + 1) * TILE_SIZE),         # 右下
    ]
    xs = np.array([c[0] for c in corners_xy])
    ys = np.array([c[1] for c in corners_xy])
    lats, lons = sinu2geo(xs, ys)
    return (float(np.min(lats)), float(np.max(lats)),
            float(np.min(lons)), float(np.max(lons)))


# ===========================
# 工具函数
# ===========================

def check_overlap(bbox1, bbox2) -> bool:
    """判断两个 bbox 是否有重叠。

    参数:
        bbox1, bbox2: (xmin, xmax, ymin, ymax) 或 (lat_min, lat_max, lon_min, lon_max)

    返回:
        True 如果有重叠
    """
    xmin1, xmax1, ymin1, ymax1 = bbox1
    xmin2, xmax2, ymin2, ymax2 = bbox2
    return not (xmax1 < xmin2 or xmin1 > xmax2 or ymax1 < ymin2 or ymin1 > ymax2)


def parse_tile_from_filename(filename: str):
    """从 MODIS HDF 文件名中解析 (h, v) 瓦片编号。

    示例:
        "MCD12Q1.A2000001.h23v03.061.xxx.hdf" → (23, 3)

    参数:
        filename: MODIS HDF 文件名

    返回:
        (h, v): 整数元组
    """
    import os
    base = os.path.basename(filename)
    try:
        h = int(base.split(".h")[1][:2])
        v = int(base.split("v")[1][:2])
    except (IndexError, ValueError):
        raise ValueError(f"无法从文件名解析 h/v: {base}")
    return h, v


# ===========================
# LC_Type5 分类名称常量
# ===========================

CLASS_NAMES = [
    "Water",
    "Evergreen Needleleaf trees",
    "Evergreen Broadleaf trees",
    "Deciduous Needleleaf trees",
    "Deciduous Broadleaf trees",
    "Shrub",
    "Grass",
    "Cereal crops",
    "Broad-leaf crops",
    "Urban and built-up",
    "Snow and ice",
    "Barren or sparse vegetation",
]

# Growth Form 聚合映射
TREE_COLS = [
    "Evergreen Needleleaf trees",
    "Evergreen Broadleaf trees",
    "Deciduous Needleleaf trees",
    "Deciduous Broadleaf trees",
]
CROP_COLS = ["Cereal crops", "Broad-leaf crops"]
SHRUB_COL = "Shrub"
HERB_COL = "Grass"
