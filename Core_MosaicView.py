"""
Core_MosaicView.py - MODIS 瓦片拼接可视化
=============================================
MODIS 瓦片拼接 + 可视化：PFT 分类图 + LAI 拼接图，带瓦片边界标注。

所有路径通过函数签名传入，无硬编码。

提供函数:
  - mosaic_pft_tiles() -> dict
      拼接 MCD12Q1 PFT 瓦片并绘图
  - mosaic_lai_tiles() -> dict
      拼接 MCD15A2H LAI 瓦片（可选保存 HDF）
  - run_mosaic_view_pipeline() -> dict
      管道函数：拼图 + 保存
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from typing import List, Optional

from Core_MODIS_IO import read_single_mcd12q1_pft, read_single_mcd15a2h


# PFT 分类颜色映射 (LC_Type5 12 类)
PFT_COLORS = [
    '#0000FF', '#006400', '#228B22', '#7FFF00', '#ADFF2F', '#FFD700',
    '#DAA520', '#FF8C00', '#CD5C5C', '#8B4513', '#D3D3D3', '#FFFFFF'
]

PIXEL_PER_DEGREE = 2400 / 10.0  # MCD12Q1 / MCD15A2H 分辨率


def mosaic_pft_tiles(
    file_paths: List[str],
    hlines: List[int],
    vlines: List[int],
) -> Optional[dict]:
    """拼接 MCD12Q1 PFT 瓦片。

    参数:
        file_paths: HDF 文件路径列表
        hlines: 横向瓦片编号
        vlines: 纵向瓦片编号

    返回:
        dict: {"mosaic_data", "lat_range", "lon_range", "tile_boundaries"}
    """
    tile_data_map = {}
    for path in file_paths:
        tile_data = read_single_mcd12q1_pft(path)
        if tile_data:
            tile_data_map[tile_data['tile']] = tile_data

    if not tile_data_map:
        print("No valid tile data found.")
        return None

    # 按 v→h 顺序排列
    ordered_tiles = []
    for v in vlines:
        for h in hlines:
            tile_id = f"h{h:02d}v{v:02d}"
            if tile_id in tile_data_map:
                ordered_tiles.append(tile_data_map[tile_id])
            else:
                print(f"Missing {tile_id}")

    all_lat_min = min(d['lat_range'][0] for d in ordered_tiles)
    all_lat_max = max(d['lat_range'][1] for d in ordered_tiles)
    all_lon_min = min(d['lon_range'][0] for d in ordered_tiles)
    all_lon_max = max(d['lon_range'][1] for d in ordered_tiles)

    total_lat_pixels = int(round((all_lat_max - all_lat_min) * PIXEL_PER_DEGREE))
    total_lon_pixels = int(round((all_lon_max - all_lon_min) * PIXEL_PER_DEGREE))
    mosaic_data = np.full((total_lat_pixels, total_lon_pixels), np.nan, dtype=np.float32)

    tile_boundaries = []
    for tile in ordered_tiles:
        pft = tile['data']
        lat_min, lat_max = tile['lat_range']
        lon_min, lon_max = tile['lon_range']

        start_lat = int(round((all_lat_max - lat_max) * PIXEL_PER_DEGREE))
        end_lat = start_lat + 2400
        start_lon = int(round((lon_min - all_lon_min) * PIXEL_PER_DEGREE))
        end_lon = start_lon + 2400

        mosaic_data[start_lat:end_lat, start_lon:end_lon] = pft
        tile_boundaries.append({
            'tile': tile['tile'],
            'lat': (lat_min, lat_max),
            'lon': (lon_min, lon_max),
        })

    return {
        'mosaic_data': mosaic_data,
        'lat_range': (all_lat_min, all_lat_max),
        'lon_range': (all_lon_min, all_lon_max),
        'tile_boundaries': tile_boundaries,
    }


def mosaic_lai_tiles(
    file_paths: List[str],
    hlines: List[int],
    vlines: List[int],
    save_hdf: bool = False,
    output_hdf: str = "",
) -> Optional[dict]:
    """拼接 MCD15A2H LAI + QC 瓦片。

    参数:
        file_paths: HDF 文件路径列表
        hlines, vlines: 瓦片范围
        save_hdf: 是否保存为 HDF
        output_hdf: HDF 输出路径

    返回:
        dict: {"lai_data", "qc_data", ...}
    """
    tile_data_map = {}
    valid_count = 0

    for path in file_paths:
        tile_data = read_single_mcd15a2h(path)
        if tile_data and tile_data.get('original_shape') == (2400, 2400):
            tile_data_map[tile_data['tile']] = tile_data
            valid_count += 1

    if not tile_data_map:
        print("No valid LAI tile data found.")
        return None

    ordered_tiles = []
    for v in vlines:
        for h in hlines:
            tile_id = f"h{h:02d}v{v:02d}"
            if tile_id in tile_data_map:
                ordered_tiles.append(tile_data_map[tile_id])

    if not ordered_tiles:
        return None

    all_lat_min = min(d['lat_range'][0] for d in ordered_tiles)
    all_lat_max = max(d['lat_range'][1] for d in ordered_tiles)
    all_lon_min = min(d['lon_range'][0] for d in ordered_tiles)
    all_lon_max = max(d['lon_range'][1] for d in ordered_tiles)

    lat_span = all_lat_max - all_lat_min
    lon_span = all_lon_max - all_lon_min
    total_lat_pixels = int(round(lat_span * PIXEL_PER_DEGREE))
    total_lon_pixels = int(round(lon_span * PIXEL_PER_DEGREE))

    mosaic_lai = np.full((total_lat_pixels, total_lon_pixels), np.nan, dtype=np.float32)
    mosaic_qc = np.full((total_lat_pixels, total_lon_pixels), 255, dtype=np.uint8)

    for tile in ordered_tiles:
        lai = tile['lai_data']
        qc = tile['qc_data']
        lat_min_t, lat_max_t = tile['lat_range']
        lon_min_t, lon_max_t = tile['lon_range']

        start_lat = int(round((all_lat_max - lat_max_t) * PIXEL_PER_DEGREE))
        end_lat = start_lat + 2400
        start_lon = int(round((lon_min_t - all_lon_min) * PIXEL_PER_DEGREE))
        end_lon = start_lon + 2400

        s_lat = max(0, start_lat)
        e_lat = min(total_lat_pixels, end_lat)
        s_lon = max(0, start_lon)
        e_lon = min(total_lon_pixels, end_lon)

        d_slat = s_lat - start_lat
        d_elat = d_slat + (e_lat - s_lat)
        d_slon = s_lon - start_lon
        d_elon = d_slon + (e_lon - s_lon)

        mosaic_lai[s_lat:e_lat, s_lon:e_lon] = lai[d_slat:d_elat, d_slon:d_elon]
        mosaic_qc[s_lat:e_lat, s_lon:e_lon] = qc[d_slat:d_elat, d_slon:d_elon]

    result = {
        'lai_data': mosaic_lai,
        'qc_data': mosaic_qc,
        'lat_range': (all_lat_min, all_lat_max),
        'lon_range': (all_lon_min, all_lon_max),
    }

    if save_hdf and output_hdf:
        _save_mosaic_hdf(result, ordered_tiles[0], output_hdf)

    return result


def _save_mosaic_hdf(mosaic: dict, reference: dict, output_path: str):
    """将拼接 LAI 数据保存为 HDF 格式。"""
    from pyhdf.SD import SD, SDC
    try:
        hdf = SD(output_path, SDC.WRITE | SDC.CREATE | SDC.TRUNC)

        for attr_name, attr_value in reference['global_attrs'].items():
            if attr_name not in ['StructMetadata.0', 'CoreMetadata.0']:
                try:
                    hdf.setattr(attr_name, attr_value)
                except Exception:
                    pass

        lai_shape = mosaic['lai_data'].shape
        lai_sds = hdf.create('Lai_500m', SDC.FLOAT32, lai_shape)
        lai_data = mosaic['lai_data'].copy()
        fill_value = reference['lai_attrs'].get('_FillValue', 255)
        scale_factor = reference['lai_attrs'].get('scale_factor', 0.1)
        lai_data[np.isnan(lai_data)] = fill_value * scale_factor
        lai_data = lai_data / scale_factor
        lai_sds[:] = lai_data.astype(np.int16)
        lai_sds.endaccess()

        qc_sds = hdf.create('FparLai_QC', SDC.UINT8, lai_shape)
        qc_sds[:] = mosaic['qc_data']
        qc_sds.endaccess()
        hdf.end()
        print(f"Mosaic HDF saved: {output_path}")
    except Exception as e:
        print(f"Save HDF failed: {e}")


def plot_pft_mosaic(
    mosaic: dict,
    date_str: str,
    output_path: str,
    china_bbox: tuple = (73, 135, 18, 54),
):
    """绘制 PFT 拼接图，标注瓦片边界和中国范围。

    参数:
        mosaic: mosaic_pft_tiles() 返回的 dict
        date_str: 日期标签 (e.g. "A2000001")
        output_path: PNG 输出路径
        china_bbox: (lon_min, lon_max, lat_min, lat_max) 中国范围
    """
    data = mosaic['mosaic_data']
    lat_min, lat_max = mosaic['lat_range']
    lon_min, lon_max = mosaic['lon_range']
    tile_boundaries = mosaic['tile_boundaries']

    cmap = ListedColormap(PFT_COLORS)

    plt.figure(figsize=(24, 18))
    plt.imshow(data, cmap=cmap, extent=[lon_min, lon_max, lat_min, lat_max],
               origin='upper', interpolation='none')
    plt.colorbar(label='PFT Category')
    plt.title(f"MODIS MCD12Q1 LC_Type5 Mosaic - {date_str}", fontsize=16)
    plt.xlabel("Longitude", fontsize=14)
    plt.ylabel("Latitude", fontsize=14)

    china_lon = (china_bbox[0], china_bbox[1])
    china_lat = (china_bbox[2], china_bbox[3])

    for boundary in tile_boundaries:
        lat_min_t, lat_max_t = boundary['lat']
        lon_min_t, lon_max_t = boundary['lon']
        tile_name = boundary['tile']

        is_china = (
            lon_min_t >= china_lon[0] and lon_max_t <= china_lon[1] and
            lat_min_t >= china_lat[0] and lat_max_t <= china_lat[1]
        )
        line_color = 'blue' if is_china else 'red'

        plt.axvline(x=lon_min_t, color=line_color, lw=0.8)
        plt.axvline(x=lon_max_t, color=line_color, lw=0.8)
        plt.axhline(y=lat_min_t, color=line_color, lw=0.8)
        plt.axhline(y=lat_max_t, color=line_color, lw=0.8)

        cx = (lon_min_t + lon_max_t) / 2
        cy = (lat_min_t + lat_max_t) / 2
        plt.text(cx, cy, tile_name, color='red', fontsize=6,
                 ha='center', va='center', weight='bold')

    plt.grid(ls='--', alpha=0.3)
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Mosaic plot saved: {output_path}")


def run_mosaic_view_pipeline(
    year: int,
    input_dir: str,
    output_dir: str,
    hlines: List[int],
    vlines: List[int],
    data_type: str = "PFT",
    doy: int = 1,
) -> dict:
    """瓦片拼接可视化管道。

    参数:
        year: 年份
        input_dir: HDF 输入目录
        output_dir: 图片输出目录
        hlines, vlines: 瓦片范围
        data_type: "PFT" (MCD12Q1) 或 "LAI" (MCD15A2H)
        doy: 年积日 (仅 LAI 类型使用)

    返回:
        dict: {"png": output_path, "mosaic": mosaic_dict}
    """
    os.makedirs(output_dir, exist_ok=True)

    if data_type.upper() == "PFT":
        pattern = f"MCD12Q1.A{year}001.*.061.*.hdf"
        date_tag = f"A{year}001"
        out_name = f"MCD12Q1_PFT_Tiles_{date_tag}.png"
    else:
        pattern = f"MCD15A2H.A*{doy:03d}.*.hdf"
        date_tag = f"A{year}{doy:03d}"
        out_name = f"MCD15A2H_LAI_Mosaic_{date_tag}.png"

    file_paths = []
    for v in vlines:
        for h in hlines:
            p = os.path.join(input_dir, pattern.replace("*", f"h{h:02d}v{v:02d}"))
            file_paths.extend(glob.glob(p))

    if not file_paths:
        raise FileNotFoundError(f"No matching files for pattern: {pattern}")

    print(f"Found {len(file_paths)} tiles, starting mosaic...")

    if data_type.upper() == "PFT":
        mosaic = mosaic_pft_tiles(file_paths, hlines, vlines)
        if mosaic:
            output_path = os.path.join(output_dir, out_name)
            plot_pft_mosaic(mosaic, date_tag, output_path)
    else:
        mosaic = mosaic_lai_tiles(file_paths, hlines, vlines)
        output_path = os.path.join(output_dir, out_name)

    return {"png": output_path, "mosaic": mosaic}
