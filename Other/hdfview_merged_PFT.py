from pyhdf.SD import SD, SDC
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
import os
from glob import glob

def read_single_mcd12q1_pft(file_path):
    """读取单个 MCD12Q1 文件的 LC_Type5 (PFT) 数据及地理范围"""
    try:
        hdf = SD(file_path, SDC.READ)
    except Exception as e:
        print(f"❌ 无法打开文件: {file_path} ({e})")
        return None

    try:
        ds = hdf.select('LC_Type5')
        data = ds[:].astype(np.float32)
        fill_value = ds.attributes().get('_FillValue', 255)
        data[data == fill_value] = np.nan
        data[data > 250] = np.nan
    except Exception as e:
        print(f"❌ 读取 LC_Type5 失败: {e}")
        hdf.end()
        return None

    try:
        meta = hdf.attr('StructMetadata.0').get().decode('utf-8', errors='ignore')
        lat_min = float(meta.split('SOUTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lat_max = float(meta.split('NORTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lon_min = float(meta.split('WESTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lon_max = float(meta.split('EASTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
    except:
        fname = os.path.basename(file_path)
        tile = fname.split('.')[2]
        h = int(tile[1:3])
        v = int(tile[4:6])
        lon_min = -180 + h * 10
        lon_max = lon_min + 10
        lat_max = 90 - v * 10
        lat_min = lat_max - 10

    hdf.end()
    tile = os.path.basename(file_path).split('.')[2]
    print(f"✅ 读取完成: {tile}")
    return {
        'tile': tile,
        'data': data,
        'lat_range': (lat_min, lat_max),
        'lon_range': (lon_min, lon_max)
    }

def mosaic_modis_files(file_paths, hlines, vlines):
    """拼接多个MODIS分块文件，严格按照NCL顺序（先v后h）"""
    if not file_paths:
        print("未提供文件路径")
        return None

    # 构建分块映射字典
    tile_data_map = {}
    for path in file_paths:
        tile_data = read_single_mcd12q1_pft(path)
        if tile_data:
            tile_data_map[tile_data['tile']] = tile_data

    if not tile_data_map:
        print("没有有效数据可拼接")
        return None

    # 按NCL顺序（先v后h）排列
    ordered_tiles = []
    for v in vlines:
        for h in hlines:
            tile_id = f"h{h:02d}v{v:02d}"
            if tile_id in tile_data_map:
                ordered_tiles.append(tile_data_map[tile_id])
            else:
                print(f"⚠️ 缺失 {tile_id}")

    # 全局经纬度范围
    all_lat_min = min([d['lat_range'][0] for d in ordered_tiles])
    all_lat_max = max([d['lat_range'][1] for d in ordered_tiles])
    all_lon_min = min([d['lon_range'][0] for d in ordered_tiles])
    all_lon_max = max([d['lon_range'][1] for d in ordered_tiles])

    pixel_per_degree = 2400 / 10.0  # 240 px per degree
    total_lat_pixels = int(round((all_lat_max - all_lat_min) * pixel_per_degree))
    total_lon_pixels = int(round((all_lon_max - all_lon_min) * pixel_per_degree))
    mosaic_data = np.full((total_lat_pixels, total_lon_pixels), np.nan, dtype=np.float32)
    tile_boundaries = []

    print(f"开始拼接，总尺寸 {total_lat_pixels}×{total_lon_pixels}")
    for tile in ordered_tiles:
        pft = tile['data']
        lat_min, lat_max = tile['lat_range']
        lon_min, lon_max = tile['lon_range']

        start_lat = int(round((all_lat_max - lat_max) * pixel_per_degree))
        end_lat = start_lat + 2400
        start_lon = int(round((lon_min - all_lon_min) * pixel_per_degree))
        end_lon = start_lon + 2400

        mosaic_data[start_lat:end_lat, start_lon:end_lon] = pft
        tile_boundaries.append({
            'tile': tile['tile'],
            'lat': (lat_min, lat_max),
            'lon': (lon_min, lon_max)
        })
        print(f"✔ 已拼接 {tile['tile']}")

    return {
        'mosaic_data': mosaic_data,
        'lat_range': (all_lat_min, all_lat_max),
        'lon_range': (all_lon_min, all_lon_max),
        'tile_boundaries': tile_boundaries
    }

def plot_mosaic_with_tiles(mosaic, date_str):
    """绘制拼接结果并标注瓦片边界"""
    if not mosaic:
        return

    data = mosaic['mosaic_data']
    lat_min, lat_max = mosaic['lat_range']
    lon_min, lon_max = mosaic['lon_range']
    tile_boundaries = mosaic['tile_boundaries']

    # 分类颜色（LC_Type5）
    colors = [
        '#0000FF', '#006400', '#228B22', '#7FFF00', '#ADFF2F', '#FFD700',
        '#DAA520', '#FF8C00', '#CD5C5C', '#8B4513', '#D3D3D3', '#FFFFFF'
    ]
    cmap = ListedColormap(colors)

    plt.figure(figsize=(24, 18))
    plt.imshow(data, cmap=cmap, extent=[lon_min, lon_max, lat_min, lat_max],
               origin='upper', interpolation='none')
    plt.colorbar(label='PFT 类别编号')
    plt.title(f"MODIS MCD12Q1 LC_Type5 全球拼接图（带瓦片边界） - {date_str}", fontsize=16)
    plt.xlabel("经度", fontsize=14)
    plt.ylabel("纬度", fontsize=14)

    # 定义中国范围
    china_lon = (73, 135)
    china_lat = (18, 54)

    for boundary in tile_boundaries:
        lat_min_tile, lat_max_tile = boundary['lat']
        lon_min_tile, lon_max_tile = boundary['lon']
        tile_name = boundary['tile']

        # 判断是否属于中国范围
        is_china = (
            lon_min_tile >= china_lon[0] and lon_max_tile <= china_lon[1] and
            lat_min_tile >= china_lat[0] and lat_max_tile <= china_lat[1]
        )
        line_color = 'blue' if is_china else 'red'

        # 绘制边框
        plt.axvline(x=lon_min_tile, color=line_color, lw=0.8)
        plt.axvline(x=lon_max_tile, color=line_color, lw=0.8)
        plt.axhline(y=lat_min_tile, color=line_color, lw=0.8)
        plt.axhline(y=lat_max_tile, color=line_color, lw=0.8)

        # 标注名称
        cx = (lon_min_tile + lon_max_tile) / 2
        cy = (lat_min_tile + lat_max_tile) / 2
        plt.text(cx, cy, tile_name, color='red', fontsize=6, ha='center', va='center', weight='bold')

    plt.grid(ls='--', alpha=0.3)
    plt.tight_layout()

    save_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             f'MCD12Q1_PFT_Tiles_{date_str}.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✅ 拼接图已保存: {save_path}")
    plt.show()

if __name__ == "__main__":
    # 遍历全部瓦片范围
    hlines = list(range(1, 36))
    vlines = list(range(1, 11))
    date_tag = "A2000001"
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Input")

    file_paths = []
    for v in vlines:
        for h in hlines:
            h_str = f"h{h:02d}"
            v_str = f"v{v:02d}"
            pattern = f"MCD12Q1.{date_tag}.{h_str}{v_str}.061.*.hdf"
            file_paths.extend(glob(os.path.join(data_dir, pattern)))

    if not file_paths:
        print("⚠️ 没找到任何匹配文件")
    else:
        print(f"✅ 找到 {len(file_paths)} 个瓦片，开始拼接...")
        mosaic = mosaic_modis_files(file_paths, hlines, vlines)
        if mosaic:
            plot_mosaic_with_tiles(mosaic, date_tag)
