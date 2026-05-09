from pyhdf.SD import SD, SDC
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import os
from glob import glob

def read_single_mcd15a2h(file_path):
    """读取单个MCD15A2H文件的数据、地理范围及原始形状"""
    try:
        hdf = SD(file_path, SDC.READ)
    except:
        print(f"无法打开文件: {file_path}")
        return None
    
    # 读取LAI数据
    try:
        lai = hdf.select('Lai_500m')
        lai_data = lai[:].astype(np.float32)
        scale_factor = lai.attributes().get('scale_factor', 0.1)
        fill_value = lai.attributes().get('_FillValue', 255)
        lai_data = lai_data * scale_factor
        lai_data[lai_data == fill_value * scale_factor] = np.nan  # 替换填充值
        
        # 读取QC数据并筛选有效值
        qc = hdf.select('FparLai_QC')
        qc_data = qc[:].astype(np.uint8)
        valid_mask = (qc_data & 0b11) <= 0b11  # 0-3级为有效数据
        lai_data = np.where(valid_mask, lai_data, np.nan)
    except Exception as e:
        print(f"读取 {file_path} 数据失败: {e}")
        hdf.end()
        return None
    
    # 获取地理范围
    try:
        # 从元数据解析（优先，更准确）
        struct_meta_bytes = hdf.attr('StructMetadata.0').get()
        struct_meta = struct_meta_bytes.decode('utf-8', errors='ignore')
        lat_min = float(struct_meta.split('SOUTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lat_max = float(struct_meta.split('NORTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lon_min = float(struct_meta.split('WESTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        lon_max = float(struct_meta.split('EASTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
    except:
        # 从分块编号计算（备选）
        file_name = os.path.basename(file_path)
        tile = file_name.split('.')[2]  # 提取hXXvYY
        h = int(tile[1:3])
        v = int(tile[4:6])
        lon_min = -180 + h * 10.0
        lon_max = lon_min + 10.0
        lat_max = 90 - v * 10.0
        lat_min = lat_max - 10.0
    
    # 记录原始数据形状（MCD15A2H固定为2400x2400）
    original_shape = lai_data.shape
    hdf.end()
    return {
        'lai_data': lai_data,
        'lat_range': (lat_min, lat_max),
        'lon_range': (lon_min, lon_max),
        'tile': tile,
        'original_shape': original_shape
    }

def mosaic_modis_files(file_paths, hlines, vlines):
    """拼接多个MODIS分块文件，严格按照NCL顺序（先v后h）处理"""
    if not file_paths:
        print("未提供文件路径")
        return None
    
    # 构建分块映射字典
    tile_data_map = {}
    for path in file_paths:
        tile_data = read_single_mcd15a2h(path)
        if tile_data and tile_data['original_shape'] == (2400, 2400):
            tile_data_map[tile_data['tile']] = tile_data
    
    if not tile_data_map:
        print("没有有效数据可拼接")
        return None
    
    # 严格按照NCL顺序处理：先v后h
    ordered_tiles = []
    for v in vlines:
        for h in hlines:
            v_str = f"v{v:02d}" if v < 10 else f"v{v}"
            tile_id = f"h{h:02d}{v_str}" if h < 10 else f"h{h}{v_str}"
            if tile_id in tile_data_map:
                ordered_tiles.append(tile_data_map[tile_id])
            else:
                print(f"警告：分块 {tile_id} 未找到，已跳过")
    
    # 计算全局地理范围
    all_lat_min = min([d['lat_range'][0] for d in ordered_tiles])
    all_lat_max = max([d['lat_range'][1] for d in ordered_tiles])
    all_lon_min = min([d['lon_range'][0] for d in ordered_tiles])
    all_lon_max = max([d['lon_range'][1] for d in ordered_tiles])
    
    # MCD15A2H固定参数：每个分块10°×10°，2400×2400像素
    pixel_per_degree = 2400 / 10.0  # 每度240像素
    
    # 计算拼接后的数据形状
    lat_span = all_lat_max - all_lat_min
    total_lat_pixels = int(round(lat_span * pixel_per_degree))
    lon_span = all_lon_max - all_lon_min
    total_lon_pixels = int(round(lon_span * pixel_per_degree))
    
    # 初始化拼接数组
    mosaic_data = np.full((total_lat_pixels, total_lon_pixels), np.nan, dtype=np.float32)
    
    # 记录分块边界
    tile_boundaries = []
    
    # 严格按照NCL顺序逐个分块填充（先v后h）
    for tile in ordered_tiles:
        lai = tile['lai_data']
        lat_min_tile, lat_max_tile = tile['lat_range']
        lon_min_tile, lon_max_tile = tile['lon_range']
        
        # 计算分块在全局数组中的索引
        start_lat = int(round((all_lat_max - lat_max_tile) * pixel_per_degree))
        end_lat = start_lat + 2400
        start_lon = int(round((lon_min_tile - all_lon_min) * pixel_per_degree))
        end_lon = start_lon + 2400
        
        # 强制裁剪到有效范围
        start_lat_clamp = max(0, start_lat)
        end_lat_clamp = min(total_lat_pixels, end_lat)
        start_lon_clamp = max(0, start_lon)
        end_lon_clamp = min(total_lon_pixels, end_lon)
        
        # 计算数据裁剪范围
        data_start_lat = start_lat_clamp - start_lat
        data_end_lat = data_start_lat + (end_lat_clamp - start_lat_clamp)
        data_start_lon = start_lon_clamp - start_lon
        data_end_lon = data_start_lon + (end_lon_clamp - start_lon_clamp)
        
        # 填充数据（注意：不需要flipud，因为MODIS数据已经是正确的方向）
        mosaic_data[start_lat_clamp:end_lat_clamp, start_lon_clamp:end_lon_clamp] = \
            lai[data_start_lat:data_end_lat, data_start_lon:data_end_lon]
        
        # 记录边界
        tile_boundaries.append({
            'tile': tile['tile'],
            'lat': (lat_min_tile, lat_max_tile),
            'lon': (lon_min_tile, lon_max_tile)
        })
    
    return {
        'mosaic_data': mosaic_data,
        'lat_range': (all_lat_min, all_lat_max),
        'lon_range': (all_lon_min, all_lon_max),
        'tile_boundaries': tile_boundaries
    }

def plot_mosaic(mosaic):
    """绘制拼接后的LAI图像，添加分块边界"""
    if not mosaic:
        return
    
    # 字体设置
    plt.rcParams["font.family"] = ["sans-serif"]
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial Unicode MS", "SimSun"]
    plt.rcParams["axes.unicode_minus"] = False
    
    # 数据和范围
    lai_mosaic = mosaic['mosaic_data']
    lat_min, lat_max = mosaic['lat_range']
    lon_min, lon_max = mosaic['lon_range']
    tile_boundaries = mosaic['tile_boundaries']
    
    # 自定义颜色映射
    colors = ['#f7f7f7', '#d9f0a3', '#addd8e', '#78c679', '#41ab5d', '#238443', '#005a32']
    cmap = LinearSegmentedColormap.from_list('veg_cmap', colors, N=100)
    
    # 绘图
    plt.figure(figsize=(24, 20))
    im = plt.imshow(
        lai_mosaic,
        cmap=cmap,
        vmin=0, vmax=8,
        extent=[lon_min, lon_max, lat_min, lat_max],
        origin='upper'
    )
    
    # 定义中国区域的大致经纬度范围（可根据需要调整）
    china_lon_range = (73, 135)
    china_lat_range = (18, 54)
    
    # 绘制分块边界并标注分块名称
    for boundary in tile_boundaries:
        lat_min_tile, lat_max_tile = boundary['lat']
        lon_min_tile, lon_max_tile = boundary['lon']
        tile_name = boundary['tile']
        
        # 判断是否在中国区域
        is_china_tile = (
            lon_min_tile >= china_lon_range[0] and lon_max_tile <= china_lon_range[1] and
            lat_min_tile >= china_lat_range[0] and lat_max_tile <= china_lat_range[1]
        )
        
        # 设置线条颜色和粗细
        line_color = 'blue' if is_china_tile else 'red'
        line_width = 0.8  # 减小线条粗细
        
        # 绘制边界
        plt.axvline(x=lon_min_tile, color=line_color, linewidth=line_width, linestyle='-')
        plt.axvline(x=lon_max_tile, color=line_color, linewidth=line_width, linestyle='-')
        plt.axhline(y=lat_min_tile, color=line_color, linewidth=line_width, linestyle='-')
        plt.axhline(y=lat_max_tile, color=line_color, linewidth=line_width, linestyle='-')
        
        # 标注分块名称（不加框，红色，字号4）
        center_lon = (lon_min_tile + lon_max_tile) / 2
        center_lat = (lat_min_tile + lat_max_tile) / 2
        plt.text(center_lon, center_lat, tile_name, 
                ha='center', va='center', fontsize=6, color='red',
                weight='bold')  # 加粗使文字更清晰
    
    plt.title('MCD15A2H 叶面积指数(LAI)拼接图（NCL顺序：先v后h）', fontsize=16)
    plt.xlabel('经度', fontsize=14)
    plt.ylabel('纬度', fontsize=14)
    cbar = plt.colorbar(im, shrink=0.8)
    cbar.set_label('叶面积指数 (m²/m²)', fontsize=12)
    plt.grid(linestyle='--', alpha=0.3)
    plt.tight_layout()
    
    # 保存拼接图
    save_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'LAI_mosaic_NCL_order.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"拼接图已保存至: {save_path}")
    plt.show()

if __name__ == "__main__":
    # ------------ 自定义参数 ------------
    # hlines = [23, 24, 25, 26, 27, 28, 29, 30]  # 经度块编号
    # vlines = [3, 4, 5, 6, 7, 8]                # 纬度块编号
    # hlines = [27, 28, 29]  # 经度块编号
    # vlines = [5, 6]                # 纬度块编号
    hlines = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]  # 经度块编号
    vlines = [1,2,3,4,5, 6,7,8,9,10]                # 纬度块编号
    target_doy = 361                           # 目标年积日
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Input")
    # ------------------------------------
    
    # 严格按照NCL顺序生成文件匹配模式（先v后h）
    file_patterns = []
    for v in vlines:          # 先遍历v（纬度块）
        for h in hlines:      # 再遍历h（经度块）
            v_str = f"v{v:02d}" if v < 10 else f"v{v}"
            h_str = f"h{h:02d}" if h < 10 else f"h{h}"
            pattern = f"MCD15A2H.A*{target_doy}.{h_str}{v_str}.*.hdf"
            file_patterns.append(pattern)
    
    # 收集文件路径
    file_paths = []
    for pattern in file_patterns:
        file_paths.extend(glob(os.path.join(data_dir, pattern)))
    
    if file_paths:
        print(f"找到 {len(file_paths)} 个{target_doy}天的分块文件，开始按NCL顺序拼接...")
        print("处理顺序：先v后h")
        mosaic = mosaic_modis_files(file_paths, hlines, vlines)
        if mosaic:
            plot_mosaic(mosaic)
    else:
        print(f"在 {data_dir} 中未找到目标文件")