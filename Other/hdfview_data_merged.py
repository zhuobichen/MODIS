from pyhdf.SD import SD, SDC
import numpy as np
import os
from glob import glob
import shutil

def read_single_mcd15a2h(file_path):
    """读取单个MCD15A2H文件的数据、地理范围及所有属性"""
    try:
        hdf = SD(file_path, SDC.READ)
    except:
        print(f"无法打开文件: {file_path}")
        return None
    
    try:
        # 读取LAI数据
        lai = hdf.select('Lai_500m')
        lai_data = lai[:].astype(np.float32)
        
        # 获取LAI属性
        lai_attrs = lai.attributes()
        scale_factor = lai_attrs.get('scale_factor', 0.1)
        fill_value = lai_attrs.get('_FillValue', 255)
        lai_data = lai_data * scale_factor
        lai_data[lai_data == fill_value * scale_factor] = np.nan
        
        # 读取QC数据并筛选有效值
        qc = hdf.select('FparLai_QC')
        qc_data = qc[:].astype(np.uint8)
        qc_attrs = qc.attributes()
        valid_mask = (qc_data & 0b11) <= 0b11
        lai_data = np.where(valid_mask, lai_data, np.nan)
        
        # 获取全局属性
        global_attrs = hdf.attributes()
        
        # 初始化结构元数据
        struct_meta = ""
        
        # 获取地理范围
        try:
            struct_meta_bytes = hdf.attr('StructMetadata.0').get()
            struct_meta = struct_meta_bytes.decode('utf-8', errors='ignore')
            lat_min = float(struct_meta.split('SOUTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
            lat_max = float(struct_meta.split('NORTHBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
            lon_min = float(struct_meta.split('WESTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
            lon_max = float(struct_meta.split('EASTBOUNDINGCOORDINATE=')[1].split('\n')[0].strip())
        except Exception as e:
            print(f"从元数据解析地理范围失败 {file_path}: {e}，使用分块编号计算")
            # 从分块编号计算（备选）
            file_name = os.path.basename(file_path)
            tile = file_name.split('.')[2]
            h = int(tile[1:3])
            v = int(tile[4:6])
            lon_min = -180 + h * 10.0
            lon_max = lon_min + 10.0
            lat_max = 90 - v * 10.0
            lat_min = lat_max - 10.0
        
        hdf.end()
        
        return {
            'lai_data': lai_data,
            'lai_attrs': lai_attrs,
            'qc_data': qc_data,
            'qc_attrs': qc_attrs,
            'global_attrs': global_attrs,
            'struct_metadata': struct_meta,
            'lat_range': (lat_min, lat_max),
            'lon_range': (lon_min, lon_max),
            'tile': tile,
            'original_shape': lai_data.shape
        }
        
    except Exception as e:
        print(f"读取 {file_path} 数据失败: {e}")
        try:
            hdf.end()
        except:
            pass
        return None

def mosaic_modis_files(file_paths, hlines, vlines):
    """拼接多个MODIS分块文件，严格按照NCL顺序（先v后h）处理"""
    if not file_paths:
        print("未提供文件路径")
        return None
    
    # 构建分块映射字典
    tile_data_map = {}
    valid_files = 0
    
    for path in file_paths:
        print(f"正在读取: {os.path.basename(path)}")
        tile_data = read_single_mcd15a2h(path)
        if tile_data and tile_data['original_shape'] == (2400, 2400):
            tile_data_map[tile_data['tile']] = tile_data
            valid_files += 1
            print(f"  成功读取分块 {tile_data['tile']}")
        else:
            print(f"  读取失败或数据形状不正确")
    
    print(f"成功读取 {valid_files}/{len(file_paths)} 个文件")
    
    if not tile_data_map:
        print("没有有效数据可拼接")
        return None
    
    # 严格按照NCL顺序处理：先v后h
    ordered_tiles = []
    missing_tiles = []
    
    for v in vlines:
        for h in hlines:
            v_str = f"v{v:02d}" if v < 10 else f"v{v}"
            h_str = f"h{h:02d}" if h < 10 else f"h{h}"
            tile_id = f"{h_str}{v_str}"
            if tile_id in tile_data_map:
                ordered_tiles.append(tile_data_map[tile_id])
                print(f"添加分块: {tile_id}")
            else:
                missing_tiles.append(tile_id)
                print(f"警告：分块 {tile_id} 未找到")
    
    if not ordered_tiles:
        print("没有找到任何有效的分块数据")
        return None
    
    print(f"开始拼接 {len(ordered_tiles)} 个分块...")
    
    # 计算全局地理范围
    all_lat_min = min([d['lat_range'][0] for d in ordered_tiles])
    all_lat_max = max([d['lat_range'][1] for d in ordered_tiles])
    all_lon_min = min([d['lon_range'][0] for d in ordered_tiles])
    all_lon_max = max([d['lon_range'][1] for d in ordered_tiles])
    
    print(f"全局地理范围: 经度 {all_lon_min:.2f}~{all_lon_max:.2f}, 纬度 {all_lat_min:.2f}~{all_lat_max:.2f}")
    
    # MCD15A2H固定参数
    pixel_per_degree = 2400 / 10.0
    
    # 计算拼接后的数据形状
    lat_span = all_lat_max - all_lat_min
    total_lat_pixels = int(round(lat_span * pixel_per_degree))
    lon_span = all_lon_max - all_lon_min
    total_lon_pixels = int(round(lon_span * pixel_per_degree))
    
    print(f"拼接后数据形状: {total_lat_pixels} x {total_lon_pixels}")
    
    # 初始化拼接数组
    mosaic_lai = np.full((total_lat_pixels, total_lon_pixels), np.nan, dtype=np.float32)
    mosaic_qc = np.full((total_lat_pixels, total_lon_pixels), 255, dtype=np.uint8)  # 默认填充值
    
    # 严格按照NCL顺序逐个分块填充
    for i, tile in enumerate(ordered_tiles):
        print(f"处理分块 {i+1}/{len(ordered_tiles)}: {tile['tile']}")
        
        lai = tile['lai_data']
        qc = tile['qc_data']
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
        
        # 填充数据
        mosaic_lai[start_lat_clamp:end_lat_clamp, start_lon_clamp:end_lon_clamp] = \
            lai[data_start_lat:data_end_lat, data_start_lon:data_end_lon]
        mosaic_qc[start_lat_clamp:end_lat_clamp, start_lon_clamp:end_lon_clamp] = \
            qc[data_start_lat:data_end_lat, data_start_lon:data_end_lon]
        
        print(f"  分块 {tile['tile']} 已放置到位置 [{start_lat_clamp}:{end_lat_clamp}, {start_lon_clamp}:{end_lon_clamp}]")
    
    # 选择一个参考文件的属性（使用第一个有效文件）
    reference_data = ordered_tiles[0]
    
    return {
        'lai_data': mosaic_lai,
        'qc_data': mosaic_qc,
        'lai_attrs': reference_data['lai_attrs'],
        'qc_attrs': reference_data['qc_attrs'],
        'global_attrs': reference_data['global_attrs'],
        'struct_metadata': reference_data['struct_metadata'],
        'lat_range': (all_lat_min, all_lat_max),
        'lon_range': (all_lon_min, all_lon_max),
        'tile_boundaries': [(tile['tile'], tile['lat_range'], tile['lon_range']) for tile in ordered_tiles],
        'original_tiles': [tile['tile'] for tile in ordered_tiles]
    }

def create_new_struct_metadata(mosaic, original_meta):
    """创建新的结构元数据"""
    if not original_meta:
        # 如果没有原始元数据，创建一个基本的
        return f"""GROUP=SwathStructure
END_GROUP=SwathStructure
GROUP=GridStructure
  GROUP=GRID_1
    GridName="MOD_Grid_MOD15A2H"
    XDim={mosaic['lai_data'].shape[1]}
    YDim={mosaic['lai_data'].shape[0]}
    UpperLeftPointMtrs=(-20015109.354000,10007554.677000)
    LowerRightMtrs=(-18903158.572594,8895604.049787)
    Projection=GCTP_SNSOID
    ProjectionCode=0
    GridOrigin=UL
    GROUP=Point
      PointName="MOD_Grid_MOD15A2H"
      WESTBOUNDINGCOORDINATE={mosaic['lon_range'][0]:.6f}
      EASTBOUNDINGCOORDINATE={mosaic['lon_range'][1]:.6f}
      NORTHBOUNDINGCOORDINATE={mosaic['lat_range'][1]:.6f}
      SOUTHBOUNDINGCOORDINATE={mosaic['lat_range'][0]:.6f}
    END_GROUP=Point
  END_GROUP=GRID_1
END_GROUP=GridStructure"""
    
    # 否则更新现有的元数据
    meta = original_meta
    meta = meta.replace(
        f"WESTBOUNDINGCOORDINATE={mosaic['tile_boundaries'][0][2][0]:.6f}",
        f"WESTBOUNDINGCOORDINATE={mosaic['lon_range'][0]:.6f}"
    )
    meta = meta.replace(
        f"EASTBOUNDINGCOORDINATE={mosaic['tile_boundaries'][0][2][1]:.6f}",
        f"EASTBOUNDINGCOORDINATE={mosaic['lon_range'][1]:.6f}"
    )
    meta = meta.replace(
        f"NORTHBOUNDINGCOORDINATE={mosaic['tile_boundaries'][0][1][1]:.6f}",
        f"NORTHBOUNDINGCOORDINATE={mosaic['lat_range'][1]:.6f}"
    )
    meta = meta.replace(
        f"SOUTHBOUNDINGCOORDINATE={mosaic['tile_boundaries'][0][1][0]:.6f}",
        f"SOUTHBOUNDINGCOORDINATE={mosaic['lat_range'][0]:.6f}"
    )
    return meta

def save_as_hdf(mosaic, output_path):
    """将拼接数据保存为HDF格式"""
    try:
        # 创建新的HDF文件
        hdf = SD(output_path, SDC.WRITE | SDC.CREATE | SDC.TRUNC)
        
        # 设置全局属性
        for attr_name, attr_value in mosaic['global_attrs'].items():
            if attr_name not in ['StructMetadata.0', 'CoreMetadata.0']:
                try:
                    hdf.setattr(attr_name, attr_value)
                except:
                    print(f"跳过属性 {attr_name}")
        
        # 创建新的结构元数据
        new_struct_meta = create_new_struct_metadata(mosaic, mosaic['struct_metadata'])
        hdf.setattr('StructMetadata.0', new_struct_meta)
        
        # 创建LAI数据集
        lai_shape = mosaic['lai_data'].shape
        lai_sds = hdf.create('Lai_500m', SDC.FLOAT32, lai_shape)
        
        # 设置LAI属性
        for attr_name, attr_value in mosaic['lai_attrs'].items():
            try:
                lai_sds.setattr(attr_name, attr_value)
            except:
                print(f"跳过LAI属性 {attr_name}")
        
        # 写入LAI数据（注意：需要将NaN转换回填充值）
        lai_data_to_write = mosaic['lai_data'].copy()
        fill_value = mosaic['lai_attrs'].get('_FillValue', 255)
        scale_factor = mosaic['lai_attrs'].get('scale_factor', 0.1)
        lai_data_to_write[np.isnan(lai_data_to_write)] = fill_value * scale_factor
        lai_data_to_write = lai_data_to_write / scale_factor  # 反向缩放
        lai_sds[:] = lai_data_to_write.astype(np.int16)  # 转回原始数据类型
        
        # 创建QC数据集
        qc_sds = hdf.create('FparLai_QC', SDC.UINT8, lai_shape)
        
        # 设置QC属性
        for attr_name, attr_value in mosaic['qc_attrs'].items():
            try:
                qc_sds.setattr(attr_name, attr_value)
            except:
                print(f"跳过QC属性 {attr_name}")
        
        # 写入QC数据
        qc_sds[:] = mosaic['qc_data']
        
        # 结束写入
        lai_sds.endaccess()
        qc_sds.endaccess()
        hdf.end()
        
        print(f"拼接数据已保存为HDF格式: {output_path}")
        return True
        
    except Exception as e:
        print(f"保存HDF文件失败: {e}")
        try:
            hdf.end()
        except:
            pass
        return False

def get_mosaic_info(mosaic):
    """获取拼接数据的详细信息"""
    if not mosaic:
        return
    
    print("\n=== 拼接数据详细信息 ===")
    print(f"数据形状: {mosaic['lai_data'].shape}")
    print(f"纬度范围: {mosaic['lat_range']}")
    print(f"经度范围: {mosaic['lon_range']}")
    print(f"包含分块数量: {len(mosaic['original_tiles'])}")
    print(f"分块列表: {mosaic['original_tiles']}")
    
    # 统计有效数据
    valid_data = np.sum(~np.isnan(mosaic['lai_data']))
    total_data = mosaic['lai_data'].size
    valid_ratio = valid_data / total_data * 100
    print(f"有效数据比例: {valid_ratio:.2f}% ({valid_data}/{total_data})")
    
    # 数据统计
    valid_values = mosaic['lai_data'][~np.isnan(mosaic['lai_data'])]
    if len(valid_values) > 0:
        print(f"LAI值范围: [{valid_values.min():.3f}, {valid_values.max():.3f}]")
        print(f"LAI平均值: {valid_values.mean():.3f}")
        print(f"LAI标准差: {valid_values.std():.3f}")

if __name__ == "__main__":
    # ------------ 自定义参数 ------------
    hlines = [23, 24, 25, 26, 27, 28, 29, 30]  # 经度块编号
    vlines = [3, 4, 5, 6, 7, 8]                # 纬度块编号
    target_doy = 361                           # 目标年积日
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Input")
    output_dir = os.path.dirname(os.path.abspath(__file__))
    # ------------------------------------
    
    print(f"开始处理 {target_doy} 天的数据...")
    print(f"经度分块: {hlines}")
    print(f"纬度分块: {vlines}")
    
    # 严格按照NCL顺序生成文件匹配模式（先v后h）
    file_patterns = []
    for v in vlines:
        for h in hlines:
            v_str = f"v{v:02d}" if v < 10 else f"v{v}"
            h_str = f"h{h:02d}" if h < 10 else f"h{h}"
            pattern = f"MCD15A2H.A*{target_doy}.{h_str}{v_str}.*.hdf"
            file_patterns.append(pattern)
    
    # 收集文件路径
    file_paths = []
    for pattern in file_patterns:
        found_files = glob(os.path.join(data_dir, pattern))
        file_paths.extend(found_files)
        if found_files:
            print(f"找到匹配 {pattern} 的文件: {len(found_files)} 个")
    
    if file_paths:
        print(f"\n总共找到 {len(file_paths)} 个{target_doy}天的分块文件")
        print("开始按NCL顺序拼接...")
        print("处理顺序：先v后h")
        mosaic = mosaic_modis_files(file_paths, hlines, vlines)
        if mosaic:
            # 显示详细信息
            get_mosaic_info(mosaic)
            
            # 保存为HDF格式
            output_filename = f"MCD15A2H.A{target_doy}.mosaic.hdf"
            output_path = os.path.join(output_dir, output_filename)
            success = save_as_hdf(mosaic, output_path)
            
            if success:
                print(f"\n成功生成拼接HDF文件!")
                print(f"输出文件: {output_path}")
                print(f"文件大小: {os.path.getsize(output_path) / (1024*1024):.2f} MB")
        else:
            print("拼接失败")
    else:
        print(f"在 {data_dir} 中未找到目标文件")