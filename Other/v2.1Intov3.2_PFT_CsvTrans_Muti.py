import pandas as pd
import os

# ===== 0. 定义转换函数（保持不变）=====
def lcB_to_growth_form(df_b):
    tree_cols = [
        "Evergreen Needleleaf trees",
        "Evergreen Broadleaf trees",
        "Deciduous Needleleaf trees",
        "Deciduous Broadleaf trees",
    ]

    crop_cols = [
        "Cereal crops",
        "Broad-leaf crops",
    ]

    df_a = pd.DataFrame()
    df_a["gridID"]   = df_b["CELLID"]
    df_a["TreeFrac"] = df_b[tree_cols].sum(axis=1)
    df_a["CropFrac"] = df_b[crop_cols].sum(axis=1)
    df_a["ShrubFrac"] = df_b["Shrub"]
    df_a["HerbFrac"]  = df_b["Grass"]

    return df_a

# ===== 1. 手动指定输入输出文件映射 =====
# 格式：{"输入文件路径": "输出文件路径"}
file_mapping = {
    "PFT_frac_2000_27km_square.csv": "grid_growth_form_cn27.csv",
    "PFT_frac_2000_9km_square.csv": "grid_growth_form_cn09.csv",
    "PFT_frac_2000_3km_square.csv": "grid_growth_form_cn03.csv",
    # 可根据需要添加更多文件对
    # "landcover_B_XXX.csv": "grid_growth_form.XXX.csv",
}

# ===== 2. 检查文件是否存在 =====
missing_files = [f for f in file_mapping.keys() if not os.path.exists(f)]
if missing_files:
    print("以下输入文件不存在，请检查路径：")
    for f in missing_files:
        print(f"  {f}")
else:
    # ===== 3. 按映射关系处理文件 =====
    for input_path, output_path in file_mapping.items():
        print(f"\n处理文件：{input_path}")
        df_b = pd.read_csv(input_path)
        
        # 转换计算
        df_a = lcB_to_growth_form(df_b)
        
        # 保存输出
        df_a.to_csv(output_path, index=False)
        print(f"已生成输出文件：{output_path}")
    
    print("\n全部处理完成。")