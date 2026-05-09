import pandas as pd
import glob

# 匹配所有目标CSV文件
csv_files = glob.glob("grid_growth_form.GD_cn*.csv")

for file in csv_files:
    print(f"[处理中] {file}")
    # 读取CSV
    df = pd.read_csv(file)
    
    # 转换列：百分数 → 小数（除以100）
    frac_cols = ["TreeFrac", "CropFrac", "ShrubFrac", "HerbFrac"]
    df[frac_cols] = df[frac_cols] / 100.0  # 核心转换：例如100 → 1.0，50 → 0.5
    
    # 保存覆盖原文件（或改为新文件名，如加"_converted"后缀）
    df.to_csv(file, index=False, encoding="utf-8-sig")
    print(f"[完成] {file} 已转换小数格式")

print("\n所有文件处理完成！")