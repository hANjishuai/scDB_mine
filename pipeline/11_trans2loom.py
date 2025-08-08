# 08_trans2loom.py（optional）

import os, sys
import loompy as lp;
import numpy as np;
import scanpy as sc;

def csv_to_loom(input_csv, output_loom):
    # 读取 CSV 文件到 DataFrame
    x = sc.read_csv(input_csv);#R中导出的表达矩阵
    
    # 将 DataFrame 转换为 Loom 文件
    row_attrs = {"Gene": np.array(x.var_names),};
    col_attrs = {"CellID": np.array(x.obs_names)};
    lp.create(output_loom,x.X.transpose(),row_attrs,col_attrs)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python 08_trans2loom.py <input_csv> <output_loom>")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_loom = sys.argv[2]
    csv_to_loom(input_csv, output_loom)