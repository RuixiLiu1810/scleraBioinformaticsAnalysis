import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# 1. 准备工作
# ==========================================
# 请将此处替换为您实际的 .h5ad 文件路径
adata = sc.read_h5ad(r'E:\20260120\scanpy_out\data\GSE235684\processed.h5ad') 

# 假设您的聚类列名为 'leiden' 或 'louvain'，请根据实际情况修改
cluster_key = 'leiden' 

# ==========================================
# 2. 定义细胞类型字典 (基于刚才的分析)
# ==========================================
# 这是根据您上传的 CSV 文件整理的注释
cell_type_map = {
    '0': 'vSMC/Pericytes',      # 血管平滑肌/周细胞
    '1': 'T Cells',             # T细胞
    '2': 'Macrophages',         # 巨噬细胞
    '3': 'Endothelial',         # 内皮细胞
    '4': 'Fibroblasts',         # 成纤维细胞
    '5': 'CD8+ Effector T',     # CD8+ 效应T细胞
    '6': 'B Cells',             # B细胞
    '7': 'Mast Cells',          # 肥大细胞
    '8': 'NK/Cytotoxic T',      # NK/细胞毒性T细胞
    '9': 'Cycling Cells',       # 增殖细胞
    '10': 'Plasma Cells',       # 浆细胞
    '11': 'Monocytes',          # 单核细胞
    '12': 'Arterial ECs',       # 动脉内皮
    '13': 'Venous/Capillary ECs', # 静脉/毛细血管内皮
    '14': 'Lymphatic ECs',      # 淋巴管内皮
    '15': 'Tregs',              # 调节性T细胞
    '16': 'Epithelial/Tumor',   # 上皮/肿瘤细胞
    '17': 'mReg DCs',           # 成熟树突状细胞
    '18': 'pDCs',               # 浆细胞样树突状细胞
    '19': 'Neutrophils',        # 中性粒细胞
    '20': 'Stromal/Adipogenic'  # 基质/脂肪前体细胞
}

# 检查数据中是否存在所有 cluster，防止报错
# (这一步是为了确保字典里的 key 和数据里的类别能对应上)
if 'adata' in locals():
    # 将 cluster 列转换为字符串以进行映射
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)
    
    # 创建一个新的列 'cell_type'
    adata.obs['cell_type'] = adata.obs[cluster_key].map(cell_type_map).fillna(adata.obs[cluster_key])

# ==========================================
# 3. 绘制聚类图 (UMAP)
# ==========================================
if 'adata' in locals():
    plt.rcParams['figure.figsize'] = (6, 6)
    
    # 画图：左边是原始数字编号，右边是注释后的名称
    sc.pl.umap(
        adata, 
        color=[cluster_key, 'cell_type'], 
        legend_loc='on data',  # 标签直接显示在图上
        legend_fontsize=8, 
        title=['Cluster Numbers', 'Cell Types'],
        frameon=False,
        show=True
    )

# ==========================================
# 4. 绘制 Marker 气泡图 (DotPlot) - 验证注释
# ==========================================
# 定义我们要检查的核心基因列表
marker_genes_dict = {
    'vSMC/Peri': ['RGS5', 'MYH11'],
    'T Cells': ['CD3D', 'IL32'],
    'Myeloid': ['CD68', 'C1QB', 'S100A8', 'TPSAB1'], # 巨噬, 单核, 肥大
    'Endo': ['PECAM1', 'VWF', 'PROX1'], # 血管, 淋巴管
    'Fibro': ['DCN', 'COL1A1'],
    'B/Plasma': ['MS4A1', 'JCHAIN'],
    'NK/Cyto': ['GNLY', 'NKG7'],
    'Epi/Tumor': ['EPCAM', 'KRT19']
}

if 'adata' in locals():
    sc.pl.dotplot(
        adata, 
        marker_genes_dict, 
        groupby='cell_type', 
        standard_scale='var', # 标准化表达量，让对比更明显
        show=True
    )