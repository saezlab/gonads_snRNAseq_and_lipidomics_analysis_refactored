import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
# import squidpy as sq
import os 
import random

random.seed(1996)

from tqdm import tqdm


genes_df = pd.read_excel("data_raw/XeniumTargetList_Mouse-Panel-1.xlsx")
type_markers = genes_df[(genes_df['reason'] != 'gonad-dependent gene') & (genes_df['reason'] != 'XCI')]['gene name'].values


adata = sc.read_10x_h5(
    filename="data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cell_feature_matrix.h5"
)

df = pd.read_csv(
    "data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cells.csv.gz"
)

df.set_index('cell_id')
adata.obs = df
adata.obs = adata.obs.set_index('cell_id')
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()



adata_types = adata[:, type_markers].copy()
adata_types.layers["counts"] = adata_types.X.copy()
sc.pp.filter_genes(adata_types, min_cells=5)
sc.pp.filter_cells(adata_types, min_counts=10)
sc.pp.normalize_total(adata_types, inplace=True)
sc.pp.log1p(adata_types)
sc.pp.pca(adata_types)
sc.pp.neighbors(adata_types)
sc.tl.umap(adata_types)
sc.tl.leiden(adata_types, 1)

sc.pl.umap(
    adata_types,
    color=[
        # "total_counts",
        # "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)

# sc.pl.umap(
#     adata_types,
#     color=[
#         # "total_counts",
#         # "n_genes_by_counts",
#         "leiden",
#     ],
#     wspace=0.4,
# )
# adata_types.obs = adata_types.obs.set_index('cell_id')
adata = adata[adata_types.obs_names, :]
adata.obs = adata_types.obs

sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, save='_markers_type.pdf')
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5, save='_markers_type.pdf')



cell_type_dict = {
    'Central': ['4', '9', '13', '33', '34', '35', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47'],
    'Portal' : ['2', '5', '7', '10', '16', '19', '21', '29', '30', ],
    'MidLobular': ['3'],
    'Immune': ['22', '24'],
    'Endothelial': ['6', '12', '14', '23'],
    'Cholangiocyte': ['18'],
    'Stellate': ['11', '20'],
    # 'Ar+': [''],
    'Other_hep': ['8', '15', '17', '25', '26', '27', '28', '31', '32', '36', '37'],
    'Kuppfer': ['0', '1']
}




celltype_colors = {
        'Central' : '#1e90ff',
        'Immune' : '#532C8A',
        'Endothelial' : '#C72228',
        'Stellate' : '#9e6762',
        'MidLobular' : '#ffe012',
        'Portal' : '#3cb44b',
        'Other_hep': '#7BAFDE',
        'Kuppfer': '#d4b7b7',
        'Cholangiocyte' : '#E7298A'
}


celltype_dict = {y: x for x, y in cell_type_dict.items()}

transformed_dict = {v: [k for k in cell_type_dict if v in cell_type_dict[k]][0] for v in set(sum(cell_type_dict.values(), []))}

adata.obs['cell_type'] = [transformed_dict[x] for x in adata.obs['leiden'].values.astype('str')]


sc.pl.spatial(adata, color='cell_type', palette= celltype_colors, spot_size=15, save='_spatial_celltypes2.pdf')


sc.tl.rank_genes_groups(adata, "cell_type", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, save='_celltype.pdf')
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5, save='_celltype.pdf')


hep_types = ['Central', 'Portal', 'MidLobular', 'Other_hep']


import json
with open('inferred_components.json', 'r') as f:
    processed_fusion_list = [set(component) for component in json.load(f)]


cell_ids = set(adata.obs_names) 
processed_fusion_list = [component.intersection(cell_ids) for component in processed_fusion_list]
# Filter out empty sets (optional)
processed_fusion_list = [component for component in processed_fusion_list if component]

# Create a DataFrame for the fusion groups
fusion_groups = pd.DataFrame({'fusion_id': range(len(processed_fusion_list)), 'cell_id': processed_fusion_list})

# Flatten the DataFrame
fusion_groups = fusion_groups.explode('cell_id').reset_index(drop=True)

# Merge with the original observations to get cell types and ensure correct ordering
fusion_groups = fusion_groups.merge(adata.obs[['cell_type']], left_on='cell_id', right_index=True)

# Group by fusion_id and cell_type to identify valid fusions
valid_fusion_groups = fusion_groups.groupby(['fusion_id', 'cell_type']).filter(lambda x: len(x) == 2 and len(x['cell_type'].unique()) == 1 and x['cell_type'].unique()[0] in hep_types)

# # Calculate mean expressions using transform and groupby
# valid_fusion_groups['sum_expression'] = adata.X[valid_fusion_groups['cell_id'],:].T @ valid_fusion_groups.groupby('fusion_id')['cell_id'].transform(len)


# Compute mean expressions for valid groups
fused_data = []
fused_obs = []

grouped = valid_fusion_groups.groupby('fusion_id')

for group_id, group in tqdm(grouped):
    cell_ids = group['cell_id'].values
    mean_expression = adata[cell_ids].X.sum(axis=0)
    fused_data.append(mean_expression)
    obs_to = adata.obs.loc[cell_ids[0]]
    # obs_to['is_doublet'] = True
    fused_obs.append(obs_to)
    

# Convert to array and DataFrame
fused_data = np.array(fused_data)
fused_data = np.squeeze(fused_data, axis=1)
fused_data = fused_data[::2]
fused_obs = pd.DataFrame(fused_obs)
fused_obs[['is_doublet']] = True
fused_obs = fused_obs.drop_duplicates()
# Identify non-fused cells
fused_cells = valid_fusion_groups['cell_id'].unique()

# fused_data = fused_adata
mask = ~adata.obs.index.isin(fused_cells)

# Append non-fused cells
non_fused_data = adata.X[mask]
non_fused_obs = adata.obs[mask]
non_fused_obs[['is_doublet']] = False
# fused_data = np.squeeze(fused_data, axis=1)

# Combine fused and non-fused data
final_data = np.vstack((fused_data, non_fused_data.toarray()))
final_obs = pd.concat([fused_obs, non_fused_obs])


import anndata as ad
# Create new AnnData object
fused_adata = ad.AnnData(X=final_data, obs=final_obs)
fused_adata.obs = final_obs

sc.pp.normalize_total(fused_adata, inplace=True)
sc.pp.log1p(fused_adata)
sc.pp.pca(fused_adata)
sc.pp.neighbors(fused_adata)
sc.tl.umap(fused_adata)

sc.pl.umap(fused_adata, color='cell_type', palette=celltype_colors, save='_fused_celltypes.pdf')
sc.pl.umap(fused_adata, color='is_doublet', save='_fused_doublet.pdf')
sc.pl.spatial(fused_adata, color='is_doublet', save='_fused_doublet.pdf',spot_size=25)


fused_adata.obsm["spatial"] = fused_adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()



x_lines = [0, 5000, 11000]
y_lines = [0, 3500, 6500, 11000, 16000, 21000]

k = 1
for x in range(1, len(x_lines)):
    for y in range(1, len(y_lines)):
        fused_adata.obs.loc[
            (fused_adata.obs['x_centroid'] > x_lines[x - 1]) &
            (fused_adata.obs['x_centroid'] < x_lines[x]) &
            (fused_adata.obs['y_centroid'] > y_lines[y - 1]) &
            (fused_adata.obs['y_centroid'] < y_lines[y]), "tma_id"
        ] = k
        if x == 1:
            gender = 'female'
        else:
            gender = 'male'
        fused_adata.obs.loc[
            (fused_adata.obs['x_centroid'] > x_lines[x - 1]) &
            (fused_adata.obs['x_centroid'] < x_lines[x]) &
            (fused_adata.obs['y_centroid'] > y_lines[y - 1]) &
            (fused_adata.obs['y_centroid'] < y_lines[y]), "sex"
        ] = gender
        print(x_lines[x-1], x_lines[x], y_lines[y-1], y_lines[y])
        k += 1



fused_adata.obs.tma_id = fused_adata.obs.tma_id.astype("str")




fused_adata.obs.groupby(by=['cell_type', 'is_doublet']).size()
# Save the new AnnData object
fused_adata.write('fused_anndata_annotated.h5ad')


fused_adata.obs[['group']] = fused_adata.obs[['is_doublet']]
fused_adata.obs['cell_id'] = fused_adata.obs_names 

fused_adata.obs[['cell_id', 'group']].to_csv('data_tidy/fused_doublet.csv')


fused_adata.obs[['group']] = fused_adata.obs[['cell_type']]
fused_adata.obs['cell_id'] = fused_adata.obs_names 

fused_adata.obs[['cell_id', 'group']].to_csv('data_tidy/fused_celltype.csv')