import scanpy as sc 
# from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap as lc
 
adata = sc.read('fused_anndata_annotated.h5ad')


celltype_colors = {
    'Central' : '#1e90ff',
    'Cholangiocyte' : '#E7298A',
    'Endothelial' : '#C72228',
    'Immune' : '#532C8A',
    'Kuppfer': '#d4b7b7',    
    'MidLobular' : '#ffe012',
    'Other_hep': '#7BAFDE',
    'Portal' : '#3cb44b',
    'Stellate' : '#9e6762',
}

celltype_order    = list(celltype_colors.keys())
colors_celltypes  = lc([celltype_colors[c] for c in celltype_colors.keys()])



adata_celltype = dat.from_scanpy(
    adata,
    cell_type_identifier="cell_type",
    sample_identifier="tma_id"
    # ,
    # covariate_df=cov_df
)

viz.stacked_barplot(adata_celltype, feature_name="samples", cmap=colors_celltypes)
# plt.show()
plt.savefig('figures/celltype_freq_cores.pdf', bbox_inches='tight',pad_inches=2)

adata_celltype = dat.from_scanpy(
    adata,
    cell_type_identifier="cell_type",
    sample_identifier="sex"
    # ,
    # covariate_df=cov_df
)

viz.stacked_barplot(adata_celltype, feature_name="samples", cmap=colors_celltypes)
# plt.show()
plt.savefig('figures/celltype_freq_sex.pdf', bbox_inches='tight',pad_inches=2)

adata_celltype = dat.from_scanpy(
    adata,
    cell_type_identifier="is_doublet",
    sample_identifier="sex"
    # ,
    # covariate_df=cov_df
)

viz.stacked_barplot(adata_celltype, feature_name="samples")
# plt.show()
plt.savefig('figures/doublet_freq_core.pdf', bbox_inches='tight',pad_inches=2)

adata_celltype = dat.from_scanpy(
    adata,
    cell_type_identifier="is_doublet",
    sample_identifier="tma_id"
    # ,
    # covariate_df=cov_df
)

viz.stacked_barplot(adata_celltype, feature_name="samples")
# plt.show()
plt.savefig('figures/doublet_freq_core.pdf', bbox_inches='tight',pad_inches=2)