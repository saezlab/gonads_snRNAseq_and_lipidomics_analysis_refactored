import scanpy as sc 
# from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap as lc

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

adata_annot = sc.read('fused_anndata_annotated.h5ad')
adata_annot.var_names = adata.var_names
adata = sc.read_10x_h5(
    filename="data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cell_feature_matrix.h5"
)

adata_annot.var_names = adata.var_names

adata_annot.obs['is_doublet'] = adata_annot.obs['is_doublet'].astype('category')

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


sex_colors = {
    'female': '#E72961',
    'male': '#1e5fff'
}

sex_order    = list(sex_colors.keys())
colors_sex  = lc([sex_colors[c] for c in sex_colors.keys()])





# List of unique cell types
cell_types = adata_annot.obs['cell_type'].unique()

def differential_expression_analysis(adata, cell_type):
    # Subset the data for the specific cell type
    adata_cell = adata[adata.obs['cell_type'] == cell_type]
    # Perform differential expression analysis between sexes
    sc.tl.rank_genes_groups(adata_cell, groupby='sex', method='wilcoxon')
    # Convert results to a DataFrame
    de_results = sc.get.rank_genes_groups_df(adata_cell, group='male')
    # Add -log10 adjusted p-values
    de_results['-log10(pvals_adj)'] = -np.log10(de_results['pvals_adj'])
    return de_results

def create_volcano_plot(de_results, cell_type, ax):
    # Set significance thresholds
    pval_thresh = 1e-4
    logfc_thresh = 1.0
    # Color the points based on significance and fold change
    colors = ['red' if (row['logfoldchanges'] < -logfc_thresh and row['pvals_adj'] < pval_thresh)
              else 'blue' if (row['logfoldchanges'] > logfc_thresh and row['pvals_adj'] < pval_thresh)
              else 'grey' for idx, row in de_results.iterrows()]
    ax.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals_adj']), c=colors, alpha=0.7)
    # Add labels for significant genes
    for idx, row in de_results.iterrows():
        if row['pvals_adj'] < pval_thresh and abs(row['logfoldchanges']) > logfc_thresh:
            ax.text(row['logfoldchanges'], -np.log10(row['pvals_adj']), row['names'], fontsize=8)
    ax.set_title(f'{cell_type}', color=celltype_colors.get(cell_type, 'black'))
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value')
    ax.axhline(-np.log10(pval_thresh), linestyle='--', color='grey', linewidth=0.8)
    ax.axvline(logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
    ax.axvline(-logfc_thresh, linestyle='--', color='grey', linewidth=0.8)

# Create a PDF to save the plots
with PdfPages('figures/DEG_volcano_plots.pdf') as pdf:
    # Determine the number of rows needed for the subplots
    nrows = int(np.ceil(len(cell_types) / 2))
    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(14, nrows * 6))
    axes = axes.flatten()
    for idx, cell_type in enumerate(cell_types):
        de_results = differential_expression_analysis(adata_annot, cell_type)
        create_volcano_plot(de_results, cell_type, axes[idx])
    # Remove any unused subplots
    for i in range(len(cell_types), len(axes)):
        fig.delaxes(axes[i])
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()
    
    
with PdfPages('figures/DEG_volcano_plots.pdf') as pdf:
    # Determine the number of rows needed for the subplots
    nrows = int(np.ceil(len(cell_types) / 2))
    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(14, nrows * 6))
    axes = axes.flatten()
    # Prepare to save DEG results to an Excel file
    writer = pd.ExcelWriter('data_tidy/DEG_results_sex_per_celltype.xlsx', engine='openpyxl')
    for idx, cell_type in enumerate(cell_types):
        de_results = differential_expression_analysis(adata_annot, cell_type)
        # Save DEG results to an Excel sheet
        de_results.to_excel(writer, sheet_name=cell_type, index=False)
        # Create the volcano plot
        create_volcano_plot(de_results, cell_type, axes[idx])
    # Save the Excel file
    writer.save()
    # Remove any unused subplots
    for i in range(len(cell_types), len(axes)):
        fig.delaxes(axes[i])
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()



def differential_expression_analysis(adata, cell_type):
    # Subset the data for the specific cell type
    adata_cell = adata[adata.obs['cell_type'] == cell_type]
    # Perform differential expression analysis between sexes
    sc.tl.rank_genes_groups(adata_cell, groupby='sex', method='wilcoxon')
    # Convert results to a DataFrame
    de_results = sc.get.rank_genes_groups_df(adata_cell, group='male')
    # Add -log10 adjusted p-values
    de_results['-log10(pvals_adj)'] = -np.log10(de_results['pvals_adj'])
    # Find maximum and minimum -log10(pvals_adj) values for dynamic clipping
    max_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).max()
    min_log_pval = de_results['-log10(pvals_adj)'].replace(-np.inf, np.nan).min()
    # Clip values for visualization
    de_results['-log10(pvals_adj)'] = de_results['-log10(pvals_adj)'].replace(np.inf, max_log_pval + 10).replace(-np.inf, min_log_pval - 10).clip(lower=min_log_pval - 10, upper=max_log_pval + 10)
    return de_results

# def create_volcano_plot(de_results, cell_type, ax):
#     # Set significance thresholds
#     pval_thresh = 0.05
#     logfc_thresh = 1.0
#     # Color the points based on significance and fold change
#     colors = ['red' if (row['logfoldchanges'] < -logfc_thresh and row['pvals_adj'] < pval_thresh)
#               else 'blue' if (row['logfoldchanges'] > logfc_thresh and row['pvals_adj'] < pval_thresh)
#               else 'grey' for idx, row in de_results.iterrows()]
#     max_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).max()
#     min_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).min()
#     ax.scatter(de_results['logfoldchanges'], de_results['-log10(pvals_adj)'], c=colors, alpha=0.7)
#     # Add labels for significant genes
#     for idx, row in de_results.iterrows():
#         if row['pvals_adj'] < pval_thresh and abs(row['logfoldchanges']) > logfc_thresh:
#             ax.text(row['logfoldchanges'], row['-log10(pvals_adj)'], row['names'], fontsize=8)
#     ax.set_title(f'Volcano Plot: {cell_type}', color=celltype_colors.get(cell_type, 'black'))
#     ax.set_xlabel('Log2 Fold Change')
#     ax.set_ylabel('-Log10 Adjusted P-value (truncated)')
#     ax.axhline(-np.log10(pval_thresh), linestyle='--', color='grey', linewidth=0.8)
#     ax.axvline(logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
#     ax.axvline(-logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
#     ax.set_ylim(bottom=min_log_pval - 10, top=max_log_pval + 10)  # Set y-axis limits to indicate truncation

def create_volcano_plot(de_results, cell_type, ax):
    # Set significance thresholds
    pval_thresh = 1e-5
    logfc_thresh = 1.0
    max_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).max()
    min_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).min()
    # Color the points based on significance and fold change
    colors = ['red' if (row['logfoldchanges'] < -logfc_thresh and row['pvals_adj'] < pval_thresh)
              else 'blue' if (row['logfoldchanges'] > logfc_thresh and row['pvals_adj'] < pval_thresh)
              else 'grey' for idx, row in de_results.iterrows()]
    ax.scatter(de_results['logfoldchanges'], de_results['-log10(pvals_adj)'], c=colors, alpha=0.7)
    # Add labels for significant genes
    texts = []
    for idx, row in de_results.iterrows():
        if row['pvals_adj'] < pval_thresh and abs(row['logfoldchanges']) > logfc_thresh:
            texts.append(ax.text(row['logfoldchanges'], row['-log10(pvals_adj)'], row['names'], fontsize=8))
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
    ax.set_title(f'Volcano Plot: {cell_type}', color=celltype_colors.get(cell_type, 'black'))
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value (truncated)')
    ax.axhline(-np.log10(pval_thresh), linestyle='--', color='grey', linewidth=0.8)
    ax.axvline(logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
    ax.axvline(-logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
    ax.set_ylim(bottom=min_log_pval - 10, top=max_log_pval + 10)  # Set y-axis limits to indicate truncation



def create_volcano_plot(de_results, cell_type, ax):
    # Set significance thresholds
    pval_thresh = 1e-5
    logfc_thresh = 0.5
    max_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).max()
    min_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).min()
    # Color the points based on significance and fold change
    # Color the points based on significance and fold change
    colors = ['red' if (row['logfoldchanges'] < -logfc_thresh and row['pvals_adj'] < pval_thresh)
              else 'blue' if (row['logfoldchanges'] > logfc_thresh and row['pvals_adj'] < pval_thresh)
              else 'grey' for idx, row in de_results.iterrows()]
    ax.scatter(de_results['logfoldchanges'], de_results['-log10(pvals_adj)'], c=colors, alpha=0.7)
    # Add labels for significant genes
    for idx, row in de_results.iterrows():
        if row['pvals_adj'] < pval_thresh and abs(row['logfoldchanges']) > logfc_thresh:
            ax.text(row['logfoldchanges'], row['-log10(pvals_adj)'], row['names'], fontsize=8, rotation=45)
    ax.set_title(f'{cell_type}', color=celltype_colors.get(cell_type, 'black'))
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value (truncated)')
    ax.axhline(-np.log10(pval_thresh), linestyle='--', color='grey', linewidth=0.8)
    ax.axvline(logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
    ax.axvline(-logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
    ax.set_ylim(bottom=min_log_pval - 10, top=max_log_pval + 10)  # Set y-axis limits to indicate truncation

# Create a PDF to save the plots
with PdfPages('DEG_volcano_sex_per_celltyp.pdf') as pdf:
    # Determine the number of rows needed for the subplots
    nrows = int(np.ceil(len(cell_types) / 2))
    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(14, nrows * 6))
    axes = axes.flatten()
    # Prepare to save DEG results to an Excel file
    writer = pd.ExcelWriter('DEG_results.xlsx', engine='openpyxl')
    for idx, cell_type in enumerate(cell_types):
        de_results = differential_expression_analysis(adata_annot, cell_type)
        # Save DEG results to an Excel sheet
        de_results.to_excel(writer, sheet_name=cell_type, index=False)
        # Create the volcano plot
        create_volcano_plot(de_results, cell_type, axes[idx])
    # Save the Excel file
    writer.save()
    # Remove any unused subplots
    for i in range(len(cell_types), len(axes)):
        fig.delaxes(axes[i])
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()
    
    
def differential_expression_analysis_pairwise(adata, cell_type1, cell_type2):
    # Subset the data for the specific cell types
    adata_pair = adata[adata.obs['cell_type'].isin([cell_type1, cell_type2])]
    # Perform differential expression analysis between cell types
    sc.tl.rank_genes_groups(adata_pair, groupby='cell_type', groups=[cell_type1, cell_type2], method='wilcoxon')
    # Convert results to a DataFrame
    de_results = sc.get.rank_genes_groups_df(adata_pair, group=cell_type1)
    # Add -log10 adjusted p-values
    de_results['-log10(pvals_adj)'] = -np.log10(de_results['pvals_adj'])
    # Find maximum and minimum -log10(pvals_adj) values for dynamic clipping
    max_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).max()
    min_log_pval = de_results['-log10(pvals_adj)'].replace(-np.inf, np.nan).min()
    # Clip values for visualization
    de_results['-log10(pvals_adj)'] = de_results['-log10(pvals_adj)'].replace(np.inf, max_log_pval + 10).replace(-np.inf, min_log_pval - 10).clip(lower=min_log_pval - 10, upper=max_log_pval + 10)
    return de_results
    
    
with PdfPages('Pairwise_DEG_volcano_hepatocytes.pdf') as pdf:
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(21, 7))
    comparisons = [('Central', 'MidLobular'), ('MidLobular', 'Portal'), ('Portal', 'Central')]
    # Prepare to save DEG results to an Excel file
    writer = pd.ExcelWriter('Pairwise_DEG_hepatocytes.xlsx', engine='openpyxl')
    for idx, (cell_type1, cell_type2) in enumerate(comparisons):
        de_results = differential_expression_analysis_pairwise(adata_annot, cell_type1, cell_type2)
        # Save DEG results to an Excel sheet
        de_results.to_excel(writer, sheet_name=f'{cell_type1}_vs_{cell_type2}', index=False)
        # Create the volcano plot
        create_volcano_plot(de_results, f'{cell_type1} vs {cell_type2}', axes[idx])
    # Save the Excel file
    writer.save()
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()



def differential_expression_analysis_doublet(adata, cell_type=None):
    if cell_type:
        # Subset the data for the specific cell type
        adata_subset = adata[adata.obs['cell_type'] == cell_type]
    else:
        # Use all data if cell_type is None
        adata_subset = adata
    # adata_subset.obs..astype('category')
    # Perform differential expression analysis between doublets and non-doublets
    sc.tl.rank_genes_groups(adata_subset, groupby='is_doublet', method='wilcoxon')
    # Convert results to a DataFrame
    de_results = sc.get.rank_genes_groups_df(adata_subset, group='True')
    # Add -log10 adjusted p-values
    de_results['-log10(pvals_adj)'] = -np.log10(de_results['pvals_adj'])
    # Find maximum and minimum -log10(pvals_adj) values for dynamic clipping
    max_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).max()
    min_log_pval = de_results['-log10(pvals_adj)'].replace(-np.inf, np.nan).min()
    # Clip values for visualization
    de_results['-log10(pvals_adj)'] = de_results['-log10(pvals_adj)'].replace(np.inf, max_log_pval + 10).replace(-np.inf, min_log_pval - 10).clip(lower=min_log_pval - 10, upper=max_log_pval + 10)
    return de_results


with PdfPages('Doublets_DEG_volcano_plots.pdf') as pdf:
    # Define the cell types of interest
    cell_types = ['Central', 'MidLobular', 'Portal']
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 14))
    axes = axes.flatten()
    # Prepare to save DEG results to an Excel file
    writer = pd.ExcelWriter('Doublets_DEG_results.xlsx', engine='openpyxl')
    for idx, cell_type in enumerate(cell_types):
        de_results = differential_expression_analysis_doublet(adata_annot, cell_type)
        # Save DEG results to an Excel sheet
        de_results.to_excel(writer, sheet_name=f'{cell_type}_Doublets', index=False)
        # Create the volcano plot
        create_volcano_plot(de_results, f'{cell_type} (Doublets)', axes[idx])
    # Analyze and plot all cells together
    de_results_all = differential_expression_analysis_doublet(adata_annot)
    # Save DEG results for all cells
    de_results_all.to_excel(writer, sheet_name='All_Cells_Doublets', index=False)
    # Create the volcano plot for all cells
    create_volcano_plot(de_results_all, 'All Cells (Doublets)', axes[-1])
    # Save the Excel file
    writer.save()
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()
    
    
    



adata_annot.obs['is_doublet'] = adata_annot.obs['is_doublet'].astype('category')
adata_annot.obs['sex'] = adata_annot.obs['sex'].astype('category')

# Lump together the specific cell types as 'all_hepatocytes'
hepatocyte_types = ['Central', 'MidLobular', 'Other_hep', 'Portal']
adata_annot.obs['lumped_cell_type'] = adata_annot.obs['cell_type'].replace(hepatocyte_types, 'all_hepatocytes')

# Define color palettes
celltype_colors = {
    'all_hepatocytes': '#1e90ff',
    'Cholangiocyte': '#E7298A',
    'Endothelial': '#C72228',
    'Immune': '#532C8A',
    'Kuppfer': '#d4b7b7',    
    'Stellate': '#9e6762',
}


def differential_expression_analysis(adata, cell_type=None, lumped=False, groupby='is_doublet'):
    if lumped:
        # Subset the data for the lumped hepatocyte types
        adata_subset = adata[adata.obs['lumped_cell_type'] == 'all_hepatocytes']
    elif cell_type:
        # Subset the data for the specific cell type
        adata_subset = adata[adata.obs['cell_type'] == cell_type]
    else:
        # Use all data if cell_type is None
        adata_subset = adata
    # Perform differential expression analysis between the specified groups
    sc.tl.rank_genes_groups(adata_subset, groupby=groupby, groups=adata_subset.obs[groupby].cat.categories.tolist(), method='wilcoxon')
    # Convert results to a DataFrame
    de_results = sc.get.rank_genes_groups_df(adata_subset, group=adata_subset.obs[groupby].cat.categories[0])
    # Add -log10 adjusted p-values
    de_results['-log10(pvals_adj)'] = -np.log10(de_results['pvals_adj'])
    # Find maximum and minimum -log10(pvals_adj) values for dynamic clipping
    max_log_pval = de_results['-log10(pvals_adj)'].replace(np.inf, np.nan).max()
    min_log_pval = de_results['-log10(pvals_adj)'].replace(-np.inf, np.nan).min()
    # Clip values for visualization
    de_results['-log10(pvals_adj)'] = de_results['-log10(pvals_adj)'].replace(np.inf, max_log_pval + 10).replace(-np.inf, min_log_pval - 10).clip(lower=min_log_pval - 10, upper=max_log_pval + 10)
    return de_results

def create_volcano_plot(de_results, title, ax):
    # Set significance thresholds
    pval_thresh = 0.05
    logfc_thresh = 1.0
    
    # Color the points based on significance and fold change
    colors = ['red' if (row['logfoldchanges'] < -logfc_thresh and row['pvals_adj'] < pval_thresh)
              else 'blue' if (row['logfoldchanges'] > logfc_thresh and row['pvals_adj'] < pval_thresh)
              else 'grey' for idx, row in de_results.iterrows()]
    
    ax.scatter(de_results['logfoldchanges'], de_results['-log10(pvals_adj)'], c=colors, alpha=0.7)
    
    # Add labels for significant genes
    for idx, row in de_results.iterrows():
        if row['pvals_adj'] < pval_thresh and abs(row['logfoldchanges']) > logfc_thresh:
            ax.text(row['logfoldchanges'], row['-log10(pvals_adj)'], row['names'], fontsize=8, rotation=45)
    
    ax.set_title(title, color='black')
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value (truncated)')
    ax.axhline(-np.log10(pval_thresh), linestyle='--', color='grey', linewidth=0.8)
    ax.axvline(logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
    ax.axvline(-logfc_thresh, linestyle='--', color='grey', linewidth=0.8)
    ax.set_ylim(bottom=min_log_pval - 10, top=max_log_pval + 10)  # Set y-axis limits to indicate truncation

# Create a PDF to save the plots
with PdfPages('Doublets_and_Sex_DEG_volcano_plots_all_hepatocytes.pdf') as pdf:
    # Define the cell types of interest
    cell_types = celltype_order 
    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(14, 28))
    axes = axes.flatten()
    # Prepare to save DEG results to an Excel file
    writer = pd.ExcelWriter('Doublets_and_Sex_DEG_results_all_hepatocytes.xlsx', engine='openpyxl')
    for idx, cell_type in enumerate(cell_types):
        # Analyze and plot doublet status for each cell type
        de_results_doublet = differential_expression_analysis(adata_annot, cell_type, groupby='is_doublet')
        de_results_doublet.to_excel(writer, sheet_name=f'{cell_type}_Doublets', index=False)
        # create_volcano_plot(de_results_doublet, f'{cell_type} (Doublets)', axes[idx * 2])
        # Analyze and plot sex differences for each cell type
        de_results_sex = differential_expression_analysis(adata_annot, cell_type, groupby='sex')
        de_results_sex.to_excel(writer, sheet_name=f'{cell_type}_Sex', index=False)
        # create_volcano_plot(de_results_sex, f'{cell_type} (Sex)', axes[idx * 2 + 1])
    # Analyze and plot all hepatocytes together for doublet status
    de_results_all_hepatocytes_doublet = differential_expression_analysis(adata_annot, lumped=True, groupby='is_doublet')
    de_results_all_hepatocytes_doublet.to_excel(writer, sheet_name='All_Hepatocytes_Doublets', index=False)
    # create_volcano_plot(de_results_all_hepatocytes_doublet, 'All Hepatocytes (Doublets)', axes[-2])
    # Analyze and plot all hepatocytes together for sex differences
    de_results_all_hepatocytes_sex = differential_expression_analysis(adata_annot, lumped=True, groupby='sex')
    de_results_all_hepatocytes_sex.to_excel(writer, sheet_name='All_Hepatocytes_Sex', index=False)
    # create_volcano_plot(de_results_all_hepatocytes_sex, 'All Hepatocytes (Sex)', axes[-1])
    # Save the Excel file
    writer.save()
    
    # plt.tight_layout()
    # pdf.savefig(fig)
    # plt.close()