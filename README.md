# Gonadal Hormone Coordination of Sex-Specific Adult Liver Physiology

This repository contains the refactored analysis code for a multi-omic Four
Core Genotypes (FCG) mouse liver study. The goal is to separate the effects of
gonadal hormones from sex chromosome complement and to connect sex-biased liver
gene expression with lipid metabolism.

## Project Overview

- **snRNA-seq / pseudobulk RNA:** 466,986 high-quality nuclei from 48 liver
  samples (43 FCG + 5 wild-type males), summarized into pseudobulk profiles by
  sample and cell type.
- **Untargeted lipidomics:** lipid species measured across FCG liver samples.
  The current MOFA path uses the RNA-shared lipidomic sample set.
- **Multi-omic integration:** MOFA2-based integration of 5 RNA cell-type views
  and 1 lipidomic view.
- **Validation:** GDX RNA/lipidomic validation, Xenium spatial validation,
  proteomics validation, siRNA knockdown lipidomics, and collaborator Xenium
  spatial-analysis scripts.
- **Human liver comparison:** optional Python workflow comparing the mouse
  hepatocyte gonadal t-value signature with human bulk liver sex differences
  from GEO accession GSE179250.

## Key Findings

- Gonadal hormones, rather than sex chromosome complement, are the dominant
  driver of sex-biased adult liver physiology.
- Hepatocytes show the strongest gonad-dependent RNA effects, with roughly
  4,000 gonad-dependent genes in the main analysis.
- Lipid abundance is strongly gonad-dependent, while sex chromosome effects on
  lipids are minimal to absent in the current analysis narrative.
- Spatial zonation modulates gonad-dependent hepatocyte expression.
- Key validated regulatory genes include `Cyp2b9`, `Cyp2b13`, `Cyp2c38`,
  `Cyp2d9`, `Sult3a1`, and `Mtmr11`.

For exact counts, prefer the generated result files over prose summaries. In
particular, workflow 03 currently contains two different printed headline lipid
counts in comments/output text, so result tables are the quantitative source.

## Four Core Genotypes Model

| Group | Sex Chromosomes | Gonads | Description |
|-------|-----------------|--------|-------------|
| XXO | XX | Ovaries | Chr3-Sry- |
| XXT | XX | Testes | Chr3-Sry+ |
| XYO | XY | Ovaries | Chr3-Sry- |
| XYT | XY | Testes | Chr3-Sry+ |

Main comparisons:

- **Gonadal effect:** Testes vs Ovaries, across chromosome backgrounds.
- **Chromosomal effect:** XY vs XX, across gonadal backgrounds.
- **Group comparisons:** XXT vs XXO, XYT vs XYO, XYO vs XXO, XYT vs XXT.

## Agent-Oriented Data Exploration

This repository is documented for direct agent use without an MCP server.

Start with:

- `AGENTS.md` - repository rules and biological context for agents.
- `docs/AGENT_GUIDE.md` - static query map for extracting results.
- `docs/workflow_plans/_meta_cross_workflow_overview.md` - cross-workflow data
  flow and shared fragility points.
- `docs/workflow_plans/NN_*.md` - per-workflow structure, dependencies,
  outputs, and known vulnerabilities.

When modifying any `scripts/workflows/NN_*` workflow file, read the matching
plan first and update it after the code change.

## Repository Structure

```text
.
|-- data/
|   |-- RNA/                         # Pseudobulk RNA CSVs and Seurat RDS inputs
|   |-- metabolomics/                # Lipidomics and structural input CSVs
|   |-- proteomics/                  # Bulk liver proteomics validation data
|   `-- Suppl_Table_1_OverviewData.xlsx
|-- scripts/
|   |-- workflows/                   # Main workflows 00-10
|   |-- core/                        # Data loading, DE, enrichment, MOFA, structural helpers
|   |-- config/                      # Dataset paths, contrasts, plotting config
|   |-- visualization/               # Plot helpers
|   |-- support_functions.R          # Compatibility wrapper to config/support_functions.R
|   `-- seurat_compatibility.R
|-- support/                         # Gene sets, regulons, hormone targets, LINEX2 network
|-- docs/
|   |-- AGENT_GUIDE.md
|   `-- workflow_plans/
|-- results/
|   |-- qc/
|   |-- pseudobulk/
|   |-- lipidomic/
|   |-- mofa/
|   |-- zonation/
|   |-- validation/
|   `-- xenium/
`-- tests/
```

## Workflow Order

Most users should start at workflow 01. Workflow 00 is a reference/HPC raw-data
preprocessing workflow and is inactive by default.

```r
# Optional reference workflow; requires original raw alignment environment
# source("scripts/workflows/00_snrnaseq_preprocessing.R")

source("scripts/workflows/01_data_preparation_qc.R")
source("scripts/workflows/02_fcg_pseudobulk_analysis.R")
source("scripts/workflows/03_fcg_lipidomic_analysis.R")
source("scripts/workflows/04_mofa_integration.R")
source("scripts/workflows/05_zonation_analysis.R")
source("scripts/workflows/06_structural_analysis.R")
source("scripts/workflows/07_validation_experiments.R")
source("scripts/workflows/08_variance_analysis.R")

# Optional external-data workflow; requires large GSE179250 inputs outside repo
# python scripts/workflows/09_human_analysis.py

# Optional collaborator Xenium scripts; require raw Xenium data under the
# workflow folder and are run from scripts/workflows/10_xenium_analysis/
# python 00_infer_doublets.py
# python 01_merge_doublets.py
# python 02_compositional_analysis.py
# python 03_DGE_analysis.py
```

Run workflows 00-08 from the repository root. Workflows 09 and 10 are optional
Python analyses with their own external or folder-local path assumptions, as
documented below and in their plan files.

## Analysis Workflows

### Workflow 00: snRNA-seq Pre-processing

**Script:** `scripts/workflows/00_snrnaseq_preprocessing.R`

Reference pipeline for raw sci-RNA-seq import, UMI QC, Seurat clustering,
scDblFinder doublet removal, cell-type annotation, pseudobulking, and
hepatocyte zonation annotation. Both FCG and GDX processing blocks are wrapped
in `if (FALSE)` and depend on original DKFZ HPC paths plus a placeholder import
function.

Key intended outputs include `pseudobulks_scc.rds`,
`heps_zonationsubsetted_pseudobulks_gonad.rds`, and GDX pseudobulk RDS files.

### Workflow 01: Data Preparation and QC

**Script:** `scripts/workflows/01_data_preparation_qc.R`

Converts preprocessed pseudobulk Seurat objects into CSV matrices and metadata.
The main RNA CSVs are split by sample columns, with all genes retained in each
file.

Key outputs:

- `data/RNA/scaled_RNA_1.csv`
- `data/RNA/scaled_RNA_2.csv`
- `data/RNA/scaled_RNA_3.csv`
- `data/RNA/meta_data_df.csv`
- `data/RNA/scaled_RNA_zonation.csv`
- `data/RNA/meta_data_df_zonation.csv`
- `data/RNA/scaled_RNA_gdxhr_{1,2,3}.csv`
- `data/RNA/meta_data_df_gdxhr.csv`
- `results/qc/`

### Workflow 02: FCG Pseudobulk RNA Analysis

**Script:** `scripts/workflows/02_fcg_pseudobulk_analysis.R`

Runs RNA differential expression for 6 FCG contrasts across 5 cell types, then
adds GDX RNA validation, GO pathway enrichment, TF activity inference, Xenium
panel validation, hepatocyte PCA coordinates, and bulk proteomics validation.

Key outputs:

- `results/pseudobulk/ttop_full_list.RData`
- `results/pseudobulk/tval_full_df.RData`
- `results/pseudobulk/tval_full_df.csv`
- `results/pseudobulk/gdxrh/gdhxr_groups_ttop_list.RData`
- `results/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv`
- `results/pseudobulk/volcanos/`
- `results/pseudobulk/pca/`
- `results/pseudobulk/pathways/pathway_scores_df_all.RData`
- `results/pseudobulk/pathways/tf_scores_df_all.RData`
- `results/pseudobulk/gdx_validation_comparison.pdf`
- `results/pseudobulk/xenium_validation_volcano.pdf`
- `results/pseudobulk/proteo_validation_comparison.pdf`

### Workflow 03: FCG Lipidomic Analysis

**Script:** `scripts/workflows/03_fcg_lipidomic_analysis.R`

Runs lipid differential abundance for 6 FCG contrasts, optional GDX lipid
validation, lipid-class enrichment, LINEX2 network extraction, and per-class
chain-length/double-bond bubble plots.

Key outputs:

- `results/lipidomic/ttop_full_list.RData`
- `results/lipidomic/gdxhr_ttop_list.RData`
- `results/lipidomic/gdxhr_t_table.csv`
- `results/lipidomic/fcg_gdx_correlation.pdf`
- `results/lipidomic/volcanos/`
- `results/lipidomic/pca/`
- `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData`
- `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv`
- `results/lipidomic/enrichment/lipid_class_enrichment_plot.pdf`
- `results/lipidomic/network/significant_lipids_network_SIF.csv`
- `results/lipidomic/network/significant_lipids_network_ATT.csv`
- `results/lipidomic/bubbles_DBC/`

### Workflow 04: MOFA Integration

**Script:** `scripts/workflows/04_mofa_integration.R`

Builds a 6-view MOFA input (5 RNA cell types + lipidomics), trains or loads
5-15 factor MOFA models, selects the 12-factor model, extracts weights, enriches
factors for gonadal type, and builds gene-lipid and TF-gene-lipid networks.

Key outputs:

- `results/mofa/MOFA_data.RData`
- `results/mofa/RNA/mofa_res_12factor.hdf5`
- `results/mofa/main_model_RNA_lipid_weight.RData`
- `results/mofa/factor_analysis/factor_enrichment_gonadal.RData`
- `results/mofa/factor_analysis/*.pdf`
- `results/mofa/networks/gene_lipid_correlations_factor1.csv`
- `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_full_{SIF,ATT}.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_{SIF,ATT}.csv`

### Workflow 05: Hepatocyte Zonation Analysis

**Script:** `scripts/workflows/05_zonation_analysis.R`

Uses workflow 01 zonation CSVs to run duplicate-correlation limma models across
periportal/pericentral and gonadal contrasts. It selects sex-specific
zonation-biased genes through an entropy/stability/coverage/specificity
threshold sweep, builds zonation-specific gene-lipid networks, and maps lipid
classes to zonation categories using workflow 03 outputs.

Key outputs:

- `results/zonation/zonation_ttop_list.RData`
- `results/zonation/zonation_limma_results.RData`
- `results/zonation/t_table_zonation.csv`
- `results/zonation/ttop_full_long.csv`
- `results/zonation/sex_specific_zonation_genes_summary.csv`
- `results/zonation/sex_specific_zonation_panel_AB.pdf`
- `results/zonation/sex_specific_zonation_heatmap_faceted.pdf`
- `results/zonation/lipid_species_zonation_association.csv`
- `results/zonation/lipid_class_association_summary.csv`
- `results/zonation/threshold_sensitivity_total.csv`
- `results/zonation/optimal_threshold_metrics.csv`
- `results/zonation/networks/cormat_*_{SIF,ATT}.csv`
- `results/zonation/networks/cormat_*_full_{SIF,ATT}.csv`
- `results/zonation/networks/compressed_cormat_*_{SIF,ATT}.csv`

### Workflow 06: Lipid Structural Analysis

**Script:** `scripts/workflows/06_structural_analysis.R`

Current active workflow is lipid-only. It parses lipid structures from workflow
03 gonadal lipid DE results and tests whether lipid chain length, double-bond
count, and their interaction explain gonadal logFC within lipid classes.

Active model:

```r
logfc ~ chain_length + db_position + chain_length:db_position
```

Key outputs:

- `results/lipidomic/structural_model/lipids_with_structure.csv`
- `results/lipidomic/structural_model/lipid_only_class_coefficients.csv`
- `results/lipidomic/structural_model/lipid_only_class_summaries.csv`
- `results/lipidomic/structural_model/lipid_only_chain_length_effects.csv`
- `results/lipidomic/structural_model/lipid_only_db_position_effects.csv`
- `results/lipidomic/structural_model/lipid_only_interaction_effects.csv`
- `results/lipidomic/structural_model/lipid_only_summary_heatmap.pdf`
- `results/lipidomic/structural_model/lipid_only_plots/`

The older gene-lipid structural-regression sections are commented out in the
current script.

### Workflow 07: Validation Experiments

**Script:** `scripts/workflows/07_validation_experiments.R`

Validates MOFA and zonation network predictions using siRNA knockdown
lipidomics from primary hepatocytes, then evaluates whether knockdown-response
consistency varies by lipid structure.

Key outputs:

- `results/validation/knockdown_PCA.pdf`
- `results/validation/knockdown_ttop_list.RData`
- `results/validation/knockdown_t_table.csv`
- `results/validation/knockdown_heatmap.pdf`
- `results/validation/global_network_validation_scores.csv`
- `results/validation/[gene]_validation_scatter.pdf`
- `results/validation/zonation_network_validation_scores.csv`
- `results/validation/[gene]_zonation_validation_scatter.pdf`
- `results/validation/validation_summary.csv`
- `results/validation/structural_consistency/[gene]_structural_consistency.pdf`
- `results/validation/structural_consistency/consistency_structural_regression_summary.csv`

### Workflow 08: Lipidomic Variance Analysis

**Script:** `scripts/workflows/08_variance_analysis.R`

Compares lipid variability between gonadal groups using CVs and Brown-Forsythe
Levene tests at both lipid and lipid-class levels.

Key outputs:

- `results/lipidomic/variance_analysis/variance_analysis_results.csv`
- `results/lipidomic/variance_analysis/significant_differential_variance_lipids.csv`
- `results/lipidomic/variance_analysis/variance_analysis_by_class.csv`
- `results/lipidomic/variance_analysis/variance_analysis_summary.csv`
- `results/lipidomic/variance_analysis/cv_difference_volcano.pdf`
- `results/lipidomic/variance_analysis/cv_difference_volcano_by_class.pdf`

### Workflow 09: Human Liver Comparison

**Script:** `scripts/workflows/09_human_analysis.py`

Optional Python workflow that compares the mouse FCG hepatocyte gonadal
t-value signature with sex-associated gene expression in the Zhou et al. human
bulk liver cohort ([`GSE179250`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179250)).
It translates mouse signature genes to human ortholog symbols, scores the
signature with decoupler ULM, runs PyDESeq2 male-vs-female differential
expression in human liver, and generates concordant mouse-human ranked-gene
plots.

The raw human count/TPM matrices and GEO SOFT metadata are too large for this
repository and are expected under the script's external `DATA_DIR`. The checked
in signature seed is:

- `results/pseudobulk/hepa_gonad_signature.csv`

The script writes outputs to its external `RESULTS_DIR`, including:

- `ulm_sex_signature_estimates.csv`
- `adata_processed.h5ad`
- `deg_zhou.csv`
- `full_results_deg_zhou_malevsfemale.csv`
- `12_49_malevsfemale_deg_zhou.csv`
- `12_49_malevsfemale_full_results_deg_zhou.csv`
- `olderthan11female_oldvsyoung_full_results_deg_zhou.csv`
- `olderthan11female_oldvsyoung_deg_statmatrix_zhou.csv`
- `pca_by_sex.pdf`
- `ulm_enrichment_12_49.pdf`
- `volcano_male_vs_female_highlight.pdf`
- `volcano_male_vs_female_highlight_1249.pdf`
- `top_ranked_genes_female_dot.pdf`
- `top_ranked_genes_male_dot.pdf`
- `ulm_enrichment_female__old_vs_young.pdf`

### Workflow 10: Xenium Spatial Analysis

**Folder:** `scripts/workflows/10_xenium_analysis/`

Optional collaborator-provided Python workflow for Xenium spatial validation.
Unlike workflows 00-09, this is a folder of notebook-style scripts rather than
a single canonical workflow script. It infers neighboring-nucleus doublets from
Xenium nucleus boundaries, annotates and fuses likely hepatocyte doublets,
creates cell-type/doublet composition plots, and runs Scanpy Wilcoxon
differential expression for sex, hepatocyte-zone, and doublet comparisons.

The scripts use relative paths and should be treated as folder-local: run them
from `scripts/workflows/10_xenium_analysis/` with raw Xenium files under
`data_raw/` and output folders such as `data_tidy/` and `figures/` available.
The raw Xenium bundle is not committed.

Main scripts:

- `00_infer_doublets.py`
- `01_merge_doublets.py`
- `02_compositional_analysis.py`
- `03_DGE_analysis.py`

Key folder-local inputs include:

- `data_raw/output-XETG00078__0018188__Region_1__20240301__143248/nucleus_boundaries.csv.gz`
- `data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cell_feature_matrix.h5`
- `data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cells.csv.gz`
- `data_raw/XeniumTargetList_Mouse-Panel-1.xlsx`

Key folder-local outputs include:

- `data_tidy/inferred_components.csv`
- `inferred_components.json`
- `fused_anndata_annotated.h5ad`
- `data_tidy/fused_doublet.csv`
- `data_tidy/fused_celltype.csv`
- `figures/celltype_freq_cores.pdf`
- `figures/celltype_freq_sex.pdf`
- `figures/doublet_freq_core.pdf`
- `data_tidy/DEG_results_sex_per_celltype.xlsx`
- `DEG_results.xlsx`
- `Pairwise_DEG_hepatocytes.xlsx`
- `Doublets_DEG_results.xlsx`
- `Doublets_and_Sex_DEG_results_all_hepatocytes.xlsx`

## Requirements

Workflows 00-08 are R scripts and expect packages from CRAN and Bioconductor.
Core R packages include:

- Data manipulation: `dplyr`, `tidyr`, `tibble`, `readr`, `readxl`, `reshape2`
- Visualization: `ggplot2`, `ggrepel`, `pheatmap`, `cowplot`, `patchwork`,
  `gridExtra`, `ggridges`, `scales`
- Statistics: `limma`, `decoupleR`, `car`, `vsn`
- Single-cell/preprocessing: `Seurat`, `glmGamPoi`, optional `scDblFinder`
- Multi-omic integration: `MOFA2` and Python `mofapy2` for retraining models
- Gene sets: `GSEABase`

Workflow 09 is Python and expects packages including `anndata`, `scanpy`,
`decoupler`, `pydeseq2`, `pypath`, `pandas`, `numpy`, `matplotlib`, `seaborn`,
`scipy`, and `adjustText`.

Workflow 10 is Python and expects packages including `scanpy`, `anndata`,
`sccoda`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `tifffile`, `shapely`,
`rtree`, `networkx`, `tqdm`, and `openpyxl`.

## Key Results Files

### RNA Differential Expression

- `results/pseudobulk/ttop_full_list.RData`
- `results/pseudobulk/tval_full_df.csv`
- `results/pseudobulk/hepa_gonad_signature.csv`

### Lipid Differential Abundance

- `results/lipidomic/ttop_full_list.RData`
- `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv`

### MOFA and Networks

- `results/mofa/main_model_RNA_lipid_weight.RData`
- `results/mofa/networks/gene_lipid_correlations_factor1.csv`
- `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv`

### Zonation

- `results/zonation/ttop_full_long.csv`
- `results/zonation/sex_specific_zonation_genes_summary.csv`
- `results/zonation/lipid_class_association_summary.csv`
- `results/zonation/networks/`

### Structural and Variance Analysis

- `results/lipidomic/structural_model/lipids_with_structure.csv`
- `results/lipidomic/structural_model/lipid_only_class_coefficients.csv`
- `results/lipidomic/variance_analysis/variance_analysis_results.csv`
- `results/lipidomic/variance_analysis/variance_analysis_by_class.csv`

### Validation

- `results/validation/global_network_validation_scores.csv`
- `results/validation/zonation_network_validation_scores.csv`
- `results/validation/structural_consistency/consistency_structural_regression_summary.csv`

### Human Liver Comparison

- `results/pseudobulk/hepa_gonad_signature.csv`
- Workflow 09 outputs are written outside this repository by default; see
  `scripts/workflows/09_human_analysis.py` and
  `docs/workflow_plans/09_human_analysis.md`.

### Xenium Spatial Analysis

- `scripts/workflows/10_xenium_analysis/`
- Workflow 10 outputs are folder-local by default; see
  `docs/workflow_plans/10_xenium_analysis.md`.

## Data Availability

The checked-in CSV and RData outputs support result lookup and most downstream
analysis. Some raw/preprocessed inputs may not be present in all checkouts:

- `data/RNA/fcg_pseudobulks.rds`
- `data/RNA/fcg_zonation_pseudobulks.rds`
- `data/RNA/gdx_pseudobulks.rds`
- MOFA `.hdf5` models if retraining/loading is needed
- Workflow 09 human liver inputs from `GSE179250`, including the count matrix,
  TPM matrix, and SOFT metadata file.
- Workflow 10 raw Xenium files under `scripts/workflows/10_xenium_analysis/data_raw/`,
  including the 10x output folder and Xenium panel workbook.

Workflow 00 requires the original raw alignment environment and is not expected
to run in a standard checkout.

## Citation

Manuscript citation to be added upon publication.

## License

GNU General Public License v3.0. See `LICENSE`.
