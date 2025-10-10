# Analysis Pipeline Execution Order

This document describes the correct order to run analysis scripts in this repository.

## Overview

The pipeline consists of 5 main steps:
1. RNA data preparation (Seurat to dataframes)
2. RNA differential analysis
3. Lipidomics differential analysis
4. Pathway enrichment analysis
5. Multi-omic integration (MOFA)

## Detailed Pipeline

### Step 1: RNA Data Preparation
**Script:** `scripts/seurat_to_df.R`

**Description:** Converts Seurat object containing pseudobulk RNA data to CSV files for downstream analysis.

**Inputs:**
- `data/RNA/fcgl_pseudobulk_gonad_20241023.rds` - Seurat object with pseudobulked scRNA-seq data

**Outputs:**
- `data/RNA/scaled_RNA_1.csv` - Scaled RNA expression (samples 1-79)
- `data/RNA/scaled_RNA_2.csv` - Scaled RNA expression (samples 80-159)
- `data/RNA/scaled_RNA_3.csv` - Scaled RNA expression (samples 160-215)
- `data/RNA/meta_data_df.csv` - Sample metadata

**Notes:** Data is split into 3 CSV files to stay under GitHub's file size limits.

---

### Step 2: RNA Differential Analysis
**Script:** `scripts/pseudobulk_DA.R`

**Description:** Performs differential expression analysis using limma on pseudobulk RNA data across different grouping variables (gonadal type, allosome type, full groups).

**Inputs:**
- `data/RNA/scaled_RNA_1.csv`
- `data/RNA/scaled_RNA_2.csv`
- `data/RNA/scaled_RNA_3.csv`
- `data/RNA/meta_data_df.csv`
- `scripts/support_functions.R`

**Outputs:**
- `results/pseudobulk/ttop_full_list.RData` - Differential expression results for all cell types and contrasts
- `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` - PCA coordinates for hepatocytes
- Various volcano plots and heatmaps in `results/pseudobulk/volcanos/`

**Key comparisons:**
- **Gonadal type:** Testes vs. Ovaries
- **Allosome type:** XX vs. XY
- **Full groups:** XX-Ovaries, XX-Testes, XY-Ovaries, XY-Testes

---

### Step 3: Lipidomics Differential Analysis
**Script:** `scripts/lipidomic_DA.R`

**Description:** Performs differential abundance analysis on lipidomics data and enrichment analysis for lipid classes.

**Inputs:**
- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` - Raw lipidomics data
- `data/Suppl_Table_1_OverviewData.xlsx` - Sample metadata
- `support/lipid_classification.csv` - Lipid class annotations
- `support/linex2_adjacency_matrix_OT.csv` - LINEX2 lipid-enzyme network
- `scripts/support_functions.R`

**Outputs:**
- `data/metabolomics/lipidomic.RData` - Processed lipidomics matrix (log2-transformed)
- `results/lipidomic/gonade_type_ttop.csv` - Full differential abundance results (testes vs. ovaries)
- `results/lipidomic/gonade_type_ttop_20.csv` - Top 20 most significant lipids
- `results/lipidomic/ttop_full_list.RData` - DA results for all contrasts
- `results/lipidomic/PCA_sample_coord_lipid.csv` - PCA coordinates
- `results/lipidomic/factor_enrichment_long.csv` - PC enrichment by sample metadata
- `results/lipidomic/significant_lipids_network_SIF.csv` - Lipid network (SIF format)
- `results/lipidomic/significant_lipids_network_ATT.csv` - Lipid attributes
- Various volcano plots in `results/lipidomic/volcanos/`

**Key findings:**
- 733/1,386 lipids are gonad-dependent (adj. p-value < 0.05)
- Triglycerides (TG) and diglycerides (DG) are testes-biased
- Sphingolipids are ovaries-biased

---

### Step 4: Pathway Enrichment Analysis
**Script:** `scripts/pseudobulk_pathways.R`

**Description:** Performs pathway enrichment analysis on differential expression results using MSigDB Gene Ontology gene sets.

**Inputs:**
- `results/pseudobulk/ttop_full_list.RData` - Differential expression results
- `support/m5.all.v2023.2.Mm.symbols.gmt` - MSigDB GO gene sets (mouse)
- `scripts/support_functions.R`

**Outputs:**
- `results/pseudobulk/tval_full_df.RData` - t-values for all genes and contrasts
- `results/pseudobulk/tval_full_df.csv` - CSV version of t-values
- `results/pseudobulk/pathways/pathway_scores_df_all.RData` - Pathway enrichment scores
- `results/pseudobulk/pathways/pathway_scores_df_gonade_top.pdf` - Heatmap of top gonad-dependent pathways
- `results/pseudobulk/pathways/pathway_scores_df_allosome_top.pdf` - Heatmap of top allosome-dependent pathways
- `results/pseudobulk/pathways/pathway_scores_df_group_top.pdf` - Heatmap of top group-specific pathways

**Key findings:**
- Testes-biased pathways: ribosomal subunit assembly, translation
- Ovaries-biased pathways: sex hormone metabolism, lipid biosynthesis

---

### Step 5: Multi-Omic Integration (MOFA)
**Script:** `scripts/MOFA_lipidcell_RNA.R`

**Description:** Integrates RNA expression (by cell type) and lipidomics data using Multi-Omics Factor Analysis (MOFA2). **Note: MOFA model training is computationally expensive and requires external dependencies. The script is designed to skip training and load pre-trained models.**

**Inputs:**
- `data/RNA/scaled_RNA_1.csv`
- `data/RNA/scaled_RNA_2.csv`
- `data/RNA/scaled_RNA_3.csv`
- `data/metabolomics/lipidomic.RData` - Processed lipidomics matrix
- `support/TF_regulons.RData` - Transcription factor regulons for enrichment
- `support/Gonadal_Hormone_Genes.csv` - Hormone-responsive gene sets (AR, ESR1)

**MOFA Training (SKIP IN NORMAL WORKFLOW):**
```r
# Lines 63-67 in MOFA_lipidcell_RNA.R (commented out by default)
# Saves MOFA input data, trains models with 5-15 factors
# Training is done via external Python script
```

**MOFA Input (saved for testing):**
- `data/mofa/MOFA_RNA_data.csv` - Input data in MOFA format (sample, feature, view, value)
  - RNA views: Hepatocytes, Endothelial, Kupffer, LSECs, Stellate
  - Lipidomic view: Liver lipidomic
  - 215 samples (after filtering to samples with both RNA and lipid data)

**Pre-trained MOFA Models (DO NOT RE-RUN):**
- `results/mofa/RNA/mofa_res_5factor.hdf5` through `mofa_res_15factor.hdf5`
- Model with 12 factors is selected as optimal (`model12`)

**Outputs:**
- `results/mofa/main_model_RNA_lipid_weight.RData` - MOFA weights with TF enrichment
- Various plots showing variance explained and factor enrichment

**Key findings:**
- Factor 1: Strongest gonadal type association (35% variance in hepatocytes)
- Factor 2: Strong lipidomic variance (40% variance)
- Transcription factor enrichment reveals hormone receptor activity (AR, ESR1) in multiple factors

**IMPORTANT FOR REFACTORING:**
- The MOFA input (`data/mofa/MOFA_RNA_data.csv`) should remain unchanged during refactoring
- If it changes, MOFA models need to be retrained (very time-consuming)
- Tests verify MOFA input consistency without re-running training

---

## Running the Complete Pipeline

### Option 1: Run individually (recommended for understanding)
```r
# Step 1
source("scripts/seurat_to_df.R")

# Step 2
source("scripts/pseudobulk_DA.R")

# Step 3
source("scripts/lipidomic_DA.R")

# Step 4
source("scripts/pseudobulk_pathways.R")

# Step 5 (loads pre-trained MOFA models)
source("scripts/MOFA_lipidcell_RNA.R")
```

### Option 2: Run with test validation
```bash
# Save reference outputs before refactoring
Rscript tests/save_reference_outputs.R

# Run the pipeline (or make your changes)
# ...

# Test that outputs haven't changed
Rscript tests/test_output_consistency.R
```

---

## Parallel / Independent Scripts

The following scripts can be run independently or in parallel as they analyze different subsets:

### Zonation Analysis
- `scripts/seurat_to_df_zonation.R` - Prepares zonation-specific RNA data
- `scripts/pseudobulk_DA_zonation.R` - DA for periportal vs. pericentral hepatocytes
- `scripts/pseudobulk_pathways_zonation.R` - Pathway analysis for zonation

### GDX-HR (Gonadectomy + Hormone Replacement) Validation
- `scripts/gdxhr_seurat_to_df.R` - Prepares GDX-HR RNA data
- `scripts/gdxhr_pseudobulk_DA.R` - DA for GDX-HR samples
- `scripts/gdxhr_lipidomic_DA.R` - Lipidomics for GDX-HR
- `scripts/comparison_fcg_gdhxr.R` - Compares FCG vs. GDX-HR results

### Wild-Type Comparisons
- `scripts/seurat_to_df_withWT.R` - Includes WT samples
- `scripts/pseudobulk_DA_WT.R` - Compares FCG to WT

### Validation Experiments
- `scripts/lipidomic_validation_DA.R` - Validation cohort lipidomics
- `scripts/validation_plot_first_experiment.R` - Plots for validation
- `scripts/proteo_analysis.R` - Proteomics validation

### Downstream Analysis
- `scripts/MOFA_lipidcell_downstream.R` - Downstream analysis of MOFA results
- `scripts/pseudobulk_TFs.R` - TF enrichment analysis
- `scripts/pseudobulk_LRs.R` - Ligand-receptor analysis
- `scripts/lipidomic_bubbles.R` - Bubble plots for lipid classes

### Figure Generation
- `scripts/figure_1_UMAP.R` - UMAP visualization
- `scripts/Fig2_Doublets_and_Sex_DEG_results_all_hepatocytes.R` - Figure 2 panels
- `scripts/Fig5_Xenium_jpc.R` - Xenium spatial validation figures
- `scripts/misc_volcanos.R` - Additional volcano plots

---

## Testing Strategy

### Before Refactoring
1. Run `tests/save_reference_outputs.R` to save checksums and statistics of current outputs
2. This creates `tests/reference_outputs/reference_checksums.rds`

### After Refactoring
1. Run your refactored code
2. Run `tests/test_output_consistency.R` to verify outputs are unchanged
3. Tests compare:
   - File checksums (to detect any changes)
   - Content checksums (to detect data changes even if file metadata differs)
   - MOFA input consistency (critical - must not change!)
   - Statistical properties (dimensions, column names, summary stats)

### Test Results
- **PASS:** Outputs are identical to reference
- **FAIL:** Outputs have changed (review differences)
- **WARNING:** File was missing in reference or new file created

If tests fail, review the changes carefully. If changes are expected (e.g., due to bug fixes or intentional improvements), you may need to update the reference outputs.

---

## Dependencies

### R Packages Required
```r
# Data manipulation
library(readr)
library(readxl)
library(reshape2)
library(dplyr)

# Visualization
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(cowplot)
library(gridExtra)
library(circlize)

# Single-cell analysis
library(Seurat)

# Differential analysis
library(limma)

# Pathway/TF enrichment
library(decoupleR)
library(GSEABase)

# Multi-omics integration
library(MOFA2)

# Testing
library(digest)
library(jsonlite)
```

### External Data Files
- MSigDB gene sets: `support/m5.all.v2023.2.Mm.symbols.gmt`
- TF regulons: `support/TF_regulons.RData`
- Lipid classifications: `support/lipid_classification.csv`
- LINEX2 network: `support/linex2_adjacency_matrix_OT.csv`
- Hormone gene sets: `support/Gonadal_Hormone_Genes.csv`

---

## Notes for Refactoring

### What Can Be Safely Changed
- Code organization and structure
- Function naming (as long as outputs remain the same)
- Adding intermediate outputs for debugging
- Improving code readability
- Adding error handling

### What Should NOT Change
- **MOFA input data** (`data/mofa/MOFA_RNA_data.csv`) - must remain identical
- Final output values (differential analysis results, pathway scores, etc.)
- File formats and column names in output files
- Random seeds (if any stochastic processes exist)

### Recommended Refactoring Steps
1. Save reference outputs
2. Refactor one script at a time
3. Test after each refactoring
4. If tests fail, investigate differences
5. Only proceed to next script if tests pass

### Common Refactoring Goals
- Modularize scripts into reusable functions
- Reduce code duplication
- Improve parameter handling
- Add logging and progress indicators
- Standardize plotting code
- Improve error messages
- Add input validation
