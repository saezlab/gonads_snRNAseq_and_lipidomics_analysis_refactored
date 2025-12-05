# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains analysis code for **"Gonadal hormone coordination of sex-specific adult liver physiology"**, a multi-omic study using the Four Core Genotypes (FCG) mouse model to quantify the relative contributions of gonadal hormones versus sex chromosome complement to hepatic sexual dimorphism.

**Main Findings:**
- **Gonadal hormones (not sex chromosomes) are the primary driver** of sex differences in liver physiology
- ~4,000 gonad-dependent genes identified in hepatocytes (>4000 DEGs at FDR ≤ 0.05)
- 733 of 1,386 lipids (53.9%) show significant gonad-dependent abundance
- **Strong coordination** between gonad-dependent gene expression and lipid metabolism
- **Spatial zonation** (periportal vs. pericentral hepatocytes) modulates gonad-dependent effects
- Key regulatory genes identified and validated: **Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9, Mtmr11**

**Dataset Scale:**
- 466,986 high-quality nuclei from 43 FCG livers + 5 wild-type males (sci-RNA-seq)
- 1,386 lipid species across 42 lipid classes (untargeted lipidomics)
- 5 major liver cell types: Hepatocytes, Immune cells, Endothelial cells, Stellate cells, Cholangiocytes

**Key Research Questions:**
- What are the relative contributions of gonadal hormones vs. sex chromosome complement?
- How do gonadal hormones coordinate cell type-specific gene expression and tissue lipid abundance?
- What is the role of liver spatial zonation in sex-specific physiology?

## Development Best Practices for Claude Code

**⚠️ CRITICAL: Always test scripts after modifications! ⚠️**

When modifying any R script in this repository, you MUST:

1. **Run the script immediately after making changes** to catch bugs early
   - Use: `Rscript workflows/[script_name].R` or `source("workflows/[script_name].R")`
   - Don't wait until multiple changes accumulate
   - Small bugs compound quickly in data analysis pipelines

2. **Use explicit namespace qualification for tidyverse functions** to avoid conflicts:
   - `dplyr::select()`, `dplyr::filter()`, `dplyr::arrange()`, `dplyr::mutate()`
   - `tidyr::pivot_longer()`, `tidyr::pivot_wider()`
   - Common conflicts: `MASS::select()`, `stats::filter()`, Bioconductor packages

3. **Add debug output liberally** when modifying functions:
   - Print object classes: `cat("Class:", class(obj), "\n")`
   - Print dimensions: `cat("Dimensions:", nrow(obj), "x", ncol(obj), "\n")`
   - Print previews: `print(head(obj, 3))`
   - Print sample names/IDs to catch matching issues
   - These can be removed later if needed, but are invaluable during development

4. **Check data format assumptions**:
   - Matrix vs. data frame
   - Wide vs. long format
   - Row names vs. columns
   - Sample name formatting (with/without suffixes like "_Testes", "_Hepatocytes")

5. **Common pitfalls in this codebase**:
   - MOFA functions return data frames, not lists (don't use `[[1]]` blindly)
   - decoupleR expects matrices in specific orientations (features × samples)
   - Sample names may have gonad or cell type suffixes that need cleaning
   - CSV files split into multiple parts need careful merging
   - Large matrices can cause memory issues - use `rm()` to free memory

6. **When errors occur**:
   - Read the full error message carefully
   - Check sample name matching issues first (very common)
   - Verify data dimensions match expectations
   - Use debug output to trace where data format changes

**Remember:** Testing immediately saves debugging time later. A small bug caught early is trivial; the same bug discovered after 10 more changes is a nightmare.

## Repository Structure

```
.
├── data/                      # Raw and processed data files
│   ├── RNA/                   # snRNA-seq data (scaled expression matrices split into 3 parts)
│   ├── metabolomics/          # Lipidomics data
│   ├── proteomics/            # Proteomics data
│   └── mofa/                  # MOFA integration data
├── scripts/                   # R analysis scripts
│   ├── support_functions.R    # Core utility functions (PCA, volcano plots, limma wrappers, etc.)
│   └── [various analysis scripts]
├── results/                   # Output files organized by analysis type
│   ├── pseudobulk/           # RNA-seq differential analysis results
│   ├── lipidomic/            # Lipidomic differential analysis results
│   ├── mofa/                 # MOFA integration results
│   └── multi-omic/           # Cross-omics integration results
│       └── structural_model/ # Lipid structural analysis results
└── support/                  # Reference files (gene sets, regulons, networks)
```

## Key Analysis Workflows

### 1. RNA-seq Pseudobulk Analysis

**Main script:** `scripts/pseudobulk_DA.R`

**Workflow:**
1. Load scaled RNA data from three CSV files in `data/RNA/` (files are split due to GitHub size constraints)
2. Perform PCA on all cell types and assess variance explained by experimental factors
3. Run differential expression analysis using limma for three contrasts:
   - **`gonadal_type`**: Testes vs. Ovaries (primary comparison showing ~4000 DEGs in hepatocytes)
   - **`allosome_type`**: XY vs. XX (minimal effect - few DEGs, mainly Y-linked or X-escape genes)
   - **`group`**: Specific group comparisons (e.g., XXT vs. XXO)
4. Generate volcano plots and save results in `results/pseudobulk/ttop_full_list.RData`
5. Identify periportal vs. pericentral hepatocytes using zonation markers (Cyp2e1, Cyp2f2)

**Key Finding:** Hepatocytes show the strongest gonad-dependent gene expression among all liver cell types.

**Related scripts:**
- `seurat_to_df.R`: Converts Seurat objects to dataframes for analysis
- `pseudobulk_pathways.R`: Pathway enrichment using decoupleR and MSigDB gene sets (GO biological processes)
- `pseudobulk_TFs.R`: Transcription factor activity inference (ESR1, AR, STAT5a/5b)
- `pseudobulk_LRs.R`: Ligand-receptor interaction analysis
- `gdxhr_*`: Gonadectomy validation experiments

### 2. Lipidomic Analysis

**Main script:** `scripts/lipidomic_DA.R`

**Workflow:**
1. Load lipidomic data from `data/metabolomics/Lipid_data_file_JPC_cleaned.csv`
2. Log2 transform and perform PCA (samples separate by gonad identity, not sex chromosome)
3. Run differential analysis using limma for multiple contrasts
4. **Lipid class enrichment** using decoupleR to identify enriched lipid categories:
   - **Male-biased (testes):** TG, DG, and glycerophospholipids (PI, PE, LPE, PA, LPC, PG, PC)
   - **Female-biased (ovaries):** Sphingolipids (SM, Cer, HexCer)
5. Network analysis using LINEX2 adjacency matrix to identify lipid-lipid conversion relationships
6. Generate bubble plots and volcano plots in `results/lipidomic/`

**Key Finding:** 733/1386 lipids (53.9%) are gonad-dependent, with NO lipids showing sex chromosome effects.

**Key concepts:**
- Lipid classes are extracted from lipid names (e.g., "PC(16:0/18:1)" → "PC")
- LINEX2 provides lipid-enzyme and lipid-lipid reaction networks from the literature
- Top 5 male-biased and female-biased lipids form an interconnected network

### 3. Multi-Omic Integration (MOFA)

**Main script:** `scripts/MOFA_lipidcell_RNA.R`

**Workflow:**
1. Combine RNA expression data (by cell type) and lipidomics into MOFA format
2. Train MOFA models with varying numbers of factors (5-15)
3. Compare models to select optimal number (model12 with 12 factors selected)
4. Extract factor scores (Z matrix) and weights (W matrix)
5. **Enrich factors for gonadal type** to identify hormone-associated variation (Factor 1 shows strongest gonad association)
6. **Perform transcription factor enrichment** on RNA weights (Nfic, Dbp, Esr1 identified)
7. **Build correlation networks** linking:
   - Top 2.5% MOFA features from hepatocyte RNA + lipidomics
   - Lipid metabolism genes (filtered by GO term "Cellular lipid metabolic process")
   - Top correlated TFs, genes, and lipids (Pearson r > 0.65)
8. **In vitro validation:** siRNA knockdowns of Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9 in primary hepatocytes

**Key Finding:** Cyp2b9, Cyp2d9, Cyp2b13, Cyp2c38 show highest correlation with multiple lipid classes and are regulated by Nfic, Dbp, Esr1.

**Key output:** `results/mofa/main_model_RNA_lipid_weight.RData`

### 4. Zonation-Specific Analysis

**Related scripts:**
- `pseudobulk_DA_zonation.R`
- `gdxhr_pseudobulk_DA_zonation.R`
- `gdxhr_pseudobulk_pathways_zonation.R`

**Key Finding:**
- 43% of gonad-dependent genes are **pericentral-biased** (Cyp2e1+, Glul+): includes Cyp2c38, Cyp2b13, Cyp2d9
- 18% of gonad-dependent genes are **periportal-biased** (Cyp2f2+, Gls2+, Hal+, Sds+): includes Cyp2b9
- Seven genes show high correlation with both TFs and lipids in zonation-specific manner:
  - **Female-biased pericentral:** Sult3a1, Abcb1a, Cyp2a4
  - **Male-biased periportal:** Mtmr11, Mindy3, Esco2, Kitl

**Validation:** Mtmr11 knockdown validated with lipidomics in primary hepatocytes

### 5. Lipid Structural Analysis

**Main scripts:**
- `lipid_rna_correaltion_for_structural_analysis.R` (note: typo in filename)
- `lipid_structure_analysis.R`
- `lipid_structure_visualisation.R`

**Workflow:**

**Step 1: Prepare Gene-Lipid Correlations** (`lipid_rna_correaltion_for_structural_analysis.R`)
1. Load MOFA model weights from `results/mofa/main_model_RNA_lipid_weight.RData`
2. Filter to top 50% of features by weight (both positive and negative ends of distribution)
3. Extract Factor 1 features (strongest gonad-associated factor)
4. Load scaled RNA expression and lipidomics data
5. Compute Pearson correlations between hepatocyte gene expression and lipid abundances
6. Filter genes to lipid metabolism genes only (GO term "GOBP_CELLULAR_LIPID_METABOLIC_PROCESS")
7. Save correlation matrix: `results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata`

**Step 2: Regression Analysis** (`lipid_structure_analysis.R`)
1. Load lipid structural data from `data/metabolomics/for_bubbles_with_names/` (one CSV per lipid class)
   - Each lipid has: chain length (C), double bond position (DB), fold change (LogFC)
2. Load gene-lipid correlation matrix from Step 1
3. For each gene, run linear regression models testing how lipid structural features predict gene-lipid correlations:
   - **Main model:** `gene_correlation ~ lipid_class + lipid_abundance + chain_length + db_position + lipid_abundance:chain_length + lipid_abundance:db_position`
   - **Class-specific models:** Same formula but run separately within each lipid class (requires ≥10 lipids per class)
4. Identify significant associations (p < 0.05) and categorize by:
   - Chain length effects
   - Double bond position effects
   - Interaction effects (abundance × structure)
5. Generate summary tables of genes with structural associations

**Step 3: Visualization** (`lipid_structure_visualisation.R`)
1. Load regression results and correlation data
2. For each gene with significant structural associations, create 5-panel plots per lipid class:
   - **Panel 1:** Gene-lipid correlation vs. Fold Change
   - **Panel 2:** Gene-lipid correlation vs. Chain Length
   - **Panel 3:** Gene-lipid correlation vs. Double Bond Position
   - **Panel 4:** Gene-lipid correlation vs. FC × Chain Length (interaction)
   - **Panel 5:** Gene-lipid correlation vs. FC × DB Position (interaction)
3. Color-code regression lines: green = significant (p < 0.05), black = non-significant
4. Save plots to `results/multi-omic/structural_model/gene_plots/`

**Key Findings:**
- Identifies genes whose association with lipids depends on lipid structural features
- Reveals substrate specificity: which genes/enzymes preferentially process lipids with specific chain lengths or saturation levels
- Example: Cyp2c38 shows structural preferences in its lipid associations
- Genes with multiple structural associations may be key regulators of specific lipid subclasses

**Biological Interpretation:**
- **Chain length associations:** Genes that preferentially act on short-chain vs. long-chain lipids
- **Double bond position associations:** Genes that distinguish saturated vs. unsaturated lipids
- **Interaction effects:** Genes whose response to lipid abundance depends on structural features

**Key Output Files:**
- `results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv`: Full correlation matrix
- `results/multi-omic/structural_model/lipids_with_structure.csv`: Lipid structural features
- `results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv`: All significant associations
- `results/multi-omic/structural_model/genes_associated_with_chain_length.csv`: Chain length-specific genes
- `results/multi-omic/structural_model/genes_associated_with_db_position.csv`: Saturation-specific genes
- `results/multi-omic/structural_model/genes_with_significant_interactions.csv`: Interaction effects
- `results/multi-omic/structural_model/genes_multiple_structural_associations.csv`: Multi-feature genes
- `results/multi-omic/structural_model/gene_lipid_structure_[gene].pdf`: Visualization for specific genes

### 6. Support Functions

**File:** `scripts/support_functions.R`

This file contains essential utility functions used across all analyses:

- `nicePCA()`: Generate comprehensive PCA plots with density distributions
- `volcano_nice()`: Create publication-quality volcano plots with asymptotic curves
- `runLimma()`: Wrapper for limma differential analysis with contrast matrix support
- `ttop_list_to_t_table()`: Convert limma toptables to t-statistic matrices
- `magicPlotMaker()`: Generate diagnostic plots (density, PCA, heatmaps, correlations)
- `import_gmt()`: Import GMT gene set files (e.g., MSigDB)
- `df_to_viper_regulon()`: Convert dataframes to VIPER regulon format
- `createLinearColors()`: Generate diverging color palettes for heatmaps
- `format_lipid_names()`: Standardize lipid nomenclature
- `clean_lipid_data()`: Filter and format lipidomic datasets

## Common Development Commands

### Running Individual Analysis Scripts

All scripts are designed to be run from the project root directory:

```r
# From R console or RStudio
source("scripts/pseudobulk_DA.R")
source("scripts/lipidomic_DA.R")
source("scripts/MOFA_lipidcell_RNA.R")
source("scripts/pseudobulk_pathways.R")

# Structural analysis pipeline (run in order)
source("scripts/lipid_rna_correaltion_for_structural_analysis.R")
source("scripts/lipid_structure_analysis.R")
source("scripts/lipid_structure_visualisation.R")
```

### Key R Packages Required

- **Data manipulation:** `dplyr`, `tidyr`, `reshape2`, `readr`, `readxl`
- **Visualization:** `ggplot2`, `ggrepel`, `pheatmap`, `ComplexHeatmap`, `cowplot`, `patchwork`, `gridExtra`, `circlize`
- **Statistics:** `limma`, `decoupleR`
- **Single-cell:** `Seurat`
- **Integration:** `MOFA2`
- **Gene sets:** `GSEABase`, `org.Hs.eg.db`, `AnnotationDbi`

## Important Data Files

### Input Data

- `data/Suppl_Table_1_OverviewData.xlsx`: Sample metadata (group assignments, gonads, chromosomes)
- `data/RNA/scaled_RNA_[1-3].csv`: Pre-scaled RNA expression data split into 3 files (466,986 nuclei)
- `data/RNA/meta_data_df.csv`: Single-cell metadata (cell types, samples, groups)
- `data/RNA/fcgl_pseudobulk_gonad_20241023.rds`: Original Seurat object with scaled data
- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv`: Cleaned lipidomic measurements (1,386 lipids)
- `data/metabolomics/for_bubbles_with_names/`: Lipid structural data (chain length, double bond position) organized by lipid class

### Support Files

- `support/TF_regulons.RData`: Transcription factor target gene sets (CollecTRI regulons)
- `support/Gonadal_Hormone_Genes.csv`: Curated hormone receptor target genes (ESR1, AR from Zheng et al. 2018)
- `support/m5.all.v2023.2.Mm.symbols.gmt`: MSigDB Gene Ontology biological process gene sets (mouse)
- `support/linex2_adjacency_matrix_OT.csv`: LINEX2 lipid-enzyme-lipid conversion network
- `support/lipid_classification.csv`: Lipid class hierarchy
- `support/Sry_Target_GeneCard.csv`: Sry transcription factor target genes
- `support/FCGL_Xenium_Panel.csv`: Custom 100-gene panel for Xenium spatial validation

### Output Data

- `results/pseudobulk/ttop_full_list.RData`: All differential expression results (list of lists: contrasts x cell types)
- `results/pseudobulk/tval_full_df.RData`: t-statistic matrix for all contrasts and cell types combined
- `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv`: PCA coordinates for hepatocytes
- `results/lipidomic/ttop_full_list.RData`: All lipidomic differential analysis results
- `results/lipidomic/gonade_type_ttop.csv`: Top gonad-dependent lipids
- `results/lipidomic/significant_lipids_network_SIF.csv`: LINEX2 network edges for significant lipids
- `results/lipidomic/significant_lipids_network_ATT.csv`: Attributes for significant lipids
- `results/mofa/main_model_RNA_lipid_weight.RData`: MOFA model weights (W matrix) with TF scores
- `results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata`: Gene-lipid correlation matrix for structural analysis
- `results/pseudobulk/pathways/pathway_scores_df_all.RData`: Pathway enrichment scores across all contrasts
- `results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv`: All significant gene-lipid structural associations
- `results/multi-omic/structural_model/genes_associated_with_chain_length.csv`: Genes with chain length preferences
- `results/multi-omic/structural_model/genes_associated_with_db_position.csv`: Genes with saturation preferences

## Analysis Patterns and Conventions

### Contrast Definitions

The project uses numeric vectors to define contrasts in limma:
- `c(-1, 2)`: Condition 2 - Condition 1 (e.g., Testes - Ovaries)
- `c(1, -4)`: Condition 1 - Condition 4 (e.g., XXT - XXO)
- Complex contrasts possible: `c(2, 2, -1, -3)` averages multiple conditions

**Example from paper:**
```r
comparisons <- list("testes_vs_ova" = c(-1,2))  # Comparing gonadal types
comparisons <- list("XY_vs_XX" = c(-1,2))  # Comparing sex chromosomes
comparisons <- list("XXT_vs_XXO" = c(1,-4))  # Comparing specific groups
```

### Color Schemes

Consistent colors used across analyses:
- **Gonadal type:** Ovaries = "red", Testes = "royalblue"
- **Allosome type:** XX = "red", XY = "royalblue" (or brown/cyan for specific comparisons)
- **Cell types:**
  - Hepatocytes = "#ca3659"
  - Endothelial cells = "#8f7b2f"
  - Immune cells = "#0d934d"
  - Stellate cells = "#1892b1"
  - Cholangiocytes = "#834f9f"
- **Heatmaps:** Red-white-blue or red-grey-blue diverging palettes

### Naming Conventions

- **Samples:** Format varies but typically `[MouseID]_[CellType]_[GonadType]`
- **Genes:** Mouse gene symbols (capitalized first letter: Esr1, Ar, Cyp2c29, Cyp2b9)
- **Lipids:** Standard notation with class and fatty acid chains: `PC(16:0/18:1)`, `TG(16:0_18:1_22:6)`
- **Cell types:** "Hepatocytes", "Endothelial cells", "Immune cells", "Stellate cells", "Cholangiocytes"
- **Groups:** XXO, XXT, XYO, XYT (or XYT-WT for wild-type controls)

## Working with Large Data Files

The RNA expression data is split into 3 CSV files due to GitHub size limits. Scripts handle this by:
1. Reading all three files separately (`scaled_RNA_1.csv`, `scaled_RNA_2.csv`, `scaled_RNA_3.csv`)
2. Merging by the "gene" column
3. Setting row names to gene symbols
4. Removing intermediate objects with `rm()` to free memory

When modifying data loading, maintain this pattern to ensure compatibility.

## Integration with decoupleR

Many analyses use decoupleR for enrichment (footprint-based activity estimation):
- **Input format:** Matrix of statistics (rows = features, columns = samples/contrasts)
- **Network format:** 3 columns (source, target, mor) where mor = mode of regulation (+1 or -1)
- **Method:** `run_ulm()` for univariate linear model (most common method used in this study)
- **Output:** Wide format with `dcast()` for heatmap visualization

**Common patterns:**
1. **Pathway enrichment:** Use MSigDB GMT files → decoupleR::run_ulm() → filter by score threshold
2. **TF activity:** Use CollecTRI regulons + custom ESR1/AR targets → decoupleR::run_ulm()
3. **Lipid class enrichment:** Create source-target network from lipid names → decoupleR::run_ulm()

## Pathway and TF Analysis

Pathway and transcription factor analyses follow this pattern:
1. Extract t-statistics from limma results (`ttop_list$[[i]]$t` column)
2. Format as matrix with genes as rows, contrasts as columns
3. Load gene sets (GMT format) or regulons (3-column dataframe)
4. Run `decoupleR::run_ulm()` with minsize parameter (default 10-20)
5. Filter top enriched pathways/TFs by score threshold (e.g., |score| > 11-12)
6. Visualize with pheatmap using custom color palettes

**Key pathways identified:**
- **Testes-biased:** Ribosomal complexes
- **Ovaries-biased:** Sex hormone metabolism, estrogen/testosterone 16/18-alpha hydroxylase pathways

**Key TFs identified:**
- ESR1, AR (directly from hormone signaling)
- STAT5a, STAT5b (growth hormone signaling)
- Nfic, Dbp, Esr1 (correlated with Cyp genes)

## Figure Generation

Key figures are organized in `Fig_panels/` and generated by specific scripts:
- `figure_1_UMAP.R`: UMAP visualization of cell types
- `Fig2_Doublets_and_Sex_DEG_results_all_hepatocytes.R`: Differential expression results
- `Fig5_Xenium_jpc.R`: Spatial transcriptomics validation with Xenium platform

Volcano plots are saved to `results/[analysis_type]/volcanos/` with consistent naming.

## Notes on Mouse Models

### Four Core Genotypes (FCG) Model

The FCG model allows separation of:
- **Gonadal effects:** Compare mice with testes vs. ovaries (regardless of chromosomes)
- **Chromosomal effects:** Compare XX vs. XY (regardless of gonads)
- **Groups:**
  - **XXO**: XX chromosomes + Ovaries (Chr3^Sry-^)
  - **XXT**: XX chromosomes + Testes (Chr3^Sry+^)
  - **XYO**: XY chromosomes + Ovaries (Chr3^Sry-^)
  - **XYT**: XY chromosomes + Testes (Chr3^Sry+^)

**Key principle:** Sry gene is translocated from Y chromosome to chromosome 3, allowing creation of XY mice with ovaries.

**Validation approaches used:**
- Wild-type controls (XYT-WT) compared to FCG XYT to control for model-specific effects
- Gonadectomy (GDX) experiments in wild-type mice (correlation: r = -0.48 with FCG)
- Xenium spatial transcriptomics validation (82 genes tested)

### Sample Sizes

- FCG samples: n=43 (approximately 10-11 per group)
- Wild-type males: n=5
- GDX experiments: separate cohorts of males and females

## Biological Context

### Key Biological Insights from the Paper

1. **Gonadal hormone dominance:** Sex differences in liver are ~95% due to gonadal hormones, ~5% due to sex chromosomes
2. **Lipid metabolism coordination:** Glycerolipids (TG, DG) male-biased; Sphingolipids (Cer, SM) female-biased
3. **Zonation matters:** Pericentral hepatocytes show different gonad-dependent patterns than periportal
4. **Cyp genes are key regulators:** Multiple Cyp2 family members directly regulate lipid abundance
5. **Hormone receptor activity:** ESR1 and AR are expressed across liver cell types and show activity
6. **GH signaling involved:** STAT5a/STAT5b activity suggests gonadal hormones work through growth hormone axis

### Experimental Validations Performed

1. **siRNA knockdowns:** Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9, Sult3a1, Mtmr11, Abcb1a in primary hepatocytes
2. **Xenium spatial validation:** 82 gonad-dependent genes validated in wild-type male vs. female tissue
3. **RNA FISH:** Cyp2d9 (male-biased), Cyp2b9 (female-biased) confirmed
4. **Proteomics:** Bulk liver tissue proteomics (33% of gonad-dependent genes show protein-level changes)
5. **GDX experiments:** Gonadectomy in wild-type mice confirms activational hormone effects

## Tips for Working with This Codebase

1. **⚠️ ALWAYS TEST IMMEDIATELY after modifying any script** - Run it to catch bugs before they compound (see "Development Best Practices" section above)
2. **Always use explicit namespacing** - `dplyr::select()`, `dplyr::filter()`, etc. to avoid package conflicts
3. **Always source `support_functions.R` first** - it contains essential utilities
4. **Check sample metadata** in `Suppl_Table_1_OverviewData.xlsx` to understand group assignments
5. **PCA plots are diagnostic** - samples should separate by gonadal type (not sex chromosome)
6. **Limma contrasts follow order** in unique(targets$condition) - verify with `print(unique(targets$condition))`
7. **MOFA factor selection** - model12 (12 factors) was chosen after comparing models 5-15
8. **Correlation networks** - Pearson r > 0.65 used as threshold for significant associations
9. **Zonation markers** - Cyp2e1 = pericentral, Cyp2f2 = periportal (Halpern et al. 2017)
10. **View() calls in code** - Remove or comment out `View()` calls when running non-interactively
