# Gonadal Hormone Coordination of Sex-Specific Adult Liver Physiology

Multi-omic analysis using the Four Core Genotypes (FCG) mouse model to quantify the relative contributions of gonadal hormones versus sex chromosome complement to hepatic sexual dimorphism.

## Project Overview

This repository contains a complete, refactored analysis pipeline for investigating sex-specific liver physiology using:
- **Single-nucleus RNA-seq (snRNA-seq)**: 466,986 high-quality nuclei from 48 liver samples (43 FCG + 5 WT males)
- **Untargeted lipidomics**: 1,386 lipid species across 42 lipid classes from 42 samples
- **Multi-omic integration**: MOFA2-based integration revealing hormone-dependent coordination

## Key Findings

- **Gonadal hormones (not sex chromosomes) are the primary driver** of sex differences in liver physiology
- ~**4,000 gonad-dependent genes** identified in hepatocytes (FDR ≤ 0.05)
- **733 of 1,386 lipids (53%)** show significant gonad-dependent abundance
- **Strong coordination** between gonad-dependent gene expression and lipid metabolism
- **Spatial zonation** (periportal vs. pericentral) modulates gonad-dependent effects
- Key regulatory genes validated: **Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9, Mtmr11**

## Biological Context: Four Core Genotypes (FCG) Model

The FCG model allows separation of gonadal hormone effects from sex chromosome effects:

| Group | Sex Chromosomes | Gonads | Description |
|-------|----------------|--------|-------------|
| XXO | XX | Ovaries | Chr3^Sry-^ females |
| XXT | XX | Testes | Chr3^Sry+^ females with testes |
| XYO | XY | Ovaries | Chr3^Sry-^ males with ovaries |
| XYT | XY | Testes | Chr3^Sry+^ males |

The Sry gene (normally on Y chromosome) is translocated to chromosome 3, enabling creation of XY mice with ovaries and XX mice with testes.

**Key comparisons:**
- **Gonadal effect**: Testes vs. Ovaries (across all chromosome backgrounds) - **PRIMARY COMPARISON**
- **Chromosomal effect**: XY vs. XX (across both gonadal backgrounds)

## AI-Assisted Data Exploration (MCP Server)

This repository includes an MCP (Model Context Protocol) server that allows Claude or other AI assistants to dynamically query the analysis results. This enables natural language questions about the data without writing code.

### Features
- Query differential expression for specific genes/lipids
- Explore pathway and TF enrichment
- Query MOFA integration results and gene-lipid networks
- Get summary statistics and compare across cell types

### Setup for Claude Code

```bash
# Install Python dependencies
pip install mcp pandas pyreadr

# Add the MCP server to Claude Code
claude mcp add --transport stdio fcgl-analysis -- python mcp/fcgl_server.py
```

### Setup for Claude Desktop

Add to your Claude Desktop configuration file:

**Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "fcgl-analysis": {
      "command": "python",
      "args": ["C:/path/to/gonads_snRNAseq_and_lipidomics_analysis_refactored/mcp/fcgl_server.py"],
      "env": {
        "FCGL_DATA_PATH": "C:/path/to/gonads_snRNAseq_and_lipidomics_analysis_refactored"
      }
    }
  }
}
```

**Note**: Replace `C:/path/to/` with the actual path to your repository.

### Example Questions

Once connected, you can ask Claude questions like:
- *"What is the differential expression of Cyp2b9 in hepatocytes?"*
- *"Show me the top 20 upregulated genes for the gonadal effect"*
- *"Which lipid classes are enriched in males vs females?"*
- *"What genes are correlated with triglycerides in the MOFA network?"*
- *"Compare Esr1 expression across all liver cell types"*

See `mcp/README.md` for full documentation of available tools.

## Repository Structure

```
.
├── data/                           # Raw and processed data files
│   ├── RNA/                        # snRNA-seq data (pseudobulk, split into 3 parts)
│   ├── metabolomics/               # Lipidomics data
│   └── Suppl_Table_1_OverviewData.xlsx  # Sample metadata
├── scripts/
│   ├── workflows/                  # Main analysis workflows (00-07)
│   ├── core/                       # Core analysis functions
│   │   ├── data_loading.R          # Data loading utilities
│   │   ├── differential_analysis.R # Differential expression/abundance
│   │   ├── enrichment_analysis.R   # Pathway and TF enrichment
│   │   ├── multiomics.R            # MOFA integration functions
│   │   └── structural_analysis.R   # Lipid structural analysis
│   ├── visualization/              # Plotting functions
│   │   ├── volcano_plots.R         # Volcano plot utilities
│   │   ├── dimred_plots.R          # PCA and dimensionality reduction
│   │   └── heatmaps.R              # Heatmap generation
│   ├── config/                     # Configuration files
│   │   ├── analysis_config.R       # Global analysis settings
│   │   ├── contrasts.R             # Experimental contrasts
│   │   ├── datasets.R              # Dataset paths
│   │   └── plotting_themes.R       # Visualization themes
│   ├── support_functions.R         # Core utility functions
│   └── legacy/                     # Original scripts (for reference)
├── support/                        # Reference files
│   ├── m5.all.v2023.2.Mm.symbols.gmt  # MSigDB GO gene sets (mouse)
│   ├── TF_regulons.RData           # CollecTRI TF target gene sets
│   ├── Gonadal_Hormone_Genes.csv   # ESR1/AR target genes
│   ├── linex2_adjacency_matrix_OT.csv  # LINEX2 lipid network
│   └── FCGL_Xenium_Panel.csv       # Spatial validation gene panel (100 genes)
├── results/                        # Analysis outputs
│   ├── qc/                         # Quality control
│   ├── pseudobulk/                 # RNA-seq differential analysis
│   │   ├── ttop_full_list.RData    # All DE results
│   │   ├── volcanos/               # Volcano plots
│   │   ├── pca/                    # PCA visualizations
│   │   └── pathways/               # Pathway and TF enrichment
│   ├── lipidomic/                  # Lipidomic differential analysis
│   │   ├── ttop_full_list.RData    # All lipid DE results
│   │   ├── volcanos/               # Lipid volcano plots
│   │   ├── enrichment/             # Lipid class enrichment
│   │   ├── network/                # LINEX2 lipid networks
│   │   └── bubbles_DBC/            # Structural bubble plots
│   ├── mofa/                       # Multi-omic integration
│   │   ├── MOFA_data.RData         # MOFA input data
│   │   ├── RNA/                    # Trained MOFA models
│   │   ├── factor_analysis/        # Factor enrichment analysis
│   │   └── networks/               # Gene-lipid correlation networks
│   ├── zonation/                   # Hepatocyte zonation analysis
│   └── multi-omic/structural_model/  # Lipid structural preferences
└── tests/                          # Validation and testing
    ├── TESTING_GUIDE.md            # Testing documentation
    └── test_*.R                    # Test scripts
```

## Analysis Workflows

The analysis is organized into **8 sequential workflows** (00-07) that can be run independently:

### Workflow 0: snRNA-seq Pre-processing Pipeline
**Script**: `scripts/workflows/00_snrnaseq_preprocessing.R`

⚠️ **NOTE**: This workflow requires access to raw sci-RNA-seq alignment data. Most users should start with Workflow 1.

Documents the raw data processing pipeline that produces the pseudobulk Seurat objects:

**Datasets processed**:
- **FCG (Four Core Genotypes)**: 43 samples + 5 wild-type males
- **GDX (Gonadectomy Hormone Replacement)**: Validation cohort

**Pipeline steps**:
1. Import raw sci-RNA-seq data from alignment output
2. Quality control: UMI filtering (remove outliers)
3. Standard Seurat workflow (normalize, scale, PCA, cluster, UMAP)
4. Doublet detection and removal (scDblFinder)
5. Cell type annotation using known markers
6. Pseudobulk aggregation by sample and cell type
7. Hepatocyte zonation annotation (periportal vs pericentral)

**Key outputs**:
- `annotated.rds` - Cell type annotated Seurat object
- `pseudobulks_scc.rds` - Main pseudobulk data
- `heps_zonationsubsetted_pseudobulks_gonad.rds` - Zonation-annotated hepatocyte pseudobulk

**Cell type markers used**:
- Hepatocytes: Mlxipl, Apoa2, Hnf4a
- Immune cells: Ptprc, Vsig4, Adgre1
- Endothelial cells: Btnl9, Clec4g, Pecam1
- Stellate cells: Lrat, Des, Reln
- Cholangiocytes: Krt19, Ankrd1, Epcam

**Zonation markers**:
- Pericentral: Cyp2e1 > 1 AND Cyp2f2 < 1
- Periportal: Cyp2e1 < 1 AND Cyp2f2 > 1

### Workflow 1: Data Preparation and Quality Control
**Script**: `scripts/workflows/01_data_preparation_qc.R`

- Loads Seurat object with snRNA-seq data (466,986 nuclei)
- Extracts and validates pseudobulk data (already aggregated)
- Performs comprehensive QC checks
- Generates initial PCA visualizations
- Outputs data files for downstream analyses

**Key outputs**:
- `data/RNA/scaled_RNA_[1-3].csv` - Expression data (split into 3 files)
- `data/RNA/meta_data_df.csv` - Sample metadata
- `results/qc/` - QC plots and reports

### Workflow 2: FCG Pseudobulk RNA-seq Analysis
**Script**: `scripts/workflows/02_fcg_pseudobulk_analysis.R`

Comprehensive differential expression analysis across 5 cell types and 6 contrasts:

**Main contrasts**:
- **Gonadal effect**: Testes vs. Ovaries (PRIMARY - reveals ~4,000 hepatocyte DEGs)
- **Chromosomal effect**: XY vs. XX (minimal effects - Y-linked + X-escapees only)
- **FCG group comparisons**: XXT vs. XXO, XYT vs. XYO, XYO vs. XXO, XYT vs. XXT

**Additional analyses**:
- Gonadectomy (GDX) validation experiment (6 contrasts)
- PCA visualization with ridgeline plots
- Pathway enrichment (MSigDB GO Biological Process)
- Transcription factor activity (CollecTRI regulons + ESR1/AR targets)
- Xenium spatial validation (100-gene panel)
- GDX correlation validation (r = -0.48)

**Key outputs**:
- `results/pseudobulk/ttop_full_list.RData` - All FCG DE results
- `results/pseudobulk/volcanos/` - Volcano plots for all contrasts
- `results/pseudobulk/pathways/` - Pathway and TF enrichment heatmaps
- `results/pseudobulk/gdx_validation_comparison.pdf` - FCG vs. GDX correlation

**Key findings**:
- Hepatocytes show strongest gonad effects (~4,000 DEGs)
- PC1 separates by gonadal type (~40% variance)
- Minimal sex chromosome effects across all cell types

### Workflow 3: FCG Lipidomic Analysis
**Script**: `scripts/workflows/03_fcg_lipidomic_analysis.R`

Untargeted lipidomics differential analysis and enrichment:

**Analyses**:
- Differential abundance testing (6 FCG contrasts)
- Gonadectomy validation experiment
- Lipid class enrichment analysis (decoupleR)
- LINEX2 network analysis (lipid-lipid conversions)
- Structural bubble plots (chain length vs. double bonds)

**Key outputs**:
- `results/lipidomic/ttop_full_list.RData` - All lipid DE results
- `results/lipidomic/enrichment/lipid_class_enrichment_plot.pdf` - Paper figure
- `results/lipidomic/network/` - LINEX2 network files (Cytoscape import)
- `results/lipidomic/bubbles_DBC/` - Structural bubble plots

**Key findings**:
- 733/1,386 lipids (53%) are gonad-dependent
- NO lipids show sex chromosome effects
- **Male-biased (Testes)**: TG, DG, glycerophospholipids (PI, PE, LPE, PA, LPC, PG, PC)
- **Female-biased (Ovaries)**: Sphingolipids (SM, Cer, HexCer)

### Workflow 4: Multi-Omic Integration (MOFA)
**Script**: `scripts/workflows/04_mofa_integration.R`

⚠️ **CRITICAL WORKFLOW** - Most innovative analysis in the manuscript

Multi-Omics Factor Analysis (MOFA2) integrating RNA expression (5 cell types) and lipidomics (1,254 lipids):

**Analyses**:
- Prepare MOFA input (6 views: 5 RNA cell types + 1 lipidomic)
- Train MOFA models (5-15 factors) or load pre-trained
- Select optimal model (model12 with 12 factors)
- Extract factor weights and compute TF enrichment
- Enrich factors for gonadal type (Factor 1: R² ~ 0.8)
- Build gene-lipid correlation networks

**Key outputs**:
- `results/mofa/MOFA_data.RData` - MOFA input (CRITICAL file)
- `results/mofa/RNA/mofa_res_12factor.hdf5` - Optimal model
- `results/mofa/networks/` - Gene-lipid correlation networks
- `results/mofa/factor_analysis/` - Variance explained and enrichment plots

**Key findings**:
- Factor 1 shows strongest gonad association (R² ~ 0.8)
- **Cyp2b9, Cyp2d9, Cyp2b13, Cyp2c38** highly correlated with multiple lipid classes
- TFs **Nfic, Dbp, Esr1** regulate Cyp genes
- Network reveals direct hormone-lipid coordination

### Workflow 5: Hepatocyte Zonation Analysis
**Script**: `scripts/workflows/05_zonation_analysis.R`

Spatial zonation analysis comparing periportal vs. pericentral hepatocytes:

**Analyses**:
- Differential expression between periportal and pericentral zones
- Zonation marker identification (Cyp2e1: pericentral, Cyp2f2: periportal)
- Interaction between zonation and gonadal type
- Pathway enrichment by zone

**Key outputs**:
- `results/zonation/` - Zonation-specific DE results and PCA plots

**Key findings**:
- **43% of gonad-dependent genes are pericentral-biased** (Cyp2e1+): includes Cyp2c38, Cyp2b13, Cyp2d9
- **18% of gonad-dependent genes are periportal-biased** (Cyp2f2+): includes Cyp2b9
- Zonation modulates hormone responsiveness

### Workflow 6: Lipid Structural Analysis
**Script**: `scripts/workflows/06_structural_analysis.R`

Identifies genes whose association with lipids depends on lipid structural features:

**Regression model**:
```r
gene_correlation ~ lipid_class + lipid_abundance + chain_length + db_position +
                   lipid_abundance:chain_length + lipid_abundance:db_position
```

**Analyses**:
- Compute gene-lipid correlations from MOFA Factor 1
- Regression analysis for structural preferences
- Identify substrate specificity patterns
- Generate 5-panel structural plots per gene

**Key outputs**:
- `results/multi-omic/structural_model/` - Regression results and summary tables
- `results/multi-omic/structural_model/gene_plots/` - Structural preference plots
- `results/multi-omic/structural_model/summary_figure_panel[A-C].pdf` - Publication-ready summary figures

**Key findings**:
- Identifies substrate specificity (chain length and saturation preferences)
- Example: Cyp2c38 shows structural preferences in lipid associations
- Reveals coordination between enzyme expression and lipid composition
- Generates integrated multi-panel summary figures for manuscript

### Workflow 7: Validation of Gene-Lipid Correlations
**Script**: `scripts/workflows/07_validation_experiments.R`

Validates gene-lipid correlations from MOFA integration using in vitro knockdown experiments.

**Analysis steps**:
- Load knockdown lipidomic data from primary hepatocytes
- Differential analysis comparing knockdown vs control
- Score consistency between predictions (correlations) and observations (knockdown effects)
- Validate both global and zonation-specific networks

**Genes validated**:
- **Global hepatocyte network**: Cyp2d9, Cyp2b9, Cyp2c38
- **Zonation-specific network**: Sult3a1, Mtmr11 (spatially biased genes)

**Key outputs**:
- `results/validation/knockdown_heatmap.pdf` - Overview of all knockdown effects
- `results/validation/[gene]_validation_scatter.pdf` - Consistency plots per gene
- `results/validation/global_network_validation_scores.csv` - Validation scores
- `results/validation/zonation_network_validation_scores.csv` - Zonation-specific scores

**Key findings**:
- Confirms that predicted gene-lipid correlations match experimental knockdown effects
- Quantifies consistency between FCGL mouse observations and in vitro experiments
- Validates zonation-specific effects for spatially biased genes

## Installation and Requirements

### R Version
- R >= 4.1.0

### Required R Packages

**Core analysis**:
```r
# Data manipulation
install.packages(c("dplyr", "tidyr", "tibble", "readr", "readxl"))

# Visualization
install.packages(c("ggplot2", "ggrepel", "pheatmap", "ComplexHeatmap",
                   "cowplot", "patchwork", "gridExtra", "circlize", "ggridges"))

# Statistics
BiocManager::install("limma")
install.packages("decoupleR")

# Single-cell
install.packages("Seurat")

# Multi-omic integration
BiocManager::install("MOFA2")

# Gene sets
BiocManager::install(c("GSEABase", "org.Mm.eg.db", "AnnotationDbi"))
```

## Running the Analyses

### Quick Start

Run all workflows sequentially:

```r
# Set working directory to repository root
setwd("/path/to/gonads_snRNAseq_and_lipidomics_analysis_refactored")

# Workflow 00 requires raw alignment data (HPC environment)
# Most users should start with Workflow 01
# source("scripts/workflows/00_snrnaseq_preprocessing.R")  # Optional: raw data processing

# Run workflows in order (starting from pre-processed data)
source("scripts/workflows/01_data_preparation_qc.R")
source("scripts/workflows/02_fcg_pseudobulk_analysis.R")
source("scripts/workflows/03_fcg_lipidomic_analysis.R")
source("scripts/workflows/04_mofa_integration.R")
source("scripts/workflows/05_zonation_analysis.R")
source("scripts/workflows/06_structural_analysis.R")
source("scripts/workflows/07_validation_experiments.R")
```

### Running Individual Workflows

Each workflow is self-contained and can be run independently (after completing prerequisite workflows):

```r
# Example: Run only RNA-seq analysis (requires Workflow 1 first)
source("scripts/workflows/01_data_preparation_qc.R")
source("scripts/workflows/02_fcg_pseudobulk_analysis.R")
```

### Expected Runtime

- **Workflow 1** (Data Prep): ~5-10 minutes
- **Workflow 2** (RNA-seq): ~15-20 minutes
- **Workflow 3** (Lipidomics): ~10-15 minutes
- **Workflow 4** (MOFA): ~2 hours (if training models) or ~10 minutes (pre-trained)
- **Workflow 5** (Zonation): ~10-15 minutes
- **Workflow 6** (Structural + Summary Figures): ~25-35 minutes
- **Workflow 7** (Validation): ~5-10 minutes

**Total runtime**: ~3 hours (with MOFA training) or ~1.5 hours (pre-trained models)

## Key Results Files

### Differential Expression/Abundance
- `results/pseudobulk/ttop_full_list.RData` - RNA-seq DE results (6 contrasts × 5 cell types)
- `results/lipidomic/ttop_full_list.RData` - Lipidomic DE results (6 contrasts)

### Multi-Omic Integration
- `results/mofa/MOFA_data.RData` - MOFA input (4.4M data points)
- `results/mofa/main_model_RNA_lipid_weight.RData` - Factor weights + TF enrichment
- `results/mofa/networks/` - Gene-lipid correlation networks

### Enrichment Analysis
- `results/pseudobulk/pathways/pathway_enrichment_heatmap.pdf` - GO pathway enrichment
- `results/pseudobulk/pathways/tf_activity_heatmap.pdf` - TF activity scores
- `results/lipidomic/enrichment/lipid_class_enrichment_plot.pdf` - Lipid class enrichment

### Validation
- `results/pseudobulk/gdx_validation_comparison.pdf` - FCG vs. gonadectomy correlation (r = -0.48)
- `results/pseudobulk/xenium_validation_volcano.pdf` - Spatial transcriptomics validation (100 genes)

## Data Availability

Raw data files required to run these analyses:

1. **snRNA-seq data**: Available in `data/RNA/`
   - `fcgl_pseudobulk_gonad_20241023.rds` - Seurat object (not in repo - contact authors)
   - `scaled_RNA_[1-3].csv` - Pseudobulk expression data (generated by Workflow 1)

2. **Lipidomics data**: Available in `data/metabolomics/`
   - `Lipid_data_file_JPC_cleaned.csv` - Cleaned lipidomic measurements

3. **Metadata**: Available in `data/`
   - `Suppl_Table_1_OverviewData.xlsx` - Sample metadata with group assignments

## Citation

If you use this code or data, please cite:

> [Manuscript citation to be added upon publication]

## Contact

For questions or issues:
- GitHub Issues: [Repository issues page]
- Email: [Contact information to be added]

## License

GNU General Public License v3.0 - see details in the Licence file.

## Acknowledgments

This refactored codebase was developed to improve reproducibility, maintainability, and code clarity while producing **identical results** to the original analysis.

**Key improvements**:
- 75-85% code reduction through modularization
- Comprehensive testing framework
- Detailed inline documentation
- Self-contained workflows
- Validated against legacy code (all results identical)
