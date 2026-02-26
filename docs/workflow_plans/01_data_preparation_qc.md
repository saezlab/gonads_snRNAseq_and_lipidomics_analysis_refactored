# Workflow 01: Data Preparation & QC

**Script:** `scripts/workflows/01_data_preparation_qc.R`
**Lines:** ~870
**Purpose:** Load Seurat pseudobulk object, extract metadata and scaled expression, perform QC, split data into 3 CSVs for downstream use. Also converts zonation Seurat object to CSV.

## Source Dependencies

```r
source("scripts/config/analysis_config.R")
source("scripts/config/datasets.R")
source("scripts/support_functions.R")
source("scripts/seurat_compatibility.R")
source("scripts/visualization/dimred_plots.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 1-30 | Load config, set output dirs |
| 2 | Load Seurat Object | 31-80 | Read `fcgl_pseudobulk_gonad_20241023.rds` |
| 3 | Extract Metadata | 81-150 | Pull cell type, gonad, chromosome info from Seurat |
| 4 | QC Metrics | 151-250 | Compute per-sample QC (gene counts, library size) |
| 5 | Extract Scaled Data | 251-350 | Get scaled expression matrix from Seurat assay |
| 6 | Verify Structure | 351-420 | Check dimensions, sample names, gene counts |
| 7 | Split & Save CSVs | 421-520 | Split into `scaled_RNA_{1,2,3}.csv` + `meta_data_df.csv` |
| 8 | Exploratory PCA | 521-650 | PCA on all cell types, variance by experimental factors |
| 9 | Data Validation | 651-750 | Cross-check metadata vs expression, sample matching |
| 10 | QC Report | 751-824 | Print summary statistics, save QC plots |
| 11 | Zonation CSV Conversion | 825-870 | Convert `fcg_zonation_pseudobulks.rds` to `scaled_RNA_zonation.csv` + `meta_data_df_zonation.csv` |
| 12 | GDX CSV Conversion | 875-960 | Convert `gdx_pseudobulks.rds` to `scaled_RNA_gdxhr_{1,2,3}.csv` + `meta_data_df_gdxhr.csv` |

## Input Files

- `data/RNA/fcg_pseudobulks.rds` - Seurat pseudobulk object (main FCG)
- `data/RNA/fcg_zonation_pseudobulks.rds` - Zonation Seurat object (hepatocyte periportal/pericentral)
- `data/RNA/gdx_pseudobulks.rds` - GDX Seurat pseudobulk object (gonadectomy validation)

## Output Files

- `data/RNA/scaled_RNA_1.csv`, `scaled_RNA_2.csv`, `scaled_RNA_3.csv` - Split expression matrices
- `data/RNA/meta_data_df.csv` - Sample metadata
- `data/RNA/scaled_RNA_zonation.csv` - Zonation scaled expression matrix (20,123 genes x 86 samples)
- `data/RNA/meta_data_df_zonation.csv` - Zonation metadata (sample, sample_name, zonation_type, gonadal_type)
- `data/RNA/scaled_RNA_gdxhr_{1,2,3}.csv` - GDX split expression matrices (20,123 genes x 110 samples)
- `data/RNA/meta_data_df_gdxhr.csv` - GDX metadata (sample, sample_name, cell_type, group)
- `results/qc/` - QC diagnostic plots

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `sobj` | Seurat | Main pseudobulk Seurat object |
| `scaled_RNA` | data.frame | Full scaled expression matrix (genes x samples) |
| `meta_data_df` | data.frame | Sample metadata with cell type, gonad, chromosome |

## Coherence Assessment

**Status: GOOD** - Linear ETL pipeline, no ordering issues.

- Flow is straightforward: Load -> Extract -> Validate -> Split -> Save
- No conditional logic or skip-ahead dependencies
- All outputs are self-contained

## Notes for Context Restoration

When returning to this workflow:
1. This is the entry point - requires raw Seurat `.rds` file
2. Downstream workflows (02-08) depend on the 3 CSV files this produces
3. The gene split is arbitrary (by column index), not biological - just for file size management
4. Sample names in the Seurat object determine naming convention for all downstream analyses
5. Section 11 converts the zonation Seurat object (`fcg_zonation_pseudobulks.rds`) to CSVs - used by workflow 05.
6. Section 12 converts the GDX Seurat object (`gdx_pseudobulks.rds`) to split CSVs - used by workflows 02 (step 3B) and 03.
