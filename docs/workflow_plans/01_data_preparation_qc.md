# Workflow 01: Data Preparation and QC

## CSV transport update - 18 November 2026

Section11 now reads `current_fcg_zonation_transport` through the plain-matrix
`pb_components()` reader. `DATASETS$fcg_zonation$rna$source_transport` is the public
source; `source_rds_local` records the excluded local original only. Dataset
source validation verifies CSV payload/content hashes without Seurat. The isolated
section11 conversion was executed and matches the approved current CSV values and
order. Full whole-cell preparation still remains an optional local preprocessing
step; downstream scientific defaults start from supplied pseudobulk CSVs.

---

## Current input routing - 16 November 2026

The optional FCG zonation conversion in section11 reads `current_fcg_zonation_rds`
and writes `current_fcg_zonation_csv` / `current_fcg_zonation_metadata` through
`scripts/config/analysis_release.R`. These resolve to the expanded-marker inputs
under `data/RNA/zonation/current/`. Whole-hepatocyte and the distinct GDXHR
conversions retain their original inputs and biological designs. The targeted
section11 conversion was executed in an isolated directory and exactly reproduced
the current FCG matrix and metadata; no canonical input was overwritten.

The section descriptions below retain historical line numbers and paths where
useful for provenance; the registry takes precedence for current zonation routing.

---

**Script:** `scripts/workflows/01_data_preparation_qc.R`  
**Lines:** ~1021  
**Purpose:** Convert preprocessed pseudobulk Seurat objects into CSV matrices/metadata for downstream workflows, run QC checks, and generate exploratory PCA plots.

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
| 1 | Configuration | 51-80 | Load config/helpers and create `results/qc` and `data/RNA`. |
| 2 | Load Seurat Object | 83-148 | Read `data/RNA/fcg_pseudobulks.rds` and inspect assay layers. |
| 3 | Extract Metadata | 151-238 | Extract pseudobulk metadata and write `meta_data_df.csv`. |
| 4 | QC Metrics | 241-331 | Count pseudobulks per mouse/cell type, design balance, and write QC tables/plots. |
| 5 | Extract Scaled Data | 334-378 | Extract `RNA` `scale.data` via Seurat compatibility helper. |
| 6 | Verify Structure | 381-413 | Validate pseudobulk sample IDs and dimensions. |
| 7 | Split and Save Data Files | 416-477 | Split by sample columns into three CSV files, keeping all genes in each file. |
| 8 | Initial Exploratory Analysis | 480-626 | Reload split files and metadata, then create PCA diagnostics. |
| 9 | Data Validation | 630-696 | Check missing/infinite values, variance, range, duplicates, and sample completeness. |
| 10 | QC Report | 699-792 | Write `results/qc/data_quality_report.txt`. |
| 11 | Zonation CSV Conversion | 796-870 | Convert `data/RNA/fcg_zonation_pseudobulks.rds` to zonation CSV/metadata if present. |
| 12 | GDX CSV Conversion | 876-988 | Convert `data/RNA/gdx_pseudobulks.rds` to split GDX CSV/metadata if present. |
| - | Complete | 990-1021 | Print output summary. |

## Input Files

- `data/RNA/fcg_pseudobulks.rds` - Main FCG pseudobulk Seurat object.
- `data/RNA/fcg_zonation_pseudobulks.rds` - Optional zonation pseudobulk Seurat object.
- `data/RNA/gdx_pseudobulks.rds` - Optional GDX pseudobulk Seurat object.
- `data/Suppl_Table_1_OverviewData.xlsx` - Used during PCA metadata merge.

## Output Files

- `data/RNA/scaled_RNA_1.csv`
- `data/RNA/scaled_RNA_2.csv`
- `data/RNA/scaled_RNA_3.csv`
- `data/RNA/meta_data_df.csv`
- `data/RNA/scaled_RNA_zonation.csv`
- `data/RNA/meta_data_df_zonation.csv`
- `data/RNA/scaled_RNA_gdxhr_1.csv`
- `data/RNA/scaled_RNA_gdxhr_2.csv`
- `data/RNA/scaled_RNA_gdxhr_3.csv`
- `data/RNA/meta_data_df_gdxhr.csv`
- `results/qc/data_quality_report.txt`
- `results/qc/pseudobulk_per_biosample.pdf`
- `results/qc/pseudobulk_per_celltype.pdf`
- `results/qc/sample_celltype_design.csv`
- `results/qc/samples_per_group_celltype.csv`
- `results/qc/PCA_all_samples_by_celltype.pdf`
- `results/qc/PCA_[CellType].pdf`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `seurat_path` | character | `data/RNA/fcg_pseudobulks.rds` |
| `seurat_obj` | Seurat | Main pseudobulk object |
| `metadata` | data.frame | Pseudobulk metadata written to `meta_data_df.csv` |
| `scaled_RNA` | data.frame | Full genes x pseudobulk sample expression matrix |
| `data_rna_dir` | character | `data/RNA` |
| `zonation_rds_path` | character | `data/RNA/fcg_zonation_pseudobulks.rds` |
| `gdx_rds_path` | character | `data/RNA/gdx_pseudobulks.rds` |

## Coherence Assessment

**Status: GOOD.** This is the main ETL entry point after raw preprocessing. Optional zonation/GDX conversions are guarded by file-existence checks and skip with warnings if the RDS files are absent.

## Notes for Context Restoration

1. The split RNA CSVs are split by **sample columns**, not genes. Each file keeps the `gene` column and all rows, which allows `load_rna_data()` to merge by `gene`.
2. Main downstream consumers are workflows 02 and 04 for main RNA data, workflow 05 for zonation CSVs, and workflow 02 for GDX RNA CSVs.
3. The Seurat objects are already pseudobulk objects; columns represent pseudobulk samples, not single cells.
4. Metadata from the raw Seurat object initially lacks `group` and `allosome_type`; downstream loaders merge those from the overview workbook.
