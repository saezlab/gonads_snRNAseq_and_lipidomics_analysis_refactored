# Workflow 00: snRNA-seq Pre-processing Pipeline

**Script:** `scripts/workflows/00_snrnaseq_preprocessing.R`  
**Lines:** ~864  
**Purpose:** Reference workflow for raw sci-RNA-seq processing: import, UMI QC, Seurat clustering, doublet removal, cell type annotation, pseudobulk aggregation, and hepatocyte zonation annotation.

**Runtime status:** The two processing sections are wrapped in `if (FALSE)` blocks. The script documents the original HPC workflow and will not process data unless those guards are changed and the placeholder importer/HPC paths are supplied.

## Source Dependencies

```r
# No source() calls; helper functions are defined inline.
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
library(glmGamPoi)

# Optional:
scDblFinder
SCpubr
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 0 | Configuration | 42-121 | Libraries, optional packages, DKFZ HPC paths, markers, and QC parameters. |
| - | Helper Functions | 125-380 | Import placeholder, UMI filtering, clustering, doublet removal, annotation, pseudobulking, zonation. |
| 1 | FCG Pre-processing | 384-594 | Main FCG cohort; inactive unless `if (FALSE)` is changed. |
| 2 | GDX Pre-processing | 599-823 | GDX validation cohort; inactive unless `if (FALSE)` is changed. |
| - | Complete | 827-864 | Print reference outputs and next step. |

## Helper Functions

| Function | Description |
|----------|-------------|
| `import_samples_seurat_raw()` | Placeholder that currently calls `stop()`; users must implement alignment-output import. |
| `filter_umi_outliers()` | Keep cells within UMI thresholds. |
| `run_clustering_workflow()` | Seurat normalize, variable features, scale, PCA, neighbors, clusters, UMAP. |
| `remove_doublets()` | Convert to SingleCellExperiment, run scDblFinder if available, and keep singlets. |
| `annotate_cell_types()` | Map cluster IDs to cell type labels. |
| `create_pseudobulk()` | `AggregateExpression()` wrapper. |
| `annotate_zonation()` | Classify hepatocytes using Cyp2e1/Cyp2f2 thresholds. |

## Input Files

- Raw sci-RNA-seq alignment output on the original DKFZ `/omics/` filesystem.
- `sci_metadata.xlsx`, sheet `metadata`, for FCG and GDX cohorts.
- GTF annotation files referenced in `config$fcg$gtf_file` and `config$gdx$gtf_file`.

## Output Files

**FCG**
- `imported.rds`
- `filtered.rds`
- `clustered.rds`
- `doubletless.rds`
- `doubletless_clustermarkers.rds`
- `annotated.rds`
- `pseudobulks_scc.rds`
- `heps.rds`
- `heps_zonationsubsetted.rds`
- `heps_zonationsubsetted_pseudobulks_gonad.rds`

**GDX**
- `imported.rds`
- `filtered.rds`
- `clustered.rds`
- `doubletless.rds`
- `annotated.rds`
- `pseudobulks_withoutpumpstatussubsetting.rds`
- `pseudobulks_withpumpstatussubsetting.rds`
- `heps.rds`
- `heps_zonationsubsetted.rds`
- `heps_zonationsubsetted_pseudobulks_withoutpumpstatussubsetting.rds`
- `heps_zonationsubsetted_pseudobulks_withpumpstatussubsetting.rds`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `config` | list | FCG and GDX base directories, raw folders, metadata files, and GTF paths. |
| `cell_type_markers` | named list | Marker genes for Hepatocytes, Immune cells, Endothelial cells, Stellate cells, and Cholangiocytes. |
| `zonation_markers` | named list | `Cyp2e1` for pericentral and `Cyp2f2` for periportal. |
| `qc_params` | list | UMI lower bound 500, upper bound mean + 3 SD, doublet rate 0.02, PCA dims, cluster resolution 0.03. |
| `fcg_cluster_mapping` | character | FCG clusters 0-8 to 5 major cell types. |
| `gdx_cluster_mapping` | character | GDX clusters 0-5 to 5 major cell types. |

## Coherence Assessment

**Status: GOOD as a reference workflow.** It is intentionally inactive and HPC-specific. It is not a portable rerun entry point.

## Vulnerabilities to Keep in Mind When Modifying

1. **Importer is a placeholder:** `import_samples_seurat_raw()` stops immediately. Implement import logic before enabling either processing block.
2. **Hardcoded HPC paths:** All active paths in `config` point to the original DKFZ filesystem.
3. **Hardcoded cluster mappings:** Cell type labels assume the original cluster identities and resolution.
4. **Optional doublet removal:** Missing `scDblFinder` skips doublet removal with a warning.
5. **`setwd()` inside inactive blocks:** Enabling blocks changes the working directory; reruns should restore or avoid side effects.

## Cluster-to-Cell-Type Mappings

**FCG**

| Cluster | Cell Type |
|---------|-----------|
| 0, 1, 7 | Hepatocytes |
| 2 | Endothelial cells |
| 3, 5, 6 | Immune cells |
| 4 | Stellate cells |
| 8 | Cholangiocytes |

**GDX**

| Cluster | Cell Type |
|---------|-----------|
| 0, 4 | Hepatocytes |
| 1 | Endothelial cells |
| 2 | Immune cells |
| 3 | Stellate cells |
| 5 | Cholangiocytes |

## Notes for Context Restoration

1. Most reproducible reruns should start with workflow 01.
2. Zonation classification is `Cyp2e1 > 1 & Cyp2f2 < 1` for pericentral, and `Cyp2e1 < 1 & Cyp2f2 > 1` for periportal.
3. FCG wild-type controls (`XYT-WT`) are excluded before creating `pseudobulks_scc.rds`.
4. GDX pseudobulks are written both with and without pump-status subsetting.
