# Workflow 00: snRNA-seq Pre-processing Pipeline

**Script:** `scripts/workflows/00_snrnaseq_preprocessing.R`
**Lines:** ~865
**Purpose:** Documents the pre-processing pipeline for raw sci-RNA-seq data: import, QC, clustering, doublet removal, cell type annotation, pseudobulking, and hepatocyte zonation annotation. Produces the Seurat pseudobulk objects consumed by workflow 01.

**NOTE:** Both processing sections are wrapped in `if (FALSE) { ... }` blocks because they require raw alignment data on the original HPC system. This script serves as documentation and a reproducibility reference. For users starting with pre-processed data, begin with workflow 01.

## Source Dependencies

```r
# No source() calls - self-contained with inline helper functions
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
library(glmGamPoi)
# Optional: scDblFinder, SCpubr
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 0 | Configuration | 1-121 | Libraries, HPC paths, markers, QC parameters |
| - | Helper Functions | 124-380 | 7 reusable functions (see below) |
| 1 | FCG Pre-processing | 384-592 | Main study cohort (43 + 5 WT samples) — `if (FALSE)` |
| 1.1 | Import raw data | 401-429 | Load metadata, import sci-RNA-seq alignment |
| 1.2 | UMI filtering | 432-445 | Lower (500) and upper (mean + 3*SD) thresholds |
| 1.3 | Clustering | 449-473 | Normalize, PCA, FindClusters (res=0.03), UMAP |
| 1.4 | Doublet removal | 476-502 | scDblFinder (2% expected rate) |
| 1.5 | Cell type annotation | 505-539 | 9 clusters -> 5 cell types |
| 1.6 | Pseudobulking | 542-558 | AggregateExpression by sample x cell type |
| 1.7 | Zonation annotation | 561-592 | Cyp2e1/Cyp2f2 markers -> pericentral/periportal |
| 2 | GDX Pre-processing | 598-820 | Validation cohort — `if (FALSE)` |
| 2.1-2.7 | (Same steps as FCG) | 610-819 | Different cluster mapping (6 clusters -> 5 types) |
| - | Workflow Complete | 826-865 | Summary printout |

## Helper Functions (defined inline)

| Function | Lines | Description |
|----------|-------|-------------|
| `import_samples_seurat_raw()` | 139-149 | **Placeholder** — users must implement for their pipeline |
| `filter_umi_outliers()` | 158-184 | Remove cells outside UMI thresholds |
| `run_clustering_workflow()` | 193-234 | Full Seurat workflow: normalize -> PCA -> cluster -> UMAP |
| `remove_doublets()` | 242-278 | scDblFinder-based doublet detection and removal |
| `annotate_cell_types()` | 286-307 | Map cluster IDs to cell type labels |
| `create_pseudobulk()` | 315-332 | AggregateExpression wrapper |
| `annotate_zonation()` | 345-380 | Periportal/pericentral classification via marker thresholds |

## Input Files

- Raw sci-RNA-seq alignment output (HPC-specific paths)
- `sci_metadata.xlsx` — Sample metadata with group, gonadal_type, allosome_type

## Output Files

**FCG:**
- `imported.rds` -> `filtered.rds` -> `clustered.rds` -> `doubletless.rds` -> `annotated.rds`
- `pseudobulks_scc.rds` — Main pseudobulk (consumed by workflow 01)
- `heps_zonationsubsetted_pseudobulks_gonad.rds` — Zonation pseudobulk (consumed by workflow 05)

**GDX:**
- Same progression as FCG
- `pseudobulks_withoutpumpstatussubsetting.rds` / `pseudobulks_withpumpstatussubsetting.rds`
- `heps_zonationsubsetted_pseudobulks_*.rds` (with/without pump subsetting)

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `config` | list | HPC paths for FCG and GDX datasets |
| `cell_type_markers` | named list | 5 cell types with 3 markers each |
| `zonation_markers` | named list | Cyp2e1 (pericentral), Cyp2f2 (periportal) |
| `qc_params` | list | UMI thresholds, doublet rate, PCA dims, cluster resolution |
| `fcg_cluster_mapping` | character(9) | FCG: clusters 0-8 -> 5 cell types |
| `gdx_cluster_mapping` | character(6) | GDX: clusters 0-5 -> 5 cell types |

## Coherence Assessment

**Status: GOOD** — Clean, well-documented reference pipeline. No runtime issues (code is in `if (FALSE)` blocks).

## Vulnerabilities to Keep in Mind When Modifying

1. **`import_samples_seurat_raw()` is a placeholder (line 139):** It calls `stop()` with instructions. Users must implement their own import logic. If someone sets `if (FALSE)` to `if (TRUE)` without implementing this function, section 1.1 will halt immediately.

2. **Cluster-to-celltype mappings are hardcoded (lines 513-523, 721-728):** The FCG mapping assumes 9 clusters (0-8) and the GDX mapping assumes 6 clusters (0-5). If clustering resolution or data changes, these mappings become invalid. The `annotate_cell_types()` function does check length but not biological correctness.

3. **HPC paths are hardcoded (lines 79-94):** All paths point to the DKFZ HPC `/omics/` filesystem. These must be updated for any other environment.

4. **Optional package handling (lines 58-71):** `scDblFinder` and `SCpubr` are checked with `requireNamespace()`. If `scDblFinder` is missing, doublet detection is silently skipped (no error, just a warning). This could produce different results compared to the original analysis.

5. **`setwd()` calls inside processing blocks (lines 408, 619):** Working directory is changed to dataset-specific paths. If the script is partially run or interrupted, the working directory may be left in an unexpected state.

## Cluster-to-Cell-Type Mappings

**FCG (9 clusters):**
| Cluster | Cell Type |
|---------|-----------|
| 0, 1, 7 | Hepatocytes |
| 2 | Endothelial cells |
| 3, 5, 6 | Immune cells |
| 4 | Stellate cells |
| 8 | Cholangiocytes |

**GDX (6 clusters):**
| Cluster | Cell Type |
|---------|-----------|
| 0, 4 | Hepatocytes |
| 1 | Endothelial cells |
| 2 | Immune cells |
| 3 | Stellate cells |
| 5 | Cholangiocytes |

## Notes for Context Restoration

When returning to this workflow:
1. This is a **reference/documentation** workflow — the processing code does NOT run by default
2. The actual pre-processed `.rds` files already exist in the `data/` directory
3. Both FCG and GDX follow the same 7-step pipeline (import -> filter -> cluster -> doublets -> annotate -> pseudobulk -> zonation)
4. Zonation logic: Cyp2e1 > 1 AND Cyp2f2 < 1 = pericentral; Cyp2e1 < 1 AND Cyp2f2 > 1 = periportal
5. Wild-type controls (XYT-WT) are excluded before pseudobulking (line 549)
6. GDX has an extra pump_status dimension not present in FCG
7. The `annotate_zonation()` function discards cells that don't clearly fall into either zone
