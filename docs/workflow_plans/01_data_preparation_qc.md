# Workflow 01: Data Preparation and QC

## Current execution contract (audited 18 November 2026)

The full entry is a local preparation script: it immediately reads the ignored
`data/RNA/fcg_pseudobulks.rds`, then rewrites the main RNA split CSVs, metadata and
QC outputs. It is not a public-clone prerequisite: supplied main CSVs already feed
02/03/04, and Workflow17 independently stages current expanded-marker transports.
The optional GDX conversion uses `data/RNA/gdx_pseudobulks.rds`; these are original
whole-cell-type pseudobulks, not Workflow17's paired zonal GDX input.

Section11 now reads the public plain-matrix CSV transport selected by
`current_fcg_zonation_transport`; no Seurat object is needed for that subsection.
It writes `current_fcg_zonation_csv` and `current_fcg_zonation_metadata` with the
established sample-name conventions. `DATASETS$fcg_zonation$rna$source_transport`
selects the same source; `source_rds_local` records an excluded local original,
not a default fallback. The isolated subsection was tested previously against
approved CSV values/order. Its availability does not bypass the earlier mandatory
main RDS when the whole entry is executed.

There is no CLI parser, source guard, validate-only or fresh-output option.
`source()` and `Rscript` execute the body and can overwrite fixed data/QC paths;
use an isolated checkout/copy to reproduce preparation without changing curated
inputs. Config `DATASETS` is loaded but does not redirect the hardcoded main/GDX
RDS filenames.00's HPC outputs require an explicit external import/name handoff.

---

**Script:** `scripts/workflows/01_data_preparation_qc.R`  
**Lines:** ~1023
**Purpose:** Convert preprocessed pseudobulk Seurat objects into CSV matrices/metadata for downstream workflows, run QC checks, and generate exploratory PCA plots.

## Source Dependencies

```r
source("scripts/config/analysis_config.R")
source("scripts/config/analysis_release.R")
source("scripts/config/datasets.R")
source("scripts/support_functions.R")
source("scripts/seurat_compatibility.R")
source("scripts/visualization/dimred_plots.R")
# Zonation section11 only:
source("scripts/core/pseudobulk_csv.R")
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
| 11 | Current zonation CSV conversion | 797-875 | Read `current_fcg_zonation_transport` and write registered FCG compatibility CSV/metadata if the manifest is present. |
| 12 | GDX CSV Conversion | 876-988 | Convert `data/RNA/gdx_pseudobulks.rds` to split GDX CSV/metadata if present. |
| - | Complete | 990-1021 | Print output summary. |

## Input Files

- `data/RNA/fcg_pseudobulks.rds` - Main FCG pseudobulk Seurat object.
- `data/RNA/pseudobulk_csv/current_fcg_zonation/manifest.json` and all listed payloads - Optional section11 source, selected by `current_fcg_zonation_transport`.
- `data/RNA/gdx_pseudobulks.rds` - Optional GDX pseudobulk Seurat object.
- `data/Suppl_Table_1_OverviewData.xlsx` - Used during PCA metadata merge.

## Output Files

- `data/RNA/scaled_RNA_1.csv`
- `data/RNA/scaled_RNA_2.csv`
- `data/RNA/scaled_RNA_3.csv`
- `data/RNA/meta_data_df.csv`
- `data/RNA/zonation/current/fcg/scaled_RNA_zonation.csv` (`current_fcg_zonation_csv`)
- `data/RNA/zonation/current/fcg/meta_data_df_zonation.csv` (`current_fcg_zonation_metadata`)
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
| `zonation_rds_path` | character | Despite its retained variable name, the registered `current_fcg_zonation_transport` JSON manifest. |
| `gdx_rds_path` | character | `data/RNA/gdx_pseudobulks.rds` |

## Coherence Assessment

**Status: optional local ETL, with a public CSV-backed zonation subsection.** Full execution requires the local main FCG RDS before reaching optional section11/12; the optional guards do not make the complete entry runnable without it. The section11 guard checks a transport manifest; section12 checks the separate local GDX RDS.

## Notes for Context Restoration

1. The split RNA CSVs are split by **sample columns**, not genes. Each file keeps the `gene` column and all rows, which allows `load_rna_data()` to merge by `gene`.
2. Main downstream consumers are02/04 for main RNA,03 for the main RNA metadata sample roster,05 for current zonation compatibility CSVs, and02 for the GDX RNA CSVs.03 reads its own GDX lipid input and does not consume these GDX RNA CSVs.
3. The Seurat objects are already pseudobulk objects; columns represent pseudobulk samples, not single cells.
4. Metadata from the raw Seurat object initially lacks `group` and `allosome_type`; downstream loaders merge those from the overview workbook.
