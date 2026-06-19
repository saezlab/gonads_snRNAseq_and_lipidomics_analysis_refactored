# Workflow 10: Xenium Spatial Analysis

**Folder:** `scripts/workflows/10_xenium_analysis/`
**Language:** Python
**Purpose:** Collaborator-provided Xenium spatial validation analysis for doublet inference, fused hepatocyte annotations, cell-type/doublet composition plots, and Scanpy differential expression.

## Source Dependencies

This is a folder of Python scripts, not a single canonical workflow file. The scripts use folder-local relative paths and should be run from:

```bash
cd scripts/workflows/10_xenium_analysis
```

Python packages imported by the scripts:

- `scanpy`
- `anndata`
- `sccoda`
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `tifffile`
- `shapely`
- `rtree`
- `networkx`
- `tqdm`
- `openpyxl`

## Script Order

| # | Script | Description |
|---|--------|-------------|
| 00 | `00_infer_doublets.py` | Reads Xenium nucleus boundaries, builds spatial polygon overlaps with `shapely`/`rtree`, identifies connected components of adjacent nuclei, and writes doublet/fusion candidates. |
| 01 | `01_merge_doublets.py` | Loads Xenium expression and cell metadata, clusters marker genes for cell-type annotation, merges valid same-hepatocyte-type doublets, assigns TMA/sex labels from spatial grid coordinates, and writes fused AnnData plus CellCharter-style grouping CSVs. |
| 02 | `02_compositional_analysis.py` | Reads the fused AnnData object and uses `sccoda` visualization helpers to plot cell-type and doublet compositions by core/sex. |
| 03 | `03_DGE_analysis.py` | Reads fused AnnData and raw Xenium expression, runs Scanpy Wilcoxon DGE for sex, hepatocyte-zone pairwise, doublet, and combined doublet/sex comparisons, then writes Excel tables and volcano PDFs. |

## Input Files

Folder-local raw inputs expected under `scripts/workflows/10_xenium_analysis/data_raw/`:

- `output-XETG00078__0018188__Region_1__20240301__143248/nucleus_boundaries.csv.gz`
- `output-XETG00078__0018188__Region_1__20240301__143248/cell_feature_matrix.h5`
- `output-XETG00078__0018188__Region_1__20240301__143248/cells.csv.gz`
- `XeniumTargetList_Mouse-Panel-1.xlsx`

Intermediate inputs produced by earlier scripts:

- `inferred_components.json` from `00_infer_doublets.py`
- `fused_anndata_annotated.h5ad` from `01_merge_doublets.py`

The raw Xenium files are large and are not committed to the repository.

## Output Files

Folder-local outputs written by the scripts include:

- `data_tidy/inferred_components.csv`
- `inferred_components.json`
- `fused_anndata_annotated.h5ad`
- `data_tidy/fused_doublet.csv`
- `data_tidy/fused_celltype.csv`
- Scanpy-generated marker/cell-type PDFs under the active figure directory.
- `figures/celltype_freq_cores.pdf`
- `figures/celltype_freq_sex.pdf`
- `figures/doublet_freq_core.pdf`
- `figures/DEG_volcano_plots.pdf`
- `data_tidy/DEG_results_sex_per_celltype.xlsx`
- `DEG_volcano_sex_per_celltyp.pdf`
- `DEG_results.xlsx`
- `Pairwise_DEG_volcano_hepatocytes.pdf`
- `Pairwise_DEG_hepatocytes.xlsx`
- `Doublets_DEG_volcano_plots.pdf`
- `Doublets_DEG_results.xlsx`
- `Doublets_and_Sex_DEG_volcano_plots_all_hepatocytes.pdf`
- `Doublets_and_Sex_DEG_results_all_hepatocytes.xlsx`

## Key Variables

| Variable | Script | Description |
|----------|--------|-------------|
| `sample_path` | `00_infer_doublets.py` | Relative path to the raw 10x Xenium output folder. |
| `x_lines`, `y_lines` | `00_infer_doublets.py`, `01_merge_doublets.py` | Hardcoded spatial grid boundaries used for TMA/sex assignment. |
| `type_markers` | `01_merge_doublets.py` | Xenium panel genes used for marker-based cell-type clustering after excluding gonad-dependent and XCI genes. |
| `cell_type_dict` | `01_merge_doublets.py` | Manual Leiden-cluster to cell-type mapping. |
| `hep_types` | `01_merge_doublets.py` | Hepatocyte classes eligible for same-type doublet fusion. |
| `fused_adata` | `01_merge_doublets.py` | Final fused AnnData object with cell-type, doublet, TMA, sex, and spatial metadata. |
| `celltype_colors` | `01_merge_doublets.py`, `02_compositional_analysis.py`, `03_DGE_analysis.py` | Shared plotting color palette for Xenium cell types. |
| `cell_types` | `03_DGE_analysis.py` | Cell types used for sex/doublet DGE loops. |

## Coherence Assessment

**Status: DOCUMENTED, collaborator-imported.** The folder is intentionally documented without code changes. It is not yet refactored into the main workflow style and should not be assumed to run from the repository root.

## Known Vulnerabilities

1. **Folder-local paths:** The scripts assume `data_raw/`, `data_tidy/`, and `figures/` relative to the current working directory.
2. **External raw data:** The raw 10x Xenium output folder and panel workbook are not committed.
3. **Notebook-style script state:** `03_DGE_analysis.py` contains repeated helper definitions and top-level code blocks inherited from collaborator analysis. It should be treated as imported analysis provenance rather than a fully normalized pipeline.
4. **Hardcoded spatial assignment:** TMA and sex labels are assigned from hardcoded `x_lines`/`y_lines` coordinate thresholds.
5. **Manual cell-type mapping:** Leiden cluster labels are mapped manually in `cell_type_dict`; changes to clustering resolution or raw inputs require revisiting this mapping.

## Notes for Context Restoration

1. Workflow 10 is independent of workflows 00-09 and does not consume their generated files.
2. Workflow 02 already includes a separate Xenium panel validation step using `support/FCGL_Xenium_Panel.csv`; workflow 10 is the collaborator spatial-analysis folder for raw Xenium output handling and DGE.
3. Keep code edits out of this folder unless explicitly requested; document assumptions and outputs instead.
4. The Xenium files provided in the 2026-06-19 missing-material zip were
   compared against this folder's `00`-`03` scripts and matched after ignoring
   CRLF line-ending differences, so no workflow 10 script files were recopied.
