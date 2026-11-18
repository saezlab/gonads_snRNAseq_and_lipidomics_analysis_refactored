# Workflow 03: FCG Lipidomic Analysis

**Script:** `scripts/workflows/03_fcg_lipidomic_analysis.R`  
**Lines:** ~953  
**Purpose:** Run lipid differential abundance across FCG contrasts, optional GDX lipid validation, lipid-class enrichment, LINEX2 network extraction, and structural bubble plots.

## Current execution contract (audited 18 November 2026)

Both `source()` and `Rscript` execute the complete body. There is no CLI parser,
validate-only, source guard or output-root option; named outputs under
`results/lipidomic/` can be overwritten. Use an isolated checkout/copy for a fresh
reproduction. The script uses literal/default loader paths, not `DATASETS` routing.

The main branch needs the supplied cleaned lipid CSV, overview workbook and main
RNA metadata roster, but no Workflow02 result. It explicitly intersects lipid
samples with the RNA mouse roster before fitting. The supplied roster retains43
FCG mice; Workflow08 instead uses all48 measured lipid mice and includes the five
WT males. Treat regenerated tables as authoritative if inputs change.

The optional GDX lipid branch reads `data/metabolomics/lipidomics_gdxhr.csv`
directly; it does not need Workflow01's GDX RNA CSVs or Workflow02 GDX results.
LINEX2 adjacency and supplied bubble-plot CSVs are separate curated resources.
Outputs supply lipid DE/enrichment to05 and lipid DE to06. Public-clone execution
requires these files and the R packages, without any source RNA/Seurat object.

## Source Dependencies

```r
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/core/differential_analysis.R")
source("scripts/core/enrichment_analysis.R")
source("scripts/config/contrasts.R")
source("scripts/config/plotting_themes.R")
source("scripts/visualization/volcano_plots.R")
source("scripts/visualization/dimred_plots.R")
source("scripts/visualization/heatmaps.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Load and Prepare Data | 77-114 | Load log2 lipidomics and metadata, then filter lipidomics to RNA-shared samples using `meta_data_df.csv`. |
| 2 | PCA Visualization | 116-151 | Lipid PCA colored by gonad/allosome. |
| 3 | FCG Differential Analysis | 154-235 | Six FCG contrasts on lipids. |
| 3.5 | GDX Lipid Differential Analysis | 238-356 | Optional GDX lipidomic validation from `lipidomics_gdxhr.csv`. |
| 4 | Volcano Plots | 358-410 | Volcano plots for lipid contrasts. |
| 4.5 | FCG-GDX Correlation | 412-505 | Optional correlation between FCG and GDX lipid t-values. |
| 5 | Lipid Class Enrichment | 507-568 | Build t-stat matrix and run class enrichment with `minsize = 5`. |
| 6 | Enrichment Plot | 570-678 | Jitter/crossbar lipid-class plot for the gonadal contrast. |
| 7 | LINEX2 Network | 681-778 | Convert LINEX2 adjacency to edge list and save significant-lipid subnetwork. |
| 7.5 | Bubble Plots | 781-864 | Plot chain length (C) vs DB for every CSV in `data/metabolomics/for_bubbles/`. |
| 8 | Summary Statistics | 867-953 | Print key lipidomic findings and output list. |

## Input Files

- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv`
- `data/Suppl_Table_1_OverviewData.xlsx` via `load_metadata()`
- `data/RNA/meta_data_df.csv` - Used to keep lipid samples shared with RNA.
- `data/metabolomics/lipidomics_gdxhr.csv` - Optional GDX lipid validation input.
- `support/linex2_adjacency_matrix_OT.csv`
- `data/metabolomics/for_bubbles/*.csv`

## Output Files

- `results/lipidomic/ttop_full_list.RData`
- `results/lipidomic/gdxhr_t_table.csv`
- `results/lipidomic/gdxhr_ttop_list.RData`
- `results/lipidomic/volcanos/*.pdf`
- `results/lipidomic/pca/*.pdf`
- `results/lipidomic/fcg_gdx_correlation.pdf`
- `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData`
- `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv`
- `results/lipidomic/enrichment/lipid_class_enrichment_plot.pdf`
- `results/lipidomic/network/significant_lipids_network_SIF.csv`
- `results/lipidomic/network/significant_lipids_network_ATT.csv`
- `results/lipidomic/bubbles_DBC/*.pdf`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `lipidomic` | data.frame | Log2 lipid abundance matrix, filtered to RNA-shared samples |
| `metadata_lipid` | data.frame | Metadata for lipidomic samples retained in the analysis |
| `ttop_full_list` | list | Nested lipid limma result list by contrast |
| `tval_matrix` | matrix | Lipid t-statistic matrix used for enrichment |
| `lipid_class_results` | data.frame | decoupleR lipid-class enrichment scores |
| `lipid_enrichment_gonade` | data.frame | Gonadal enrichment scores used for plotting |

## Coherence Assessment

**Status: GOOD.** The workflow is linear after optional GDX input checks. It feeds workflow 05 with lipid DE/enrichment outputs and workflow 06 with `ttop_full_list.RData`.

## Key Constants and Caveats

- Lipid class enrichment uses `minsize = 5`.
- The enrichment plot marks classes with `abs(enrichment) > 1.8`.
- The loader reads the supplied cleaned lipid matrix, then the entry filters samples to the RNA-shared cohort; use generated results for feature counts.
- The script's header and final printed count differ (`733/1,386` vs `773/1,386`); generated result files should be treated as the quantitative source for exact counts.

## Notes for Context Restoration

1. WF03's `ttop_full_list.RData` is nested (`contrast -> comparison`) and is handled tolerantly in workflow 05.
2. LINEX2 outputs are in `results/lipidomic/network/`; older root-level network CSVs may exist but are not written by the current script.
3. Bubble plots use `data/metabolomics/for_bubbles/`, while workflow 06 parses lipid names directly from limma results.
