# Workflow 02: FCG Pseudobulk Analysis

**Script:** `scripts/workflows/02_fcg_pseudobulk_analysis.R`
**Lines:** ~1370
**Purpose:** Differential expression analysis across FCG groups and cell types using limma, pathway/TF enrichment, GDX and Xenium validation.

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
| 1 | Data Loading | 1-80 | Load 3 CSVs, merge by gene, load metadata |
| 2 | PCA | 81-200 | PCA per cell type, variance by gonad/chromosome |
| 3 | DE Analysis (FCG) | 201-400 | limma for 6 contrasts x 5 cell types |
| 3B | DE Analysis (GDX) | 401-500 | GDX validation DE analysis |
| 4 | Volcano Plots | 501-600 | Volcano plots for all contrasts x cell types |
| 5 | T-statistic Matrix | 601-700 | Convert ttop lists to t-value matrix, save |
| 6 | Pathway Enrichment | 701-850 | decoupleR ULM with MSigDB GO-BP gene sets |
| 6B | Pathway Heatmaps | 851-950 | Full + focused heatmaps of pathway scores |
| 7C | GDX Validation | 951-1050 | Correlation between FCG and GDX t-values |
| 7D | Xenium Validation | 1051-1100 | Spatial transcriptomics validation (82 genes) |
| 8 | TF Activity | 1101-1180 | CollecTRI + ESR1/AR regulons, TF heatmap |
| 9 | Save PCA Coords | 1181-1210 | Save hepatocyte PCA coordinates for MOFA |
| 10 | Proteomics Validation | 1211-1340 | WT bulk liver proteomics (male vs female), VSN normalization, limma DE, correlation with FCG RNA, concordance analysis |

## Input Files

- `data/RNA/scaled_RNA_{1,2,3}.csv` - Expression matrices (from workflow 01)
- `data/RNA/meta_data_df.csv` - Sample metadata
- `support/m5.all.v2023.2.Mm.symbols.gmt` - MSigDB gene sets
- `support/TF_regulons.RData` - CollecTRI TF regulons
- `support/Gonadal_Hormone_Genes.csv` - ESR1/AR target genes
- `support/FCGL_Xenium_Panel.csv` - Xenium validation gene panel
- `data/proteomics/P3225_fractionated_protein.tsv` - TMT proteomics (WT male vs female)
- `data/proteomics/metadata.csv` - Proteomics sample metadata (12 TMT channels)

## Output Files

- `results/pseudobulk/ttop_full_list.RData` - All DE results (list of lists)
- `results/pseudobulk/tval_full_df.RData` - T-statistic matrix (binary)
- `results/pseudobulk/tval_full_df.csv` - T-statistic matrix (CSV, with ID column; used by Step 7C)
- `results/pseudobulk/volcanos/*.pdf` - Volcano plots (30 files: 6 contrasts x 5 cell types)
- `results/pseudobulk/pca/*.pdf` - PCA plots per cell type
- `results/pseudobulk/pathways/*.pdf` - Pathway heatmaps
- `results/pseudobulk/gdx_validation_comparison.pdf` - GDX vs FCG correlation
- `results/pseudobulk/xenium_validation_volcano.pdf` - Xenium validation
- `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` - Hepatocyte PCA coords
- `results/pseudobulk/proteo_validation_comparison.pdf` - Proteomics vs RNA correlation

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `scaled_RNA` | data.frame | Merged expression matrix (genes x samples) |
| `meta_data_df` | data.frame | Sample metadata |
| `targets` | data.frame | Design matrix for limma |
| `ttop_full_list` | list | Nested list: contrast -> cell_type -> toptable |
| `tval_full_df` | data.frame | T-statistics matrix |
| `comparisons` | named list | Contrast vectors (from `contrasts.R`) |
| `pathway_scores_df` | data.frame | decoupleR enrichment scores |

## Coherence Assessment

**Status: GOOD** — No open issues. The tval_full_df CSV/RData mismatch was fixed (Step 5 now exports both formats).

## Contrasts (from `scripts/config/contrasts.R`)

6 contrasts applied to each of 5 cell types:
1. `testes_vs_ova` (Testes vs Ovaries) - primary gonadal comparison
2. `XY_vs_XX` - sex chromosome comparison
3. `XYT_vs_XYO` - within XY: testes vs ovaries
4. `XYT_vs_XXT` - within Testes: XY vs XX
5. `XXT_vs_XXO` - within XX: testes vs ovaries
6. `XYO_vs_XXO` - within Ovaries: XY vs XX

## Notes for Context Restoration

When returning to this workflow:
1. This is the main DE analysis pipeline - most downstream workflows depend on its outputs
2. `ttop_full_list.RData` is the central result used by workflows 04, 05, 06, 07
3. Pathway enrichment uses `run_ulm()` from decoupleR with `minsize = 20`
4. Focused pathway heatmap uses score threshold |score| > 11-12 (may need tuning)
5. GDX validation expects correlation r ~ -0.48 with FCG (negative due to direction convention)
6. Cell types: Hepatocytes, Endothelial_cells, Immune_cells, Stellate_cells, Cholangiocytes
7. Proteomics validation uses WT mice (not FCG): male=testes, female=ovaries. VSN-normalized TMT data. ~33% of RNA DEGs confirmed at protein level.
