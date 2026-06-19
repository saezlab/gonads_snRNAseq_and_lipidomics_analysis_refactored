# Workflow 02: FCG Pseudobulk Analysis

**Script:** `scripts/workflows/02_fcg_pseudobulk_analysis.R`  
**Lines:** ~1452  
**Purpose:** Run FCG pseudobulk RNA differential expression across five liver cell types, enrich pathways/TFs, and perform GDX, Xenium, and proteomics validation.

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
| 1 | Load and Prepare Data | 77-97 | Load workflow 01 RNA CSVs and merged metadata. |
| 2 | PCA Visualization | 100-220 | PCA per cell type plus PC1 ridgeline plots by gonad/allosome. |
| 3 | FCG Differential Expression | 223-326 | Six FCG contrasts across five cell types, saved as nested ttop list. |
| 3B | GDX Differential Analysis | 329-475 | Optional GDX RNA validation if workflow 01 GDX CSVs exist. |
| 4 | Volcano Plots | 477-542 | Volcano plots for each FCG contrast and cell type. |
| 5 | T-statistic Matrix | 545-606 | Build and save `tval_full_df.RData` and `tval_full_df.csv`. |
| 6 | Pathway Enrichment | 609-636 | GO-BP ULM enrichment with `minsize = 5`. |
| 7 | Pathway Heatmap | 639-692 | Full pathway heatmap. |
| 7B | Focused Pathway Heatmap | 694-846 | Key pathway heatmap for gonadal contrast. |
| 7C | GDX Validation Comparison | 848-982 | Correlate FCG and GDX t-values. |
| 7D | Xenium Validation | 984-1078 | Optional Xenium panel volcano plot. |
| 8 | TF Activity | 1090-1164 | CollecTRI/custom hormone receptor TF activity with `minsize = 20`. |
| 9 | PCA Coordinates | 1166-1198 | Save hepatocyte PCA coordinates. |
| 10 | Proteomics Validation | 1200-1399 | Optional WT bulk liver proteomics validation with VSN and limma. |
| - | Complete | 1401-1452 | Print output summary. |

## Input Files

- `data/RNA/scaled_RNA_{1,2,3}.csv` - Main RNA expression from workflow 01.
- `data/RNA/meta_data_df.csv` - Main RNA metadata from workflow 01.
- `data/RNA/scaled_RNA_gdxhr_{1,2,3}.csv` and `data/RNA/meta_data_df_gdxhr.csv` - Optional GDX RNA validation inputs.
- `support/m5.all.v2023.2.Mm.symbols.gmt` - GO-BP gene sets.
- `support/TF_regulons.RData` - CollecTRI regulons.
- `support/Gonadal_Hormone_Genes.csv` - Custom ESR1/AR targets, used by helper functions.
- `support/FCGL_Xenium_Panel.csv` - Xenium validation panel.
- `data/proteomics/P3225_fractionated_protein.tsv`
- `data/proteomics/metadata.csv`

## Output Files

- `results/pseudobulk/ttop_full_list.RData`
- `results/pseudobulk/tval_full_df.RData`
- `results/pseudobulk/tval_full_df.csv`
- `results/pseudobulk/gdxrh/gdhxr_groups_ttop_list.RData`
- `results/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv`
- `results/pseudobulk/volcanos/volcano_*.pdf`
- `results/pseudobulk/pca/*.pdf`
- `results/pseudobulk/pathways/pathway_scores_df_all.RData`
- `results/pseudobulk/pathways/pathway_enrichment_heatmap.pdf`
- `results/pseudobulk/pathways/pathway_enrichment_heatmap_focused.pdf`
- `results/pseudobulk/pathways/tf_scores_df_all.RData`
- `results/pseudobulk/pathways/tf_activity_heatmap.pdf`
- `results/pseudobulk/gdx_validation_comparison.pdf`
- `results/pseudobulk/xenium_validation_volcano.pdf`
- `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv`
- `results/pseudobulk/proteo_validation_comparison.pdf`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `scaled_RNA` | data.frame | Main genes x pseudobulk samples matrix |
| `metadata` | data.frame | RNA metadata merged with overview workbook |
| `cell_types` | character | Hepatocytes, Endothelial cells, Immune cells, Stellate cells, Cholangiocytes |
| `contrasts` | character | Six FCG contrast config names |
| `ttop_full_list` | list | Nested list: contrast display name -> cell type -> topTable |
| `tval_full_df` | data.frame | Wide t-statistic matrix |
| `pathway_results` | data.frame | decoupleR GO enrichment scores |
| `tf_results` | data.frame | decoupleR TF activity scores |

## Coherence Assessment

**Status: GOOD.** The workflow's products are used mainly for RNA result lookup and validation figures. Current workflow 04 and 05 scripts do not consume workflow 02 DE outputs.

## Contrasts

1. `gonadal_type` -> `Testes vs Ovaries`
2. `allosome_type` -> `XY vs XX`
3. `XXT_vs_XXO`
4. `XYT_vs_XYO`
5. `XYO_vs_XXO`
6. `XYT_vs_XXT`

## Notes for Context Restoration

1. `results/pseudobulk/ttop_full_list.RData` is the primary RNA DE archive. The current script does not write `ttop_full_long.csv`.
2. Pathway enrichment intentionally uses `minsize = 5` to match the legacy/default decoupleR behavior.
3. TF activity uses `minsize = 20`.
4. GDX validation expects the workflow 01 GDX CSVs and writes into `results/pseudobulk/gdxrh/` (spelled `gdxrh` in the script).
5. Proteomics validation is optional and skipped if either proteomics input file is missing.
6. The original collaborator P3225 fractionated TMT proteomics R Markdown file
   is retained as provenance at
   `scripts/james/proteomics/JamesCleland_P3225_V1.Rmd`. Workflow 02 does not
   source it; it uses the curated checked-in `data/proteomics/` files listed
   above.
