# Workflow 03: FCG Lipidomic Analysis

**Script:** `scripts/workflows/03_fcg_lipidomic_analysis.R`
**Lines:** ~938
**Purpose:** Differential abundance analysis of 1,386 lipids across FCG groups, lipid class enrichment, LINEX2 network analysis, bubble plots.

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
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Data Loading | 1-80 | Load lipidomic CSV, log2 transform, load metadata |
| 2 | PCA | 81-180 | PCA colored by gonad type and sex chromosome |
| 3 | DE Analysis (FCG) | 181-350 | limma for 6 contrasts on lipids |
| 3.5 | DE Analysis (GDX) | 351-450 | GDX lipidomic DE for validation |
| 4 | Volcano Plots | 451-530 | Volcano plots for all contrasts |
| 4.5 | GDX Correlation | 531-600 | FCG vs GDX t-value correlation plot |
| 5 | Lipid Class Enrichment | 601-700 | decoupleR ULM on lipid class membership network |
| 5B | Enrichment Plot | 701-750 | Jitter plot of lipid t-values by class |
| 6 | LINEX2 Network | 751-850 | Load adjacency matrix, extract significant lipid subnetwork |
| 7 | Bubble Plots | 851-920 | Bubble plots showing lipid class x structural features |
| 8 | Summary | 921-938 | Print summary statistics |

## Input Files

- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` - Raw lipidomic data (1,386 lipids)
- `data/Suppl_Table_1_OverviewData.xlsx` - Sample metadata
- `support/linex2_adjacency_matrix_OT.csv` - LINEX2 lipid reaction network
- `support/lipid_classification.csv` - Lipid class hierarchy

## Output Files

- `results/lipidomic/ttop_full_list.RData` - All lipidomic DE results
- `results/lipidomic/gdxhr_ttop_list.RData` - GDX DE results
- `results/lipidomic/volcanos/*.pdf` - Volcano plots (6 files)
- `results/lipidomic/pca/*.pdf` - PCA plots
- `results/lipidomic/fcg_gdx_correlation.pdf` - FCG vs GDX correlation
- `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData` - Enrichment scores
- `results/lipidomic/enrichment/lipid_class_enrichment_plot.pdf` - Jitter enrichment plot
- `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv` - Enrichment table

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `lipidomic` | data.frame | Log2-transformed lipid abundances (lipids x samples) |
| `ttop_full_list` | list | DE results per contrast |
| `lipid_class_results` | data.frame | decoupleR enrichment scores per lipid class |
| `lipid_enrichment_gonade` | data.frame | Gonadal enrichment with direction |

## Coherence Assessment

**Status: GOOD** - Linear pipeline with clear data flow.

- Enrichment threshold of 1.8 is used for lipid class significance (same threshold reused in workflow 05)
- LINEX2 network section is self-contained
- `lipid_class_enrichment_scores.RData` is consumed by workflow 05 (section 14)

## Key Constants

- Enrichment threshold: `|score| > 1.8` for lipid class significance
- DE threshold: FDR < 0.05 for significant lipids
- 733/1386 lipids significant for gonadal type
- 0 lipids significant for sex chromosome complement

## Notes for Context Restoration

When returning to this workflow:
1. `ttop_full_list.RData` naming differs from workflow 02 - here it's a flat list (contrasts only, no cell types)
2. `lipid_class_enrichment_scores.RData` is critical for workflow 05's lipid-zonation analysis
3. Lipid class extraction: `gsub("[(].*", "", lipid_name)` strips everything after first parenthesis
4. Male-biased: TG, DG, glycerophospholipids; Female-biased: sphingolipids (SM, Cer, HexCer)
5. GDX correlation direction may be negative due to contrast definition conventions
