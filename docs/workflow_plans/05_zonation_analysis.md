# Workflow 05: Zonation Analysis

**Script:** `scripts/workflows/05_zonation_analysis.R`
**Lines:** ~2523 (LARGEST workflow)
**Purpose:** Hepatocyte zonation (periportal vs pericentral) analysis: DE by zone x gonad, sex-specific zonation gene identification via entropy-optimized threshold, correlation networks per zone x gonad category, lipid-zone association heatmaps.

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

Headers and console output are now aligned (1-16).

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 1-40 | Load config, set output dirs |
| 2 | Load Zonation Data | 41-93 | Load zonation pseudobulk, metadata |
| 3 | PCA | 95-132 | PCA by condition, gonad, zonation |
| 4 | PCA Enrichment | 134-183 | Test variance explained by metadata factors |
| 5 | DE Analysis | 185-237 | limma for 7 contrasts (incl. interaction) |
| 6 | Volcano Plots | 239-288 | Volcano plots for all contrasts |
| 7 | Pathway Enrichment | 290-374 | **RUN COMMENTED OUT** - loads from `zonation_limma_results.RData` |
| 8 | Lipidomics + TF | 376-439 | Load lipidomic data, compute TF scores, define sample groups |
| 9 | Zonation vs Sex Bias | 441-533 | Scatter: zonation t-value vs gonad t-value |
| 10 | Zone-Specific Sex Bias | 535-620 | PP vs PC t-values scatter (collaborator style) |
| 11 | Optimal Threshold | 622-774 | Entropy x stability x coverage x specificity optimization. Produces `sensitivity_df`, `proportion_total`, `optimal_metrics`, `optimal_threshold` |
| 12 | Identify Genes | 776-1151 | Sex-specific zonation gene selection + quadrant scatter + bar + heatmap |
| 13 | Correlation Networks | 1153-1407 | 4 networks (TPC, TPP, OPC, OPP) x 3 versions (sparse, full, compressed) |
| 14 | Lipid Association | 1409-2455 | Binary heatmap, class proportions, butterfly chart, PG example, jitter. Uses `sensitivity_df`, `proportion_total`, `optimal_metrics` from section 11 |
| 15 | Save Results | 2457-2483 | Save ttop_list, limmaRes, metadata, summary stats |
| 16 | Workflow Complete | 2485-2523 | Summary printout |

## Input Files

- `data/RNA/heps_zonationsubsetted_pseudobulks_gonad.rds` - Zonation-annotated hepatocyte pseudobulk (from workflow 00)
- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` - Lipidomic data
- `results/lipidomic/ttop_full_list.RData` - Lipidomic DE results (from workflow 03)
- `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData` - Enrichment scores (from workflow 03)
- `support/TF_regulons.RData` - TF regulons
- `support/Gonadal_Hormone_Genes.csv` - ESR1/AR targets
- `results/zonation/zonation_limma_results.RData` - Pre-computed pathway results (loaded in section 7)

## Output Files

**DE Results:**
- `results/zonation/zonation_ttop_list.RData` - All zonation DE results
- `results/zonation/zonation_limma_results.RData` - Full limma results
- `results/zonation/t_table_zonation.csv` - T-statistic table

**Figures:**
- `results/zonation/PCA_by_*.pdf` - PCA plots (3 files)
- `results/zonation/volcanos/*.pdf` - Volcano plots (7 files)
- `results/zonation/zonation_vs_sex_bias_scatter.pdf` - Scatter: zonation vs gonad
- `results/zonation/zone_specific_sex_bias_scatter.pdf` - PP vs PC scatter
- `results/zonation/sex_specific_zonation_panel_AB.pdf` - **MAIN FIGURE**: quadrant + bar
- `results/zonation/sex_specific_zonation_quadrant.pdf` - Full quadrant scatter
- `results/zonation/sex_specific_zonation_heatmap_faceted.pdf` - 137 genes heatmap
- `results/zonation/sex_specific_zonation_summary_bar.pdf` - Bar chart
- `results/zonation/optimal_threshold_detection.pdf` - Threshold optimization curves
- `results/zonation/lipid_zonation_threshold_sensitivity_total.pdf` - Sensitivity area plot
- `results/zonation/lipid_class_zonation_proportions.pdf` - Stacked bar proportions
- `results/zonation/lipid_class_zonation_proportions_horizontal.pdf` - Butterfly chart
- `results/zonation/lipid_class_zonation_proportions_PG_mini.pdf` - PG mini
- `results/zonation/gonad_to_zonation_explanation_PG.pdf` - PG explanatory figure
- `results/zonation/lipid_class_zonation_jitter_vertical.pdf` - Jitter plot

**Data:**
- `results/zonation/sex_specific_zonation_genes_summary.csv` - Gene list with stats (137 genes)
- `results/zonation/lipid_species_zonation_association.csv` - Lipid-level associations
- `results/zonation/lipid_category_binary_matrix.csv` - Binary lipid x category matrix
- `results/zonation/lipid_class_association_summary.csv` - Class-level summary
- `results/multi-omic/cormat_*_{SIF,ATT}.csv` - 4 x 3 = 12 network files
- `results/multi-omic/compressed_cormat_*_{SIF,ATT}.csv` - 4 compressed networks

## Key Variables

| Variable | Type | Defined in | Used in | Description |
|----------|------|-----------|---------|-------------|
| `t_table` | data.frame | Sec 5 | Sec 9-14 | T-statistics for all contrasts (genes x contrasts) |
| `optimal_threshold` | numeric | Sec 11 | Sec 12, 14 | Entropy-optimized threshold (~3.2) |
| `sensitivity_df` | data.frame | Sec 11 | Sec 14 | Lipid zonation at each threshold |
| `proportion_total` | data.frame | Sec 11 | Sec 14 | Total proportions by threshold |
| `optimal_metrics` | data.frame | Sec 11 | Sec 14 | Metrics for threshold selection |
| `testes_pericentral_genes` | data.frame | Sec 12 | Sec 13 | Selected TPC genes |
| `testes_periportal_genes` | data.frame | Sec 12 | Sec 13 | Selected TPP genes |
| `ovaries_pericentral_genes` | data.frame | Sec 12 | Sec 13 | Selected OPC genes |
| `ovaries_periportal_genes` | data.frame | Sec 12 | Sec 13 | Selected OPP genes |
| `scaled_RNA` | data.frame | Sec 2 | Sec 13-14 | Zonation-specific expression |
| `lipidomic` | data.frame | Sec 8 | Sec 13-14 | Lipidomic data |
| `TF_scores_df` | data.frame | Sec 8 | Sec 13 | TF activity scores |
| `jitter_data` | data.frame | Sec 14 | Sec 14 | Lipid jitter plot data |
| `classes_with_enrichment` | character | Sec 14 | Sec 14 | Ordered lipid classes |
| `lipid_class_results` | data.frame | Sec 14 (loaded) | Sec 14 | From workflow 03 enrichment |
| Sample group vectors | character | Sec 8 | Sec 13-14 | `testes_central_samples`, `ovaries_portal_samples`, etc. |

## Coherence Assessment

**Status: GOOD** — No open bugs. Section numbering and stale references have been fixed.

## Vulnerabilities to Keep in Mind When Modifying

1. **Pathway enrichment caching (section 7):** The `run_ulm()` call is commented out; section 7 loads pre-computed results from `pathway_enrichment.RData`. This is inside a `tryCatch` so it fails gracefully. If you need to re-run pathway enrichment (e.g. after changing contrasts), uncomment lines 334-347 and comment out line 350.

2. **Long-distance variable dependencies (section 11 -> 14):** `sensitivity_df`, `proportion_total`, and `optimal_metrics` are computed in section 11 but used ~700 lines later in section 14. The flow is correct but if you modify or skip section 11, section 14 will fail. Always run sections 11-14 together.

3. **Sample name regex in network generation (section 13):** The `gsub("_Peri.*_Testes", "", ...)` pattern strips zonation/gonad suffixes to match lipidomic sample names. If sample naming conventions change upstream, these regexes need updating too.

4. **`get_correlated_lipids()` is defined locally (section 11, line 633):** Not in an external file — it's a helper function within this workflow. If you refactor it out, update both the definition at line 633 and the calls in sections 11 and 14 (lines 670-676, 2251-2270).

## 7 Contrasts

1. `portal_vs_Central` - Overall zonation effect
2. `testes_vs_ova` - Overall gonad effect
3. `portalova_vs_centralova` - Zonation within ovaries
4. `portaltestes_vs_centraltestes` - Zonation within testes
5. `centraltestes_vs_centralova` - Gonad effect in pericentral
6. `portaltestes_vs_portalova` - Gonad effect in periportal
7. `testes_effect_zonation` - **INTERACTION**: (Testes PC-PP) vs (Ovaries PC-PP)

## Gene Selection Logic (Section 12)

```
optimal_threshold ≈ 3.2 (from entropy optimization in section 11)

Testes-specific zonation genes:
  |interaction_t| > optimal_threshold AND |zonation_ovaries_t| < 1
  Split by: portaltestes_vs_centraltestes (neg = pericentral, pos = periportal)

Ovaries-specific zonation genes:
  |interaction_t| > optimal_threshold AND |zonation_testes_t| < 1
  Split by: portalova_vs_centralova (neg = pericentral, pos = periportal)
```

**Result:** 62 TPC + 38 TPP + 18 OPC + 19 OPP = 137 genes total

## Notes for Context Restoration

When returning to this workflow:
1. This is the LARGEST and most complex workflow - 2523 lines with many interdependencies
2. The optimal threshold (~3.2) is the key parameter driving gene selection
3. Section 11 computes the threshold; section 12 uses it for gene selection; section 13 builds networks from those genes; section 14 maps lipids to network categories and visualizes threshold metrics
4. The 4 gene lists (TPC/TPP/OPC/OPP) are the biological foundation for the correlation networks
5. `make_cor_network()` and `expand_attributes()` come from `support_functions.R`
6. Lipid enrichment scores (`lipid_class_enrichment_scores.RData`) come from workflow 03
7. The butterfly chart and jitter plots reuse the enrichment-based lipid class ordering from workflow 03
8. Key genes: Cyp2b13, Sult3a1, Cyp2a22 (OPC); Stambpl1, Cyp2b9, Cyp2a4 (OPP); Esco2, Kitl, Mup17 (TPC); Grm8, Mtmr11, Mindy3 (TPP)
