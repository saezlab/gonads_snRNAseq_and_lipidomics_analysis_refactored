# Workflow 05: Zonation Analysis

**Script:** `scripts/workflows/05_zonation_analysis.R`  
**Lines:** ~2720  
**Purpose:** Analyze hepatocyte periportal/pericentral zonation, gonad-by-zonation interactions, sex-specific zonation-biased genes, zonation-specific gene-lipid networks, and lipid-class zonation associations.

## Source Dependencies

```r
source("scripts/config/analysis_config.R")
source("scripts/config/plotting_themes.R")
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
```

The script loads additional RData/CSV resources inside later sections.

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 41-68 | Load config/helpers and create `results/zonation/`, `volcanos/`, and `networks/`. |
| 2 | Load Zonation Pseudobulk | 70-99 | Load workflow 01 zonation CSVs and report sample format. |
| 3 | PCA Visualization | 102-138 | PCA by zonation, gonad, and combined condition. |
| 4 | PCA Enrichment | 141-189 | Correlate/test PCA coordinates against metadata and plot a heatmap. |
| 5 | Differential Analysis | 192-289 | Limma with `duplicateCorrelation()` blocked on mouse ID; 7 contrasts. |
| 6 | Volcano Plots | 292-340 | Volcano plot for each contrast. |
| 7 | Pathway Enrichment | 343-426 | Load GO gene sets, convert t-values, but active ULM call is commented out; loads cached `pathway_enrichment.RData`. |
| 8 | Lipidomics + TF Activity | 429-487 | Load lipidomics, TF regulons, ESR1/AR targets, compute TF activity scores, and write `t_table_zonation.csv`. |
| 9 | Zonation Bias vs Sex Bias | 490-610 | Scatter comparing overall zonation t-values and gonad t-values. |
| 10 | Zone-Specific Sex Bias | 613-697 | Scatter comparing ovary-zone and testes-zone t-values. |
| 10b | Zonation-Specific Sex Difference Quadrant | 700-814 | Additional quadrant plot for zone-specific gonad effects. |
| 11 | Optimal Threshold Detection | 817-968 | Sweep interaction thresholds 0-10, compute entropy/stability/coverage/specificity, and choose `optimal_threshold`. |
| 12 | Identify Zonation-Biased Genes | 971-1336 | Select sex-specific zonation genes, make quadrant/bar/panel/heatmap, and write gene summary. |
| 13 | Correlation Networks | 1339-1586 | Build four zonation x gonad networks, each with standard, full, and compressed CSV outputs. |
| 14 | Lipid-Gene Association Visualizations | 1589-2634 | Load WF03 lipid DE/enrichment and WF05 networks, map lipids to categories, export association/sensitivity tables, and plot class-level figures. |
| 15 | Save Results | 2637-2680 | Save limma objects, metadata, summary stats, and long-format zonation DE table. |
| 16 | Complete | 2683-2720 | Print outputs and expected biological findings. |

## Input Files

- `data/RNA/scaled_RNA_zonation.csv` - Zonation pseudobulk expression from workflow 01.
- `data/RNA/meta_data_df_zonation.csv` - Zonation metadata from workflow 01.
- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` - Loaded through `load_lipidomic_data()`.
- `support/m5.all.v2023.2.Mm.symbols.gmt` - GO biological process gene sets.
- `support/TF_regulons.RData` - CollecTRI regulons.
- `support/Gonadal_Hormone_Genes.csv` - Custom ESR1/AR targets.
- `results/zonation/pathway_enrichment.RData` - Cached pathway scores loaded in section 7 because active ULM run is commented out.
- `results/lipidomic/ttop_full_list.RData` - Workflow 03 lipid DE results loaded in section 14.
- `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData` - Workflow 03 lipid-class enrichment scores loaded in section 14.

## Output Files

**DE and tables**
- `results/zonation/zonation_ttop_list.RData`
- `results/zonation/zonation_limma_results.RData`
- `results/zonation/t_table_zonation.csv`
- `results/zonation/ttop_full_long.csv`
- `results/zonation/metadata.csv`
- `results/zonation/summary_statistics.csv`

**Gene-selection outputs**
- `results/zonation/zonation_specific_sex_differences_quadrant.pdf`
- `results/zonation/zonation_biased_genes_heatmap.pdf`
- `results/zonation/sex_specific_zonation_quadrant.pdf`
- `results/zonation/sex_specific_zonation_summary_bar.pdf`
- `results/zonation/sex_specific_zonation_panel_AB.pdf`
- `results/zonation/sex_specific_zonation_heatmap_faceted.pdf`
- `results/zonation/sex_specific_zonation_genes_summary.csv`

**Network outputs**
- `results/zonation/networks/cormat_testes_portal_{SIF,ATT}.csv`
- `results/zonation/networks/cormat_testes_central_{SIF,ATT}.csv`
- `results/zonation/networks/cormat_Ovaries_portal_{SIF,ATT}.csv`
- `results/zonation/networks/cormat_ova_central_{SIF,ATT}.csv`
- `results/zonation/networks/cormat_*_full_{SIF,ATT}.csv` for the same four contexts.
- `results/zonation/networks/compressed_cormat_*_{SIF,ATT}.csv` for the same four contexts.

**Lipid association and threshold outputs**
- `results/zonation/lipid_species_zonation_association.csv`
- `results/zonation/lipid_category_binary_matrix.csv`
- `results/zonation/lipid_class_association_summary.csv`
- `results/zonation/threshold_sensitivity_by_class.csv`
- `results/zonation/threshold_sensitivity_total.csv`
- `results/zonation/optimal_threshold_metrics.csv`
- `results/zonation/optimal_threshold_detection.pdf`
- `results/zonation/lipid_zonation_threshold_sensitivity_total.pdf`
- `results/zonation/lipid_class_zonation_proportions.pdf`
- `results/zonation/lipid_class_zonation_proportions_horizontal.pdf`
- `results/zonation/lipid_class_zonation_proportions_PG_mini.pdf`
- `results/zonation/gonad_to_zonation_explanation_PG.pdf`
- `results/zonation/lipid_class_zonation_jitter_vertical.pdf`

**Other figures**
- `results/zonation/PCA_by_zonation.pdf`
- `results/zonation/PCA_by_gonad.pdf`
- `results/zonation/PCA_by_condition.pdf`
- `results/zonation/PCA_metadata_enrichment.pdf`
- `results/zonation/zonation_vs_sex_bias_scatter.pdf`
- `results/zonation/zone_specific_sex_bias_scatter.pdf`
- `results/zonation/volcanos/volcano_*.pdf`

## Key Variables

| Variable | Type | Defined in | Used in | Description |
|----------|------|------------|---------|-------------|
| `scaled_RNA` | data.frame | Sec 2 | Sec 3-14 | Zonation pseudobulk matrix |
| `meta_data_df` | data.frame | Sec 2 | Sec 3-15 | Zonation metadata |
| `comparisons` | list | Sec 5 | Sec 5-15 | Seven zonation/gonad contrasts |
| `mouse_id` | character | Sec 5 | Sec 5 | Blocking factor for `duplicateCorrelation()` |
| `ttop_list` | list | Sec 5 | Sec 6-15 | Limma toptables by contrast |
| `t_table` | data.frame | Sec 8 | Sec 9-14 | Zonation t-statistic table with `gene` column |
| `TF_scores_df` | data.frame | Sec 8 | Sec 13 | TF activity matrix |
| `lipidomic` | data.frame | Sec 8 | Sec 11-14 | Lipidomic matrix |
| `optimal_threshold` | numeric | Sec 11 | Sec 12-14 | Interaction t-statistic cutoff selected by combined score |
| `sensitivity_df` | data.frame | Sec 11 | Sec 14 | Lipid zonation assignment at each threshold |
| `proportion_total` | data.frame | Sec 11 | Sec 14 | Threshold-level zonation proportions |
| `optimal_metrics` | data.frame | Sec 11 | Sec 14 | Entropy/stability/coverage/specificity metrics |
| `testes_pericentral_genes` / `testes_periportal_genes` | data.frame | Sec 12 | Sec 13-14 | Testes-specific zonation gene sets |
| `ovaries_pericentral_genes` / `ovaries_periportal_genes` | data.frame | Sec 12 | Sec 13-14 | Ovaries-specific zonation gene sets |
| `classes_with_enrichment` | character | Sec 14 | Sec 14 | Lipid classes ordered by workflow 03 enrichment score |

## Coherence Assessment

**Status: GOOD.** The current script is coherent and the plan now reflects the CSV-based zonation input, duplicate-correlation DE model, threshold-selection section, and workflow 03 lipid dependencies.

## Vulnerabilities to Keep in Mind When Modifying

1. **Pathway enrichment cache:** Section 7 has the `decoupleR::run_ulm()` call and save line commented out, then loads `results/zonation/pathway_enrichment.RData`. If that file is absent, the `tryCatch` skips pathway reporting.
2. **Long-distance threshold dependency:** `sensitivity_df`, `proportion_total`, and `optimal_metrics` are computed in section 11 and consumed in section 14 more than 900 lines later.
3. **Sample-name regex:** Network sections strip `_Peri.*_Testes` and `_Peri.*_Ovaries` suffixes to match lipidomic sample names. Upstream sample-format changes will break RNA-lipid matching.
4. **WF05 -> WF07 filenames:** Workflow 07 expects the exact `cormat_testes_portal_full_SIF.csv`, `cormat_testes_central_full_SIF.csv`, `cormat_Ovaries_portal_full_SIF.csv`, and `cormat_ova_central_full_SIF.csv` names.
5. **Locally defined helper:** `get_correlated_lipids()` is defined inside section 11 and reused in section 14; it is not imported from a helper module.

## Seven Contrasts

1. `portal_vs_Central` - Overall zonation effect.
2. `testes_vs_ova` - Overall gonad effect.
3. `portalova_vs_centralova` - Zonation effect within ovaries.
4. `portaltestes_vs_centraltestes` - Zonation effect within testes.
5. `centraltestes_vs_centralova` - Gonad effect in pericentral hepatocytes.
6. `portaltestes_vs_portalova` - Gonad effect in periportal hepatocytes.
7. `testes_effect_zonation` - Interaction: whether zonation differs by gonad.

## Gene Selection Logic

```r
testes_zonation_bias_genes:
  abs(testes_effect_zonation) > optimal_threshold
  abs(portalova_vs_centralova) < 1

ovaries_zonation_bias_genes:
  abs(testes_effect_zonation) > optimal_threshold
  abs(portaltestes_vs_centraltestes) < 1

general_zonation_bias_genes:
  abs(testes_effect_zonation) < 1
  abs(portal_vs_Central) > optimal_threshold
```

Zone assignment is then based on the sign of `portaltestes_vs_centraltestes` or `portalova_vs_centralova`: negative means pericentral and positive means periportal in the current script.

## Notes for Context Restoration

1. Workflow 05 is the largest workflow and should be read in sections, not skimmed.
2. The active input is the workflow 01 CSV pair, not a Seurat RDS.
3. Section 5 deliberately blocks on mouse ID because each mouse contributes both zonation pseudobulks.
4. Section 14 depends on workflow 03's lipid DE and lipid-class enrichment outputs.
