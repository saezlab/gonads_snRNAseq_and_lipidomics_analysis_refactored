# Cross-Workflow Dependencies and Coherence Report

**Last audited:** 2026-05-27
**Ground truth:** `scripts/workflows/*.R` plus `scripts/workflows/09_human_analysis.py`

This meta overview is synchronized with the current workflow scripts and the
per-workflow plan files in this directory.

## Dependency Graph

```text
00_snrnaseq_preprocessing
  Reference/HPC workflow: raw sci-RNA-seq -> Seurat RDS files
  |
  v
01_data_preparation_qc
  Seurat RDS -> main RNA CSVs, metadata, zonation CSVs, GDX RNA CSVs, QC
  |
  +---> 02_fcg_pseudobulk_analysis
  |       RNA DE, GDX RNA validation, pathway/TF scores, Xenium, proteomics
  |
  +---> 03_fcg_lipidomic_analysis
  |       lipid DE, GDX lipid validation, lipid-class enrichment, LINEX2 network
  |       |
  |       +---> 05_zonation_analysis
  |       |       uses lipid DE and lipid-class enrichment
  |       |
  |       +---> 06_structural_analysis
  |               lipid-only structural analysis from lipid DE
  |
  +---> 04_mofa_integration
  |       RNA + lipidomics -> MOFA weights and gene-lipid/TF-gene-lipid networks
  |       |
  |       +---> 07_validation_experiments
  |               validates global MOFA network predictions
  |
  +---> 05_zonation_analysis
          uses zonation RNA CSVs and writes zonation-specific networks
          |
          +---> 07_validation_experiments
                  validates zonation network predictions

06_structural_analysis
  |
  +---> 07_validation_experiments
          uses structural lipid feature table for KO consistency plots

08_variance_analysis
  standalone; uses raw lipidomic data plus metadata only

09_human_analysis
  optional external-data Python workflow; uses the mouse hepatocyte gonad
  signature plus external GSE179250 human liver count/TPM/SOFT files
```

Important current-script clarifications:

- WF02 DE outputs are not currently consumed by WF04 or WF05. Optional WF09
  uses a checked-in hepatocyte gonad signature derived from the WF02-style
  Testes vs Ovaries hepatocyte t-values.
- WF05 does not load WF04 outputs; it uses WF01 zonation CSVs, WF03 lipid
  DE/enrichment outputs, support regulons, and raw lipidomics.
- WF06 currently runs a lipid-only structural analysis. Older gene-lipid
  structural sections are commented out.
- WF07 section 11 depends on WF06 `lipids_with_structure.csv`.

## Key File Dependencies

| Producing workflow | File | Consuming workflow(s) |
| --- | --- | --- |
| 00 | `data/RNA/fcg_pseudobulks.rds` | 01 |
| 00 | `data/RNA/fcg_zonation_pseudobulks.rds` | 01 |
| 00 | `data/RNA/gdx_pseudobulks.rds` | 01 |
| 01 | `data/RNA/scaled_RNA_{1,2,3}.csv` | 02, 04 |
| 01 | `data/RNA/meta_data_df.csv` | 02, 03 |
| 01 | `data/RNA/scaled_RNA_zonation.csv` | 05 |
| 01 | `data/RNA/meta_data_df_zonation.csv` | 05 |
| 01 | `data/RNA/scaled_RNA_gdxhr_{1,2,3}.csv` | 02 |
| 01 | `data/RNA/meta_data_df_gdxhr.csv` | 02 |
| 02 | `results/pseudobulk/ttop_full_list.RData` | Result lookup; no current scripted downstream workflow |
| 02 | `results/pseudobulk/tval_full_df.csv` | 02 internal GDX/proteomics validation and result lookup |
| 02-derived/manual export | `results/pseudobulk/hepa_gonad_signature.csv` | 09, as the mouse hepatocyte gonad signature seed; script currently expects a copy in external `DATA_DIR` |
| 02 | `results/pseudobulk/pathways/pathway_scores_df_all.RData` | Result lookup |
| 02 | `results/pseudobulk/pathways/tf_scores_df_all.RData` | Result lookup |
| 03 | `results/lipidomic/ttop_full_list.RData` | 05, 06 |
| 03 | `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData` | 05 |
| 03 | `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv` | Result lookup |
| 03 | `results/lipidomic/network/significant_lipids_network_{SIF,ATT}.csv` | Cytoscape/result lookup |
| 04 | `results/mofa/main_model_RNA_lipid_weight.RData` | Result lookup |
| 04 | `results/mofa/networks/gene_lipid_correlations_factor1.csv` | Result lookup |
| 04 | `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | 07 |
| 04 | `results/mofa/networks/hepatocyte_TF_gene_lipid_network_*_{SIF,ATT}.csv` | Cytoscape/result lookup |
| 05 | `results/zonation/pathway_enrichment.RData` | 05 internal cached pathway load |
| 05 | `results/zonation/zonation_ttop_list.RData` | Zonation result archive |
| 05 | `results/zonation/t_table_zonation.csv` | Result lookup |
| 05 | `results/zonation/ttop_full_long.csv` | Result lookup |
| 05 | `results/zonation/sex_specific_zonation_genes_summary.csv` | Result lookup |
| 05 | `results/zonation/lipid_class_association_summary.csv` | Result lookup |
| 05 | `results/zonation/lipid_species_zonation_association.csv` | Result lookup |
| 05 | `results/zonation/threshold_sensitivity_total.csv` | Result lookup |
| 05 | `results/zonation/optimal_threshold_metrics.csv` | Result lookup |
| 05 | `results/zonation/networks/cormat_*_full_SIF.csv` | 07 |
| 06 | `results/lipidomic/structural_model/lipids_with_structure.csv` | 07 |
| 06 | `results/lipidomic/structural_model/lipid_only_*effects.csv` | Result lookup, when files exist |
| 06 | `results/lipidomic/structural_model/lipid_only_class_{coefficients,summaries}.csv` | Result lookup |
| 07 | `results/validation/validation_summary.csv` | Result lookup |
| 07 | `results/validation/global_network_validation_scores.csv` | Result lookup |
| 07 | `results/validation/zonation_network_validation_scores.csv` | Result lookup |
| 07 | `results/validation/structural_consistency/consistency_structural_regression_summary.csv` | Result lookup |
| 08 | `results/lipidomic/variance_analysis/variance_analysis_summary.csv` | Result lookup |
| 08 | `results/lipidomic/variance_analysis/variance_analysis_results.csv` | Result lookup |
| 08 | `results/lipidomic/variance_analysis/variance_analysis_by_class.csv` | Result lookup |
| External GEO `GSE179250` | Human liver count matrix, TPM matrix, and SOFT metadata | 09 |
| External copied mouse table | `DATA_DIR/male_vs_female_DA_res.csv` | 09 concordant ranked-gene dot plots |

## Coherence Status

| Workflow | Script status | Plan-file status | Notes |
| --- | --- | --- | --- |
| 00 | GOOD | GOOD | Reference workflow in `if (FALSE)` blocks; requires original HPC context. |
| 01 | GOOD | GOOD | Produces main RNA, zonation, and GDX RNA CSVs. |
| 02 | GOOD | GOOD | Main RNA DE workflow; current script does not write RNA `ttop_full_long.csv`. |
| 03 | GOOD | GOOD | Central lipid DE workflow; feeds WF05 and current WF06. |
| 04 | GOOD | GOOD | Independent of WF02 DE outputs; feeds WF07 global network validation. |
| 05 | GOOD | GOOD | Uses workflow 01 zonation CSVs, workflow 03 lipid outputs, and cached pathway enrichment. |
| 06 | GOOD | GOOD | Active path is lipid-only structural analysis under `results/lipidomic/structural_model/`. |
| 07 | GOOD | GOOD | Includes structural consistency analysis dependent on WF06. |
| 08 | GOOD | GOOD | Standalone lipid variance workflow with current `variance_analysis_*` outputs. |
| 09 | GOOD | GOOD | Optional Python human comparison; depends on external GSE179250 files and writes outside the repo by default. |

## Cross-Workflow Vulnerabilities

1. **WF05 pathway cache:** Section 7 loads `results/zonation/pathway_enrichment.RData`; the active ULM run is commented out. A fresh environment without this file skips pathway reporting.
2. **WF05 threshold variables:** `sensitivity_df`, `proportion_total`, and `optimal_metrics` are created in section 11 and consumed much later in section 14.
3. **WF05 -> WF07 network names:** WF07 expects exact mixed names from WF05: `cormat_testes_portal_full_SIF.csv`, `cormat_testes_central_full_SIF.csv`, `cormat_Ovaries_portal_full_SIF.csv`, and `cormat_ova_central_full_SIF.csv`.
4. **WF06 -> WF07 structural feature table:** WF07 section 11 expects `results/lipidomic/structural_model/lipids_with_structure.csv`.
5. **WF03 quantitative count drift:** The WF03 script header and final printed text contain different headline lipid counts. Use generated result files for exact quantitative lookups.
6. **WF09 external-data paths:** The script uses absolute collaborator/HPC `DATA_DIR` and `RESULTS_DIR` values. The large human liver inputs are not committed, and `results/pseudobulk/hepa_gonad_signature.csv` must be copied to the expected external `DATA_DIR` path before running without code changes.

## Shared Utility Functions

`scripts/support_functions.R` is a compatibility entry point that sources
`scripts/config/support_functions.R`. Workflows and plan files use the documented
`scripts/support_functions.R` path.

| Function | Used in | Purpose |
| --- | --- | --- |
| `runLimma()` | 02, 03, 05, 07 | Wrapper for limma differential analysis |
| `limma_res_to_ttop_list()` | 02, 03, 05, 07 | Convert limma results to topTable lists |
| `ttop_list_to_t_table()` | 02, 03, 05, 07 | Convert topTable list to a t-value matrix/table |
| `volcano_nice()` | 02, 03, 05 | Publication-style volcano plots |
| `nicePCA()` | 01, 02, 03, 05, 07 | PCA plotting helper |
| `make_cor_network()` | 04, 05 | Build correlation network from two matrices |
| `expand_attributes()` | 04, 05 | Add DE statistics to network nodes |
| `get_correlated_lipids()` | 05 local | Get lipids correlated with a gene set |
| `generate_lipid_structural_data()` | 06 | Parse lipid names and append logFC |
| `parse_lipid_names_to_df()` | 07 | Parse lipid structures for KO consistency analysis |

## Recommended Pre-Run Checks

Before running any workflow, verify:

1. Working directory is the project root.
2. Required upstream files from the dependency table exist.
3. R packages are installed and loadable.
4. No `View()` calls remain in the code path you will run.
5. Edited tidyverse calls use explicit namespaces where conflicts are likely.
6. The relevant workflow plan still matches the script if you have changed code.

## Agent Documentation Entry Point

Use `docs/AGENT_GUIDE.md` as the static replacement for the discontinued MCP
query layer. It maps common biological/data questions to result files that
answer them directly.

Keep this meta overview and `docs/AGENT_GUIDE.md` synchronized whenever a
workflow starts writing a new machine-readable result or changes a dependency.
