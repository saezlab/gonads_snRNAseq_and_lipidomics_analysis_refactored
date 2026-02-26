# Cross-Workflow Dependencies & Coherence Report

## Dependency Graph

```
00_preprocessing (raw data -> Seurat)  [if(FALSE) — requires HPC]
  |
  v
01_data_preparation_qc (Seurat -> 3 CSVs + metadata + zonation CSVs)
  |
  +---> 02_fcg_pseudobulk_analysis (RNA DE, pathways, TFs)
  |       |
  |       +---> 05_zonation_analysis (uses ttop patterns, TF regulons)
  |       +---> 07_validation_experiments (uses DE results pattern)
  |
  +---> 03_fcg_lipidomic_analysis (lipid DE, enrichment)
  |       |
  |       +---> 05_zonation_analysis (uses lipid ttop, enrichment scores)
  |       +---> 08_variance_analysis (independent, same raw data)
  |
  +---> 04_mofa_integration (RNA + lipid integration)
          |
          +---> 05_zonation_analysis (uses MOFA weights for TF scores)
          +---> 06_structural_analysis (uses MOFA weights for gene-lipid cors)
          +---> 07_validation_experiments (uses MOFA network files)
```

## Key File Dependencies

| Producing Workflow | File | Consuming Workflow(s) |
|-------------------|------|----------------------|
| 00 | `data/RNA/fcgl_pseudobulk_gonad_*.rds` | 01 |
| 00 | `data/RNA/fcg_zonation_pseudobulks.rds` | 01 (section 11) |
| 00 | `data/RNA/gdx_pseudobulks.rds` | 01 (section 12) |
| 01 | `data/RNA/scaled_RNA_{1,2,3}.csv` | 02, 04, 05, 06 |
| 01 | `data/RNA/meta_data_df.csv` | 02, 04, 05 |
| 01 | `data/RNA/scaled_RNA_zonation.csv` | 05 |
| 01 | `data/RNA/meta_data_df_zonation.csv` | 05 |
| 01 | `data/RNA/scaled_RNA_gdxhr_{1,2,3}.csv` | 02 |
| 01 | `data/RNA/meta_data_df_gdxhr.csv` | 02, 03 |
| 02 | `results/pseudobulk/ttop_full_list.RData` | 04, 05 |
| 02 | `results/pseudobulk/tval_full_df.csv` | 02 (internal: step 7C) |
| 03 | `results/lipidomic/ttop_full_list.RData` | 05 |
| 03 | `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData` | 05 |
| 04 | `results/mofa/main_model_RNA_lipid_weight.RData` | 05, 06 |
| 04 | `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | 07 |
| 05 | `results/multi-omic/cormat_*_full_SIF.csv` (4 files) | 07 |
| 05 | `results/zonation/zonation_limma_results.RData` | 05 (internal: pathway load) |

## Coherence Status

All workflows have been audited. **No open bugs.** See individual plan files for per-workflow vulnerabilities.

| Workflow | Status | Vulnerabilities | Plan File |
|----------|--------|-----------------|-----------|
| 00 | GOOD | 5 (HPC paths, placeholder import, cluster mappings, setwd, optional packages) | `00_snrnaseq_preprocessing.md` |
| 01 | GOOD | 0 | `01_data_preparation_qc.md` |
| 02 | GOOD | 0 | `02_fcg_pseudobulk_analysis.md` |
| 03 | GOOD | 0 | `03_fcg_lipidomic_analysis.md` |
| 04 | GOOD | 0 | `04_mofa_integration.md` |
| 05 | GOOD | 4 (pathway caching, long-distance deps, sample regex, local helper) | `05_zonation_analysis.md` |
| 06 | GOOD | 2 (hardcoded focus genes, structural file matching) | `06_structural_analysis.md` |
| 07 | GOOD | 1 (inconsistent but matched network filenames) | `07_validation_experiments.md` |
| 08 | GOOD | 0 | `08_variance_analysis.md` |

## Cross-Workflow Vulnerabilities

These are fragility points that span multiple workflows:

1. **WF05 pathway caching (section 7):** `run_ulm()` is commented out; loads from `zonation_limma_results.RData`. On a fresh run without this file, fails gracefully via `tryCatch` but pathway enrichment will be missing.

2. **WF05 long-distance variable dependency (section 11 -> 14):** `sensitivity_df`, `proportion_total`, `optimal_metrics` computed in section 11, consumed ~700 lines later in section 14. If section 11 is skipped, section 14 will fail.

3. **WF07 network file naming:** Filenames use mixed conventions (`testes_portal`, `Ovaries_portal`, `ova_central`) that must exactly match what WF05 writes. Currently matched but fragile if either side is renamed.

## Shared Utility Functions

These functions from `scripts/support_functions.R` are used across multiple workflows:

| Function | Used In | Purpose |
|----------|---------|---------|
| `runLimma()` | 02, 03, 05, 07 | Wrapper for limma DE analysis |
| `volcano_nice()` | 02, 03, 05 | Publication-quality volcano plots |
| `nicePCA()` | 01, 02, 03, 05 | PCA with density distributions |
| `ttop_list_to_t_table()` | 02, 03, 05 | Convert ttop list to t-value matrix |
| `make_cor_network()` | 04, 05 | Build correlation network from two matrices |
| `expand_attributes()` | 04, 05 | Add DE statistics to network nodes |
| `get_correlated_lipids()` | 05 (local) | Get lipids correlated with gene set |
| `magicPlotMaker()` | 01 | Generate diagnostic plots |

## Recommended Pre-Run Checks

Before running any workflow, verify:

1. **Working directory** is project root
2. **Required `.RData` files** from upstream workflows exist
3. **R packages** are loaded (check `support_functions.R` for full list)
4. **No `View()` calls** in the code (breaks non-interactive execution)
5. **Explicit namespacing** for `dplyr::select()`, `dplyr::filter()`, etc.
