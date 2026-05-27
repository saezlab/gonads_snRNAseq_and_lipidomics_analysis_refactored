# Agent Guide

This repository does not require an MCP server. Agents should use the checked-in
documentation, workflow scripts, and machine-readable result files as the stable
interface for reproducing analyses and extracting results.

## First Files To Read

1. `AGENTS.md` - operating rules for this repository, including testing and
   workflow-plan requirements.
2. `README.md` - biological overview, workflow order, installation notes, and
   major result locations.
3. `docs/workflow_plans/_meta_cross_workflow_overview.md` - cross-workflow
   dependency graph and shared fragility points.
4. `docs/workflow_plans/NN_*.md` - the plan file for the workflow you are about
   to inspect or modify.
5. `docs/pipeline_rerun_plan.md` - expected headline values for full reruns.

When modifying `scripts/workflows/NN_*`, read the matching workflow plan first
and update it after the code change.

## Reproducibility Path

Workflow 00 documents raw sci-RNA-seq processing and depends on the original HPC
alignment environment. Most reproducible reruns should start at workflow 01:

```r
source("scripts/workflows/01_data_preparation_qc.R")
source("scripts/workflows/02_fcg_pseudobulk_analysis.R")
source("scripts/workflows/03_fcg_lipidomic_analysis.R")
source("scripts/workflows/04_mofa_integration.R")
source("scripts/workflows/05_zonation_analysis.R")
source("scripts/workflows/06_structural_analysis.R")
source("scripts/workflows/07_validation_experiments.R")
source("scripts/workflows/08_variance_analysis.R")

# Optional human comparison; requires external GSE179250 files and Python deps
# python scripts/workflows/09_human_analysis.py
```

Run from the repository root for workflows 00-08. After editing an R script,
run that script immediately or record why it could not be run. Workflow 09 uses
absolute external data/output paths in the script and is not part of the default
R rerun path.

## Static Query Map

Use these files directly for common questions that the discontinued MCP server
previously answered.

| Question | Primary file(s) | Notes |
| --- | --- | --- |
| Dataset summary and workflow order | `README.md`; `docs/workflow_plans/_meta_cross_workflow_overview.md` | Start here for context restoration. |
| RNA differential expression for a gene | `results/pseudobulk/ttop_full_list.RData`; `results/pseudobulk/tval_full_df.csv` | Current workflow 02 writes RData plus t-value CSV, not a long-format RNA CSV. |
| Top RNA DE genes by contrast/cell type | `results/pseudobulk/ttop_full_list.RData` | Nested list: contrast display name -> cell type -> limma topTable. |
| RNA t-statistic matrix | `results/pseudobulk/tval_full_df.csv`; `results/pseudobulk/tval_full_df.RData` | CSV has `ID` plus combined contrast/cell-type columns. |
| Mouse hepatocyte gonad signature | `results/pseudobulk/hepa_gonad_signature.csv` | Two-column signature seed (`ID`, `t_Hepatocytes`) used by workflow 09 for human-liver scoring. It is a near-identical earlier export of the current hepatocyte gonad t-values. |
| Lipid differential abundance | `results/lipidomic/ttop_full_list.RData`; `results/lipidomic/gonade_type_ttop.csv` | `logFC > 0` is testes-biased for the main gonadal contrast. `gonade_type_ttop.csv` is a generated lookup table present in results, not written by the current workflow 03 script. |
| Lipid class enrichment | `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv`; `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData` | Produced by workflow 03 using decoupleR. |
| Lipid conversion network | `results/lipidomic/network/significant_lipids_network_SIF.csv`; `results/lipidomic/network/significant_lipids_network_ATT.csv` | Current workflow 03 writes these under `results/lipidomic/network/`. |
| MOFA factor and feature weights | `results/mofa/main_model_RNA_lipid_weight.RData`; `results/mofa/factor_analysis/` | Factor 1 is the main gonad-associated factor. |
| Gene-lipid correlations | `results/mofa/networks/gene_lipid_correlations_factor1.csv`; `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | Correlation matrix rows are genes; columns are lipids. |
| TF-gene-lipid network | `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv`; `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_ATT.csv` | Cytoscape-ready compressed regulatory network. |
| Zonation DE and summaries | `results/zonation/ttop_full_long.csv`; `results/zonation/t_table_zonation.csv`; `results/zonation/summary_statistics.csv` | Workflow 05 writes long-format zonation DE. |
| Sex-specific zonation genes | `results/zonation/sex_specific_zonation_genes_summary.csv` | Main selected spatially biased gene table. |
| Lipid-zonation associations | `results/zonation/lipid_class_association_summary.csv`; `results/zonation/lipid_species_zonation_association.csv`; `results/zonation/threshold_sensitivity_total.csv` | Connects lipid species/classes to zonation categories and threshold behavior. |
| Structural lipid associations | `results/lipidomic/structural_model/lipids_with_structure.csv`; `results/lipidomic/structural_model/lipid_only_class_coefficients.csv`; `results/lipidomic/structural_model/lipid_only_class_summaries.csv`; `results/lipidomic/structural_model/lipid_only_chain_length_effects.csv`; `results/lipidomic/structural_model/lipid_only_db_position_effects.csv`; `results/lipidomic/structural_model/lipid_only_interaction_effects.csv` | Workflow 06 currently tests lipid-only structural effects on gonadal logFC. Older gene-lipid structural sections are commented out. |
| Validation experiments | `results/validation/validation_summary.csv`; `results/validation/global_network_validation_scores.csv`; `results/validation/zonation_network_validation_scores.csv`; `results/validation/structural_consistency/consistency_structural_regression_summary.csv` | Workflow 07 compares predicted correlations with knockdown lipidomics and structural consistency. |
| Lipid variance analysis | `results/lipidomic/variance_analysis/variance_analysis_summary.csv`; `results/lipidomic/variance_analysis/variance_analysis_results.csv`; `results/lipidomic/variance_analysis/variance_analysis_by_class.csv` | Workflow 08 compares variability between gonadal groups. |
| Human liver comparison | `scripts/workflows/09_human_analysis.py`; `docs/workflow_plans/09_human_analysis.md` | Optional Python workflow using external `GSE179250` human liver inputs. Outputs are written to the script's external `RESULTS_DIR`, not this repository, by default. |

## Minimal Lookup Examples

```r
# RNA DE lookup for one gene from current workflow 02 output
load("results/pseudobulk/ttop_full_list.RData")
tt <- ttop_full_list[["Testes vs Ovaries"]][["Hepatocytes"]]
subset(tt, ID == "Cyp2b9")

# Main gonad-dependent lipids from generated lookup table
lipids <- read.csv("results/lipidomic/gonade_type_ttop.csv")
head(lipids[order(lipids$adj.P.Val), ], 20)

# Zonation lookup for one gene
zonation <- read.csv("results/zonation/ttop_full_long.csv")
subset(zonation, gene == "Cyp2e1")

# MOFA factor 1 gene-lipid correlations
cors <- read.csv("results/mofa/networks/gene_lipid_correlations_factor1.csv",
                 row.names = 1, check.names = FALSE)
cors["Cyp2d9", , drop = FALSE]
```

## Interpretation Conventions

- The main gonadal contrast is Testes vs Ovaries.
- For RNA and lipid limma tables, positive `logFC` generally means higher in
  testes for Testes vs Ovaries; negative `logFC` means higher in ovaries.
- Mouse gene symbols use standard capitalization, for example `Esr1`, `Ar`,
  `Cyp2b9`, and `Mtmr11`.
- Lipid names use class plus chain notation, for example `PC(16:0/18:1)` or
  `TG(16:0_18:1_22:6)`.
- Many result files are regenerated by workflows. Prefer changing scripts, not
  hand-editing generated results.

## Agent Checklist

Before changing code:

1. Run `git status --short` and preserve unrelated user changes.
2. Read the relevant workflow plan in `docs/workflow_plans/`.
3. Check upstream/downstream dependencies in
   `docs/workflow_plans/_meta_cross_workflow_overview.md`.

After changing code:

1. Use explicit namespaces for tidyverse functions in modified R code.
2. Run the modified script immediately, or document the blocker.
3. Update the relevant workflow plan if sections, outputs, variables, or
   vulnerabilities changed.
4. Re-check `git status --short` and summarize only the files you changed.
