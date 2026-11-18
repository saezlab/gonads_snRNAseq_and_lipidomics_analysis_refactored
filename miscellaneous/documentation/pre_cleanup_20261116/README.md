# Gonadal Hormone Coordination of Sex-Specific Adult Liver Physiology

This repository contains the refactored analysis code for a multi-omic Four
Core Genotypes (FCG) mouse liver study. The goal is to separate the effects of
gonadal hormones from sex chromosome complement and to connect sex-biased liver
gene expression with lipid metabolism.

## Project Overview

- **snRNA-seq / pseudobulk RNA:** 466,986 high-quality nuclei from 48 liver
  samples (43 FCG + 5 wild-type males), summarized into pseudobulk profiles by
  sample and cell type.
- **Untargeted lipidomics:** lipid species measured across FCG liver samples.
  The current MOFA path uses the RNA-shared lipidomic sample set.
- **Multi-omic integration:** MOFA2-based integration of 5 RNA cell-type views
  and 1 lipidomic view.
- **Validation and provenance:** GDX RNA/lipidomic validation, Xenium spatial
  validation, proteomics validation, siRNA knockdown lipidomics, and
  collaborator Xenium/MALDI/proteomics source materials.
- **Human liver comparison:** optional Python workflow comparing the mouse
  hepatocyte gonadal t-value signature with human bulk liver sex differences
  from GEO accession GSE179250.

## Key Findings

- Gonadal hormones, rather than sex chromosome complement, are the dominant
  driver of sex-biased adult liver physiology.
- Hepatocytes show the strongest gonad-dependent RNA effects, with roughly
  4,000 gonad-dependent genes in the main analysis.
- Lipid abundance is strongly gonad-dependent, while sex chromosome effects on
  lipids are minimal to absent in the current analysis narrative.
- Spatial zonation modulates gonad-dependent hepatocyte expression.
- Key validated regulatory genes include `Cyp2b9`, `Cyp2b13`, `Cyp2c38`,
  `Cyp2d9`, `Sult3a1`, and `Mtmr11`.

For exact counts, prefer the generated result files over prose summaries. In
particular, workflow 03 currently contains two different printed headline lipid
counts in comments/output text, so result tables are the quantitative source.

## Four Core Genotypes Model

| Group | Sex Chromosomes | Gonads | Description |
|-------|-----------------|--------|-------------|
| XXO | XX | Ovaries | Chr3-Sry- |
| XXT | XX | Testes | Chr3-Sry+ |
| XYO | XY | Ovaries | Chr3-Sry- |
| XYT | XY | Testes | Chr3-Sry+ |

Main comparisons:

- **Gonadal effect:** Testes vs Ovaries, across chromosome backgrounds.
- **Chromosomal effect:** XY vs XX, across gonadal backgrounds.
- **Group comparisons:** XXT vs XXO, XYT vs XYO, XYO vs XXO, XYT vs XXT.

## Agent-Oriented Data Exploration

This repository is documented for direct agent use without an MCP server.

Start with:

- `AGENTS.md` - repository rules and biological context for agents.
- `docs/AGENT_GUIDE.md` - static query map for extracting results.
- `docs/workflow_plans/_meta_cross_workflow_overview.md` - cross-workflow data
  flow and shared fragility points.
- `docs/workflow_plans/NN_*.md` - per-workflow structure, dependencies,
  outputs, and known vulnerabilities.

When modifying any `scripts/workflows/NN_*` workflow file, read the matching
plan first and update it after the code change.

## Repository Structure

```text
.
|-- data/
|   |-- RNA/                         # Pseudobulk RNA CSVs and Seurat RDS inputs
|   |-- metabolomics/                # Lipidomics and structural input CSVs
|   |-- proteomics/                  # Bulk liver proteomics validation data
|   `-- Suppl_Table_1_OverviewData.xlsx
|-- scripts/
|   |-- workflows/                   # Main workflows 00-10
|   |-- core/                        # Data loading, DE, enrichment, MOFA, structural helpers
|   |-- config/                      # Dataset paths, contrasts, plotting config
|   |-- visualization/               # Plot helpers
|   |-- moon/                        # Selected TF-to-lipid workflow, ancillary searches, historical archive
|   |-- pkn/                         # Optional mouse lipid-expanded prior construction
|   |-- james/                       # Collaborator source/provenance material
|   |-- support_functions.R          # Compatibility wrapper to config/support_functions.R
|   `-- seurat_compatibility.R
|-- support/                         # Gene sets, regulons, hormone targets, LINEX2 network
|-- docs/
|   |-- AGENT_GUIDE.md
|   `-- workflow_plans/
|-- results/
|   |-- qc/
|   |-- pseudobulk/
|   |-- lipidomic/
|   |-- mofa/
|   |-- causal_modelling/            # Frozen COSMOS resources and regenerated selected MOON result
|   |-- zonation/
|   |-- validation/
|   `-- xenium/
`-- tests/
```

## Workflow Order

Most users should start at workflow 01. Workflow 00 is a reference/HPC raw-data
preprocessing workflow and is inactive by default.

```r
# Optional reference workflow; requires original raw alignment environment
# source("scripts/workflows/00_snrnaseq_preprocessing.R")

source("scripts/workflows/01_data_preparation_qc.R")
source("scripts/workflows/02_fcg_pseudobulk_analysis.R")
source("scripts/workflows/03_fcg_lipidomic_analysis.R")
source("scripts/workflows/04_mofa_integration.R")
source("scripts/workflows/05_zonation_analysis.R")
source("scripts/workflows/06_structural_analysis.R")
source("scripts/workflows/07_validation_experiments.R")
source("scripts/workflows/08_variance_analysis.R")

# Optional external-data workflow; requires large GSE179250 inputs outside repo
# python scripts/workflows/09_human_analysis.py

# Optional collaborator Xenium scripts; require raw Xenium data under the
# workflow folder and are run from scripts/workflows/10_xenium_analysis/
# python 00_infer_doublets.py
# python 01_merge_doublets.py
# python 02_compositional_analysis.py
# python 03_DGE_analysis.py
```

The October 2026 validation work packages use frozen outputs from the original
analyses above and repository-local new inputs:

```sh
Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --validate-only
Rscript scripts/workflows/12_wdsd_network_comparison.R --validate-only
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --validate-only
```

Run workflows 00-08 and 11-14 from the repository root. Workflows 11-14 use the `Rscript` commands documented below. Workflows 09 and 10 are optional
Python analyses with their own external or folder-local path assumptions, as
documented below and in their plan files.

## Analysis Workflows

### Optional Factor 1 causal modelling

The selected manuscript TF-to-lipid analysis has one entry point, run from the
repository root:

```sh
Rscript scripts/moon/run_factor1_tf_to_lipid.R --validate-only
Rscript scripts/moon/run_factor1_tf_to_lipid.R --output=results/causal_modelling/factor1_tf_to_lipid_replay
```

The first command checks configuration and dependencies without writing output.
The workflow consumes the frozen lipid-expanded mouse PKN and audited Factor 1
inputs, then regenerates compressed MOON scoring, RNA-coherence filtering,
network extraction, expected-TF diagnostics and provenance/plots. Settings are
in `scripts/moon/config/factor1_tf_to_lipid.json`; results go to
`results/causal_modelling/factor1_tf_to_lipid/` by default. Existing nonempty
output folders are refused, so the example supplies a fresh replay destination.
The selected configuration uses
six pruning steps, four scoring layers and cutoffs 1.85/0.55. Its reference
comparison checks against the frozen 388-edge/153-node manuscript network.
The full refit and exports passed on 2026-09-15: graph identities, compression
mapping, score identities/levels and numerical scores within 1e-10 matched the
frozen reference, including the expected-TF exclusions. The generated
`reference_comparison.json` records the checks for each execution.

Parameter searches use `scripts/moon/ancillary/search_factor1_parameters.R`
and the same scoring/extraction core. The original pilot and expanded searches
remain in `scripts/moon/archive/2026-09/`, with unchanged results under
`results/causal_modelling/cosmos_*`. They are historical records, not the active
commands.

The separate lipid-to-TF workflow is now available:

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf.R --validate-only --include-ph73
Rscript scripts/moon/run_factor1_lipid_to_tf.R --include-ph73 \
  --output=results/causal_modelling/factor1_lipid_to_tf_replay
```

It audits signalling coverage before fitting, retains the PKN's original edge
direction, and uses lipid Factor 1 weights upstream and hepatocyte TF footprints
downstream. The original PKN's negative results and the optional acid/base-
normalized precursor sensitivity are kept separate. Results, saved compressed
fits and supporting audits are under
`results/causal_modelling/factor1_lipid_to_tf/`; see the
[lipid-to-TF analysis report](docs/cosmos_factor1_lipid_to_tf.md) for the selected
exploratory network, limitations and ancillary reproduction commands.

A separate directional-STITCH relaxed-mapping screen treats each plausible
feature-to-chemical identity as an auditable upstream candidate, without copying
the lipid weight into the ULM calculation. The active selected-analysis entry
point is `scripts/moon/run_factor1_lipid_to_tf_relaxed.R`, configured by
`scripts/moon/config/factor1_lipid_to_tf_relaxed.json`, with default output
under `results/causal_modelling/factor1_lipid_to_tf_relaxed/`. It refits the
predeclared five-edge LPS sensitivity and checks the score, compression and
strict route against the saved legacy reference. The depth and dual-threshold
search is a separate ancillary workflow under
`scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R`.
It writes independently to
`results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search/`; the
candidate audit writes under the selected analysis and the evidence audit writes
under the parameter search. Candidate-local score traces and strict
score-layer-constrained routes remain separate from broad sign-consistent
connectivity, so duplicated chemical mappings cannot be interpreted as
independent regulators. The pre-cleanup discovery topology/search under
`results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory/` is retained
only as legacy audit material. See
[the relaxed lipid-to-TF report](docs/cosmos_factor1_lipid_to_tf_relaxed.md).
The LPS→Gpr34→Gnai1→Src→{Ar, Esr1} route is a 0.5/0.5 low-threshold
sensitivity with no compact result at either higher cutoff; AR and ESR1 are
supplied terminal observations.

See [the MOON workflow guide](scripts/moon/README.md),
[the TF-to-lipid workflow plan](docs/workflow_plans/moon_factor1_tf_to_lipid.md),
[the relaxed lipid-to-TF workflow plan](docs/workflow_plans/moon_factor1_lipid_to_tf_relaxed.md) and
[the Results and Methods draft](docs/manuscript_drafts/cosmos_results_methods.md).

### Workflow 00: snRNA-seq Pre-processing

**Script:** `scripts/workflows/00_snrnaseq_preprocessing.R`

Reference pipeline for raw sci-RNA-seq import, UMI QC, Seurat clustering,
scDblFinder doublet removal, cell-type annotation, pseudobulking, and
hepatocyte zonation annotation. Both FCG and GDX processing blocks are wrapped
in `if (FALSE)` and depend on original DKFZ HPC paths plus a placeholder import
function.

Key intended outputs include `pseudobulks_scc.rds`,
`heps_zonationsubsetted_pseudobulks_gonad.rds`, and GDX pseudobulk RDS files.

### Workflow 01: Data Preparation and QC

**Script:** `scripts/workflows/01_data_preparation_qc.R`

Converts preprocessed pseudobulk Seurat objects into CSV matrices and metadata.
The main RNA CSVs are split by sample columns, with all genes retained in each
file.

Key outputs:

- `data/RNA/scaled_RNA_1.csv`
- `data/RNA/scaled_RNA_2.csv`
- `data/RNA/scaled_RNA_3.csv`
- `data/RNA/meta_data_df.csv`
- `data/RNA/scaled_RNA_zonation.csv`
- `data/RNA/meta_data_df_zonation.csv`
- `data/RNA/scaled_RNA_gdxhr_{1,2,3}.csv`
- `data/RNA/meta_data_df_gdxhr.csv`
- `results/qc/`

### Workflow 02: FCG Pseudobulk RNA Analysis

**Script:** `scripts/workflows/02_fcg_pseudobulk_analysis.R`

Runs RNA differential expression for 6 FCG contrasts across 5 cell types, then
adds GDX RNA validation, GO pathway enrichment, TF activity inference, Xenium
panel validation, hepatocyte PCA coordinates, and bulk proteomics validation.

Key outputs:

- `results/pseudobulk/ttop_full_list.RData`
- `results/pseudobulk/tval_full_df.RData`
- `results/pseudobulk/tval_full_df.csv`
- `results/pseudobulk/gdxrh/gdhxr_groups_ttop_list.RData`
- `results/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv`
- `results/pseudobulk/volcanos/`
- `results/pseudobulk/pca/`
- `results/pseudobulk/pathways/pathway_scores_df_all.RData`
- `results/pseudobulk/pathways/tf_scores_df_all.RData`
- `results/pseudobulk/gdx_validation_comparison.pdf`
- `results/pseudobulk/xenium_validation_volcano.pdf`
- `results/pseudobulk/proteo_validation_comparison.pdf`

### Workflow 03: FCG Lipidomic Analysis

**Script:** `scripts/workflows/03_fcg_lipidomic_analysis.R`

Runs lipid differential abundance for 6 FCG contrasts, optional GDX lipid
validation, lipid-class enrichment, LINEX2 network extraction, and per-class
chain-length/double-bond bubble plots.

Key outputs:

- `results/lipidomic/ttop_full_list.RData`
- `results/lipidomic/gdxhr_ttop_list.RData`
- `results/lipidomic/gdxhr_t_table.csv`
- `results/lipidomic/fcg_gdx_correlation.pdf`
- `results/lipidomic/volcanos/`
- `results/lipidomic/pca/`
- `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData`
- `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv`
- `results/lipidomic/enrichment/lipid_class_enrichment_plot.pdf`
- `results/lipidomic/network/significant_lipids_network_SIF.csv`
- `results/lipidomic/network/significant_lipids_network_ATT.csv`
- `results/lipidomic/bubbles_DBC/`

### Workflow 04: MOFA Integration

**Script:** `scripts/workflows/04_mofa_integration.R`

Builds a 6-view MOFA input (5 RNA cell types + lipidomics), trains or loads
5-15 factor MOFA models, selects the 12-factor model, extracts weights, enriches
factors for gonadal type, and builds gene-lipid and TF-gene-lipid networks.

Key outputs:

- `results/mofa/MOFA_data.RData`
- `results/mofa/RNA/mofa_res_12factor.hdf5`
- `results/mofa/main_model_RNA_lipid_weight.RData`
- `results/mofa/factor_analysis/factor_enrichment_gonadal.RData`
- `results/mofa/factor_analysis/*.pdf`
- `results/mofa/networks/gene_lipid_correlations_factor1.csv`
- `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_full_{SIF,ATT}.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_{SIF,ATT}.csv`

### Workflow 05: Hepatocyte Zonation Analysis

**Script:** `scripts/workflows/05_zonation_analysis.R`

Uses workflow 01 zonation CSVs to run duplicate-correlation limma models across
periportal/pericentral and gonadal contrasts. It selects sex-specific
zonation-biased genes through an entropy/stability/coverage/specificity
threshold sweep, builds zonation-specific gene-lipid networks, and maps lipid
classes to zonation categories using workflow 03 outputs.

Key outputs:

- `results/zonation/zonation_ttop_list.RData`
- `results/zonation/zonation_limma_results.RData`
- `results/zonation/t_table_zonation.csv`
- `results/zonation/ttop_full_long.csv`
- `results/zonation/sex_specific_zonation_genes_summary.csv`
- `results/zonation/sex_specific_zonation_panel_AB.pdf`
- `results/zonation/sex_specific_zonation_heatmap_faceted.pdf`
- `results/zonation/lipid_species_zonation_association.csv`
- `results/zonation/lipid_class_association_summary.csv`
- `results/zonation/threshold_sensitivity_total.csv`
- `results/zonation/optimal_threshold_metrics.csv`
- `results/zonation/networks/cormat_*_{SIF,ATT}.csv`
- `results/zonation/networks/cormat_*_full_{SIF,ATT}.csv`
- `results/zonation/networks/compressed_cormat_*_{SIF,ATT}.csv`

The October 2026 WDSD work package is now [Workflow 12](#workflow-12-wdsd-pooled-and-zonation-networks-wp2). The former 05b/05c entry points remain thin compatibility wrappers for its individual phases.

### Workflow 06: Lipid Structural Analysis

**Script:** `scripts/workflows/06_structural_analysis.R`

Current active workflow is lipid-only. It parses lipid structures from workflow
03 gonadal lipid DE results and tests whether lipid chain length, double-bond
count, and their interaction explain gonadal logFC within lipid classes.

Active model:

```r
logfc ~ chain_length + db_position + chain_length:db_position
```

Key outputs:

- `results/lipidomic/structural_model/lipids_with_structure.csv`
- `results/lipidomic/structural_model/lipid_only_class_coefficients.csv`
- `results/lipidomic/structural_model/lipid_only_class_summaries.csv`
- `results/lipidomic/structural_model/lipid_only_chain_length_effects.csv`
- `results/lipidomic/structural_model/lipid_only_db_position_effects.csv`
- `results/lipidomic/structural_model/lipid_only_interaction_effects.csv`
- `results/lipidomic/structural_model/lipid_only_summary_heatmap.pdf`
- `results/lipidomic/structural_model/lipid_only_plots/`

The older gene-lipid structural-regression sections are commented out in the
current script.

### Workflow 07: Validation Experiments

**Script:** `scripts/workflows/07_validation_experiments.R`

Validates MOFA and zonation network predictions using siRNA knockdown
lipidomics from primary hepatocytes, then evaluates whether knockdown-response
consistency varies by lipid structure.

Key outputs:

- `results/validation/knockdown_PCA.pdf`
- `results/validation/knockdown_ttop_list.RData`
- `results/validation/knockdown_t_table.csv`
- `results/validation/knockdown_heatmap.pdf`
- `results/validation/global_network_validation_scores.csv`
- `results/validation/[gene]_validation_scatter.pdf`
- `results/validation/zonation_network_validation_scores.csv`
- `results/validation/[gene]_zonation_validation_scatter.pdf`
- `results/validation/validation_summary.csv`
- `results/validation/structural_consistency/[gene]_structural_consistency.pdf`
- `results/validation/structural_consistency/consistency_structural_regression_summary.csv`

The independent in vivo Cyp2c38 work package is now [Workflow 11](#workflow-11-cyp2c38-in-vivo-validation-wp1). Former 07b/07c entry points remain thin compatibility wrappers; the primary-hepatocyte siRNA workflow above remains separate.

### Workflow 08: Lipidomic Variance Analysis

**Script:** `scripts/workflows/08_variance_analysis.R`

Compares lipid variability between gonadal groups using CVs and Brown-Forsythe
Levene tests at both lipid and lipid-class levels.

Key outputs:

- `results/lipidomic/variance_analysis/variance_analysis_results.csv`
- `results/lipidomic/variance_analysis/significant_differential_variance_lipids.csv`
- `results/lipidomic/variance_analysis/variance_analysis_by_class.csv`
- `results/lipidomic/variance_analysis/variance_analysis_summary.csv`
- `results/lipidomic/variance_analysis/cv_difference_volcano.pdf`
- `results/lipidomic/variance_analysis/cv_difference_volcano_by_class.pdf`

### Workflow 09: Human Liver Comparison

**Script:** `scripts/workflows/09_human_analysis.py`

Optional Python workflow that compares the mouse FCG hepatocyte gonadal
t-value signature with sex-associated gene expression in the Zhou et al. human
bulk liver cohort ([`GSE179250`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179250)).
It translates mouse signature genes to human ortholog symbols, scores the
signature with decoupler ULM, runs PyDESeq2 male-vs-female differential
expression in human liver, and generates concordant mouse-human ranked-gene
plots.

The raw human count/TPM matrices and GEO SOFT metadata are too large for this
repository and are expected under the script's external `DATA_DIR`. The checked
in signature seed is:

- `results/pseudobulk/hepa_gonad_signature.csv`

The script writes outputs to its external `RESULTS_DIR`, including:

- `ulm_sex_signature_estimates.csv`
- `adata_processed.h5ad`
- `deg_zhou.csv`
- `full_results_deg_zhou_malevsfemale.csv`
- `12_49_malevsfemale_deg_zhou.csv`
- `12_49_malevsfemale_full_results_deg_zhou.csv`
- `olderthan11female_oldvsyoung_full_results_deg_zhou.csv`
- `olderthan11female_oldvsyoung_deg_statmatrix_zhou.csv`
- `pca_by_sex.pdf`
- `ulm_enrichment_12_49.pdf`
- `volcano_male_vs_female_highlight.pdf`
- `volcano_male_vs_female_highlight_1249.pdf`
- `top_ranked_genes_female_dot.pdf`
- `top_ranked_genes_male_dot.pdf`
- `ulm_enrichment_female__old_vs_young.pdf`

### Workflow 10: Xenium Spatial Analysis

**Folder:** `scripts/workflows/10_xenium_analysis/`

Optional collaborator-provided Python workflow for Xenium spatial validation.
Unlike workflows 00-09, this is a folder of notebook-style scripts rather than
a single canonical workflow script. It infers neighboring-nucleus doublets from
Xenium nucleus boundaries, annotates and fuses likely hepatocyte doublets,
creates cell-type/doublet composition plots, and runs Scanpy Wilcoxon
differential expression for sex, hepatocyte-zone, and doublet comparisons.

The scripts use relative paths and should be treated as folder-local: run them
from `scripts/workflows/10_xenium_analysis/` with raw Xenium files under
`data_raw/` and output folders such as `data_tidy/` and `figures/` available.
The raw Xenium bundle is not committed.

Main scripts:

- `00_infer_doublets.py`
- `01_merge_doublets.py`
- `02_compositional_analysis.py`
- `03_DGE_analysis.py`

Key folder-local inputs include:

- `data_raw/output-XETG00078__0018188__Region_1__20240301__143248/nucleus_boundaries.csv.gz`
- `data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cell_feature_matrix.h5`
- `data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cells.csv.gz`
- `data_raw/XeniumTargetList_Mouse-Panel-1.xlsx`

Key folder-local outputs include:

- `data_tidy/inferred_components.csv`
- `inferred_components.json`
- `fused_anndata_annotated.h5ad`
- `data_tidy/fused_doublet.csv`
- `data_tidy/fused_celltype.csv`
- `figures/celltype_freq_cores.pdf`
- `figures/celltype_freq_sex.pdf`
- `figures/doublet_freq_core.pdf`
- `data_tidy/DEG_results_sex_per_celltype.xlsx`
- `DEG_results.xlsx`
- `Pairwise_DEG_hepatocytes.xlsx`
- `Doublets_DEG_results.xlsx`
- `Doublets_and_Sex_DEG_results_all_hepatocytes.xlsx`

### Workflow 11: Cyp2c38 In Vivo Validation (WP1)

**Script:** `scripts/workflows/11_cyp2c38_in_vivo_validation.R`.

One workflow performs KO lipidomics intake/QC, weighted limma contrasts and
fixed-reference validation, then native ULM, full ovarian/testicular networks,
contrast-null bootstrap, permutation and conditional calibration diagnostics.
Its sectioned entry point uses `scripts/core/wp1_workflow.R` and the existing
statistical helpers. The integrated export now lives in
`data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv`.

Run from the repository root:

```sh
Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --validate-only
Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --output=results/validation/cyp2c38_in_vivo_rerun
```

The defaults retain base seed 20261006, 1,999 permutations and 999 combined
bootstrap draws, followed by ULM seed 20261007, 4,999 contrast-null draws,
1,999 permutations, 120 calibration datasets/199 inner draws and 399 injected
datasets. An explicit `--bootstraps` overrides both phases unless
`--base-bootstraps` is also given; use `--base-bootstraps=999` to preserve the
historical base setting. Completed outputs are protected; choose a fresh
`--output` directory. Validation reads local raw inputs without saved WP1 results.

Outputs remain in `results/validation/cyp2c38_in_vivo/`, with the requested
ULM-focused report and tables in `ulm_followup/`. See [the plan](docs/workflow_plans/11_cyp2c38_in_vivo_validation.md)
and [the current ULM report](results/validation/cyp2c38_in_vivo/ulm_followup/README.md).

### Workflow 12: WDSD Pooled and Zonation Networks (WP2)

**Script:** `scripts/workflows/12_wdsd_network_comparison.R`.

One workflow runs the frozen pooled Figure 3 comparison first, then the four
sex/zonation reference comparisons and within-cohort WD-minus-SD changes. Its
combined report follows that order, with annotation overlap, QC, sensitivities
and reproducibility documented afterwards. The statistical models and original
feature selection are unchanged; overall pooled and sex-conditioned null tests
remain separate. `scripts/core/wp2_workflow.R` orchestrates the phases and
`scripts/core/wp2_report.R` writes the combined summary.

```sh
Rscript scripts/workflows/12_wdsd_network_comparison.R --validate-only
Rscript scripts/workflows/12_wdsd_network_comparison.R --output=results/wdsd_network_comparison_rerun
Rscript scripts/workflows/12_wdsd_network_comparison.R --report-only
```

Defaults retain seed 20261007, 19,999 pooled random pairing draws, all 5,040
exact pairings per zonal context, 4,999 null-bootstrap draws and 999 paired
stability draws. Full runs require a fresh output directory; `--report-only`
refreshes the combined Markdown from completed tables without rerunning models.
The pooled-first fresh run no longer needs prior WP2 matching ledgers.

Inputs are repository-local under `data/metabolomics/wdsd/`,
`data/metabolomics/fcg/` and `data/RNA/wdsd/`. The [input inventory and hash manifest](data/WP1_WP2_INPUTS.md)
document unchanged copies and the preserved shared archives. CF459/WP3 inputs
remain outside these completed work packages.

Outputs remain in `results/wdsd_network_comparison/`, with pooled artifacts in
`followup/`. See [the plan](docs/workflow_plans/12_wdsd_network_comparison.md)
and [the pooled-first report](results/wdsd_network_comparison/README.md).
`build_wp2_report_pdf.py` renders that report and the ordered scientific figures.

### Workflow 13: Western-Diet Cross-Cohort Robustness (WP3)

**Script:** `scripts/workflows/13_wd_crosscohort_robustness.R`.

Compares male Western-diet gene–lipid networks between CF453 (seven mice) and
CF459 (four mice), separately periportal/pericentral. Original Figure 5 genes
and molecular pairs remain frozen; all common measurable full-network pairs
are primary, with the originally displayed pairs secondary. No CF459
standard-diet arm or new MOFA fit is involved. A dated extension compares CF459
also with the original 22 testes-bearing FCG SD mice (11 XXT/11 XYT) and
computes `delta r = cor(CF459, FCG) - cor(CF459, CF453)` on identical pairs.
User-confirmed ages are 3 months for CF459/FCG and 10 months for CF453.
The original research plan itself is amended in [the plan copy](docs/research_plans/2026-10-05_Western_diet_and_Cyp2c38_analysis_plan.md).

```sh
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --validate-only
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --output=results/wd_crosscohort_robustness_reproduction
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --report-only
# Add to saved original WP3 data; target must be fresh:
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --fcg-comparison-only --comparison-output=results/wd_crosscohort_robustness/fcg_comparison_reproduction
python3 scripts/reports/build_wp3_report_pdf.py
```

Defaults: seed 20261014, 4,999 independent-marginal CF459 null bootstraps and
999 paired two-cohort stability draws, reused for the direct comparison with
an additional independently seeded FCG paired bank stratified by XXT/XYT.
The same CF459 draw enters both anchors; positive-difference fractions are
descriptive stability, not null p-values. The extension checks another 22 FCG
individual omissions. Exact CF459 (24) and complementary
CF453 (5,040) RNA–lipid pairing permutations test distinct conditional nulls;
all 28 joint and 11 single-cohort omissions describe influence. No within-draw
feature deletion is allowed. Native ULM/nominal p-values accompany correlation
magnitudes and animal-level null/stability checks.

Inputs reside under `data/metabolomics/wd/` and `data/RNA/wd/`, with byte/hash
provenance in [the WP3 input guide](data/WP3_INPUTS.md). CF459 decimal-comma
CSV import, updated Seurat layer access and mouse pairing are checked explicitly.
Analytical QC, Fill and richer chemical annotation are unavailable in CF459;
QC/Fill sensitivity restrictions therefore apply only to CF453.

Outputs are isolated in `results/wd_crosscohort_robustness/`; full runs require
a fresh output directory. See [the workflow plan](docs/workflow_plans/13_wd_crosscohort_robustness.md)
and [the technical report](results/wd_crosscohort_robustness/README.md).
The report retains n=4 uncertainty and unresolved background/exposure/processing
provenance. CF459 is relatively closer to FCG on the full pairs (r=0.052/0.080
versus -0.063/-0.019 with CF453), but absolute alignment is weak and Mtmr11
does not follow this ranking. Cohort differences cannot identify causal age or
diet effects. Extension files also record an exact source-index contingency
certificate for all-zero null correlation vectors, qualifying the preserved
historical bootstrap results without changing observed comparisons.

### Workflow 14: Standalone CF453 Figure 6

**Script:** `scripts/workflows/14_cf453_standalone_figure6.R`.

This independent analysis uses CF453 alone: 28 mice in a sex x diet design,
with paired hepatocyte zones and matched whole-liver lipids. It estimates
sex-dependent zonation and its diet modification, selects CF453 genes de novo,
and exports independent Figure 6 A-F panel files for manual assembly. It does
not use original-cohort selected networks or cross-cohort molecular overlap.

```bash
Rscript scripts/workflows/14_cf453_standalone_figure6.R --output=results/cf453_standalone_figure6_reproduction
Rscript scripts/workflows/14_cf453_standalone_figure6.R --plot-only
```

Full reruns require a destination without saved models/analysis. Outputs in
`results/cf453_standalone_figure6/` include PDF/SVG/PNG panels, source data,
factorial model tables, complete association coefficients, paired bootstrap
indices, legends and independent numerical/visual checks. The networks remain
exploratory at seven mice per sex/diet. See the
[workflow plan](docs/workflow_plans/14_cf453_standalone_figure6.md) and
[standalone results](results/cf453_standalone_figure6/README.md).

### Workflow 15: Four Direct Diet-Impact Figure 6 Panels

This earlier compact Figure 6 presentation directly models how WD versus SD
changes sex differences, zonation and sex-dependent zonation in CF453. It uses
the validated Workflow 14 model-input cache, preserves the historical panels,
and adds no network-based selection. Four independent PDF/SVG/PNG panels cover
RNA interactions, their gene-level decomposition, mouse examples and the direct
sex-specific lipid response. [Results and panel index](results/cf453_diet_impact_figure6/README.md)
and the [workflow plan](docs/workflow_plans/15_cf453_diet_impact_figure6.md) define
the new testing families and interpretation.

```bash
Rscript scripts/workflows/15_cf453_diet_impact_figure6.R --output=results/cf453_diet_impact_figure6_reproduction
Rscript scripts/workflows/15_cf453_diet_impact_figure6.R --plot-only
python3 scripts/reports/package_cf453_diet_impact.py
```

### Workflow 16: Selected Zonation Programmes and Diet Association Changes

The current Figure 5-style revision uses the standalone cohort of 28
approximately 10-month-old wild-type mice (seven per sex/diet), without an
original-cohort comparison. Three independent panel groups show reciprocal
SD/WD-selected zonation-gene trajectories, relative lipid-class association
changes for those same selected programmes, and exploratory gained/lost
gene-lipid association examples. Fixed genes, lipid features and matched
mouse resampling make the diet comparison explicit. Individual associations
and class shifts remain exploratory at n=7; no selected correlation passes
the declared family-wide BH threshold. See the
[panel index](results/cf453_diet_rewiring_figure6/README.md) and
[plan](docs/workflow_plans/16_cf453_diet_rewiring_figure6.md).

```bash
Rscript scripts/workflows/16_cf453_diet_rewiring_figure6.R --output=results/cf453_diet_rewiring_figure6_reproduction
Rscript scripts/workflows/16_cf453_diet_rewiring_figure6.R --plot-only
python3 scripts/reports/package_cf453_diet_rewiring.py
```

Workflow 16 reads validated Workflow 14 models and Workflow 15 direct-diet
evidence; historical outputs are preserved. Each panel has PDF/SVG/PNG and a
source CSV, with methods, selection/sensitivity tables and a portable ZIP for
manual figure assembly.

## Collaborator Provenance Materials

### Expanded-marker zonation update, 16 November2026

The versioned rerun uses James's supplied multi-marker zonation pseudobulks
for FCG, gonadectomy and standard/Western-diet mice. It preserves historical
analyses and produces a linked panel trace, manuscript numerical replacement
ledger, updated source tables and Figure6 layout alternatives. Figure4A and
S4B-C are explicitly deferred to James. See the
[report](results/zonation_marker_update_20261116/report/zonation_update_report.html)
and [Workflow17 plan](docs/workflow_plans/17_expanded_marker_zonation_update.md).

```bash
Rscript scripts/workflows/17_expanded_marker_zonation_update.R --validate-only
Rscript scripts/workflows/17_expanded_marker_zonation_update.R --report-only
```

Additional collaborator source files that are not part of the default workflow
order are kept under `scripts/james/`:

- `scripts/james/maldi/Cleland_analysis_expt2_recal.ipynb` - attachment 1
  from the 2026-06-19 missing-material bundle; negative-mode MALDI-MSI liver
  zonation notebook using external imaging/AnnData/SCiLS inputs.
- `scripts/james/proteomics/JamesCleland_P3225_V1.Rmd` - attachment 6 from
  the same bundle; original P3225 fractionated TMT proteomics R Markdown
  analysis. Workflow 02 does not source this file; it uses the curated
  proteomics inputs in `data/proteomics/`.
- Attachments 2-5 from the same bundle were Xenium scripts already present in
  `scripts/workflows/10_xenium_analysis/`; the archived copies matched the
  checked-in scripts after ignoring CRLF line-ending differences.

## Requirements

Workflows 00-08 are R scripts and expect packages from CRAN and Bioconductor.
Core R packages include:

- Data manipulation: `dplyr`, `tidyr`, `tibble`, `readr`, `readxl`, `reshape2`
- Visualization: `ggplot2`, `ggrepel`, `pheatmap`, `cowplot`, `patchwork`,
  `gridExtra`, `ggridges`, `scales`
- Statistics: `limma`, `decoupleR`, `car`, `vsn`
- Single-cell/preprocessing: `Seurat`, `glmGamPoi`, optional `scDblFinder`
- Multi-omic integration: `MOFA2` and Python `mofapy2` for retraining models
- Gene sets: `GSEABase`

Workflow 09 is Python and expects packages including `anndata`, `scanpy`,
`decoupler`, `pydeseq2`, `pypath`, `pandas`, `numpy`, `matplotlib`, `seaborn`,
`scipy`, and `adjustText`.

Workflow 10 is Python and expects packages including `scanpy`, `anndata`,
`sccoda`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `tifffile`, `shapely`,
`rtree`, `networkx`, `tqdm`, and `openpyxl`.

## Key Results Files

### RNA Differential Expression

- `results/pseudobulk/ttop_full_list.RData`
- `results/pseudobulk/tval_full_df.csv`
- `results/pseudobulk/hepa_gonad_signature.csv`

### Lipid Differential Abundance

- `results/lipidomic/ttop_full_list.RData`
- `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv`

### MOFA and Networks

- `results/mofa/main_model_RNA_lipid_weight.RData`
- `results/mofa/networks/gene_lipid_correlations_factor1.csv`
- `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv`

### Zonation

- `results/zonation/ttop_full_long.csv`
- `results/zonation/sex_specific_zonation_genes_summary.csv`
- `results/zonation/lipid_class_association_summary.csv`
- `results/zonation/networks/`

### Structural and Variance Analysis

- `results/lipidomic/structural_model/lipids_with_structure.csv`
- `results/lipidomic/structural_model/lipid_only_class_coefficients.csv`
- `results/lipidomic/variance_analysis/variance_analysis_results.csv`
- `results/lipidomic/variance_analysis/variance_analysis_by_class.csv`

### Validation

- `results/validation/global_network_validation_scores.csv`
- `results/validation/zonation_network_validation_scores.csv`
- `results/validation/structural_consistency/consistency_structural_regression_summary.csv`

### Human Liver Comparison

- `results/pseudobulk/hepa_gonad_signature.csv`
- Workflow 09 outputs are written outside this repository by default; see
  `scripts/workflows/09_human_analysis.py` and
  `docs/workflow_plans/09_human_analysis.md`.

### Xenium Spatial Analysis

- `scripts/workflows/10_xenium_analysis/`
- Workflow 10 outputs are folder-local by default; see
  `docs/workflow_plans/10_xenium_analysis.md`.

### Collaborator Provenance

- `scripts/james/maldi/Cleland_analysis_expt2_recal.ipynb`
- `scripts/james/proteomics/JamesCleland_P3225_V1.Rmd`
- `scripts/james/README.md`

## Data Availability

The checked-in CSV and RData outputs support result lookup and most downstream
analysis. Some raw/preprocessed inputs may not be present in all checkouts:

- `data/RNA/fcg_pseudobulks.rds`
- `data/RNA/fcg_zonation_pseudobulks.rds`
- `data/RNA/gdx_pseudobulks.rds`
- MOFA `.hdf5` models if retraining/loading is needed
- Workflow 09 human liver inputs from `GSE179250`, including the count matrix,
  TPM matrix, and SOFT metadata file.
- Workflow 10 raw Xenium files under `scripts/workflows/10_xenium_analysis/data_raw/`,
  including the 10x output folder and Xenium panel workbook.
- MALDI-MSI raw imaging/AnnData/registration/SCiLS inputs referenced by
  `scripts/james/maldi/Cleland_analysis_expt2_recal.ipynb`.
- Raw collaborator P3225 proteomics files referenced by
  `scripts/james/proteomics/JamesCleland_P3225_V1.Rmd`; the curated validation
  inputs used by workflow 02 are checked in under `data/proteomics/`.

Workflow 00 requires the original raw alignment environment and is not expected
to run in a standard checkout.

## Citation

Manuscript citation to be added upon publication.

## License

GNU General Public License v3.0. See `LICENSE`.
