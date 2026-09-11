# Cross-Workflow Dependencies and Coherence Report

**Last audited:** 2026-06-19
**Ground truth:** `scripts/workflows/*.R`,
`scripts/workflows/09_human_analysis.py`, and
`scripts/workflows/10_xenium_analysis/`

Non-workflow collaborator provenance files under `scripts/james/` are
documented here but are not treated as cross-workflow dependencies.

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

10_xenium_analysis
  optional collaborator Python workflow folder; raw Xenium output + panel
  workbook -> doublet calls, fused AnnData, composition plots, and DGE tables

Non-workflow collaborator provenance
  scripts/james/maldi/Cleland_analysis_expt2_recal.ipynb
  scripts/james/proteomics/JamesCleland_P3225_V1.Rmd
  retained for manuscript provenance; external inputs are not committed
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
- WF10 is a folder-local collaborator workflow, not a downstream dependency of
  workflows 00-09.
- `scripts/james/maldi/Cleland_analysis_expt2_recal.ipynb` and
  `scripts/james/proteomics/JamesCleland_P3225_V1.Rmd` are collaborator
  provenance materials only. They are not sourced by active workflows.

## Non-Workflow Collaborator Material

| Bundle item | Repository location | Notes |
| --- | --- | --- |
| Attachment 1, MALDI | `scripts/james/maldi/Cleland_analysis_expt2_recal.ipynb` | Negative-mode MALDI-MSI liver zonation notebook. External AnnData, imaging, registration, and SCiLS inputs are not committed. |
| Attachments 2-5, Xenium | `scripts/workflows/10_xenium_analysis/00_infer_doublets.py`; `01_merge_doublets.py`; `02_compositional_analysis.py`; `03_DGE_analysis.py` | Already present in workflow 10. The 2026-06-19 archive copies matched the checked-in scripts after ignoring CRLF line-ending differences. |
| Attachment 6, proteomics | `scripts/james/proteomics/JamesCleland_P3225_V1.Rmd` | Original P3225 fractionated TMT proteomics R Markdown analysis. Workflow 02 uses curated checked-in inputs under `data/proteomics/`, not this Rmd directly. |

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
| External Xenium output folder | `scripts/workflows/10_xenium_analysis/data_raw/output-XETG00078__0018188__Region_1__20240301__143248/` | 10 |
| External Xenium panel workbook | `scripts/workflows/10_xenium_analysis/data_raw/XeniumTargetList_Mouse-Panel-1.xlsx` | 10 |
| 10 | `scripts/workflows/10_xenium_analysis/inferred_components.json` | 10 internal, `01_merge_doublets.py` |
| 10 | `scripts/workflows/10_xenium_analysis/fused_anndata_annotated.h5ad` | 10 internal, `02_compositional_analysis.py` and `03_DGE_analysis.py` |
| 10 | `scripts/workflows/10_xenium_analysis/data_tidy/*.csv`; `scripts/workflows/10_xenium_analysis/data_tidy/*.xlsx`; `scripts/workflows/10_xenium_analysis/figures/*.pdf` | Folder-local result lookup when generated |

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
| 10 | IMPORTED | DOCUMENTED | Optional collaborator Xenium folder; notebook-style scripts with raw-data and working-directory assumptions. |

## Cross-Workflow Vulnerabilities

1. **WF05 pathway cache:** Section 7 loads `results/zonation/pathway_enrichment.RData`; the active ULM run is commented out. A fresh environment without this file skips pathway reporting.
2. **WF05 threshold variables:** `sensitivity_df`, `proportion_total`, and `optimal_metrics` are created in section 11 and consumed much later in section 14.
3. **WF05 -> WF07 network names:** WF07 expects exact mixed names from WF05: `cormat_testes_portal_full_SIF.csv`, `cormat_testes_central_full_SIF.csv`, `cormat_Ovaries_portal_full_SIF.csv`, and `cormat_ova_central_full_SIF.csv`.
4. **WF06 -> WF07 structural feature table:** WF07 section 11 expects `results/lipidomic/structural_model/lipids_with_structure.csv`.
5. **WF03 quantitative count drift:** The WF03 script header and final printed text contain different headline lipid counts. Use generated result files for exact quantitative lookups.
6. **WF09 external-data paths:** The script uses absolute collaborator/HPC `DATA_DIR` and `RESULTS_DIR` values. The large human liver inputs are not committed, and `results/pseudobulk/hepa_gonad_signature.csv` must be copied to the expected external `DATA_DIR` path before running without code changes.
7. **WF10 folder-local paths:** Scripts use relative paths such as `data_raw/`, `data_tidy/`, and `figures/`; run from `scripts/workflows/10_xenium_analysis/` or adjust paths.
8. **WF10 notebook-style code:** The DGE script contains repeated helper definitions and assumes interactive/script state in places. Treat it as collaborator-imported analysis code rather than a fully refactored pipeline.

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

### Experimental mouse COSMOS PKN (2026-09-08)

`scripts/pkn/run_build.py` is a separate, optional prior-construction pipeline,
documented in `docs/cosmos_mouse_pkn_build.md`. It reads frozen external network
resources, the current WF03 lipid DE feature IDs and WF04 Factor 1 lipid IDs for
coverage/rule matching. It writes `results/cosmos_mouse_pkn/`, including signed
R/CSV priors, rich reaction/evidence tables, mapping candidates and validation.
It neither changes the existing workflows nor uses their correlation/abundance
signs to create causal edges. The optional Factor 1 MOON pilot below consumes
the measured-species rules prior; workflows 00–10 remain unchanged.

### Optional Factor 1 MOON pilot (2026-09-08)

`scripts/moon/prepare_factor1_inputs.R`, `audit_factor1_model.py` and
`map_factor1_inputs.R` read the WF04 saved loading/TF-score bundle, selected
MOFA HDF5 model, sample metadata, existing TF regulons and PKN mappings. They
write audited, jointly oriented inputs under `results/cosmos_factor1_pilot/inputs/`.
`scripts/moon/run_factor1_moon.R` consumes these inputs and
`results/cosmos_mouse_pkn/pkn_lipid_rules.rds` for six-layer ULM scoring,
compression, iterative RNA coherence and signed TF-to-lipid extraction.
`audit_factor1_thresholds.R` reuses the final scores for extraction sensitivity;
`export_factor1_network.py` adds prior evidence, candidate paths and figures.
The main result is documented in `docs/cosmos_factor1_moon_pilot.md`: the
tutorial extraction thresholds yield no complete paths with compression, while
an independently converged uncompressed run retains 8 TFs and 4 lipids. Relaxed
and stricter settings, compression effects and bounded catalyst routes are
audited separately. These outputs are not consumed by existing
workflows or treated as validated causal mechanisms.

### Optional Factor 1 MOON parameter search (2026-09-09)

`scripts/moon/search_factor1_parameter_space.R` extends the pilot without
altering workflows 00--10, the saved MOFA fit, PKN, or input measurements. It
scores fixed-radius terminal-observation PKNs across propagation depths 4--12,
separately for compressed and uncompressed representations where appropriate,
then reuses each saved score state for a coarse and local-refined dual-threshold
grid. `audit_factor1_depth_topology.R`,
`summarize_factor1_parameter_search.py`, `audit_refined_candidates.R` and
`export_factor1_network.py` write topology coverage, all threshold rows,
independent edge/sign/connectivity checks, GraphML and reaction-evidence
exports under `results/cosmos_factor1_parameter_search/`.

The original uncompressed U5/U4 preference is superseded by the compressed
review below. These saved runs remain representation diagnostics, since higher
uncompressed scores do not establish stronger evidence than compressed scores.
`docs/cosmos_factor1_moon_parameter_search.md` preserves the historical search
and exact artifacts. These outputs are not consumed by existing workflows.

### Optional compressed MOON review and TF audit (2026-09-09)

`docs/cosmos_factor1_compressed_review.md` is the current interpretation entry
point. `search_compressed_threshold_review.R` evaluates 7,563 additional cutoff
pairs on saved baseline compressed scores and eight new compressed sensitivity
states, writing to `results/cosmos_factor1_compressed_review/`. The new states
come from `run_factor1_compressed_sensitivities.R` (radii and uniform initial
RNA coherence) and `run_factor1_inputs_only_compression.R` (stronger compression
with original-ID coherence and recompression after deletions). The primary
exploratory C03 has 388 edges at radius 6/depth 4 and cutoffs 1.85/0.55, with
Ahr/Hnf4a/Ppara/Tbx2 and three DGs. Stronger compression reproduces its exact
graph. A wider radius-12/depth-8 compressed alternative has 468 edges at
2.15/1.075 with Klf10/Tbx2 and two DGs. Cutoffs across changed scoring
backgrounds are not ranks of biological confidence.

`audit_factor1_focus_tf_retention.R` and
`validate_factor1_compressed_tf_sensitivities.R` establish that Esr1/Nfic/Dbp
survive preprocessing and compression but have weak positive inferred scores
against negative supplied TF activities; MOON removes them after scoring.
All 66 audited TF/run comparisons retain this sign conflict.
`export_compressed_score_retention.R` saves all pre-sign-filter score rows
with explicit return/exclusion and candidate-membership flags. Exact PKN
projection, coherence, sign and connectivity validation, reaction provenance
and figures accompany the networks. No workflow 00--10 dependency, accepted
input mapping, MOFA fit or raw PKN is changed. These are mechanistic hypotheses
with narrow lipid-rule coverage, not validated causal mechanisms.

### Manuscript drafting and relocated causal-model archive (2026-09-09)

The user moved the COSMOS result folders into `results/causal_modelling/`.
`docs/manuscript_drafts/cosmos_results_methods.md` drafts the Results and
Methods additions for the supplied revised paper; its companion notes map
scripts and current artifacts and distinguish curated PKN edges from species
hypotheses. `scripts/moon/audit_manuscript_causal_model.R` was executed against
the relocated archive. It verifies input hashes and the final compression
projection and re-extracts exactly the selected388-edge/153-node C03 network
from frozen scores, without refitting MOFA or MOON. Its output is
`results/causal_modelling/manuscript_audit/`. Historical generation scripts and
embedded RDS/manifest paths retain their old locations; full source replay
requires path migration and reconciliation of recorded script versions. No
workflow00--10 dependency or scientific result was changed by this drafting
and audit step.

Use `docs/AGENT_GUIDE.md` as the static replacement for the discontinued MCP
query layer. It maps common biological/data questions to result files that
answer them directly.

Keep this meta overview and `docs/AGENT_GUIDE.md` synchronized whenever a
workflow starts writing a new machine-readable result or changes a dependency.
