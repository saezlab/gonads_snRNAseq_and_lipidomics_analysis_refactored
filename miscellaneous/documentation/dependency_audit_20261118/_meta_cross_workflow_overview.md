## CSV input distribution update — 18 November 2026

Public RNA inputs now use the seven registry-selected manifests under `data/RNA/pseudobulk_csv/`; every CSV part is strictly below 50,000,000 bytes. RNA/Seurat RDS objects are local provenance only. Current scientific values and historical input roles are unchanged. Read [the CSV input contract](../PSEUDOBULK_CSV_INPUTS.md); older RDS path descriptions below identify original production sources, not files to commit.

# Current release and archived references — 16 November 2026

The expanded-marker default is selected by `scripts/config/analysis_release.json`; main zonation/diet entry is Workflow17. Current inputs are `data/RNA/zonation/current/`, and the approved outputs remain `results/zonation_marker_update_20261116/`. Workflow17 defaults to scientific execution; manuscript auditing is opt-in.

See [the maintained execution graph](../../PIPELINE.md), [the organization map](../REPOSITORY_ORGANIZATION.md) and [the causal dependency map](../CAUSAL_ANALYSIS_MAP.md). Old `results/zonation` and old zone input paths in the historical details below resolve to `miscellaneous` through the release registry; they must not be interpreted as current defaults. Former full CF45314/15/16 presentations and WF12/13 comparison outputs are archived. Current14/15 model stages feed16.

Original05/07 bodies are preserved under `scripts/core/engines/`; the maintained wrappers/adapters distinguish current scientific execution from frozen replay. Main non-zonal workflows and independent measurements retain their documented roles. Exact moves/hashes are in `docs/provenance/`.

The following older audit details are retained as design/provenance context; for current file locations the registry and updated individual plans take precedence. A pre-cleanup snapshot is also retained under `miscellaneous/documentation/pre_cleanup_20261116/`.

---

# Cross-Workflow Dependencies and Coherence Report

**Last audited:** 2026-06-19
**Cyp2c38 extensions added:** 2026-10-06 (07b), 2026-10-07 (07c).
**WDSD extension added:** 2026-10-07 (05b, meeting WP2).
**Pooled Figure 3 / lipid annotation follow-up:** 2026-10-07 (05c).
**WP1/WP2 consolidation and local-input audit:** 2026-10-09 (active 11/12; former entries are phase wrappers).
**WP3 added:** 2026-10-14 (Workflow 13, male WD cross-cohort robustness).
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

11_cyp2c38_in_vivo_validation (WP1)
  original WF01 inputs + WF02 checked groups + WF04 fixed network
  + data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv
  -> base KO intake/models/historical concordance
  -> newly written base outputs + original paired FCG inputs
  -> full/historical sex networks, native ULM, null bootstrap and calibration
  outputs: results/validation/cyp2c38_in_vivo/ and ulm_followup/
  former 07b/07c entries remain thin phase compatibility wrappers

12_wdsd_network_comparison (WP2)
  WF04 frozen Figure 3 branch + WF01 whole-hepatocyte matrices/groups
  + local WDSD whole-hepatocyte/lipid inputs + raw FCG annotation export
  -> pooled conservation, sex-conditioned tests and fresh annotation audit
  -> followup/ outputs first
  WF05 frozen Figure 5 branches + original zonation/lipid matrices/groups
  + local WDSD paired zone/lipid inputs
  -> four sex-zone SD comparisons and descriptive WD-SD/LOCO analyses
  -> combined pooled-first report, vector PDF built separately
  outputs: results/wdsd_network_comparison/ and followup/
  former 05c/05b entries remain thin phase compatibility wrappers

13_wd_crosscohort_robustness (WP3)
  WF05 frozen male Figure 5 genes/full/displayed pairs + original reference replay
  + local CF453 WDSD inputs and local CF459 WD-only inputs
  -> male WD CF453 n7 versus CF459 n4, separate periportal/pericentral
  + original paired FCG SD n22 (11 XXT/11 XYT), exact archive replay
  -> matched-edge CF459 closeness to FCG versus CF453; stratified three-cohort stability
  -> strict shared measurements; exact conditional pairing/null bootstrap
     and paired/39-omission stability; report/scientific figures
  optional WF12 numerical WD replay; no saved WP3 output needed for intake
  outputs: results/wd_crosscohort_robustness/

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
- WF11 consolidates the independent in vivo Cyp2c38 WP1. It reads local KO
  lipidomics plus original WF01 inputs, WF02 checked group metadata and WF04's
  fixed network. It passes newly created base results to its ULM/calibration
  phase; `--validate-only` passes fresh intake and needs no prior WP1 output.
  No MOFA refit or dependence on WF06 structural outputs is involved.
- WF12 runs pooled Figure 3 first, followed by the four Figure 5 sex-zone
  networks and diet comparisons. Both phases read local inputs under `data/`.
  Pooled analysis replays the frozen sparse 192-pair branch and 139-by-64
  precursor, not the permissive 33%-tail network. It builds its annotation
  audit directly from raw/cleaned inputs, with no earlier WP2 output dependency.
  Zonal analysis replays original log2 lipid coefficients and displayed selection.
- WF11/12 retain all original statistical helper calculations and separate
  output directories. Full runs require fresh destinations; WP12 `--report-only`
  intentionally refreshes combined presentation. PDF builders are separate.
  They do not rerun WF04/05/07, train MOFA, or consume CF459 WP3 inputs.
- Former 07b/07c and 05c/05b scripts are compatibility phase wrappers; their
  phase plans retain original production results and historical report provenance.
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
| 05 | `results/zonation/networks/cormat_*_full_SIF.csv` | 07; 12 zonation frozen reference |
| 05 + original matrices/groups + local WDSD inputs | Four frozen displayed/full networks; zonation RNA, log2 lipids and FCG groups | 12 sex-zone SD conservation and diet comparison |
| 12 zonation phase | `results/wdsd_network_comparison/` | Native ULM, exact pairing and null-bootstrap tests, SD/WD coefficients, LOCO/QC, Figure 6 and combined report |
| 04 + 01 + local WDSD inputs/raw annotations | Frozen Figure 3 species pairs and whole-hepatocyte/log2 lipid matrices | 12 pooled conservation, annotation audit and conditional sensitivity |
| 12 pooled phase | `results/wdsd_network_comparison/followup/` | Strict/conditional ledgers; 108 strict displayed pairs; pooled and sex-conditioned nulls, centered sensitivities and component report |
| Imported local input | `data/metabolomics/wdsd/lipidomics_cf453.csv` | 12 both phases |
| Imported local input | `data/RNA/wdsd/wdsd_{hep,zonation}_pseudobulks.rds` | 12 pooled / zonation phases respectively |
| Imported local input | `data/metabolomics/fcg/lipidomics_fcg.csv` | 12 raw annotation audit; supplements cleaned matrix |
| Imported optional reference | `docs/reference_figures/FCGL_Fig3.ai` | Figure 3 selection/provenance, not required for numerical replay |
| 06 | `results/lipidomic/structural_model/lipids_with_structure.csv` | 07 |
| 06 | `results/lipidomic/structural_model/lipid_only_*effects.csv` | Result lookup, when files exist |
| 06 | `results/lipidomic/structural_model/lipid_only_class_{coefficients,summaries}.csv` | Result lookup |
| 07 | `results/validation/validation_summary.csv` | Result lookup |
| 01 + original input | `data/RNA/scaled_RNA_{1,2,3}.csv`; RNA metadata; cleaned original lipid matrix | 11 fixed-reference replay and sex-specific recomputation |
| 02 | `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | 11 FCG group/allosome metadata, checked against overview workbook |
| 04 | `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | 11 fixed Cyp2c38 universe |
| Imported local input | `data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv` | 11 base KO phase |
| 11 base phase | `results/validation/cyp2c38_in_vivo/` | 11 ULM phase and historical result lookup; separate from original WF07 |
| 11 ULM phase | `results/validation/cyp2c38_in_vivo/ulm_followup/` | Native ULM, six reference networks, null draws and calibration/injection diagnostics |
| Imported local CF459 inputs | `data/metabolomics/wd/lipidomics_wd.csv`; `data/RNA/wd/wd_{,zonation_}pseudobulks.rds` | 13 WD-only male cross-cohort comparison; cell-type object for intake QC only |
| 05 + local CF453/CF459 inputs | Frozen male original Figure 5 genes/pairs and paired zone/lipid matrices | 13 full primary/displayed secondary networks |
| 12 optional numerical replay | `results/wdsd_network_comparison/network_edges.csv` | 13 independently reconstructed CF453 male WD coefficients checked where eligible universes overlap |
| 13 | `results/wd_crosscohort_robustness/` | Correlation/ULM and exact conditional nulls, bootstrap/omission ledgers, QC, sensitivities and report |
| Import provenance | `data/WP1_WP2_INPUTS.md`; `data/wp1_wp2_input_migration_manifest.csv` | Local input roles, archive sources and six verified byte/SHA256 copies |
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
| 05b | COMPATIBILITY | DOCUMENTED | Thin wrapper for WF12 zonation/diet phase; original scientific specification/results retained. |
| 05c | COMPATIBILITY | DOCUMENTED | Thin wrapper for WF12 pooled/annotation phase; original scientific specification/results retained. |
| 06 | GOOD | GOOD | Active path is lipid-only structural analysis under `results/lipidomic/structural_model/`. |
| 07 | GOOD | GOOD | Includes structural consistency analysis dependent on WF06. |
| 07b | COMPATIBILITY | DOCUMENTED | Thin wrapper for WF11 base phase; original scientific specification/results retained. |
| 07c | COMPATIBILITY | DOCUMENTED | Thin wrapper for WF11 ULM/calibration phase; original scientific specification/results retained. |
| 08 | GOOD | GOOD | Standalone lipid variance workflow with current `variance_analysis_*` outputs. |
| 09 | GOOD | GOOD | Optional Python human comparison; depends on external GSE179250 files and writes outside the repo by default. |
| 10 | IMPORTED | DOCUMENTED | Optional collaborator Xenium folder; notebook-style scripts with raw-data and working-directory assumptions. |
| 11 | VERIFIED | DOCUMENTED | Complete WP1; local intake, fresh validation, full two-phase smoke and observed-result replay passed; production results preserved. |
| 12 | VERIFIED | DOCUMENTED | Complete pooled-first WP2 default run and report passed; 12 tables within 1e-10 and 16 resampling/index banks identical; local inputs. |
| 13 | VERIFIED | DOCUMENTED | Complete male WD cross-cohort run; strict 610 lipids, independent 500 checks within 6.41e-14, weak observed agreement and substantial n=4 bootstrap degeneracy. |

WP3 uses strict full/displayed fixed-original universes but defines eligibility
from the two male WD cohorts only, independently of CF453 SD variability. It
shares each resampling/omission decision across zones/features. CF459 has no
analytical QC or SD arm; metadata `m-ctrl/f-ctrl` is WD-only cohort provenance.
The 24- and 5,040-pairing tests condition on different fixed observed networks;
paired bootstrap and omissions are descriptive. See the dedicated plan for
implemented tests and numerical results.

## WP3 production verification, 14 October 2026

Workflow 13 completed with all 24/5,040 pairings, 4,999 null draws, 999 paired draws,
39 omission configurations and three full-count sensitivities. All 610 shared
strict original lipids yield 15,250 periportal and 18,300 pericentral full pairs;
full coefficient-vector r = -0.06310/-0.01857. The complete result matches an
independent default-count pilot exactly. Independent actual replay passes 500
checks, maximum error 6.40e-14. Bootstrap undefined vectors are explicitly
reported, particularly pericentral genes supported by one distinguishing CF459
mouse; no dynamic deletion is allowed. All 239 existing WP1/WP2 files retain
identical SHA256. Final PDF/content and numerical QA live in the new result
folder. Read the dedicated plan/report for interpretation and limitations.

## Work package consolidation verification, 9 October 2026

Six required input/provenance files (30,965,755 bytes) were copied unchanged
into the repository and verified against source SHA256; shared archives were
preserved. Both actual `--validate-only` commands work without saved WP1/WP2
outputs and create no files. All four legacy wrappers passed actual validation.
CLI routing, phase order, source loading and nonempty-output guards passed
independent checks. The WF11 complete reduced-draw smoke reproduced observed
models and all ten native ULM rows; reduced Monte Carlo results stayed temporary.
The WF12 complete fresh default run reproduced all 12 checked numerical tables
and exactly matched 16 primary resampling/index/LOCO banks. Evidence is under
`tests/outputs/wp1_wp2_migration/`; permanent actual replay checks are in
`tests/test_work_package_workflows.R`.

Primary statistical output/provenance files remain historical records. New
runs record current local inputs/code and package phase manifests. The combined
WP2 Markdown/PDF preserves all table values and original scientific vector
pages, with pooled comparison before zonation/diet. PDF QA is recorded in
`results/wdsd_network_comparison/QA/`.

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

### Active selected TF-to-lipid MOON workflow (2026-09-15)

`scripts/moon/run_factor1_tf_to_lipid.R` is the current reproduction entry point,
with settings in `scripts/moon/config/factor1_tf_to_lipid.json`, shared operations
in `scripts/moon/lib/workflow.R`, and the unchanged safeguard implementation in
`scripts/moon/moon_safe_helpers.R`. The dedicated plan is
`docs/workflow_plans/moon_factor1_tf_to_lipid.md`.

Dependency flow: WF04 saved Factor 1 outputs and frozen mouse PKN -> previously
audited mapped input bundle -> selected compressed scoring/refit ->
`results/causal_modelling/factor1_tf_to_lipid/`. The main workflow starts from
`results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds` and
`results/causal_modelling/cosmos_mouse_pkn/pkn_lipid_rules.rds`; it neither refits
MOFA nor reconstructs the PKN or input mapping. Workflows 00–10 do not consume
these optional outputs.

The full main workflow passed on 2026-09-15: pruning 6, scoring 4 and thresholds
1.85/0.55 reproduced the frozen C03 signed edge/node sets, compression mapping,
score IDs/levels and numerical scores within 1e-10. It converged in three
coherence iterations, produced 388 edges/153 nodes, reproduced Esr1/Nfic/Dbp
sign-filter exclusions and exported evidence/GraphML/figures. Run provenance,
reference comparison and exact input/output hashes accompany the regenerated
outputs. `--validate-only` checks setup without writing files; reruns require a
new `--output` directory because nonempty destinations are refused.

`scripts/moon/ancillary/search_factor1_parameters.R` reuses the same compressed
core for explicit pruning, scoring-layer and threshold grids, recording every
combination. Default pruning is 6 and layers are 4,5,6,8. This generic ancillary
runner does not reproduce the historical adaptive 7,563-pair search by default.
The lipid-to-TF analysis has a separate entry point described next.

### Active lipid-to-TF exploratory MOON workflow (2026-09-15)

`scripts/moon/run_factor1_lipid_to_tf.R` consumes the same frozen Factor 1 input
bundle, with lipid loadings upstream and TF footprints downstream. It uses the
same compression/scoring/coherence helpers without changing the forward result.
Its configuration and dedicated plan are
`scripts/moon/config/factor1_lipid_to_tf.json` and
`docs/workflow_plans/moon_factor1_lipid_to_tf.md`.

Dependency flow: frozen Factor 1 inputs plus original/normalized mouse PKN and
mapping evidence -> coverage/shortcut audits -> direction-specific compressed
fits -> ancillary threshold refinement, score diagnostics and provenance export
under `results/causal_modelling/factor1_lipid_to_tf/`. Neither MOFA nor the PKN is
rebuilt. Workflows 00–10 and the selected TF-to-lipid result do not consume these
outputs. The original PKN's lack of measured-lipid regulatory routes is retained
as a result; the acid/base-normalized precursor graph is an explicitly optional
sensitivity. Detailed outcomes and reproduction commands are in
`docs/cosmos_factor1_lipid_to_tf.md`.

### Active relaxed regulatory lipid-to-TF candidate workflow (2026-09-16)

`scripts/moon/run_factor1_lipid_to_tf_relaxed.R` is an explicitly separate
exploratory analysis using the directionality-audited STITCH overlay and the
long-form relaxed chemical mapping ledger. Its configuration is
`scripts/moon/config/factor1_lipid_to_tf_relaxed.json`; its default output is
`results/causal_modelling/factor1_lipid_to_tf_relaxed/`. It uses only direct
chemical→protein and protein→protein edges, excluding reaction traversal.
Candidate chemistry is represented as topology seeds, while Factor 1 lipid
weights are retained after scoring as assignment-level sign-concordance evidence
because MOON's numeric ULM score is calculated from terminal TF activities.

The ancillary parameter search is
`scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R`, with
`scripts/moon/config/factor1_lipid_to_tf_relaxed_parameter_search.json` and
default output `results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search/`.
Candidate-local score traces and direct-edge fingerprints are audited by
`ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R` under the selected
analysis output. Score-layer evidence for every materialized strict search
solution is written under `<parameter-search>/evidence_audit/` by
`ancillary/audit_factor1_lipid_to_tf_relaxed_score_layer_evidence.py` and
selected display-network exports are handled by
`ancillary/export_factor1_lipid_to_tf_relaxed_network.py`.

The prior `results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory/`
folder is a legacy pre-cleanup discovery topology/search, not the active result.
The active main workflow refits and compares the low-threshold
LPS→Gpr34→Gnai1→Src→{Ar, Esr1} sensitivity with that legacy reference, but no
compact score-layer solution survived either threshold above 0.5. AR and ESR1
are supplied terminal observations. A broad sign-consistent connectivity
extraction remains topology only, not a recursive-MOON score explanation. No
workflow 00–10, manuscript TF→lipid result, frozen PKN or original lipid→TF
analysis consumes this screen. See
`docs/cosmos_factor1_lipid_to_tf_relaxed.md` for the result and caveats.

The sections below describe historical experiments. Their scripts now live in
`scripts/moon/archive/2026-09/` with a SHA256/path migration manifest and role
index; their bytes and saved results were preserved. Historical `results/cosmos_*`
paths refer to the corresponding `results/causal_modelling/cosmos_*` folders
unless quoted as original execution/provenance strings. Archived scripts retain
those original path assumptions and are not current runnable commands.

### Experimental mouse COSMOS PKN (2026-09-08)

`scripts/pkn/run_build.py` is a separate, optional prior-construction pipeline,
documented in `docs/cosmos_mouse_pkn_build.md`. It reads frozen external network
resources, the current WF03 lipid DE feature IDs and WF04 Factor 1 lipid IDs for
coverage/rule matching. It writes `results/cosmos_mouse_pkn/`, including signed
R/CSV priors, rich reaction/evidence tables, mapping candidates and validation.
It neither changes the existing workflows nor uses their correlation/abundance
signs to create causal edges. The optional Factor 1 MOON pilot below consumes
the measured-species rules prior; workflows 00–10 remain unchanged.

### Historical Factor 1 MOON pilot (2026-09-08)

`scripts/moon/archive/2026-09/prepare_factor1_inputs.R`, `scripts/moon/archive/2026-09/audit_factor1_model.py` and
`scripts/moon/archive/2026-09/map_factor1_inputs.R` read the WF04 saved loading/TF-score bundle, selected
MOFA HDF5 model, sample metadata, existing TF regulons and PKN mappings. They
write audited, jointly oriented inputs under `results/cosmos_factor1_pilot/inputs/`.
`scripts/moon/archive/2026-09/run_factor1_moon.R` consumes these inputs and
`results/cosmos_mouse_pkn/pkn_lipid_rules.rds` for six-layer ULM scoring,
compression, iterative RNA coherence and signed TF-to-lipid extraction.
`scripts/moon/archive/2026-09/audit_factor1_thresholds.R` reuses the final scores for extraction sensitivity;
`scripts/moon/archive/2026-09/export_factor1_network.py` adds prior evidence, candidate paths and figures.
The main result is documented in `docs/cosmos_factor1_moon_pilot.md`: the
tutorial extraction thresholds yield no complete paths with compression, while
an independently converged uncompressed run retains 8 TFs and 4 lipids. Relaxed
and stricter settings, compression effects and bounded catalyst routes are
audited separately. These outputs are not consumed by existing
workflows or treated as validated causal mechanisms.

### Historical Factor 1 MOON parameter search (2026-09-09)

`scripts/moon/archive/2026-09/search_factor1_parameter_space.R` extends the pilot without
altering workflows 00--10, the saved MOFA fit, PKN, or input measurements. It
scores fixed-radius terminal-observation PKNs across propagation depths 4--12,
separately for compressed and uncompressed representations where appropriate,
then reuses each saved score state for a coarse and local-refined dual-threshold
grid. `scripts/moon/archive/2026-09/audit_factor1_depth_topology.R`,
`scripts/moon/archive/2026-09/summarize_factor1_parameter_search.py`, `scripts/moon/archive/2026-09/audit_refined_candidates.R` and
`scripts/moon/archive/2026-09/export_factor1_network.py` write topology coverage, all threshold rows,
independent edge/sign/connectivity checks, GraphML and reaction-evidence
exports under `results/cosmos_factor1_parameter_search/`.

The original uncompressed U5/U4 preference is superseded by the compressed
review below. These saved runs remain representation diagnostics, since higher
uncompressed scores do not establish stronger evidence than compressed scores.
`docs/cosmos_factor1_moon_parameter_search.md` preserves the historical search
and exact artifacts. These outputs are not consumed by existing workflows.

### Historical compressed MOON review and TF audit (2026-09-09)

`docs/cosmos_factor1_compressed_review.md` records the historical selection
evidence for the current manuscript configuration. `scripts/moon/archive/2026-09/search_compressed_threshold_review.R` evaluates 7,563 additional cutoff
pairs on saved baseline compressed scores and eight new compressed sensitivity
states, writing to `results/cosmos_factor1_compressed_review/`. The new states
come from `scripts/moon/archive/2026-09/run_factor1_compressed_sensitivities.R` (radii and uniform initial
RNA coherence) and `scripts/moon/archive/2026-09/run_factor1_inputs_only_compression.R` (stronger compression
with original-ID coherence and recompression after deletions). The primary
exploratory C03 has 388 edges at radius 6/depth 4 and cutoffs 1.85/0.55, with
Ahr/Hnf4a/Ppara/Tbx2 and three DGs. Stronger compression reproduces its exact
graph. A wider radius-12/depth-8 compressed alternative has 468 edges at
2.15/1.075 with Klf10/Tbx2 and two DGs. Cutoffs across changed scoring
backgrounds are not ranks of biological confidence.

`scripts/moon/archive/2026-09/audit_factor1_focus_tf_retention.R` and
`scripts/moon/archive/2026-09/validate_factor1_compressed_tf_sensitivities.R` establish that Esr1/Nfic/Dbp
survive preprocessing and compression but have weak positive inferred scores
against negative supplied TF activities; MOON removes them after scoring.
All 66 audited TF/run comparisons retain this sign conflict.
`scripts/moon/archive/2026-09/export_compressed_score_retention.R` saves all pre-sign-filter score rows
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
hypotheses. `scripts/moon/archive/2026-09/audit_manuscript_causal_model.R` was executed against
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


## WP3 original-cohort extension, 14 October 2026

The user expanded WP3 to compare CF459 with the original male FCG SD network
as well as CF453 WD, on precisely the same frozen 610-lipid full/displayed
universe. Ages clarified by the user are CF459/FCG 3 months and CF453 10 months.
The extension uses the original paired FCG matrices (11 XXT + 11 XYT) as a
measured anchor, not only as feature provenance. No original FCG/WP1/WP2
output changes. Native ULM and conditional null checks accompany all three
pairwise correlations and the direct CF459 closeness difference. Paired
three-cohort stability reuses CF459 simultaneously for both anchors and resamples
FCG within chromosome strata. Extension data and an exact-zero null-draw audit
are isolated in `results/wd_crosscohort_robustness/fcg_comparison/`.
Read the Workflow 13 plan and current report for verification and limitations.

## Standalone CF453 Figure 6, 30 October 2026

Workflow 14 reads only the three local CF453 files through the shared WDSD
intake helper. The canonical-lipid formatting function is reused without any
original-cohort input loading. There is no data dependency on workflows 02-13
or their frozen results. Inputs -> paired sex/diet/zone RNA model and mouse-level
sex/diet lipid model -> CF453 RNA-selected gene candidates -> eight within-group
association contexts -> independent panel files and source tables under
`results/cf453_standalone_figure6/`. No downstream workflow currently consumes
the network outputs. Workflow 15 consumes its validated models.rds only.
Read `14_cf453_standalone_figure6.md` for selection, FDR families,
units, QC universes, paired resampling and interpretation limits.

## Compact direct CF453 diet effects, 2 November 2026

Workflow 15 reads Workflow 14's validated RNA/lipid model-input cache and checks
its three CF453 raw-input hashes. It refits the eight-cell paired RNA model for
direct diet x sex, diet x zone and three-way interactions, with new declared
FDR families and a covariance-aware omnibus F test. It reuses the direct lipid
sex x diet model. No network selection or association coefficient is consumed.
Four independent panel groups under `results/cf453_diet_impact_figure6/` replace
the larger set for manuscript presentation, preserving all historical outputs.
Read `15_cf453_diet_impact_figure6.md` for precise weights and verification.

## Selected zonation programmes and diet associations, 2 November 2026

Workflow 16 consumes validated Workflow 14 model inputs/zonation statistics
and Workflow 15 direct three-way evidence. Its reciprocal 20-gene anchor
selections use RNA alone. Fixed unions of those selected genes support the
relative lipid-class association and gained/lost example panels. No original
FCG input or frozen network is loaded. New coefficients, mouse bootstrap
indices, omission stability, class/rest rates, observed sensitivities and
three independent panel groups live under
`results/cf453_diet_rewiring_figure6/`. Earlier outputs are preserved and no
downstream workflow currently consumes these results. Read
`16_cf453_diet_rewiring_figure6.md` for thresholds, selection and uncertainty.

The4November2026 Figure6A layout comparison reads only the saved Workflow16
anchors and writes two isolated plot-only candidates to
`layout_comparison_20261104/` within that output root. It deduplicates40
selections to34 unique genes for display, preserving SD/WD/both membership
and the original direct three-way evidence. No downstream analysis consumes
these candidate exports; original A/B/C and all model/network results are preserved.

## Expanded marker update 16 November2026

Workflow17 stages the supplied score-based FCG/GDX/CF453 zonation pseudobulks
under `results/zonation_marker_update_20261116/inputs/`. Isolated FCG05 outputs
feed isolated zonation07 validation. Updated CF453 inputs feed14models-only,
15models-only,16 and both6A layouts; the repaired GDX branch is supporting.
The report maps panels and89 manuscript claims plus24 pictured network edges
to old/new numbers and source files. Figure4A/S4B-C are deferred to James.
All old data/results and independent measurements are preserved. WF12 zonal
references and WF13's CF459 compartments remain historical; they are not mixed
with the new definitions. Read the Workflow17 plan for adapters, QA and entry.
