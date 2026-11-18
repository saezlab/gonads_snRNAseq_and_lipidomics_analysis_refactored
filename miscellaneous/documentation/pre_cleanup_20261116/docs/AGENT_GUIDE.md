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

# Optional collaborator Xenium folder; run from scripts/workflows/10_xenium_analysis/
# python 00_infer_doublets.py
# python 01_merge_doublets.py
# python 02_compositional_analysis.py
# python 03_DGE_analysis.py
```

Run from the repository root for workflows 00-08. After editing an R script,
run that script immediately or record why it could not be run. Workflow 09 uses
absolute external data/output paths in the script and is not part of the default
R rerun path. Workflow 10 uses folder-local relative paths and raw Xenium inputs
that are not committed.

Unrefactored collaborator source material that is not part of the default rerun
path is kept under `scripts/james/`. See `scripts/james/README.md` before using
those files.

## October 2026 Work Packages

The current individual entry points are Workflow 11 (WP1, KO contrasts then
ULM/calibration) and Workflow 12 (WP2, pooled Figure 3 then zonation/diet). Workflow 13 (WP3)
compares male CF459 WD with CF453 WD and the original male FCG SD network
within each zone, including a direct matched-edge difference in agreement.
Use `Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --validate-only`
and `Rscript scripts/workflows/12_wdsd_network_comparison.R --validate-only`
for checks against local inputs without requiring saved work-package outputs.
Full reruns require fresh `--output` destinations. Workflow 12 `--report-only`
rebuilds the pooled-first summary from completed outputs. Its PDF builder is
`scripts/reports/build_wp2_report_pdf.py`.

[data/WP1_WP2_INPUTS.md](../data/WP1_WP2_INPUTS.md) and the SHA256 migration
manifest identify the imported data. Historical phase scripts 07b/07c/05b/05c
are thin compatibility wrappers. Historical recorded provenance is preserved;
new runs record the local inputs and new entry/core files.

[WP3 local inputs](../data/WP3_INPUTS.md) and its import manifest document the
CF459 files. Workflow 13 raw validation requires no saved WP3 results; the
existing WP2 WD coefficients provide an optional numerical replay check.

## Standalone CF453 Figure 6 (30 October 2026)

Workflow 14 uses CF453 alone for paired sex x diet x zone RNA inference and
within-sex/diet RNA-lipid associations. It does not transfer original selected
genes or networks. Read `docs/workflow_plans/14_cf453_standalone_figure6.md`.
Run `Rscript scripts/workflows/14_cf453_standalone_figure6.R --output=results/cf453_standalone_figure6_reproduction`
from the root. The destination must not contain saved models or analysis.
`--plot-only` refreshes independent PDF/SVG/PNG panels from saved analysis;
`--models-only` and `--finish` support two-stage execution. Outputs include
`analysis.rds`, full statistical tables, paired bootstrap indices, per-panel
source data, legends and independent QA. Association displays are exploratory
at seven mice per sex/diet, and lipids are whole-liver measurements.

## Compact direct diet-impact Figure 6 (2 November 2026)

Workflow 15 is the earlier four-panel direct-interaction presentation. It directly
tests diet changes in sex, zonation and sex-dependent zonation, plus sex-specific
lipid diet response. Read `docs/workflow_plans/15_cf453_diet_impact_figure6.md`.
`Rscript scripts/workflows/15_cf453_diet_impact_figure6.R` uses validated
Workflow 14 models.rds and writes a fresh `results/cf453_diet_impact_figure6/`.
`--plot-only` regenerates exactly four PDF/SVG/PNG groups from saved results.
The new pooled three-interaction FDR family differs from Workflow 14; raw
three-way estimates/p values replay, but adjusted discovery counts differ.
Workflow 14's full historical analysis is preserved.

## Selected-programme diet associations (2 November 2026)

Workflow 16 is the current three-panel Figure 5-style revision. Read
`docs/workflow_plans/16_cf453_diet_rewiring_figure6.md` before modifications.
It uses Workflow 14 models.rds and Workflow 15 direct-diet q annotations,
checks original CF453 input hashes, and exports to
`results/cf453_diet_rewiring_figure6/`. `rewiring_analysis.rds` and `tables/`
contain the full eligibility ledger, reciprocal anchor selections, fixed
programme unions, both-diet coefficients, class/rest denominators and mouse
bootstrap/omission stability. `panel_manifest.csv` maps exactly three
independent PDF/SVG/PNG groups to source CSVs. Gain/loss networks are
threshold-defined selected examples, not calibrated differential-correlation
discoveries. All class 95% bootstrap ranges include zero; no coefficient
passes BH across the 76,772 selected correlation tests. Small classes and
whole-liver lipid interpretation require particular care.

Figure6A layout candidates requested4November2026 are isolated under
`results/cf453_diet_rewiring_figure6/layout_comparison_20261104/`.
`scripts/reports/cf453_figure6a_layouts.R` reads saved anchors only and
exports a34-gene/four-state effect heatmap and paired SD/WD dots, with the
same genes/effects/direct-q annotations. `common_gene_ledger.csv` and the
two source CSVs document all136 values per view. The original panels and
statistical outputs remain the historical deliverables until a layout is chosen.

## Static Query Map

The16November2026 expanded-marker update is isolated in
`results/zonation_marker_update_20261116/`; read
`docs/workflow_plans/17_expanded_marker_zonation_update.md` before modifying.
Entry17 stages the three supplied score pseudobulks, reruns isolated FCG05/07,
paired GDX support and CF45314models->15models->16->both6Alayouts. The report
and89-claim/24-edge replacement ledger trace main panels, sources and manuscript
values. Figure4A/S4B-C are explicitly deferred to James. Historical results,
independent assays and the manuscript are preserved; do not mix old WF12/13
references with these new compartment definitions. --validate-only checks
setup; --report-only reads saved versioned results and regenerates the audit.

Use these files directly for common questions that the discontinued MCP server
previously answered.

| Question | Primary file(s) | Notes |
| --- | --- | --- |
| WP1/WP2 local inputs and migration hashes | `data/WP1_WP2_INPUTS.md`; `data/wp1_wp2_input_migration_manifest.csv` | Five frozen runtime inputs plus Figure 3 provenance; source archives preserved; byte-identical local copies. |
| Dataset summary and workflow order | `README.md`; `docs/workflow_plans/_meta_cross_workflow_overview.md` | Start here for context restoration. |
| RNA differential expression for a gene | `results/pseudobulk/ttop_full_list.RData`; `results/pseudobulk/tval_full_df.csv` | Current workflow 02 writes RData plus t-value CSV, not a long-format RNA CSV. |
| Top RNA DE genes by contrast/cell type | `results/pseudobulk/ttop_full_list.RData` | Nested list: contrast display name -> cell type -> limma topTable. |
| RNA t-statistic matrix | `results/pseudobulk/tval_full_df.csv`; `results/pseudobulk/tval_full_df.RData` | CSV has `ID` plus combined contrast/cell-type columns. |
| Mouse hepatocyte gonad signature | `results/pseudobulk/hepa_gonad_signature.csv` | Two-column signature seed (`ID`, `t_Hepatocytes`) used by workflow 09 for human-liver scoring. It is a near-identical earlier export of the current hepatocyte gonad t-values. |
| Lipid differential abundance | `results/lipidomic/ttop_full_list.RData`; `results/lipidomic/gonade_type_ttop.csv` | `logFC > 0` is testes-biased for the main gonadal contrast. `gonade_type_ttop.csv` is a generated lookup table present in results, not written by the current workflow 03 script. |
| Lipid class enrichment | `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv`; `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData` | Produced by workflow 03 using decoupleR. |
| Cyp2c38 in vivo KO validation (October 2026 WP1) | `results/validation/cyp2c38_in_vivo/README.md`; `scripts/workflows/11_cyp2c38_in_vivo_validation.R`; `docs/workflow_plans/11_cyp2c38_in_vivo_validation.md` | Separate from primary-hepatocyte siRNA. Pooled/female/male KO-minus-WT lipid effects, exact fixed FCG reference replay, sex-matched predictions, stratified mouse permutation and paired-animal bootstrap. Use animal-level concordance inference; native ULM p-values are descriptive. |
| Cyp2c38 native ULM reporting, separate gonadal networks and bootstrap calibration | `results/validation/cyp2c38_in_vivo/ulm_followup/README.md`; `ULM_summary.csv` and `network_manifest.csv` in that directory; `scripts/workflows/11_cyp2c38_in_vivo_validation.R`; `docs/workflow_plans/11_cyp2c38_in_vivo_validation.md` | Current requested reporting: native ULM score and nominal p, explicit KO-only and KO-plus-reference bootstrap p. Full networks have 1,254 features / 856 strict KO matches; historical 828 / 552 retained. Negative ULM denotes expected response. Bootstrap tests zero tested KO contrast for every lipid, a stronger null than zero alignment. Calibration and injection diagnostics condition on a fitted residual scenario. |
| Lipid conversion network | `results/lipidomic/network/significant_lipids_network_SIF.csv`; `results/lipidomic/network/significant_lipids_network_ATT.csv` | Current workflow 03 writes these under `results/lipidomic/network/`. |
| Experimental mouse COSMOS prior | `results/causal_modelling/cosmos_mouse_pkn/pkn_lipid_curated.rds`; `results/causal_modelling/cosmos_mouse_pkn/pkn_lipid_rules.rds`; `docs/cosmos_mouse_pkn_build.md` | Rebuilt mouse signaling/metabolism plus lipid activity regulation and separate measured-species reaction hypotheses. Evidence, mapping coverage and source hashes are in the same results folder. |
| Comparison with the mouse brain lipid atlas CORNETO resources | `docs/cosmos_lamanno_pkn_comparison.md`; `scripts/pkn/compare_lamanno_lipid_network.py`; `support/cosmos_mouse/comparison/lamanno_2026/`; `results/causal_modelling/lamanno_2026_pkn_comparison/` | Read-only comparison. Recovered Table 10 has 350 lipid pairs and 137 possible measured-label matches; a broader deposited lipid draft has 4,562 pairs and 392 possible matches. Sum projections are conditional and often fail retained-chain checks; 29/30 printed table enzyme symbols already occur in our PKN. The complete signed CORNETO graph was not located in the cited deposits. |
| Expanded chain-compatible measured lipid reaction hypotheses | `docs/cosmos_lipid_species_expansion.md`; `scripts/pkn/build_expanded_lipid_rules.py`; `scripts/pkn/expanded_lipid_matching.py`; `scripts/pkn/validate_expanded_lipid_rules.py`; `scripts/pkn/validate_expanded_lipid_rules.R`; `support/cosmos_mouse/expanded_species_templates.json`; `results/causal_modelling/cosmos_mouse_pkn_expanded_species/` | Separate additive resource; follow-up MOON sensitivity is documented in the next row. Candidate products increase from 128 to 419 (331,032 edges), or 431 with typed glycan alternatives (331,472 edges). Preserves complete directional Rhea chemistry, retained chains, joint PC/DG requirements and three-step Pemt. Independent structural/projection and COSMOS checks passed. Mapping ledger is not a finalized input vector; reachability is unfiltered and ignores activity signs and joint substrate availability. |
| Factor 1 MOON with chain-compatible reaction expansion | `docs/cosmos_factor1_expanded_moon_rerun.md`; `docs/cosmos_reaction_to_pkn_rules.md`; `scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R`; `docs/workflow_plans/moon_factor1_tf_to_lipid_expanded_species.md`; `results/causal_modelling/factor1_tf_to_lipid_expanded_species/selected_network.json` | Separate conservative sensitivity: 155 preserved + 268 new one-to-one observations; 369 survive cofactor-filtered pruning. Preferred pruning 6/scoring 4, compressed, adjacent scoring levels, cutoffs 0.75/0.75: 138 edges, eight TFs, 13 measured features. Six layers adds Hoxa10 and edges, with the same features/cutoffs. Initial stronger-cutoff Cer topology depends on water and lacks recursive scoring witnesses. Seven fits and seven extraction sets independently passed. LPC annotation overlap remains explicit. |
| Ranked KO/KD validation sets for Factor 1 TF-to-lipid hypotheses | `docs/cosmos_factor1_ko_validation_targets.md`; `docs/workflow_plans/moon_factor1_ko_validation_prioritization.md`; `scripts/moon/ancillary/audit_factor1_ko_expression.R`; `scripts/moon/ancillary/propose_factor1_ko_validation.py`; `results/causal_modelling/factor1_ko_validation_targets/` | Current alternating_prioritization: five enzyme-only/five TF-containing sets, six singles/four pairs, nine genes; add three constituent singles for 13 conditions. Every pair tests a possible alternate producer of the same output. Cebpb/Pparg are source-corrected unscored follow-ups, not selected or sign-coherent Factor 1 routes; RXRA is a contextual intermediary. Unsupported source identities/ENPP2 are screened out; TF/PLD2 compensation requires PLD1 suppression and retained PLD2. Original graphs unchanged; no MOON refit or quantitative KO/synergy prediction. |
| Directionality-audited STITCH v5 exploratory overlay | `results/causal_modelling/cosmos_mouse_pkn_stitch_directional/pkn_exploratory_stitch_directional.rds`; `stitch_directional_evidence.csv`; `stitch_directional_manifest.json`; `stitch_directional_validation.json` | Explicitly opt-in sensitivity resource. Raw chemical-actor/action records only, with score and mapping audits; excludes every orientation-unknown cached STITCH row. It does not change `pkn_lipid_rules`; it is used only by the separate relaxed lipid→TF sensitivity, not the manuscript TF→lipid or original lipid→TF fits. STITCH v5 is from 2016, so directionality alone is not biological validation. |
| Measured lipid→protein regulatory coverage audit | `docs/cosmos_lipid_regulatory_measurement_coverage.md`; `docs/cosmos_lipid_regulation_resource_audit.md`; `results/causal_modelling/lipid_regulatory_mapping_audit/`; `results/causal_modelling/lipid_regulatory_resolution_audit/` | Separates planned from implemented resources and exact, compatible-isomer and broad-class matches across all 1,254 measurements, accepted 155 inputs, historical 64 lipids and two LPC candidates. The scored conservative PKN contains GtoPdb regulation only; cached BRENDA/STITCH retain unrecoverable orientation, while the raw directional-STITCH overlay is assessed through a separate conditional relaxed-mapping ledger and used only in the relaxed lipid→TF sensitivity. MRC assumed signs were rejected and CellPhoneDB was not imported. |
| Measured lipid reaction-product coverage audit | `docs/cosmos_lipid_reaction_target_mapping_audit.md`; `results/causal_modelling/lipid_reaction_target_mapping_audit/`; `scripts/pkn/audit_measured_lipid_reaction_targets.py` | Separates any-role species-rule connectivity from actual directed product-side endpoints. All 149 species rules preserve source→reaction→target semantics; 123 rule features are product-side, and all 74 post-pruning TF-to-lipid observations are product nodes. |
| Selected Factor 1 TF-to-lipid MOON workflow | `scripts/moon/run_factor1_tf_to_lipid.R`; `scripts/moon/config/factor1_tf_to_lipid.json`; `docs/workflow_plans/moon_factor1_tf_to_lipid.md`; `results/causal_modelling/factor1_tf_to_lipid/` | Active entry point for the manuscript configuration: pruning 6, scoring 4, compressed ULM with RNA coherence, thresholds 1.85/0.55. Regenerates scores, extraction, expected-TF diagnostics, provenance and plots; checks the refit against frozen C03 (388 edges/153 nodes). Full refit/export passed on 2026-09-15; graph/compression identities, score IDs/levels and numerical scores within 1e-10 matched. `reference_comparison.json` records checks. `--validate-only` checks setup without writing output; reruns need a new `--output` destination. |
| Factor 1 lipid-to-TF signalling hypotheses | `scripts/moon/run_factor1_lipid_to_tf.R`; `scripts/moon/config/factor1_lipid_to_tf.json`; `docs/workflow_plans/moon_factor1_lipid_to_tf.md`; `docs/cosmos_factor1_lipid_to_tf.md`; `results/causal_modelling/factor1_lipid_to_tf/` | Separate direction using frozen lipid loadings upstream and TF footprints downstream; original PKN directions and compressed scoring retained. Direct coverage and lipid-only precursor policies yield no original-PKN solution. `--include-ph73` adds the existing acid/base-normalized curated precursor sensitivity. The report separates coarse-grid status, finer extraction, score/positive-control audits and selected exploratory network; reaction routes require unmeasured co-substrates and do not establish lipid signalling or flux. |
| Relaxed Factor 1 lipid-to-TF regulatory candidates | `scripts/pkn/audit_factor1_relaxed_regulatory_mapping.py`; `scripts/moon/run_factor1_lipid_to_tf_relaxed.R`; `scripts/moon/config/factor1_lipid_to_tf_relaxed.json`; `scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R`; `scripts/moon/config/factor1_lipid_to_tf_relaxed_parameter_search.json`; `scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R`; `scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_score_layer_evidence.py`; `scripts/moon/ancillary/export_factor1_lipid_to_tf_relaxed_network.py`; `docs/cosmos_factor1_lipid_to_tf_relaxed.md`; `results/causal_modelling/factor1_lipid_to_tf_relaxed/` | Active direct chemical→protein plus protein→protein screen on the directional-STITCH overlay. The mapping ledger retains conditional assignments; Factor 1 lipid values are post-score annotations, not ULM inputs. The separate parameter search writes to `results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search/` and its evidence audit to `<search>/evidence_audit/`. The older `factor1_lipid_to_tf_relaxed_regulatory/` folder is a legacy discovery topology/search. The main refit reproduces the low-threshold LPS→Gpr34→Gnai1→Src→{Ar, Esr1} sensitivity; no compact strict solution survived either cutoff above 0.5, and AR/ESR1 are supplied terminal observations. |
| Ancillary Factor 1 parameter search | `scripts/moon/ancillary/search_factor1_parameters.R`; `scripts/moon/lib/workflow.R`; `scripts/moon/README.md` | Compressed-only configurable grids reuse the same core as the selected workflow and retain every threshold combination. Defaults test pruning 6 and layers 4,5,6,8. This is a new reusable runner, not a replay of the historical 7,563-pair search. |
| Factor 1 TF-to-lipid MOON pilot (historical) | `docs/cosmos_factor1_moon_pilot.md`; `results/causal_modelling/cosmos_factor1_pilot/`; `scripts/moon/archive/2026-09/run_factor1_moon.R` | Original six-layer pilot. Compressed tutorial thresholds yield no paths; independently converged uncompressed analysis yields 8 TFs/4 lipids. The selected manuscript result uses the active depth-4 workflow above. Archived scripts retain historical paths. |
| Factor 1 MOON parameter search (historical) | `docs/cosmos_factor1_moon_parameter_search.md`; `results/causal_modelling/cosmos_factor1_parameter_search/`; `scripts/moon/archive/2026-09/search_factor1_parameter_space.R` | Original depth 4--12 and dual-threshold search. Its uncompressed U5/U4 preference is superseded by the compressed review below; retain these runs as representation diagnostics. Scores, topology audit, exact PKN validation, GraphML and evidence exports remain available. |
| Factor 1 compressed MOON review and TF retention (historical selection evidence) | `docs/cosmos_factor1_compressed_review.md`; `results/causal_modelling/cosmos_factor1_compressed_review/`; `scripts/moon/archive/2026-09/search_compressed_threshold_review.R` | 7,563 additional threshold pairs. Baseline C03 is compressed radius 6/depth 4 at 1.85/0.55: 388 edges, Ahr/Hnf4a/Ppara/Tbx2 and 3 DGs; stronger compression reproduces its exact graph. A radius-12/depth-8 compressed alternative has 468 edges at 2.15/1.075. Esr1/Nfic/Dbp survive preprocessing but are removed by MOON upstream-sign disagreement. Complete score-retention flags, 66 TF/run audits, compression validation, provenance and GraphML exports are saved. |
| COSMOS/MOON manuscript Results and Methods | `docs/manuscript_drafts/cosmos_results_methods.md`; `docs/manuscript_drafts/cosmos_results_methods_notes.md`; `scripts/moon/run_factor1_tf_to_lipid.R` | Draft additions with verified PKN scale, the selected 388-edge compressed network and expected-TF exclusions. The active workflow refits the selected configuration against frozen inputs and compares it with the archived result. Original scripts and their SHA256 migration manifest are in `scripts/moon/archive/2026-09/`; frozen results and historical source-hash records are preserved. |
| MOFA factor and feature weights | `results/mofa/main_model_RNA_lipid_weight.RData`; `results/mofa/factor_analysis/` | Factor 1 is the main gonad-associated factor. |
| Gene-lipid correlations | `results/mofa/networks/gene_lipid_correlations_factor1.csv`; `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | Correlation matrix rows are genes; columns are lipids. |
| TF-gene-lipid network | `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv`; `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_ATT.csv` | Cytoscape-ready compressed regulatory network. |
| Zonation DE and summaries | `results/zonation/ttop_full_long.csv`; `results/zonation/t_table_zonation.csv`; `results/zonation/summary_statistics.csv` | Workflow 05 writes long-format zonation DE. |
| Sex-specific zonation genes | `results/zonation/sex_specific_zonation_genes_summary.csv` | Main selected spatially biased gene table. |
| WDSD baseline conservation and diet networks (WP2) | `scripts/workflows/12_wdsd_network_comparison.R`; `docs/workflow_plans/12_wdsd_network_comparison.md`; `results/wdsd_network_comparison/README.md` | The combined report presents pooled Figure 3 first, then four frozen Figure 5 references compared with CF453 SD, then CF453 WD versus SD. `baseline_inference.csv` contains native ULM and exact pairing/null-bootstrap probabilities; `network_edges.csv` contains all paired coefficients, delta orientation and 49-omission stability. Seven mice per sex/diet, two RNA zones sharing whole-liver lipid measurements. Molecular edges are not independent animals; diet changes are descriptive. |
| Lipid overlap and pooled Figure 3 (WP2 follow-up) | `scripts/workflows/12_wdsd_network_comparison.R`; `docs/workflow_plans/12_wdsd_network_comparison.md`; `results/wdsd_network_comparison/followup/README.md` | Raw annotation audit: 752 strict matches, 959 alias-covered original rows, 854 conditional compatibility pairs. `lipid_overlap_audit/` retains every candidate and ambiguity. `pooled_inference.csv` compares 108 displayed pairs and 5,004 precursor pairs across 43 FCG/14 SD whole-hepatocyte profiles. Overall absolute-tail and within-sex upper-tail animal tests answer different questions; sex-centered statistics are descriptive. No zonation or MOFA refit. |
| Male Western-diet cross-cohort robustness (WP3) | `scripts/workflows/13_wd_crosscohort_robustness.R`; `docs/workflow_plans/13_wd_crosscohort_robustness.md`; `results/wd_crosscohort_robustness/README.md`; `crosscohort_inference.csv`; `network_edges.csv` | CF453 seven male WD mice versus CF459 four male WD mice; two zones, frozen original full/displayed pairs, strict three-way lipid matching. Coefficient-vector r/native ULM, exact conditional 24 and complementary 5,040 pairing tests, approximate marginal-null bootstrap, paired/omission stability. Original male FCG SD anchor n=22 (11 XXT/11 XYT), same-edge delta r and paired three-cohort stability in `fcg_comparison/`. User-confirmed ages: CF459/FCG 3 months, CF453 10 months. No within-CF459 SD arm; background/exposure and analytical-QC metadata remain limitations. |
| Lipid-zonation associations | `results/zonation/lipid_class_association_summary.csv`; `results/zonation/lipid_species_zonation_association.csv`; `results/zonation/threshold_sensitivity_total.csv` | Connects lipid species/classes to zonation categories and threshold behavior. |
| Structural lipid associations | `results/lipidomic/structural_model/lipids_with_structure.csv`; `results/lipidomic/structural_model/lipid_only_class_coefficients.csv`; `results/lipidomic/structural_model/lipid_only_class_summaries.csv`; `results/lipidomic/structural_model/lipid_only_chain_length_effects.csv`; `results/lipidomic/structural_model/lipid_only_db_position_effects.csv`; `results/lipidomic/structural_model/lipid_only_interaction_effects.csv` | Workflow 06 currently tests lipid-only structural effects on gonadal logFC. Older gene-lipid structural sections are commented out. |
| Validation experiments | `results/validation/validation_summary.csv`; `results/validation/global_network_validation_scores.csv`; `results/validation/zonation_network_validation_scores.csv`; `results/validation/structural_consistency/consistency_structural_regression_summary.csv` | Workflow 07 compares predicted correlations with knockdown lipidomics and structural consistency. |
| Lipid variance analysis | `results/lipidomic/variance_analysis/variance_analysis_summary.csv`; `results/lipidomic/variance_analysis/variance_analysis_results.csv`; `results/lipidomic/variance_analysis/variance_analysis_by_class.csv` | Workflow 08 compares variability between gonadal groups. |
| Human liver comparison | `scripts/workflows/09_human_analysis.py`; `docs/workflow_plans/09_human_analysis.md` | Optional Python workflow using external `GSE179250` human liver inputs. Outputs are written to the script's external `RESULTS_DIR`, not this repository, by default. |
| Xenium spatial analysis | `scripts/workflows/10_xenium_analysis/`; `docs/workflow_plans/10_xenium_analysis.md` | Optional collaborator Python workflow folder for Xenium doublet inference, cell-type annotation, composition plots, and Scanpy DGE. Raw Xenium inputs and outputs are folder-local by default. |
| MALDI-MSI collaborator notebook | `scripts/james/maldi/Cleland_analysis_expt2_recal.ipynb`; `scripts/james/README.md` | Provenance notebook for negative-mode MALDI-MSI liver zonation analysis. Raw AnnData, imaging, registration, and SCiLS inputs are external. |
| Original P3225 proteomics Rmd | `scripts/james/proteomics/JamesCleland_P3225_V1.Rmd`; `data/proteomics/`; `docs/workflow_plans/02_fcg_pseudobulk_analysis.md` | Collaborator R Markdown source for the fractionated TMT proteomics analysis. Workflow 02 uses the curated checked-in proteomics files, not the Rmd directly. |

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
