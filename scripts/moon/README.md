# Factor 1 causal modelling with MOON

The [causal analysis map](../../docs/CAUSAL_ANALYSIS_MAP.md) distinguishes current
analyses from frozen archived inputs and documents the 16 November cleanup.
Use `ancillary/validate_frozen_archives.R` for saved-artifact checks without fitting.
Full-fit software gates remain unchanged.

## Chain-compatible expanded PKN sensitivity

The 6 October 2026 follow-up is separate from the selected manuscript replay.
Read [`docs/cosmos_factor1_expanded_moon_rerun.md`](../../docs/cosmos_factor1_expanded_moon_rerun.md)
and [`docs/cosmos_reaction_to_pkn_rules.md`](../../docs/cosmos_reaction_to_pkn_rules.md).
The preferred conservative sensitivity has 138 edges, 94 nodes, eight TFs and
13 measured lipid features at 0.75/0.75, pruning 6 and scoring 4. It blocks
generic cofactor substrate traversal and requires consecutive scoring layers.
The initial stronger-cutoff Cer topology depends on water and remains diagnostic.

Ancillary preparation, fitting, extraction, currency-policy construction,
provenance/route export and independent validation are respectively in
`ancillary/prepare_expanded_factor1_mapping.py`,
`ancillary/run_expanded_factor1_tf_to_lipid.R`,
`ancillary/reextract_expanded_factor1_thresholds.R`,
`ancillary/build_currency_traversal_policy.py`,
`ancillary/export_expanded_factor1_network.py` and
`ancillary/audit_expanded_factor1_rerun.R`. The separate configuration is
`config/factor1_tf_to_lipid_expanded_species.json`. Commands and the final
four/six-layer comparison are in the report; `ancillary/plot_expanded_factor1_routes.R`
renders audited examples. Results are under
`results/causal_modelling/factor1_tf_to_lipid_expanded_species/`.

For the lipid→TF extension, read the
[regulatory measurement-coverage audit](../../docs/cosmos_lipid_regulatory_measurement_coverage.md)
before interpreting the long LPC precursor routes. It distinguishes the
implemented regulation resources from the initial proposal and documents
specific fatty-acid candidates missed by the original mapper. It does not
remap inputs or alter either fitted analysis.

## Directionality-audited STITCH sensitivity

[`results/causal_modelling/cosmos_mouse_pkn_stitch_directional/`](../../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/)
contains a separate exploratory PKN in which all cache-normalized STITCH rows
were replaced with raw chemical-actor records with explicit direction and
sign. The active relaxed lipid→TF configuration
`config/factor1_lipid_to_tf_relaxed.json` uses this overlay and its separate
conditional mapping ledger. The selected TF→lipid baseline and the exact-mapped
lipid→TF configuration retain the original PKN. The direct-coverage and no-route
statements below apply to that original `pkn_lipid_rules` resource, not to the
relaxed overlay analysis. Do not substitute either prior for the other.

The retained overlay, mapping/evidence tables and selected states ship with the
repository. Rebuilding the overlay additionally requires the two large local-only
STITCH action/link downloads recorded in its source manifest; those raw files
are not needed merely to read or fit from the retained overlay.

The main workflow reproduces the **TF activities → lipid weights** network in
the manuscript Results. It starts from the frozen, previously audited Factor 1
inputs and lipid-expanded mouse PKN, then recalculates pruning, compression,
MOON scores, RNA coherence, solution extraction and expected-TF diagnostics.
It compares the new results with the archived selected C03 network before
exporting evidence and figures. PKN construction and MOFA fitting are separate
upstream analyses; this command does not rebuild those resources.

## Run the selected analysis

From the project root:

```sh
Rscript scripts/moon/run_factor1_tf_to_lipid.R --validate-only
Rscript scripts/moon/run_factor1_tf_to_lipid.R
```

The default output is `results/causal_modelling/factor1_tf_to_lipid/`. This
directory now contains the verified refit produced on 15 September 2026. For
another run, supply a new directory:

```sh
Rscript scripts/moon/run_factor1_tf_to_lipid.R \
  --output=results/causal_modelling/factor1_tf_to_lipid_replay
```

Existing nonempty output directories are refused. Paths in the configuration
and command options resolve from the project root; the entry point can also
be invoked by its absolute path from another directory. `--validate-only`
checks dependencies and hashes without scoring or writing output.

The software gate is part of that validation, not an optional full-fit check.
As confirmed from installed package metadata on 18 November 2026, local
decoupleR is 2.9.7 while these replay configurations require 2.17.0. Thus the
three main MOON `--validate-only` commands stop before fitting in this
environment. Saved-artifact extraction checks use a separate entry point and
do not certify that a fresh fit is ready:

```sh
Rscript scripts/moon/ancillary/validate_frozen_archives.R
```

That command reads existing selected/relaxed states and their frozen references,
checks hashes and re-extracts networks without new MOON scoring. By default it
prints a report; `--output=NEW_JSON` explicitly writes one new report and refuses
an existing filename. The 13 compact archived dependencies in
[`active_frozen_dependencies.json`](../../docs/provenance/active_frozen_dependencies.json)
are deliberately shipped. Other historical search payloads remain local-only.

The checked-in configuration is
[`config/factor1_tf_to_lipid.json`](config/factor1_tf_to_lipid.json):

| Setting | Selected value |
| --- | --- |
| Upstream | Hepatocyte Factor 1 TF footprints, absolute activity > 2 |
| Downstream | Exact-mapped, unscaled Factor 1 lipid weights |
| Orientation | Jointly oriented RNA, TF and lipid weights; testes-positive factor |
| RNA signature | TF–target coherence only |
| Directed pruning distance | 6 steps in both reachability searches |
| MOON scoring layers | 4 |
| Compression | Identical signed children; protect inputs and regulon identifiers |
| Coherence | Fixed compression map; repeat until no edges change, maximum 10 passes |
| Primary / secondary thresholds | Strict absolute inferred-score cutoffs > 1.85 / > 0.55 |
| Result | 388 signed edges, 153 nodes, four TFs and three measured DG species |

Pruning distance, scoring layers and route length are distinct. Neither of the
first two is imposed as a length limit on the final extracted directed walks.
Measured observations are terminal, exempt from inferred-score cutoffs and
cannot seed the primary threshold. The selected configuration is a documented
representative from the search, not a calibrated statistical optimum.

The replay checks cosmosR 1.21.2 at commit
`b2ed195738479c4ad204cace37fd11ae2fd02bb8`, decoupleR 2.17.0,
igraph 2.3.3 and jsonlite 2.0.0, plus the PKN, mapped-input and compatibility-helper
MD5 hashes. R 4.6.1 was used for the verified refit; the current R/session details
are recorded on every run. Export uses Python 3's standard library and base R
plotting. The main result must match archived edge identities, score node IDs,
levels and compression mapping, with numeric score differences below 1e-10.
Different software versions need a deliberate validation update, not a silent
replacement of the reference result.

## Files and responsibilities

| File or directory | Role |
| --- | --- |
| `run_factor1_tf_to_lipid.R` | Main selected-result orchestration and manuscript comparison |
| `config/factor1_tf_to_lipid.json` | Accepted input paths, versions, settings and frozen reference |
| `lib/workflow.R` | Shared loading, pruning, compressed scoring, coherence, extraction checks and diagnostics |
| `moon_safe_helpers.R` | Unchanged, tested compression/decompression and dual-threshold contracts; path and bytes retained for historical provenance |
| `ancillary/search_factor1_parameters.R` | Additional compressed parameter grids using the same fitting functions |
| `ancillary/export_factor1_network.py` | Exact signed-edge provenance, GraphML and illustrative path export |
| `ancillary/plot_factor1_network.R` | Figures called by the exporter |
| `ancillary/audit_direct_lipid_regulatory_mapping.py` | Read-only per-measurement audit of the frozen crosswalk against raw and retained regulatory records |
| `ancillary/audit_regulatory_resolution_candidates.py` | Read-only ChEBI specificity audit for compatible regulatory chemicals missed by the frozen candidate list; no identity reassignment |
| `ancillary/summarize_lipid_regulatory_coverage.py` | Consolidates the two coverage audits into one row per measured lipid and deduplicated resource/panel summaries |
| `archive/2026-09/` | Unchanged historical preparation, pilot, searches and sensitivities, with a path/SHA256 migration manifest |

The main output keeps the fitted graph and its compression map together:

- `run_state.rds`, `prepared_topology.rds`, `compression.rds`,
  `pkn_final_scoring.rds` and `pkn_final_uncompressed.rds` contain the model state.
- `pruning_ledger.csv`, `coherence_ledger.csv` and
  `coherence_removed_edges.csv` record preprocessing and fitting.
- `moon_scores.csv` contains returned original-ID scores;
  `moon_scores_all.csv` additionally exposes the final upstream-sign exclusions.
- `focus_tf_retention.csv` verifies that Esr1, Nfic and Dbp survived topology
  and compression but were excluded after scoring because their inferred and
  supplied activity signs disagreed. The diagnostic pass never changes fitting.
- `network_SIF.csv`, `network_ATT.csv`, `solution.rds` and `network.graphml`
  contain the selected network. `network_evidence_network*`,
  `candidate_paths_network*` and `figures/` document reactions and routes.
- `results_summary.json`, `reference_comparison.json` and
  `extraction_audit.json` report the headline values and validation.
- `run_manifest.json`, `dependency_manifest.csv`, `output_manifest.csv`,
  `sessionInfo.txt` and `export.log` record execution and provenance.

## Supporting parameter exploration

The supporting search recalculates each requested pruning-distance/scoring-layer
combination once, then extracts every requested dual-threshold pair from its
saved scores. It uses the same compressed fitting implementation as the main
workflow. For example:

```sh
Rscript scripts/moon/ancillary/search_factor1_parameters.R \
  --pruning-steps=6,12 --layers=4,8 \
  --primary=1.85,2.15 --secondary=0.55,1.075 \
  --output=results/causal_modelling/factor1_parameter_followup
```

Use `--help` for the complete options and `--validate-only` to inspect a grid
without fitting it. All combinations, including empty solutions and failures,
are retained in the search tables. Eligibility requires a converged, audited
network with fewer than 500 edges and at least two upstream and two downstream
inputs. A frontier helps compare thresholds and coverage; it does not establish
biological confidence or select a universal optimum.

The historical 7,563-combination compressed review and its stronger-compression,
initial-RNA-coherence and TF-support audits remain under
`miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/`, with the generating
scripts documented in [`archive/2026-09/README.md`](archive/2026-09/README.md).
The new grid runner does not claim to recreate that entire historical search.
Archived scripts contain historical paths and should not be run as current
entry points without migration.

## Lipid-to-TF analysis

The second direction uses the same fitting functions with a separate biological
input decision, graph policy and configuration. It preserves the original PKN
edge direction. Run:

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf.R --validate-only --include-ph73
Rscript scripts/moon/run_factor1_lipid_to_tf.R --include-ph73 \
  --output=results/causal_modelling/factor1_lipid_to_tf_replay
Rscript scripts/moon/ancillary/refine_lipid_tf_thresholds.R \
  --results=results/causal_modelling/factor1_lipid_to_tf_replay
Rscript scripts/moon/ancillary/audit_lipid_tf_scores.R \
  --analysis=results/causal_modelling/factor1_lipid_to_tf_replay
```

The main runner audits direct regulation and lipid-only precursor traversal
before scoring. It starts from the 155 accepted lipid loadings and all 353
PKN-present TF footprints; the 49 strong TFs are a separate observation-set
sensitivity. The original `pkn_lipid_rules` PKN has no measured lipid with a direct regulatory
edge and no lipid-only route to TF observations. `--include-ph73` explicitly
adds the existing acid/base-normalized curated PKN sensitivity. Omitting the
flag preserves the original-PKN negative results, for which subsequent
refinement is unnecessary.

The completed run is under `results/causal_modelling/factor1_lipid_to_tf/`.
There are 72 initial configurations, 16 successful compressed fits, 752 coarse
extractions and 1,431 fine extractions. A documented representative from the
26 compact fine-grid candidates uses **12 pruning steps, 12 scoring layers,
primary 1.65 and secondary 1.10**, giving **455 edges/233 nodes**, one
LPC(18:0/0:0) input and 64 TF observations. It uses the all-TF panel, converges
in three coherence passes and is exported in `refinement/selected_all_tfs/`.
Both retained chemical-regulatory transitions lead from PGF₂α/PGD₂ to Ptger3.
The strong-TF sensitivity yields no compact network in the tested refinement.

This network is an exploratory lipid-remodelling/availability hypothesis.
LPC is regenerated by the acylation/hydrolysis cycle and arachidonate is supplied
by unmeasured arachidonoyl-CoA. It is not a direct LPC receptor interaction or
evidence for synthesis of arachidonate from LPC's stearoyl chain. Esr1 and Dbp
are retained **observations**, not independently predicted controls. The
LPC16/DHA/Rxra branch is absent from the selected solution and is sensitive to
which TFs are supplied as observations.

| Active component | Responsibility |
| --- | --- |
| `run_factor1_lipid_to_tf.R` | Input/software validation, coverage audit, feasibility gates and compressed fitting/coarse extraction |
| `config/factor1_lipid_to_tf.json` | Direction, graph/observation policies, grid and frozen dependencies |
| `ancillary/audit_lipid_tf_coverage.py` | Exact/ambiguous mapping distinctions, currency-shortcut diagnostics and lipid-only graph construction |
| `ancillary/refine_lipid_tf_thresholds.R` | Refine saved fits, record all threshold outcomes, materialize and export the documented representative |
| `ancillary/audit_lipid_tf_scores.R` | Recompute final-graph scores; diagnose lipid sign exclusions, scoring support and expected-TF status |
| `ancillary/export_lipid_tf_network.py` | Separate lipid/TF/MOON quantities, exact signed provenance with explicit aliases, GraphML and illustrative paths |
| `ancillary/audit_lipid_tf_selected_reactions.py` | Complete selected reaction equations, catalyst alternatives and original receptor-binding evidence |
| `ancillary/audit_normalized_rules_topology.py` and `ancillary/audit_normalized_rules_pruning.R` | Verify that restoring measured-species rules leaves every tested pruned graph, input vector and compression map unchanged |

The main `analysis_summary.json` reports the initial coarse grid; its zero
eligible count precedes the 26 finer candidates in
`refinement/refinement_summary.json`. Per-run RDS objects, compression maps,
scores and ledgers are in `runs/`. `score_audit/` independently reproduces all
16 retained score tables exactly, without fitting coherence again. Dependency
hashes and source snapshots distinguish initial grid execution from subsequent
selection/export code.

Full biological interpretation, score bottlenecks, reaction evidence and audit
entry points are in [the analysis report](../../docs/cosmos_factor1_lipid_to_tf.md)
and [workflow plan](../../docs/workflow_plans/moon_factor1_lipid_to_tf.md).

## Relaxed regulatory lipid-to-TF screen

The relaxed screen is a separate exploratory workflow for the
directionality-audited STITCH overlay. It uses an audited long-form
feature-to-chemical ledger rather than the frozen direct lipid inputs. Candidate
chemicals are topology seeds; Factor 1 lipid weights are **post-score
sign-concordance annotations**, not numeric ULM inputs. The graph contains only
direct chemical-to-protein and protein-to-protein edges, without reaction
traversal.

The active selected-analysis entry point and configuration are:

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R --validate-only
Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R
```

They use `config/factor1_lipid_to_tf_relaxed.json` and write to
`results/causal_modelling/factor1_lipid_to_tf_relaxed/` by default. The separate
parameter search uses
`ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R` and
`config/factor1_lipid_to_tf_relaxed_parameter_search.json`, with default output
under `results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search/`.
It retains threshold and depth sensitivity separately from the selected main
analysis, materializes only strict consecutive-score-layer solutions, and does
not designate an automatic optimum.

Use `ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R` to audit
candidate-local score traces and duplicate direct-edge fingerprints. Its package
score traces are flagged when they are not a complete descending
candidate-to-observation walk; such traces are not display routes. Use
`ancillary/audit_factor1_lipid_to_tf_relaxed_score_layer_evidence.py` with the
parameter-search directory to attach edge-level PKN evidence to every
materialized strict solution under `<search>/evidence_audit/`; and use
`ancillary/export_factor1_lipid_to_tf_relaxed_network.py` to export an audited,
score-layer-constrained selected display network. A broad sign-consistent
connectivity view is not a recursive-MOON score explanation.

The pre-cleanup discovery topology/search remains in
`miscellaneous/older_results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory/` as legacy
audit material, rather than an active result. The detailed assumptions,
discovery result and caveats are in
[`docs/cosmos_factor1_lipid_to_tf_relaxed.md`](../../docs/cosmos_factor1_lipid_to_tf_relaxed.md)
and [`docs/workflow_plans/moon_factor1_lipid_to_tf_relaxed.md`](../../docs/workflow_plans/moon_factor1_lipid_to_tf_relaxed.md).

The selected main workflow exactly replays the discovery's conditional
LPS→Gpr34→Gnai1→Src→{Ar, Esr1} route as a low-threshold sensitivity: every
compact score-layer result used 0.5/0.5 cutoffs, and none survived either cutoff
above 0.5. AR and ESR1 are supplied terminal observations, not independent
validation targets.

## KO and knockdown validation proposal

The 13 October 2026 source/expression audit and revised ranking are in
[`docs/cosmos_factor1_ko_validation_targets.md`](../../docs/cosmos_factor1_ko_validation_targets.md).
Ten overlapping designs alternate five enzyme-only and five TF-containing
sets. Every pair targets a plausible alternate route retaining the same lipid
output. Six singles/four pairs plus constituent Pld2, Plcb1 and Plcg1 singles
give 13 target-containing conditions; no justified triplet is currently proposed.

Use `ancillary/audit_factor1_ko_expression.R` for actual-count/context evidence
and `ancillary/propose_factor1_ko_validation.py` with
`support/cosmos_mouse/factor1_ko_validation_review.json` for reviewed ranked/path
exports. Current outputs are in
`results/causal_modelling/factor1_ko_validation_targets/alternating_prioritization/`.
Initial exports and source snapshots are retained. Unsupported source identities
are excluded from audit paths. Cebpb/Pparg regulatory prefixes are separately
labelled unscored follow-up hypotheses, not selected, sign-coherent Factor 1
MOON explanations. RXRA and all TF-to-enzyme claims require mediator/context
gates. No PKN mutation, new fit or quantitative KO/synergy simulation occurs.

Recomputing the expression audit additionally needs the local original
`data/RNA/fcg_pseudobulks.rds` whole-cell count object and the archived culture
count CSV under `miscellaneous/older_inputs/legacy_rna/`; neither is shipped.
The seven public CSV transports do not supply this whole-cell FCG counts
object. The proposal exporter can consume its retained expression-audit CSV
instead, with its reviewed graph/source inputs. Keep those two preparation
routes distinct; neither command fits MOON.
