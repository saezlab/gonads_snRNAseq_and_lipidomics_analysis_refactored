# Factor 1 TF→lipid: chain-compatible expanded PKN sensitivity

This ancillary workflow extends the frozen TF→lipid analysis without replacing
the manuscript result or its reference assertions. It uses
`scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R` and
`scripts/moon/config/factor1_tf_to_lipid_expanded_species.json`.

## Inputs and scope

- Frozen Factor 1 TF activities, RNA coherence weights and regulons:
  `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds`.
- Chain-compatible PKN:
  `results/causal_modelling/cosmos_mouse_pkn_expanded_species/pkn_expanded_chain.rds`.
- Mapping preparation: `scripts/moon/ancillary/prepare_expanded_factor1_mapping.py`.
  Preserve all 155 old observation identities and values, adding 268 previously
  unmapped measured product features once each. The resulting 423 observations
  include 418 active product endpoints and five preserved non-product mappings.
  Exclude conditional glycan alternatives; retain unresolved positional and
  stereochemical assumptions as hypotheses.
- The `frozen_control` mode uses the original 155 downstream observations
  verbatim on the expanded graph. The `expanded` mode uses the extended mapping.
  Both keep the same TF cutoff, joint factor orientation, compression and RNA
  coherence policy.

## Operations and audit contracts

Validate frozen software/input hashes and expanded graph hash. Terminalize
observed lipids, prune at six steps in each search, protect input and regulon
identifiers during same-children compression, score four layers, and repeat RNA
coherence until convergence within ten iterations. Save the exact final scoring
graph, original-ID projection, scores and compression map together.

Extract the fixed 1.85/0.55 comparison, then evaluate the recorded bounded
threshold grid on the same scores. Record every valid grid combination,
including empty networks. A compact representative requires fewer than 500
edges and at least two TFs and two lipid readouts. Rank eligible combinations by
secondary then primary cutoff, then coverage. This is a reproducible selection
policy within a grid, not a statistical optimum.

Diagnose Esr1/Nfic/Dbp from actual pre-sign-filter scores and pipeline retention;
do not assert that the baseline exclusion pattern must recur. Export exact
signed provenance and complete reaction records, checking for generic cofactor
bridges and omitted participants in the extracted graph. Full reaction records
do not make ordinary MOON edges enforce joint substrate availability or flux.

## Outputs and execution

Default root: `results/causal_modelling/factor1_tf_to_lipid_expanded_species/`.
`inputs/` contains the mapping audit. Separate `expanded_inputs/` and
`frozen_inputs_control/` directories contain fits, ledgers, fixed comparisons,
threshold tables, selected networks and expected-TF diagnostics. Existing
nonempty fit directories are refused. Commands and final validation results
will be recorded in the result report after execution.

## Executed result and revised interpretation (6 October 2026)

Read `docs/cosmos_factor1_expanded_moon_rerun.md` and the result-root
`selected_network.json`. Seven fits and seven extraction sets passed the
independent checker `scripts/moon/ancillary/audit_expanded_factor1_rerun.R`.
The initial 101-edge higher-cutoff Ar/Arnt→Cer topology depends on a water bridge
and lacks consecutive recursive scoring support; preserve it as diagnostic
material, not the preferred interpretation.

The declared conservative sensitivity removes only generic cofactor
substrate-participation edges listed in
`support/cosmos_mouse/moon_currency_substrate_edges_v2.csv`. Complete chemical
masters and chemical→protein regulatory interactions are retained. Its builder
is `scripts/moon/ancillary/build_currency_traversal_policy.py --alias-complete`.
V1 is preserved for audit; the alias-complete V2 yields exactly identical
prepared graphs, compression maps and scores after pruning in these data.

`--currency-policy=no_generic_substrate --currency-edges=CSV` selects this
explicit scoring sensitivity. `--layer-consistent` restricts selected extraction
to edges with source level exactly one above the target. Fixed-cutoff broad
topological comparators remain separately saved under `fixed_thresholds/`.

The final fits are `currency_v2_expanded_layers4/` and
`currency_v2_expanded_layers6/`. Re-extraction with
`scripts/moon/ancillary/reextract_expanded_factor1_thresholds.R` evaluates 165
valid threshold pairs per depth and preserves invalid/empty outcomes. Preferred
`currency_v2_expanded_layers4_lower_thresholds/`: cutoffs 0.75/0.75, 138 edges,
94 nodes, eight TFs and 13 measured features. Six layers: 177 edges, 117 nodes,
nine TFs and the same 13 features at the same cutoffs. Choose four layers for
parsimony within this grid; neither is a calibrated optimum. Esr1/Nfic/Dbp
survive preparation but are removed by upstream sign conflict in these fits.

`scripts/moon/ancillary/export_expanded_factor1_network.py` writes exact signed
provenance, GraphML, full masters, participant/currency audits and adjacent-layer
routes. All displayed primary edges have adjacent levels, all displayed catalyst
links are retained, and no currency bridges remain. The primary has 36 audited
TF/readout pairs. `scripts/moon/ancillary/plot_expanded_factor1_routes.R` produces
three illustrative routes in `figures/`; its PNG was visually inspected.

Count 13 measured features, not independently resolved species: two LPC
annotation pairs may overlap chemically. No weights are copied to aliases, but
assay-feature independence and joint substrate availability are not established.
The stronger original thresholds give no conservative expanded network in these
tests. Preserve this as an exploratory sensitivity alongside the manuscript fit.

## 16 November 2026 repository cleanup

Historical causal result folders now live under
`miscellaneous/older_results/causal_modelling/`. Active configurations and
consumer defaults explicitly use their frozen inputs/reference states there.
The per-file old/new path and SHA256 ledger is
`docs/provenance/causal_cleanup_20261116.csv`; the active/frozen/archive dependency
map is `docs/CAUSAL_ANALYSIS_MAP.md`. No numerical settings, graph policy, mapping,
score or fitted result changed. Saved-artifact extraction/reference replay passed
without fitting or scoring. Full-fit version gates remain unchanged; the installed
decoupleR 2.9.7 differs from the pinned 2.17.0 and blocks full-fit validation.
