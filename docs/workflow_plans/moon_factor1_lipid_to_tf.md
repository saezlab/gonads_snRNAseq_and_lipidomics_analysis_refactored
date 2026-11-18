# Optional MOON workflow: Factor 1 lipid-to-TF hypotheses

## Scope and dependencies

The entry point is `scripts/moon/run_factor1_lipid_to_tf.R`, configured by
`scripts/moon/config/factor1_lipid_to_tf.json`. Read this plan and
`scripts/moon/README.md` before modifying it. It uses the shared compressed
fitting/extraction implementation in `scripts/moon/lib/workflow.R` and the
unchanged `scripts/moon/moon_safe_helpers.R`. It does not reverse PKN edges,
retrain MOFA, rebuild the original PKN or change the accepted input identities.

Inputs are the frozen Factor 1 bundle in
`miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/` and prior/evidence files
in `results/causal_modelling/cosmos_mouse_pkn/`. The forward configuration is
read only to validate the shared input/software contract. Workflow 04 is the
upstream source of the factor analysis; workflows 00–10 do not consume the new
outputs. The TF-to-lipid manuscript result is independent and remains unchanged.

## Measurement and graph decisions

- Upstream: 155 accepted species-level lipid Factor 1 weights, comprising six
  ChEBI nodes and 149 explicit measured-species hypothesis nodes. These are
  abundance associations, not measured receptor activity or imposed perturbations.
- Downstream: all 353 finite nonzero PKN-present hepatocyte TF footprints.
  A separate sensitivity uses the 49 with absolute activity above two.
  Observations are terminal, so this changes both topology and scoring background.
- RNA loadings are reserved for the existing internal TF–target coherence loop.
  All quantities retain the joint testes-positive Factor 1 orientation.
- The primary direct-regulation policy allows chemical→protein and
  protein→protein edges from the original PKN. The primary precursor policy
  additionally allows lipid→reaction→lipid traversal, using the frozen lipid
  ontology classification. It removes nonlipid chemical shortcuts and
  protein→reaction edges; reaction catalyst annotations remain evidence.
- `--include-ph73` adds the existing acid/base-normalized **curated** PKN as an
  explicitly exploratory sensitivity. Its six accepted input IDs are unchanged;
  it does not append the measured-species rule layer. Any separate audit of that
  distinction must be labelled and must not overwrite the original PKN.

## Execution and stages

Run from the project root (absolute invocation also works):

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf.R --validate-only --include-ph73
Rscript scripts/moon/run_factor1_lipid_to_tf.R --include-ph73 \
  --output=results/causal_modelling/factor1_lipid_to_tf_replay
```

The default output is `results/causal_modelling/factor1_lipid_to_tf/`; nonempty
destinations are refused. Omitting `--include-ph73` runs the two original-PKN
policies, including their negative results.

1. Validate the frozen input/software contract and record measurement semantics.
2. Run `ancillary/audit_lipid_tf_coverage.py` on the frozen CSV resources. Record
   direct interaction coverage, exact versus ambiguous mapping opportunities,
   unrestricted/currency-filtered diagnostic paths, and lipid-only graphs.
3. Alternate directed forward/backward pruning at 6, 8 and 12 steps, with
   observed TFs terminal. Record absent topology and fewer than three available
   TF observations as explicit feasibility failures rather than fabricated scores.
4. Fit available configurations at 4, 6, 8 and 12 scoring layers. Compress only
   identical signed-child groups while protecting inputs and regulon identifiers.
   Keep the compression map fixed during up to ten RNA-coherence passes.
5. Extract each valid coarse primary/secondary pair from the final saved scores.
   Require convergence, exact final-PKN edges, sign agreement, terminal TFs and
   every node on an upstream-to-observation walk. Compact eligibility requires
   fewer than 500 edges, at least one upstream lipid and two TF observations.
6. Use the ancillary refinement and score-audit entry points described in the
   analysis report for finer extraction and independent final-score diagnostics.
   They read saved fits; they do not change the fitted PKN or upstream sign filter.
7. Export a selected sensitivity with `ancillary/export_lipid_tf_network.py`,
   explicitly supplying the normalized prior, node catalogue and chemical alias
   table. Preserve original signed evidence endpoints and separate input weights,
   inferred MOON scores, observed TF activities and RNA annotations.

## Output and validation contract

`run_summary.csv` preserves every policy/panel/pruning/layer configuration,
including skipped or failed fits. `threshold_grid.csv` preserves empty networks.
`analysis_summary.json` describes the initial coarse grid; later refinement is
reported separately, so a zero coarse-grid candidate count is not a contradiction.
Each successful `runs/<run_id>/` contains the complete fit, compression map,
pruning/coherence ledgers, final scoring and original-ID graphs, and scores.
Input decisions, dependency hashes, run manifest and session information are
saved alongside the grid. The coverage audit has its own SHA256 manifest.

The initial full run on 15 September 2026 evaluated 72 configurations, fitted
16 and extracted 752 coarse threshold pairs. The original PKN yielded neither
direct measured-lipid regulation nor lipid-only routes to observed TFs. See
`docs/cosmos_factor1_lipid_to_tf.md` for the completed refinement, selected
sensitivity, biological interpretation and exact reproduction commands.

The targeted refinement evaluated 1,431 further pairs, with 26 compact
all-TF candidates and none in the strong-TF panel. The documented representative
uses pruning 12, scoring 12 and thresholds 1.65/1.10: 455 edges, 233 nodes,
LPC(18:0/0:0) and 64 TF observations. It is materialized and exported by
`ancillary/refine_lipid_tf_thresholds.R` into `refinement/selected_all_tfs/`.
The independent score audit reproduced all 16 retained score tables exactly.
The normalized-rule topology/pruning audits reproduced identical prepared
graphs, input vectors and compression maps for all six tested configurations;
restoring the disconnected rule layer cannot change these fits.
`ancillary/audit_lipid_tf_selected_reactions.py` covers all 51 selected reaction
nodes/31 Rhea accessions and both prostanoid–Ptger3 regulatory records.

## Known limitations to preserve

Pruning distance, scoring layers and final displayed route lengths are distinct.
Input lipid magnitudes act as retained measurements; the current upstream
constraint checks sign agreement, not a dose-response model. Observed TFs are
level zero and exempt from inferred-score cutoffs. Scores across different
observation sets are not directly comparable measures of confidence.

Lipid-only filtering prevents nonlipid cofactor shortcuts but does not establish
atom transfer, mass balance, flux, rate limitation or the availability of omitted
co-substrates. Audit complete reaction equations before naming an upstream lipid
as the source of a signalling molecule. Acid/base normalization is an identifier
sensitivity, not evidence for equivalent cellular exposure. Bulk liver abundance
does not measure hepatocyte subcellular ligand availability. All resulting routes
remain mechanistic hypotheses, with narrow prior-knowledge coverage.

Run modified R scripts immediately and update this plan, the analysis report and
the static result-query map when input/output semantics change.

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
