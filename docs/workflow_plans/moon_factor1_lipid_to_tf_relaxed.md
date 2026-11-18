# Active MOON workflow: relaxed Factor 1 lipid-to-TF candidates

## Scope

This workflow tests conditional lipid-to-signalling hypotheses using the
directionality-audited STITCH overlay and a long-form relaxed chemical-mapping
ledger. It is separate from the manuscript TF-to-lipid workflow and from the
original reaction-network lipid-to-TF analysis. It neither refits MOFA nor
rebuilds either PKN.

The active entry point is
`scripts/moon/run_factor1_lipid_to_tf_relaxed.R`, configured by
`scripts/moon/config/factor1_lipid_to_tf_relaxed.json`. Its default output is
`results/causal_modelling/factor1_lipid_to_tf_relaxed/`.

The pre-cleanup discovery topology and search are preserved under
`miscellaneous/older_results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory/`. That folder
is a legacy audit record, not the active result directory or a current
reproduction target.

The selected configuration is the strong-TF sensitivity with pruning distance
4, four recursive score layers and primary/secondary thresholds of 0.5/0.5. It
expects the conditional `LPS(18:0)` lysophosphatidylserine candidate and the
supplied `Ar` and `Esr1` terminal activities. The configuration records hashes
of the legacy reference fit and strict solution so that an active replay can be
checked without treating the legacy directory as its output.

## Analysis contract

```text
Factor 1 inputs + directional STITCH overlay + relaxed mapping ledger
  -> direct chemical-to-protein + protein-to-protein graph
  -> compressed MOON fits with TF activities as terminal observations
  -> candidate score and duplicate-fingerprint audit
  -> parameter search with score-layer-constrained route extraction
  -> evidence audit and network export
```

Workflows 00–10 do not consume these outputs.

The following points must remain explicit in all code, figures and prose:

- Each mapping row is a conditional measurement-to-chemical hypothesis. It is
  not an asserted lipid identity. Broad class, ancestor and low-resolution
  mappings are annotation-only.
- Chemical candidates are topology seeds. Factor 1 lipid weights are copied
  into the mapping ledger only after scoring, to assess assignment-level sign
  concordance. They are not numeric ULM inputs, and values from multiple
  measurements are never summed or averaged.
- The graph contains only direct `Metab__chemical -> protein` and
  `protein -> protein` edges. It excludes reactions and
  metabolite-to-metabolite edges, preventing precursor-chain routes from being
  interpreted as lipid signalling.
- TF activities are supplied terminal observations. In particular, AR and ESR1
  occurring at route termini are not independently predicted validation
  targets.
- STITCH v5 is historical and exploratory. Explicit edge direction improves
  auditability but does not establish a liver mechanism.

## Active components

| Component | Responsibility | Default destination |
| --- | --- | --- |
| `run_factor1_lipid_to_tf_relaxed.R` | Validate frozen inputs, construct the direct regulatory graph, refit the selected compressed configuration, validate the strict route against the saved legacy reference, and write model, score and provenance artifacts. | `results/causal_modelling/factor1_lipid_to_tf_relaxed/` |
| `config/factor1_lipid_to_tf_relaxed.json` | Defines the selected analysis inputs, graph policy, compression/scoring settings and output contract. | — |
| `ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R` | Reuse the active fitting contract across pruning, depth and dual-threshold settings; retain every tested outcome and materialize score-layer-constrained candidates. | `results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search/` |
| `config/factor1_lipid_to_tf_relaxed_parameter_search.json` | Defines the ancillary search grid and compact-network eligibility policy. | — |
| `ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R` | Check source hashes, retain every feature-to-candidate assignment, report candidate-local MOON score traces with topology flags, and deduplicate direct-edge fingerprints. | `<selected analysis>/candidate_audit/` |
| `ancillary/audit_factor1_lipid_to_tf_relaxed_score_layer_evidence.py` | Attach PKN evidence to every materialized strict score-layer solution without changing a fit. | `<parameter search>/evidence_audit/` |
| `ancillary/export_factor1_lipid_to_tf_relaxed_network.py` | Export the selected score-layer-constrained network, attributes, provenance and graph formats. | Analysis-specific export directory |
| `scripts/pkn/audit_factor1_relaxed_regulatory_mapping.py` | Validate and regenerate the long-form mapping ledger independently of MOON fitting. | `results/causal_modelling/factor1_relaxed_lipid_regulatory_mapping/` |

The active main workflow contains the single documented analysis configuration.
The parameter search is ancillary evidence: it does not silently replace the
main configuration or label a threshold pair as statistically optimal.

The current ancillary grid tests pruning distances 4, 6 and 8; recursive
depths 2, 4 and 6; all-TF and strong-TF panels; primary thresholds 0.5–3.0; and
secondary thresholds 0–2.0. It materializes every strict score-layer network
with fewer than 500 edges, a nonzero secondary cutoff, at least one independent
direct-edge fingerprint group and at least two observed TFs, then collapses
duplicate signed-edge topologies only for the stability summary.

## Execution

Run from the repository root. Start by checking the mapping and active
configuration:

```sh
python3 scripts/pkn/audit_factor1_relaxed_regulatory_mapping.py --self-test
Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R --validate-only
Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R
```

Use a fresh explicit output directory for a replay. The main workflow refuses
to overwrite a nonempty result directory. The ancillary search, candidate
audit, evidence audit and exporter accept their own analysis/output arguments.
The evidence audit defaults to
`results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search/evidence_audit/`
and attaches PKN records only to the materialized strict search solutions. Use
`--help` before invoking a new analysis directory.

```sh
Rscript scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R --help
Rscript scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R --help
python3 scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_score_layer_evidence.py --help
python3 scripts/moon/ancillary/export_factor1_lipid_to_tf_relaxed_network.py --help
```

After changing an R component, immediately run its validation or a full run,
then rerun every downstream audit/export affected by the change. Update this
plan, `scripts/moon/README.md`, the analysis report and the static query map
when paths, outputs, selection rules or caveats change.

## Discovery record and interpretation limits

The legacy discovery folder recorded 19 conditional assignments from 12 Factor
1 features to 18 chemical source nodes. Its broad 405-edge sign-consistent
network included duplicate FA(18:2) chemical candidates with an identical
direct-edge fingerprint. It is a connectivity view only: nonconsecutive MOON
score levels mean it cannot explain a recursive ULM score.

The score-layer audit of that discovery run identified a repeatable strong-TF
sensitivity, which the active main workflow refits and checks against the saved
reference:

```text
LPS(18:0) candidate lysophosphatidylserine
  -> Gpr34 -> Gnai1 -> Src -> {Ar, Esr1}
```

It was recovered only at primary/secondary thresholds of 0.5/0.5; no compact
score-layer solution survived when either threshold was higher than 0.5. The
route is therefore a low-threshold, conditional signalling hypothesis. Its LPS
mapping is compatible-more-specific acid/base mapping and does not resolve sn
position. AR and ESR1 are supplied downstream observations, so their presence
does not validate the route independently.

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
