# Relaxed Factor 1 lipid-to-TF candidate screen

This exploratory screen asks whether a measured Factor 1 lipid can be assigned
to a plausible chemical feature with direct protein-regulatory prior knowledge,
and whether the resulting candidate is compatible with supplied hepatocyte TF
activities. It produces conditional signalling hypotheses rather than evidence
of lipid identity, hepatic exposure, receptor occupancy, temporal order or
causal regulation.

The cleaned workflow is organised around one active analysis and separate
ancillary audits:

| Purpose | Active path |
| --- | --- |
| Selected analysis | `scripts/moon/run_factor1_lipid_to_tf_relaxed.R` |
| Selected-analysis configuration | `scripts/moon/config/factor1_lipid_to_tf_relaxed.json` |
| Default selected-analysis output | `results/causal_modelling/factor1_lipid_to_tf_relaxed/` |
| Parameter search | `scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R` |
| Parameter-search configuration | `scripts/moon/config/factor1_lipid_to_tf_relaxed_parameter_search.json` |
| Default parameter-search output | `results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search/` |
| Candidate audit | `scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R` → `<selected analysis>/candidate_audit/` |
| Score-layer evidence audit | `scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_score_layer_evidence.py` → `<parameter search>/evidence_audit/` |
| Network exporter | `scripts/moon/ancillary/export_factor1_lipid_to_tf_relaxed_network.py` |

The earlier directory
`miscellaneous/older_results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory/` is retained
as a legacy pre-cleanup discovery topology/search. It should be consulted to
audit the original exploration, but it is not the active result destination.
The workflow plan in
`docs/workflow_plans/moon_factor1_lipid_to_tf_relaxed.md` defines the current
execution contract.

The selected active configuration refits the predeclared LPS sensitivity at four
pruning steps, four recursive score layers and 0.5/0.5 primary/secondary
thresholds using the strong-TF panel. It compares the refit with the saved
legacy score state, compression map and strict solution. It records the expected
LPS candidate and the supplied AR/ESR1 observations as audit checks; it does not
turn those observations into independent validation.

## Input and PKN policy

The mapping ledger records feature-to-candidate possibilities in long form. A
measurement can map to multiple candidate chemicals, but each candidate source
is scored once; its copied Factor 1 value is then evaluated separately for every
mapping row. This retains ambiguity without treating duplicated chemistry as
independent measurements.

Factor 1 lipid weights are not numeric MOON inputs. In this direction,
`cosmosR::moon()` computes recursive ULM scores from the terminal TF-activity
observations, then the lipid weights are used after scoring as sign-concordance
annotations. The TFs at route termini are therefore supplied observations, not
independent predictions.

The direct signalling graph includes only candidate chemical-to-protein and
protein-to-protein edges from the directional-STITCH overlay. It excludes
reaction and metabolite-to-metabolite edges, so a reported route cannot arise
from a long metabolic precursor chain. Directional STITCH records remain a
historical, exploratory knowledge source; their explicit direction permits an
audit, not a mechanistic conclusion.

## Selected low-threshold sensitivity and discovery record

The pre-cleanup discovery analysis contained 19 conditional assignments from 12
Factor 1 lipid features to 18 candidate chemical nodes. No assignment was a
strict direct identity. Duplicate FA(18:2) candidates had the same direct-edge
fingerprint and score, showing why they must be treated as one unresolved
candidate family rather than independent support.

A broad 405-edge sign-consistent graph from that exploration is retained only
as topology. Some of its edges skip recursive MOON score levels, so it is not a
score explanation and must not be used as a causal route display.

Applying consecutive descending score-layer, sign-consistency and
candidate-to-observed-TF walk requirements recovered a repeatable strong-TF
sensitivity:

```text
LPS(18:0) candidate lysophosphatidylserine
  -> Gpr34 -> Gnai1 -> Src -> {Ar, Esr1}
```

The direct LPS-to-GPR34 relation is carried by a GuideToPharmacology interaction
record; the protein-signalling edges are signed OmniPath records projected to
mouse. The `LPS(18:0)` feature maps conditionally to a more-specific
lysophosphatidylserine candidate through compatible-more-specific acid/base
logic. It does not resolve sn position.

This route is an informative hypothesis, but its support is limited. It
appeared only at the lowest tested nonzero primary and secondary thresholds
(0.5/0.5); no compact score-layer solution remained when either threshold was
raised above 0.5. AR and ESR1 are terminal TF activities supplied to the model,
so their recovery is not independent validation of LPS signalling in
hepatocytes.

The active selected workflow regenerates this exact strict route as a six-node,
five-edge network and verifies it against the saved discovery reference. Its
separate `score_explanation/` output is a candidate-local package trace; it is
retained for score auditing and is not substituted for the strict route display.

## Reproduction and audit sequence

Run the selected analysis from the repository root after validating the mapping
ledger and configuration:

```sh
python3 scripts/pkn/audit_factor1_relaxed_regulatory_mapping.py --self-test
Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R --validate-only
Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R \
  --output=results/causal_modelling/factor1_lipid_to_tf_relaxed_replay
Rscript scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R \
  --analysis=results/causal_modelling/factor1_lipid_to_tf_relaxed_replay
```

Run the parameter search separately when assessing depth and dual-threshold
sensitivity; give it a fresh output directory when the default already exists.
The evidence audit attaches PKN records to every materialized strict solution,
not to broad sign-consistent connectivity or candidate-local traces:

```sh
Rscript scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R \
  --output=results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search_replay
python3 scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_score_layer_evidence.py \
  --search=results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search_replay
```

Each script records its supported options through `--help`.

The active output directories above are generated by the cleaned scripts. The
legacy discovery folder is retained to audit the original exploration and must
not be silently substituted for a regenerated active analysis.
