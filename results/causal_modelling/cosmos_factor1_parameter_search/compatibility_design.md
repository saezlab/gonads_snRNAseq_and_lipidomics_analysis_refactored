# Systematic parameter-search design and compatibility review

Reviewed on 2026-09-09. This sub-audit inspects the current runner, protected
compression/extraction helpers and saved six-layer results. It does not run a
second copy of the planned depth search or modify the main runner.

## Objective and fixed quantities

Find nonempty, independently audited Factor 1 TF-to-lipid hypothesis networks
with **strictly fewer than 500 edges**, preferring higher extraction thresholds.
This is a readability/selection objective, not an optimization of causal truth
or statistical significance. Preserve every attempted setting and its outcome.

Hold the requested PKN variant, factor ID and orientation, exact lipid input
weights, footprint-derived TF activities, TF input cutoff `abs(activity)>2`, RNA
coherence evidence, ULM statistic and node identifier domain fixed. Read values
from `results/cosmos_factor1_pilot/inputs/factor1_inputs.rds`; do not roundtrip
input scores through CSV. Below-machine-precision changes around zero can alter
the package's sign filter.

Only vary the declared topology/depth policy, compression state and extraction
cutoffs. Do not silently change molecular mappings, include generic lipid-class
matches, use TF mRNA as activity or select TFs by their desired biological names
while searching for a conveniently sized graph.

## Depth and reachability are separate parameters

A practical initial depth grid is **7, 8, 10 and 12**, with the six-layer pilot
retained as a reference. `moon(n_layers=...)` is an upper limit; record the
realized maximum level and node counts at each level. It may stop before the
requested depth or encounter sparse/constant layers that yield nonfinite ULM
scores. Such outcomes are failures or early stops to report, not empty
scientifically meaningful solutions.

The current pilot uses the same number for forward pruning radius, backward
pruning radius and `n_layers`. When all three are increased together, more
observed lipids can enter the feasible subgraph. This changes ULM's measurement
background and potentially scores at **every** level. It is not a clean test of
propagation depth alone.

Two defensible designs are:

1. **Fixed search space:** prune once at the largest declared radii (for example
   12 forward/12 backward), freeze the resulting input ID sets, and compare
   `n_layers=7,8,10,12` on that initial subgraph. Each depth still receives its
   own complete RNA-coherence loop. This best separates propagation depth from
   measurement-background changes.
2. **Tutorial-style joint expansion:** use equal radii and `n_layers` at each
   tested depth, rebuilding from the same full original PKN every time. Label
   this a joint topology/depth search. Save retained TF/lipid ID sets, counts
   and pairwise observation overlap; distinguish gains due to new measurements
   from gains due to extra scored layers.

The first design is preferable for isolating the depth effect. The second is
reasonable for the user's exploratory network-construction objective if its
confounding is made explicit. Do not incrementally mutate the preceding depth's
coherence-filtered PKN: those removals belong to that run, not to all deeper
runs. Likewise, scoring the already filtered six-layer PKN to 12 layers is only
a conditional diagnostic, not a complete 12-layer analysis.

`n_layers` does not impose a maximum length on every path displayed by the
current extractor. The extractor enforces seed-to-observation reachability,
which can include score-consistent cycles and longer detours. If a maximum
display path length is desired, declare and audit it separately; do not claim
that all exported routes contain at most `n_layers` edges.

## Threshold grid grounded in the saved score distributions

Absolute inferred scores from the six-layer saved runs, after the TF sign
filter and in original node IDs:

| Quantile | Protected compression | No compression |
|---|---:|---:|
| Median | 0.672948 | 1.174858 |
| 75% | 0.982048 | 2.762994 |
| 90% | 1.640690 | 2.878271 |
| 95% | 1.976279 | 3.319783 |
| 99% | 2.574165 | 3.968904 |
| Maximum | 8.511728 | 13.429717 |

These are descriptive quantiles, not confidence levels. Original-ID quantiles
repeat shared compressed scores among merged members, so also inspect unique
compressed-node distributions when diagnosing ULM's effective background.

An economical prespecified coarse grid is:

```r
primary <- c(1.7, seq(2, 5, by = 0.25))
secondary <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 3)
grid <- expand.grid(primary = primary, secondary = secondary)
grid <- grid[grid$secondary <= grid$primary, ]
```

This includes the earlier primary/strict neighborhoods and reaches beyond the
uncompressed 99th percentile. Retain the old relaxed `(1,0)` setting as a
coverage reference, rather than ranking it as a preferred high-threshold
candidate. Every graph uses the same `moon_extract_paths()` convention:

- Strict inequalities: an inferred node passes when `abs(score) > cutoff`.
- Level-0 measured weights are exempt from the secondary threshold, but never
  become primary nodes simply because their raw weight is numerically large.
- An edge needs at least one primary inferred endpoint, while its inferred
  endpoints must both pass the secondary threshold.
- Measured lipid nodes are sinks; enforce sign agreement and unconditional
  upstream-to-level-0 reachability after thresholding.

Coarse-grid thresholds are analysis choices, not optimized significance
boundaries. If finer resolution is helpful, predeclare at most one local
refinement pass using 0.05 increments within the adjacent coarse intervals
around an eligible high-threshold frontier. Save the coarse and fine grid
results together. Do not repeatedly tune a threshold immediately below a
particular desired gene's score, or select an epsilon that gives an exact
desired edge count. A compact solution need not be as close as possible to 499
edges.

Apply all extraction settings to the **same cached final score table and PKN**
for a given depth/compression run. Changing a display threshold does not justify
re-running MOON or the coherence loop. Save exact RDS scores for reuse; CSV
exports are for inspection.

## Compression is a sensitivity axis, not just a performance switch

Run protected same-child compression and no compression from the same initial
subgraph, each with its own bounded coherence loop. The installed compressor
must receive canonical `source,target,interaction` columns. Protect complete
signed-child signature groups containing any input or TF-regulon node, use the
complete original-to-compressed map and verify its graph projection, as the
existing safe helper does.

The six-layer audit demonstrates why both representations matter: compressed
level 1 has 103 scoring nodes versus 398 uncompressed endpoint enzymes, changing
the next ULM layer's background and degrees of freedom. Later inferred signs
and score magnitudes can differ substantially. Both fully converged pilot runs
retain 14 compatible TFs, but only 11 are shared; their relaxed edge-set Jaccard
is approximately 0.197. Thus a high-threshold compact graph appearing under only
one compression choice must be presented as representation-sensitive.

Use identical numeric threshold grids in both representations for an auditable
comparison, while acknowledging that the numerical score distributions differ.
Do not treat the same numerical cutoff as the same calibrated error rate. Report
selected networks from each representation and their shared TF/lipid/edge
content. The highest-threshold candidate in one run does not automatically
supersede a lower-threshold but recurring route supported by both.

## Acceptance, ranking and avoiding selective reporting

Hard acceptance criteria:

1. The run uses finite unique input/scored node values and has reached a stable
   coherence graph. If a loop cap is reached, explicitly rescore the final PKN
   and mark the run nonconverged; exclude it from the preferred shortlist.
2. Every edge exists with the same sign in the final original-ID PKN projected
   from the scoring graph. Scores, graphs and compression maps agree.
3. Every retained edge is score-sign-consistent, every retained node can be
   reached from an intended compatible TF and can reach a measured lipid, and
   no measured lipid has an outgoing edge.
4. The final graph has `0 < n_edges < 500`, at least one intended TF and at least
   one measured lipid. Do not silently redefine “under 500” to include 500.

A reproducible ranking convention is to prefer the **higher secondary cutoff**
first, then the higher primary cutoff. The secondary cutoff controls the weakest
allowed inferred node, so this avoids calling a graph broadly high-scoring
merely because a few endpoints have high primary scores. Break ties by more
measured lipids, more TFs, then shallower depth. This is an explicit preference
rule, not a claim that it maximizes biological validity. A primary-first
preference is also possible, but should be declared before looking at which TF
names it selects.

Preserve the within-run threshold Pareto frontier: a cutoff pair is dominated
only if another eligible pair has both thresholds at least as large and one
strictly larger. Do not conceal threshold/coverage tradeoffs behind one selected
result. Flag one-TF/one-lipid results as minimal coverage; retain them as valid
under the user's stated constraint, while additionally showing candidates with
multiple TFs and multiple lipids where available.

Safeguards against selecting a misleading attractive graph:

- Save all attempted settings, including empty, oversized, failed and
  nonconverged runs, together with the declared ranking rule and selected row.
- Use parameter/coverage/structure criteria for selection, not whether Esr1,
  Dbp, Nfic or another desired gene is recovered.
- Report edge/TF/lipid overlap with adjacent grid settings and with the other
  compression representation. A single narrow feasible point warrants less
  confidence than a repeated route across neighboring settings.
- Distinguish the diagram-size threshold search from statistical testing. The
  search does not yield calibrated p-values or establish causal support, and
  exploring many settings is not independent replication.
- Preserve reaction provenance and the previous single-measured-target caveat.
  Higher thresholds do not fix sparse endpoint footprints, duplicate reaction
  representation, cofactor shortcuts or unverified lipid-species specificity.

## Independent summary checker and ranker

`scripts/moon/rank_parameter_candidates.R` provides:

```r
moon_validate_parameter_summary(summary)
moon_rank_parameter_candidates(summary, max_edges = 500L)
```

Required columns are `run_id`, `depth`, `compression`, `primary`, `secondary`,
`edges`, `upstream_tfs`, `downstream_lipids`, `converged` and `audit_passed`.
`run_id` identifies exactly one fixed final scored/coherence-filtered graph.

The validator checks that increasing both thresholds on that same graph cannot
increase retained edges, TFs or measured lipids. Such a violation indicates
mixed run identities, inconsistent extraction or a results-table bug. It never
compares monotonicity across different depth/compression runs. The ranking
helper returns the entire grid with eligibility, rank, within-run threshold
frontier and multiple-TF/multiple-lipid flags. It does not silently select or
discard biological hypotheses.

The script's immediate standalone tests passed for the strict 500-edge boundary,
empty and nonconverged exclusions, declared secondary-first ranking, threshold
frontier flags and an intentionally nonmonotone sweep. Re-run these tests with
`Rscript scripts/moon/rank_parameter_candidates.R`.
