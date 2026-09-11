# Factor 1 TF-to-lipid MOON pilot

Run date: 2026-09-08. Input prior: `results/cosmos_mouse_pkn/pkn_lipid_rules.rds`.

> The expanded depth and threshold search in
> `docs/cosmos_factor1_moon_parameter_search.md` supersedes this pilot for
> candidate selection. This document remains the detailed record of input
> construction, the initial compression diagnostic and the first six-layer run.

## Outcome

The pilot runs successfully and produces signed TF-to-lipid candidate networks,
but the result depends materially on network compression and extraction
thresholds. The compressed analysis following the requested example yields no
complete network at the example's primary/secondary thresholds of 1.7/1.0.
An independently converged **uncompressed** analysis at the same thresholds
retains **8 TFs, 4 measured lipids, 473 nodes and 1,676 edges**. These are
exploratory hypotheses; this test does not establish a robust TF-specific
mechanism.

Both analyses follow the six-layer ULM and iterative TF-target coherence
structure of the pinned [COSMOS_basic example](https://github.com/saezlab/COSMOS_basic/blob/7108c6b849c712b1282440dc122c9d783f563da8/scripts/net_compr_MOON.R),
using installed cosmosR 1.21.2 and decoupleR 2.17.0. No CARNIVAL optimization
or external solver was used.

## Inputs and coverage

- Upstream measurements are the **hepatocyte TF footprint activities inferred
  from Factor 1 RNA loadings**, not TF transcript loadings. The 367 saved TF
  activities were reproduced with maximum absolute difference 2.26e-12.
  Following the example, 51 pass `abs(activity) > 2`; 49 are present in this PKN
  and remain after topology pruning. This cutoff is exploratory, not an FDR
  criterion. Fifteen of the 367 original footprint tests have BH-adjusted
  p-value below 0.05; these p-values are not mechanism probabilities.
- Downstream measurements start from **all 1,254 Factor 1 lipid loadings**.
  We accept 155 unique feature-to-node mappings: 149 explicit measured-species
  rule nodes and six unambiguous ChEBI matches at the reported resolution.
  Class/hierarchy and merely compatible structures are excluded. The input
  audit records exclusions for every feature.
- Separate forward and reverse pruning radii of six retain **74 lipid nodes**.
  The mapped 155 include seven of the previously selected 64 correlation-network
  lipids, but this pilot is not restricted to those 64 features.
- All accepted observations have unknown compartment; none were copied across
  compartments. Bulk liver lipid measurements are paired with a hepatocyte RNA
  signature, so cell-type attribution of the lipid signal remains an assumption.
- All 20,123 hepatocyte RNA loadings are retained separately for TF-target
  coherence. Loading presence or magnitude was not treated as detectable gene
  expression. The tutorial's legacy expression filter is incompatible with this
  PKN's reaction-instance identifiers, so no expression-presence filter was run.

The fitted model started with 12 factors and retained ten. Saved Factor 1
matches native HDF5 factor 1 for both hepatocyte RNA and lipids. It is already
testes-oriented: mean sample coordinates are -0.7423 for ovaries and +0.7086 for
testes, so all inputs use a common sign multiplier of +1. Positive loading means
agreement with this axis, not a directly estimated differential effect.
Factor 1 explains **8.906% of hepatocyte RNA variance** and **29.958% of lipid
variance**. Its gonadal-group association has R² = 0.969836, a different quantity
from variance explained in the omics views. Existing fits with different initial
factor counts show strong coordinate agreement, but no bootstrap or repeated-seed
stability analysis was performed.

See `results/cosmos_factor1_pilot/inputs/README.md` for input schemas, source
hashes, model alignment, TF target counts, mapping exclusions and reproduction.

## Scoring and extraction

The 254,589-edge PKN is cleaned and observations are made terminal before
reachability pruning, matching `moon()`'s scoring convention. The pruned graph
has 44,127 edges and 6,673 nodes. Protected compression merges 2,720 original
nodes into 378 groups, producing 36,569 edges; 2,591 input/regulon identifiers
are explicitly protected. A complete original-to-compressed mapping is saved.

The compressed run converges after five coherence iterations, with 26,776
compressed edges in the final scoring graph. It scores all six layers and
retains 14 of the 49 upstream TFs with matching signs. Decompression produces
6,510 scored original nodes, including all 74 retained measurements.
The uncompressed run starts from the same initial pruned graph and performs
its **own** coherence loop. It converges after four iterations, with 34,973
edges in its final scoring graph and 6,514 scored nodes.

Extraction uses strict `>` comparisons. Each edge must touch an inferred node
above the primary cutoff; both inferred endpoints must exceed the secondary
cutoff. **Measured level-0 lipid weights are exempt from the secondary cutoff
and never promoted to primary nodes based on raw weight magnitude.** This
preserves the requested measurements without mixing loading and ULM score
scales. Every retained edge must agree with the node-score signs, and every
retained node must be reachable from an intended TF and able to reach an
observation. Outgoing measurement edges are excluded throughout.

| Analysis | Primary / secondary | TFs | Lipids | Nodes | Edges |
| --- | ---: | ---: | ---: | ---: | ---: |
| Compressed, tutorial thresholds | 1.7 / 1.0 | 0 | 0 | 0 | 0 |
| Compressed, relaxed | 1.0 / 0.0 | 14 | 10 | 1,476 | 6,320 |
| Compressed, intermediate sensitivity | 1.7 / 0.5 | 9 | 3 | 505 | 1,552 |
| Uncompressed, tutorial thresholds | 1.7 / 1.0 | 8 | 4 | 473 | 1,676 |
| Uncompressed, stricter extraction | 2.0 / 1.5 | 3 | 2 | 187 | 537 |
| Uncompressed, relaxed | 1.0 / 0.0 | 14 | 11 | 3,647 | 13,246 |

The uncompressed primary network contains **Klf10, Nfe2l2, Nr0b2, Nr2f1,
Nr4a2, Rarb, Tbx15 and Tbx2**, connected to **DG(20:0_18:1), DG(24:0_18:2),
DG(24:1_18:2) and PS(20:0_22:6)**. The stricter extraction retains Klf10,
Nfe2l2 and Tbx2, and the DG(24:0_18:2) and DG(24:1_18:2) observations. Retention in a
thresholded graph is not independent evidence that these TFs control those
species.

Only three uncompressed primary TF-to-lipid routes are at most six edges, all
following **Tbx15 → Igf1r → Pdpk1 → Prkcd → Smpd1 → species-rule reaction → DG**.
They reach the three DG species above. Their weakest inferred score is 1.4385,
versus 0.4922 for the same routes under compression. The other 29 primary pairs
need longer routes in the union graph. The stricter uncompressed network has
six TF/lipid pairs, all longer than six edges. These explicit route audits
prevent interpreting six scoring layers as six-edge explanations for every TF.

## Why the result needs caution

**Compression changes the statistical problem.** It preserves signed-child
topology with the protections implemented here, but it changes the number and
distribution of downstream variables presented to later ULM layers. It is not
only a speed improvement. A diagnostic on the same frozen final graph showed
substantial score/sign differences; the full independent comparison above
confirms that complete networks depend on this choice. The diagnostic
`uncompressed_sensitivity_scores.csv` at the root of the results folder is a
fixed-graph comparison; the independently converged analysis is in the
`uncompressed_sensitivity/` subdirectory. They should not be conflated.
The two full runs retain 14 sign-compatible TFs each, with 11 shared. Their
relaxed networks share ten of 11 lipid endpoints but have an edge Jaccard
similarity of only 0.197; endpoint recurrence is stronger than agreement on the
connecting mechanism. These comparisons are in `fullcomparison.json` and the
accompanying CSVs.

**Endpoint support is sparse.** All 398 enzyme-instance parents immediately
upstream of the 74 measured lipids have only one measured lipid target each.
ULM scores these relative to the entire input signature; the score is not
replicated evidence for enzyme substrate specificity. In the compressed run,
the strong terminal reaction scores for the large positive DG/PS loadings
still pass through weak intermediates: Smpd1 has score 0.8673, while Ptdss1
and a serine intermediate have score 0.8319. The secondary cutoff of one
therefore breaks those paths.

**Dense reachability is not TF specificity.** In the compressed relaxed graph,
all 14 TFs connect to all ten endpoints: 140 reachable pairs. Only 59 pairs
have a route of at most six edges. Restricting routes to gene/protein
intermediates followed by an annotated catalyst → reaction instance → measured
product retains seven endpoints; requiring at most six edges retains 53 pairs.
The three PS endpoints disappear under this route restriction because their
retained mixed routes depend on metabolic intermediates such as serine.
The union of those 53 shorter catalyst routes is an 81-node, 101-edge inspection
network; its scores are copied from the relaxed analysis, not refitted.

**Reaction rules remain hypotheses.** The measured-species edges transfer
generic class reaction templates to compatible reported chains. A path through
an Smpd1- or Ptdss1-annotated rule does not establish experimentally demonstrated
activity on the exact measured lipid species. Reaction traversal is qualitative:
positive substrate → reaction → product edges do not encode mass balance,
substrate depletion, cofactor sufficiency, flux or enzyme kinetics.

**The RNA coherence screen is internal consistency.** It uses the same RNA
loadings that generated the TF footprints. MOON's upstream vector filters
inferred TF signs; its magnitudes do not act as forcing strengths or require
every supplied TF to appear. The resulting ULM scores are not causal effects
or mechanism-level p-values. Six scoring layers do not guarantee that every
route in their union is at most six edges; route lengths are exported explicitly.
Esr1, Dbp and Nfic are all inferred before upstream sign filtering in both
analyses, but their inferred signs oppose their negative RNA-derived TF
activities. Their absence is therefore a sign conflict, not missing PKN coverage.

## Files to inspect

All paths below are relative to `results/cosmos_factor1_pilot/`:

- `uncompressed_sensitivity/solution_primary_{SIF,ATT}.csv`: independently
  converged uncompressed network at the tutorial thresholds; the integrated
  Cytoscape file is `uncompressed_sensitivity/solution_primary.graphml`.
- `uncompressed_sensitivity/solution_strict_{SIF,ATT}.csv`: stricter candidate
  subnetwork from the same uncompressed scores.
- `solution_primary_{SIF,ATT}.csv`: the explicitly empty compressed primary
  result; retained to avoid disguising threshold sensitivity.
- `solution_relaxed_{SIF,ATT}.csv` and `.graphml`: full compressed exploratory
  network, with measured and inferred values separately annotated.
- `solution_catalyst_bounded6.graphml` and associated SIF/ATT tables: compact
  union of shorter catalyst routes from the compressed relaxed network.
- `candidate_paths*`, `network_evidence*` and `figures/`: auditable paths,
  exact-edge provenance, reaction assumptions and visualizations. Corresponding
  uncompressed outputs are in its subdirectory.
- `moon_scores*`, `iterations/`, `pkn_final_*`, `compression_mapping.csv`,
  `coherence_ledger.csv`, `run_manifest.json`: all scores, exact scoring graphs,
  exclusions, configuration and software provenance. Each independent run has
  its own copies.
- `score_support_audit.csv`, `strong_endpoint_predecessors.csv`,
  `tf_unconstrained_compression_audit.csv`, and `score_support_README.md`: why
  endpoint scores and compression affect the result.

## Reproduce

Run from the repository root. The input README documents the Python model-audit
dependencies; the scoring stage uses installed R packages and frozen local data.

```sh
Rscript scripts/moon/prepare_factor1_inputs.R
python3 scripts/moon/audit_factor1_model.py
Rscript scripts/moon/map_factor1_inputs.R
Rscript scripts/moon/moon_safe_helpers.R
Rscript scripts/moon/run_factor1_moon.R
Rscript scripts/moon/audit_factor1_thresholds.R
Rscript scripts/moon/run_factor1_moon.R --no-compression
Rscript scripts/moon/audit_score_support.R
Rscript scripts/moon/compare_factor1_runs.R
python3 scripts/moon/export_factor1_network.py --stem solution_primary
python3 scripts/moon/export_factor1_network.py --stem solution_relaxed
python3 scripts/moon/export_factor1_network.py --stem solution_intermediate
python3 scripts/moon/export_factor1_network.py --results results/cosmos_factor1_pilot/uncompressed_sensitivity --stem solution_primary
python3 scripts/moon/export_factor1_network.py --results results/cosmos_factor1_pilot/uncompressed_sensitivity --stem solution_strict
```

The helper's synthetic tests cover protected duplicate groups, actual and zero
compression, complete decompression, separate score/weight thresholds, wrong-sign
upstream nodes, terminal observations, dead-end cycles and disconnected coherent
components. Both actual runs pass finite/unique-score, preserved-observation,
exact-final-PKN, sign-consistency and TF-to-observation reachability checks.
The evidence exporter independently checks every exported edge against the
requested PKN and its provenance. Existing workflows and analyses were not rerun
or modified.
