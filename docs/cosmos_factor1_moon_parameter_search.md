# Factor 1 TF-to-lipid MOON parameter search

Historical analysis record. Result folders were archived unchanged on 16 November 2026. See [the causal analysis map](CAUSAL_ANALYSIS_MAP.md) for current entry points and the frozen-reference dependencies. Commands using archived implementations retain historical assumptions.

> **Workflow cleanup, 2026-09-15:** The active selected-result entry point is `scripts/moon/run_factor1_tf_to_lipid.R`; see `scripts/moon/README.md`. Its refit from the frozen mapped inputs and PKN reproduced the selected 388-edge/153-node graph, scores and expected-TF exclusions, with evidence export and plots. The exploratory scripts below are preserved under `scripts/moon/archive/2026-09/`; their original path assumptions remain historical. Old `results/cosmos_*` references below identify files now under `results/causal_modelling/cosmos_*`. Historical command lists are provenance, not current reproduction commands.

> **Recommendation superseded by the [compressed-network review](cosmos_factor1_compressed_review.md) (2026-09-09).** The U5/U4 runs below genuinely bypassed same-children compression. Their higher cutoffs do not justify preferring them over compressed inference. Keep them as historical representation diagnostics. The new primary exploratory candidate uses compression, radius 6/depth 4, cutoffs 1.85/0.55 and 388 edges; a compressed radius-12/depth-8 alternative has 468 edges at 2.15/1.075. The new review also establishes that Esr1, Nfic and Dbp were removed by post-score sign filtering, not pre-score pruning.

Run date: 2026-09-09. This extends the initial Factor 1 MOON pilot using
`results/causal_modelling/cosmos_mouse_pkn/pkn_lipid_rules.rds` and the frozen, exact-mapped
Factor 1 input bundle in `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/`. It is an
extraction and score-propagation sensitivity analysis; it does not refit MOFA,
change the PKN, or introduce new measurements.

## Historical outcome of the initial search

The search produces a small, auditable family of signed TF-to-lipid hypothesis
networks with strictly fewer than 500 edges. The historical threshold-priority
representative was the uncompressed, five-layer, radius-6 network extracted at
primary/secondary thresholds 2.05/1.60:

| Role | Representation | Primary / secondary | Nodes | Edges | Retained TFs | Measured lipids |
| --- | --- | ---: | ---: | ---: | --- | --- |
| Threshold-priority representative (U5) | Uncompressed, scoring depth 5, topology radius 6 | 2.05 / 1.60 | 123 | 318 | Klf10, Nfe2l2 | DG(24:0_18:2), DG(24:1_18:2) |
| Compact sensitivity companion (U4) | Uncompressed, scoring depth 4, topology radius 6 | 1.90 / 1.60 | 83 | 208 | Nfe2l2, Nr4a2 | DG(24:0_18:2), DG(24:1_18:2) |
| Compression sensitivity (C4) | Compressed, scoring depth 4, topology radius 6 | 2.00 / 0.85 | 110 | 255 | Ppara, Tbx2 | DG(24:0_18:2), DG(24:1_18:2) |

The original ranking selected U5 by maximizing numerical thresholds across
representations. That preference is withdrawn: compression changes redundant
information and the regression background, so thresholds do not establish
comparable evidence strength across representations. U4/U5 remain diagnostic
comparisons; C4 is the compressed baseline C01 in the expanded review. The
updated broader compressed candidate is C03 at 1.85/0.55 with 388 edges.

The raw threshold-first frontier also contains a 282-edge uncompressed depth-5
extraction at 1.95/1.65, but it retains only one TF. It remains in the complete
grid and shortlist; it was not allowed to define the multi-endpoint
representative.

This is a reproducible downstream reaction hypothesis, not a robust,
TF-specific causal model. The three representative graphs share the same two
DG observations and the `Pdpk1 -> Prkcd -> Smpd1 -> reaction instance -> DG`
tail, but they have no common upstream TF endpoint. The PKN records this tail
as a Smpd1-associated, LINEX2 rule-instantiated PC-to-DG projection of the
generic `RHEA:10605` template. It does not establish Smpd1 activity on these
exact acyl-chain species, flux, substrate availability, or tissue-specific
enzyme activity.

## Search design

The upstream vector contains the 49 supplied Factor 1 hepatocyte TF
footprint activities that remained after terminal-observation topology pruning.
The downstream vector contains the 74 exact-mapped measured lipid nodes that
remained in each prepared topology. Factor 1 RNA loadings were used only for
iterative TF-target coherence. Measured lipids were terminal level-0
observations and were exempt from inferred-score thresholds.

The topology radius was held fixed within each scoring-depth sweep, so that
`n_layers` changes MOON propagation rather than the measurement background.
The search completed 22 independently converged scoring configurations:

- Radius 6: scoring depths 4, 5, 6, 7, 8, 10 and 12, each uncompressed and
  compressed.
- Radius 7: the same scoring depths uncompressed, as a minimally wider
  topology sensitivity.
- Radius 12: uncompressed depth 6, as a deliberately broad,
  cycle-dominated topology stress test.

For each run, the coarse grid used primary thresholds 1.70 and 2.00--5.00 in
0.25 increments, secondary thresholds 0.50, 0.75 and 1.00--3.50 in 0.25
increments, and only pairs where secondary was no greater than primary. Around
the strongest eligible configuration per run, a prespecified local grid sampled
both thresholds in 0.05 increments over plus/minus 0.25. All 3,792 extraction
rows are retained in `parameter_grid.csv`; repeated thresholds within a run
reuse the exact saved scores and PKN.

Eligibility requires convergence, all extraction audits, at least one TF and
one lipid, and `0 < edges < 500`. A multi-endpoint shortlist additionally
requires at least two TFs and two lipids. The ranking prioritizes the secondary
threshold, then the primary threshold, lipid coverage, TF coverage and fewer
edges. This makes the preference for stricter thresholds explicit while
preventing a single TF-to-lipid path from defining the local refinement window.

## Depth and topology findings

The terminal-observation topology audit found that topology depth 5 reaches all
49 retained TFs, all 74 retained lipids and all 3,626 reachable TF-lipid pairs.
Increasing topology radius from 5 through 14 adds no new pairwise reachability
or shorter paths; almost all added edges occur in directed cycles. Greater
topology radius is therefore not equivalent to more biological coverage here.

Only scoring depths 4 and 5 yield any nonempty sub-500-edge extraction in these
original grids; the expanded compressed review supersedes these sampled minima. The smallest nonempty network at uncompressed depths 6 or greater
has 537 edges; the compressed depth-5 minimum is 739 edges and compressed depth
6 or greater starts at 852 edges. Depths 8, 10 and 12 all saturate at realized
scoring level 7. The wide radius-12 depth-6 stress test likewise reaches a
minimum of 537 edges, so it produces no eligible solution under the strict cap.

The radius-7 sensitivity exactly reproduces the radius-6 SIF files for U4 and
U5: their SHA-256 hashes are identical at their respective thresholds. The
extra radius-7 cycle context therefore did not alter either selected solution.

## Candidate comparison and interpretation

All exported candidates passed independent checks against their exact saved
states: unique node and edge identities, exact final-PKN membership, signed
score consistency, primary and secondary threshold adherence, terminal measured
lipids, intended-TF sign agreement, convergence, and a directed TF-to-lipid
walk through every retained node.

U4 and U5 share 78 nodes and 188 signed edges (edge Jaccard 0.556), whereas U4
and C4 share 70 signed edges (edge Jaccard 0.178). The U4/U5 comparison is
particularly informative: all 2,161 inferred nodes common to their saved score
tables have exactly identical scores and levels. Depth 5 adds level-5 nodes and
additional surviving routes; it does not improve the score of a shared node.
At matched thresholds 1.90/1.60, the 208 U4 edges are a strict subset of the
363 U5 edges, with no additional measured lipid endpoint.

Klf10 is retained in U5 with an inferred score of -1.600186, only 0.000186
above the secondary cutoff. Nfe2l2 is 0.028999 above that cutoff in both U4 and
U5. This is a threshold-margin diagnostic, not a confidence interval, but it
means Klf10 must be reported as especially threshold-sensitive. The varying TF
endpoints should be treated as candidates for follow-up, not as a consensus
regulator list.

The retained shortest TF-to-DG routes are 7--10 edges in U4, 9--11 edges in U5
and nine edges in C4; no representative has a route of six or fewer edges.
`n_layers` is the depth of MOON score propagation, not a maximum path length in
the extracted network. Do not describe these as four-step or five-step
mechanisms.

## Primary artifacts

The historical uncompressed U5 network is available as:

- `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth5_uncompressed/solution_p2_05_s1_60_SIF.csv`
- `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth5_uncompressed/solution_p2_05_s1_60_ATT.csv`
- `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth5_uncompressed/solution_p2_05_s1_60.graphml`
- `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth5_uncompressed/candidate_paths_solution_p2_05_s1_60.csv`
- `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth5_uncompressed/network_evidence_solution_p2_05_s1_60_reactions.csv`

The compact U4 counterpart is in the analogous
`runs/radius6_depth4_uncompressed/solution_p1_90_s1_60_*` files. Its
illustrative route figure is
`runs/radius6_depth4_uncompressed/figures/candidate_paths_solution_p1_90_s1_60_network.pdf`.

The exhaustive grid, cross-configuration rankings and plots are in
`miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/summary/`. The independent validation
of all saved radius-6 representatives is in
`miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/refined_candidate_validation.md`; the
direct U4/U5/C4 comparison is in
`miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/refined_candidate_audit/refined_candidate_comparison.md`.

## Historical execution commands

The following records the historical experiment sequence. Archived scripts retain their original path dependencies and are not active reproduction commands:

```sh
Rscript scripts/moon/archive/2026-09/audit_factor1_depth_topology.R
Rscript scripts/moon/archive/2026-09/search_factor1_parameter_space.R --phase=coarse --topology-radius=6
Rscript scripts/moon/archive/2026-09/search_factor1_parameter_space.R --phase=refine --topology-radius=6
Rscript scripts/moon/archive/2026-09/search_factor1_parameter_space.R --phase=coarse --topology-radius=7 --compression=uncompressed
Rscript scripts/moon/archive/2026-09/search_factor1_parameter_space.R --phase=refine --topology-radius=7 --compression=uncompressed
python3 scripts/moon/archive/2026-09/summarize_factor1_parameter_search.py --grid miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/parameter_grid.csv
Rscript scripts/moon/archive/2026-09/audit_refined_candidates.R
python3 scripts/moon/archive/2026-09/export_factor1_network.py --results miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth5_uncompressed --stem solution_p2_05_s1_60
```

The historical radius-12 stress test is already in the saved grid; its original coarse-run options were `--topology-radius=12 --depths=6 --compression=uncompressed`.
