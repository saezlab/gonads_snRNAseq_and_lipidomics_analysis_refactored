# Factor 1 MOON: compressed alternatives and TF retention audit

Historical analysis record. Result folders were archived unchanged on 16 November 2026. See [the causal analysis map](CAUSAL_ANALYSIS_MAP.md) for current entry points and the frozen-reference dependencies. Commands using archived implementations retain historical assumptions.

> **Workflow cleanup, 2026-09-15:** The active selected-result entry point is `scripts/moon/run_factor1_tf_to_lipid.R`; see `scripts/moon/README.md`. Its refit from the frozen mapped inputs and PKN reproduced the selected 388-edge/153-node graph, scores and expected-TF exclusions, with evidence export and plots. The exploratory scripts below are preserved under `scripts/moon/archive/2026-09/`; their original path assumptions remain historical. Old `results/cosmos_*` references below identify files now under `results/causal_modelling/cosmos_*`. Historical command lists are provenance, not current reproduction commands.

Completed 2026-09-09. This review supersedes the representation preference in
[the previous parameter-search report](cosmos_factor1_moon_parameter_search.md).
All previous outputs remain available as historical sensitivity results.

## Decision

The previously favored U5 network really was scored **without same-children
compression**. Its higher numerical cutoffs do not justify preferring it over
a compressed result. Use compressed networks as the primary analysis and keep
the uncompressed results as representation diagnostics.

For biological inspection, start with **baseline C03: 388 edges, topology
radius 6, scoring depth 4, primary/secondary cutoffs 1.85/0.55**. It connects
Ahr, Hnf4a, Ppara and Tbx2 to three measured DG species. Stronger compression
reproduces exactly the same graph at these cutoffs. A separate compressed
radius-12/depth-8 result has 468 edges at 2.15/1.075, but changes the TF endpoints
to Klf10 and Tbx2; its higher cutoffs do not establish greater biological validity.

**Esr1, Nfic and Dbp were not pruned before scoring.** They survive topology
preparation and compression, acquire weak positive lipid-driven scores, and
are subsequently omitted by MOON's upstream-sign filter because their supplied
Factor 1 footprint activities are negative. None of the tested reasonable
compression, topology or RNA-coherence alternatives resolves that conflict.

## What was actually compressed

The previous uncompressed branch in `scripts/moon/archive/2026-09/search_factor1_parameter_space.R`
passes every node as protected to `moon_safe_compress()`. The wrapper finds no
eligible groups and never calls `cosmosR::compress_same_children()`. The scored
graph is exactly the original 44,127-edge, 6,673-node prepared PKN.

The corresponding compressed branch merges 2,720 original nodes into 378 groups,
giving 36,569 edges and 4,331 scoring nodes before RNA coherence. It protects
upstream/downstream inputs and all input-regulon sources and targets, with
whole-group protection to prevent the installed compressor renaming protected
nodes. The saved original-to-compressed projection passes validation.

The [authors' compression guidance](https://github.com/saezlab/COSMOS_basic#network-compression)
describes the reason for compression: parents with identical children can
repeat identical downstream evidence in a subsequent upstream fit.
[ULM is a regression t-statistic](https://saezlab.github.io/decoupleR/reference/run_ulm.html),
so changing these repeated features and the regression background can change
scores. Larger uncompressed values cannot be ranked as stronger evidence.
This does not prove that every individual score increase is an artifact:
propagation, background composition and coherence filtering also differ.

An export named `pkn_final_uncompressed.rds`, or a GraphML containing original
gene names, does not imply uncompressed scoring. Compressed scores are deliberately
decompressed for original-ID extraction and presentation. Inspect the run state,
compression map and `pkn_final_scoring.rds` to establish what was scored.

The conservative compressor is a safe **one-pass** implementation, not complete
elimination of all redundancy. The initial compressed graph still contains
257 identical signed-child groups involving 857 parents; 243 groups contain
protected IDs. Child merging and subsequent edge deletion can also create new
equivalences. No new bug was found in the wrapper's declared contract.
Full implementation findings and counts are in
[`compression_method_audit.md`](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/compression_method_audit.md).

## Search design and compressed solutions

The accepted input mapping was unchanged: Factor 1 hepatocyte TF footprint
activities with absolute activity above 2 as upstream inputs; exact-mapped
Factor 1 lipid loadings as terminal downstream observations; hepatocyte RNA
weights for transcriptional coherence. The baseline prepared topology has
49 supplied TFs and 74 lipid observations. These are **supplied**, not all
sign-compatible returned TF scores. The Factor 1 axis remains testes-positive.

We tested **7,563 additional state/threshold pairs**:

- 2,376 pairs on the unchanged radius-6 compressed states at depths
  4, 5, 6, 7, 8, 10 and 12. The last four states are exactly identical, leaving
  four unique score states. There are 253 nonempty results below 500 edges,
  including 249 with at least two TFs and two lipids. All occur at depth 4.
- 5,187 pairs on eight newly scored compressed sensitivity states: radii
  4, 5, 7 and 12 at depth 8; initial RNA filtering at radius 6/depths 4 and 8;
  and stronger compression at radius 6/depths 4 and 8. There are 759 nonempty
  results below 500 edges, including 749 with at least two TFs and two lipids.

The coarse grid spans primary 0.5–3.0 and secondary 0–2.0 in steps of 0.1,
with secondary no greater than primary. Refinement uses 0.025 increments
within ±0.1 of threshold/coverage Pareto points and near-cap 500–999-edge
results. This includes promising near-cap cases missed by refining only
already eligible solutions. Counts describe sampled parameter pairs, not
distinct biological networks or an exhaustive real-valued search.

Extraction uses the established local `moon_extract_paths()` implementation
of dual-threshold selection with explicit sign and connectivity checks.
Inferred nodes must exceed the secondary cutoff, every retained edge touches
a primary node, measured level-0 observations are exempt from inferred-score
cutoffs, and every node lies on a directed intended-TF-to-observation walk.
This is the guarded helper used in the pilot, rather than an unmodified call
to the package reducer. Networks were not trimmed manually to meet the edge cap.

| Role and local candidate ID | Radius / depth | Primary / secondary | Nodes / edges | TFs | Measured lipids |
|---|---|---|---|---|---|
| Broader baseline, **threshold_search/C03** | 6 / 4 | **1.85 / 0.55** | **153 / 388** | Ahr, Hnf4a, Ppara, Tbx2 | Three DGs |
| Stricter baseline, threshold_search/C01 | 6 / 4 | 2.0 / 0.85 | 110 / 255 | Ppara, Tbx2 | Two DGs |
| Wider-topology sensitivity, sensitivity_threshold_search/C01 | 12 / 8 | **2.15 / 1.075** | **177 / 468** | Klf10, Tbx2 | Two DGs |
| Maximum baseline endpoint coverage, threshold_search/C02 | 6 / 4 | 1.875 / 0.325 | 166 / 446 | Ahr, Hnf4a, Nfe2l2, Ppara, Tbx15, Tbx2 | Three DGs |

The two common DGs are **DG(24:0_18:2)** and **DG(24:1_18:2)**. The third is
**DG(20:0_18:1)**. Candidate IDs are local to their respective search directories.

C03 is a practical representative inside a fixed-score threshold plateau:
primary 1.85 or 1.875 combined with secondary 0.55 or 0.575 produces exactly
the same 388 signed edges. Its TF scores are Ahr −0.75655, Hnf4a −0.60165,
Ppara −0.98041 and Tbx2 −2.37620. Raising secondary to 0.60 leaves Hnf4a only
0.00165 above the cutoff, whereas 0.55 gives a 0.05165 margin. C03 deliberately
does not maximize every numerical cutoff. These margins measure threshold
sensitivity for fixed scores, not score uncertainty or statistical confidence.
C02 gains Nfe2l2 at only +0.33190, just 0.00690 above its secondary cutoff.

The radius-12 candidate retains Klf10 at −1.30626 and Tbx2 at −2.40575;
Klf10 is 0.23126 above its secondary cutoff. It is a viable compressed
higher-cutoff alternative. However, the wider graph changes regression context
and includes additional cyclic topology; it should remain a sensitivity
comparison rather than win solely through higher numerical scores.

All ten exported baseline/sensitivity candidates pass independent saved-state
checks of final-PKN membership, exact score identity, compression projection,
cutoffs, edge signs and directed input-to-observation connectivity. Candidate
exports use original biological IDs, with measured weights and inferred scores
kept distinct. Four main examples also have exact prior-evidence joins.

## Sensitivity of the broader compressed candidate

At the same primary/secondary cutoffs, 1.85/0.55:

| State | Edges | TFs / lipids | Signed-edge Jaccard with baseline C03 |
|---|---:|---:|---:|
| Baseline, radius 6 / depth 4 | 388 | 4 / 3 | 1 |
| Protect measured inputs only, radius 6 / depth 4 | 388 | 4 / 3 | **1** |
| Initial RNA filter, radius 6 / depth 4 | 368 | 4 / 3 | 0.913924 |
| Smaller topology, radius 4 / depth 8 | 429 | 4 / 3 | 0.783843 |

Every row retains the same four TFs and three DGs. The stronger compression
variant allows unmeasured regulon nodes to merge, then decompresses scores for
RNA coherence against original-ID edges and recompresses after deletions.
Thus reduced protection does not silently discard the regulon's identity-based
checks. Every iteration and the final original-to-scored-PKN projection were
validated. Exact reproduction of C03 under this stronger compression is useful
methodological support, although it is not independent biological validation.

## Why Esr1, Nfic and Dbp disappear

All three are present in the raw PKN, pass the input-activity filter, survive
topological pruning, retain their original IDs during compression, and remain
in the final scoring PKN. The following values are identical across the seven
baseline compressed depth settings:

| TF | Supplied Factor 1 footprint activity | Lipid-driven MOON score before sign filtering | First scoring layer | Informative children in that layer |
|---|---:|---:|---:|---|
| Esr1 | −7.08875 | **+0.13180** | 3 | Lpl |
| Nfic | −4.34923 | **+0.71859** | 4 | Sp1, Nr3c1, Src |
| Dbp | −7.66088 | **+0.34516** | 4 | Esr1 |

The two score columns describe different evidence and their magnitudes are not
interchangeable. The relevant exclusion test is their **opposing signs**.
Diagnostic recomputation without the upstream sign filter exposes the omitted
rows and reproduces all returned scores within 1e-10. These diagnostic scores
were never used to bypass the sign constraint in solution extraction.

The installed `moon()` computes all layers and then applies upstream sign
agreement. Consequently, Esr1's positive raw score can contribute to Dbp before
both rows are removed from the returned table. A missing score row therefore
does not establish that a node was absent from the scoring graph.

The sparse first-layer footprints explain the weak support:

- **Esr1:** only Lpl is informative when Esr1 is first scored. Lpl scores
  −0.10897, which is above the preceding layer's mean, −0.23681 across 50 nodes.
  ULM includes unrelated layer nodes with zero regulator coefficients; this
  activating child therefore yields a weak positive statistic despite its
  negative score. It is not a signed uncentered target average.
- **Nfic:** Sp1 (+0.53433), Nr3c1 (+0.13180) and Src (−0.05520) are all above
  the 231-node background mean, −0.23527, producing a positive statistic.
- **Dbp:** its sole informative child is Esr1, again above that background.

Their initial scores already disagree with the supplied signs, before iterative
RNA coherence. None of their own outgoing edges is removed by the standard
post-score coherence loop, because the three upstream rows have already been
omitted. Other TFs' coherence filtering can still change their layer backgrounds.
All five informative edges have local OmniPath GRN evidence. Esr1→Lpl is absent
from the original custom receptor footprint used to infer Esr1 activity; only
3 of Esr1's 302 pruned outgoing signed edges overlap that input regulon. This
is a mismatch between evidence footprints worth retaining in interpretation,
not a reason to force a different sign or silently change input construction.

Increasing scoring depth does not rescore these nodes with later children:
the installed algorithm scores each node at its first eligible layer. Lower
extraction cutoffs cannot restore rows removed by the upstream-sign filter.

## Attempts to retain the three TF controls

The audit covers 42 TF/run comparisons in prior compressed/uncompressed states
and 24 in eight new compressed sensitivities. All new runs retain Esr1, Nfic
and Dbp through preparation and score them with the same informative-child
counts, one, three and one. Their final raw scores stay positive:

| TF | Final raw-score range across audited compressed settings |
|---|---:|
| Esr1 | +0.0897 to +0.1327 |
| Nfic | +0.6281 to +0.9142 |
| Dbp | +0.3341 to +0.3690 |

The tested changes include narrower and wider topology radii, depths through
12 in the baseline, stronger compression, and a uniform **initial RNA-coherence
filter** for all 49 supplied TF activities. The latter removes an edge only if
its exact signed triple has both transcription provenance and input-regulon
support and disagrees with TF activity × RNA target weight. It removes 1,750
edges from 48 TFs, once before ordinary scoring/coherence. The complete deletion
set was independently reconstructed. Matched depth-4/depth-8 preprocessing
runs avoid confounding compression policy with scoring depth. No change restores
the three controls to the returned score table.

Their physiological relevance is compatible with a failure of this particular
Factor 1 lipid-derived network to support their direction. The audit identifies
limited informative layer coverage and sign disagreement; it cannot distinguish
incomplete tissue/species knowledge from an inadequate model of this signature.
Retain these TFs as explicit coverage controls in the report, not as inferred
members of the solution. Do not flip their activities, disable sign filtering,
or selectively delete inconvenient children merely to include them.

## Interpretation limits

The broader graph's shortest TF-to-lipid paths contain 7–9 edges; the wider
radius-12 candidate has 7–10-edge paths. Scoring depth is not an extracted
path-length bound. The selected mechanisms converge on an Smpd1-associated
PC-to-DG reaction tail instantiated from lipid rules. Evidence for a generic
reaction does not establish activity on these exact measured acyl-chain species,
flux or hepatocyte-specific substrate availability. Compression and algorithmic
checks improve the analysis; they do not validate the species-level mechanism.

The search has no independent holdout, sample resampling or parameter-selection
correction. MOON cutoffs here are analysis choices, not calibrated p-values or
FDR thresholds. The demonstrated robustness is conditional on the saved Factor 1
inputs and current PKN. The narrow DG coverage remains a substantive limitation.

## Deliverables

All outputs are under `miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/`:

- [Broader compressed GraphML](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03/network.graphml),
  [SIF](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03/network_SIF.csv),
  [node attributes](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03/network_ATT.csv),
  and [illustrative routes PDF](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03/figures/candidate_paths_network_network.pdf).
- [Higher-cutoff compressed GraphML](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/sensitivity_threshold_search/C01/network.graphml)
  and [routes PDF](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/sensitivity_threshold_search/C01/figures/candidate_paths_network_network.pdf).
- [Complete baseline score-retention table](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/score_retention/radius6_depth4_compressed_scores_all.csv):
  2,238 original-ID rows before sign filtering, 2,207 returned, 31 upstream
  sign conflicts. Explicit flags distinguish returned scores, excluded scores
  and membership of each baseline solution.
- [Detailed TF audit](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/tf_audit/README.md),
  [baseline search](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/README.md),
  [sensitivity search](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/sensitivity_threshold_search/README.md)
  and [summary figure](../miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/figures/compressed_review_summary.pdf).

Search directories contain every tested cutoff, selected candidates, endpoint
margins, shortest routes and independent validation. Scoring directories retain
the exact input vectors, compression map, scoring PKN, original-ID final PKN,
coherence ledger, scores, diagnostic pre-sign-filter scores and run manifests.
Provenance-enriched exports distinguish curated edges from measured-species rule
hypotheses. The installed cosmosR version is 1.21.2, commit
`b2ed195738479c4ad204cace37fd11ae2fd02bb8`; environment records are retained.

## Historical execution commands

The following records the historical experiment sequence at the scripts’ archived locations. These scripts were executed before archival; saved-state validation and exported-network checks passed. Their original path assumptions remain unchanged, so use the active entry point above for current reproduction.

```sh
Rscript scripts/moon/archive/2026-09/audit_factor1_focus_tf_retention.R
Rscript scripts/moon/archive/2026-09/run_factor1_compressed_sensitivities.R
Rscript scripts/moon/archive/2026-09/run_factor1_compressed_sensitivities.R --radii=6 --depth=4 --prescore-rna
Rscript scripts/moon/archive/2026-09/run_factor1_compressed_sensitivities.R --radii=6 --depth=8 --prescore-rna
Rscript scripts/moon/archive/2026-09/run_factor1_inputs_only_compression.R --depth=4
Rscript scripts/moon/archive/2026-09/run_factor1_inputs_only_compression.R --depth=8
Rscript scripts/moon/archive/2026-09/validate_factor1_compressed_tf_sensitivities.R
Rscript scripts/moon/archive/2026-09/search_compressed_threshold_review.R
Rscript scripts/moon/archive/2026-09/search_compressed_threshold_review.R --sensitivity
Rscript scripts/moon/archive/2026-09/audit_compressed_candidate_neighborhood.R
Rscript scripts/moon/archive/2026-09/export_compressed_score_retention.R
Rscript scripts/moon/archive/2026-09/plot_compressed_review.R
python3 scripts/moon/archive/2026-09/export_factor1_network.py --results miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03 --stem network
python3 scripts/moon/archive/2026-09/export_factor1_network.py --results miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/sensitivity_threshold_search/C01 --stem network
```

These optional analyses do not alter workflows 00–10, the MOFA fit, accepted
input mapping or raw PKN, and are not consumed by existing analysis workflows.
