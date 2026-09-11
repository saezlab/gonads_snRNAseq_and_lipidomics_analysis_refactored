# Compressed-score threshold review

This is an extraction-only review of the **unchanged, converged compressed** scoring states from the previous Factor 1 search. The original input mapping, MOON scores, RNA coherence filtering, and scoring PKNs are immutable inputs. Lowering extraction cutoffs cannot restore a TF that is absent from those score tables.

## Search design and validation

`scripts/moon/search_compressed_threshold_review.R` tests primary cutoffs 0.5–3.0 and secondary cutoffs 0–2.0 in steps of 0.1, with secondary ≤ primary. It refines by 0.025 within ±0.1 of the coarse threshold/coverage Pareto frontier and near-cap 500–999-edge results. Including near-cap infeasible cases avoids the earlier design's restriction to already-eligible runs. All 2,376 tested threshold pairs remain in `grid_complete.csv`, including empty results.

Depths 7, 8, 10, and 12 have exactly identical score tables, final PKNs, and input vectors; `score_state_aliases.csv` records this equivalence. Thus four unique scoring states cover the seven requested depths. The counts below describe sampled settings; the grid is not an exhaustive enumeration of all possible real-valued cutoffs.

| Compressed scoring depth | Tested cutoffs | Empty results | Smallest nonempty result | Results with 1–499 edges |
|---|---:|---:|---:|---:|
| 4 | 882 | 388 | 235 | 253 |
| 5 | 498 | 312 | 649 | 0 |
| 6 | 498 | 312 | 747 | 0 |
| 7 (also 8, 10, 12) | 498 | 312 | 753 | 0 |

Of the 253 eligible results, 249 retain at least two TFs and two measured lipid nodes. `coverage_threshold_pareto.csv` retains tradeoffs in both thresholds and TF/lipid coverage. `threshold_frontier.csv` considers thresholds alone and therefore can favor smaller endpoint sets.

All extractions use the same `moon_extract_paths()` dual-threshold function. Measured level-0 lipid weights are exempt from inferred-node cutoffs and cannot create primary nodes. No edge-budget truncation or manual TF rescue is used. Repeated cleanup of an identical immutable PKN is memoized for grid speed; each exported candidate is independently re-extracted without that cache.

Validation confirms that each final original-ID PKN projects exactly to the PKN scored in its compressed state, virtual-node groups retain identical scores after decompression, every solution edge belongs to the matched final PKN, every inferred node is strictly above the secondary cutoff, every edge touches a primary node, every edge is sign-consistent, and every node is on an intended TF-to-observation directed walk. `selected_validation.csv` records the checks. This does not imply that every walk is short or unique.

## Candidate choices

| Candidate | Primary / secondary | Nodes / edges | TFs / measured lipids | Interpretation |
|---|---|---|---|---|
| C01 | 2.0 / 0.85 | 110 / 255 | 2 / 2 | Highest secondary cutoff on this grid retaining multiple TFs and lipids |
| C03 | 1.85 / 0.55 | 153 / 388 | 4 / 3 | Practical balance between coverage, edge count, and cutoff margins |
| C02 | 1.875 / 0.325 | 166 / 446 | 6 / 3 | Maximum endpoint coverage; includes weak inferred TF scores |

**C03** retains Ahr, Hnf4a, Ppara, and Tbx2 and the measured species DG(20:0_18:1), DG(24:0_18:2), and DG(24:1_18:2). Their MOON scores are −0.75655, −0.60165, −0.98041, and −2.37620 respectively. Those are distinct from their much larger RNA-footprint input activities. A secondary cutoff of 0.55 gives Hnf4a a 0.05165 margin; the threshold-maximizing 0.60 choice would leave only a 0.00165 margin. The 1.85/0.55 candidate is deliberately not the mathematical threshold frontier.

C01 retains Ppara and Tbx2 with the two 24-carbon DG species. C02 additionally retains Nfe2l2 and Tbx15, but Nfe2l2's MOON score is only 0.33190, just 0.00690 above its cutoff. Greater endpoint coverage therefore should not be interpreted as stronger mechanistic support.

The independent `C03_local_threshold_neighborhood.csv` audit makes this selection explicit:

| Primary / secondary | Edges | Hnf4a secondary margin | Minimum inferred-node secondary margin | Minimum primary-node margin | Same edge set as C03 |
|---|---:|---:|---:|---:|---|
| 1.85 / 0.55 | 388 | 0.051645 | 0.043884 | 0.033190 | Yes |
| 1.85 / 0.575 | 388 | 0.026645 | 0.018884 | 0.033190 | Yes |
| 1.875 / 0.55 | 388 | 0.051645 | 0.043884 | 0.008190 | Yes |
| 1.875 / 0.575 | 388 | 0.026645 | 0.018884 | 0.008190 | Yes |
| 1.85 / 0.60 | 386 | 0.001645 | 0.001645 | 0.033190 | No; two edges removed |

C03 is a representative relaxed point on an unchanged-network plateau, farther from the next node-score thresholds than nearby more aggressive cutoffs. It is neither a unique optimum nor a statistically validated robustness result. These margins describe threshold sensitivity for fixed scores, not uncertainty in the scores themselves.

The shortest retained TF-to-lipid routes have 7–9 edges in C03 and 9 edges in C01. `n_layers = 4` controls scoring propagation, not a hard path-length cap in the extracted network. The selected routes converge on Smpd1 and measured-species lipid reaction hypotheses. They remain indirect, prior-dependent hypotheses; these cutoffs are analysis choices rather than universal significance levels. Esr1, Nfic, and Dbp are absent from every unchanged compressed score state and are not restored by any threshold setting here.

Each Cxx directory includes an original-ID SIF/ATT, GraphML, exact RDS extraction, and audit JSON. `selected_endpoints_and_margins.csv` and `selected_shortest_routes.csv` provide endpoint evidence and representative paths. Provenance-enriched exports may be added separately by the main review.

## Reproduce

Run from the repository root:

```sh
Rscript scripts/moon/search_compressed_threshold_review.R
```

The script reuses per-state extraction-grid caches. The independent topology and preprocessing sensitivities use `--sensitivity` and write to a separate `sensitivity_threshold_search` directory so that alternative scoring backgrounds do not enter this baseline ranking.
