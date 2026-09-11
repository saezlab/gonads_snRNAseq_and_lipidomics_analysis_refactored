# Separate topology and preprocessing sensitivity extractions

This directory applies the same audited dual-threshold search to eight newly scored, converged **compressed** states supplied by the main review. It is kept separate from `../threshold_search`, which uses the unchanged previous compressed scoring states. Every state is extracted against its own final original-ID PKN after verifying its projection to its own scored compressed PKN. Original measurements, state scores, and state PKNs are not modified by this extraction script.

The coarse grid uses primary 0.5–3.0 and secondary 0–2.0 by 0.1, with secondary ≤ primary. Local 0.025 refinement includes eligible threshold/coverage Pareto points and near-cap 500–999-edge points. The complete 5,187-row grid contains 759 nonempty results with fewer than 500 edges, including 749 with at least two TFs and two lipid observations. These are sampled cutoffs, not an exhaustive search over real numbers.

| Scoring state | Tested cutoffs | Empty | Smallest nonempty edge count | Nonempty <500 |
|---|---:|---:|---:|---:|
| Radius12, depth8 | 498 | 271 | 468 | 8 |
| Radius4, depth8 | 789 | 383 | 356 | 170 |
| Radius5, depth8 | 498 | 324 | 590 | 0 |
| Radius7, depth8 | 498 | 312 | 755 | 0 |
| Radius6, depth4, inputs-only compression protection | 882 | 375 | 232 | 280 |
| Radius6, depth8, inputs-only compression protection | 498 | 312 | 647 | 0 |
| Radius6, depth4, pre-score RNA filtering | 1,026 | 395 | 214 | 301 |
| Radius6, depth8, pre-score RNA filtering | 498 | 312 | 657 | 0 |

Scoring backgrounds differ across these states. Higher ULM cutoffs across different topologies or compression policies do not independently establish greater biological validity. The depth4/depth8 pairs separate preprocessing effects from propagation depth; judging only the depth8 preprocessing variants against the preferred baseline depth4 would confound these two choices.

## Exported examples

All candidate IDs below are local to this sensitivity directory and distinct from identically named baseline directories.

| ID | State | Primary / secondary | Nodes / edges | TFs / lipids |
|---|---|---|---|---|
| C01 | Radius12, depth8 | 2.15 / 1.075 | 177 / 468 | 2 / 2 |
| C02 | Pre-score RNA, depth4 | 2.0 / 0.85 | 102 / 237 | 2 / 2 |
| C03 | Pre-score RNA, depth4 | 1.875 / 0.325 | 159 / 419 | 7 / 3 |
| C04 | Inputs-only protection, depth4 | 1.975 / 0.85 | 112 / 262 | 2 / 2 |
| C05 | Inputs-only protection, depth4 | 1.875 / 0.30 | 162 / 431 | 6 / 3 |
| C06 | Radius4, depth8 | 1.925 / 0.775 | 143 / 356 | 2 / 2 |
| C07 | Radius4, depth8 | 1.85 / 0.425 | 177 / 457 | 5 / 3 |

C01 is a compressed alternative with higher cutoffs, retaining Klf10 (MOON score −1.30626) and Tbx2 (−2.40575) and DG(24:0_18:2) plus DG(24:1_18:2). Klf10 is 0.23126 above the secondary cutoff. Its shortest TF-to-lipid paths are 7–10 edges, and the enlarged radius primarily changes available network context. This result should be compared with the baseline and smaller-radius sensitivities rather than described as uniquely preferred because it has the highest numerical cutoffs.

The depth4 strict pre-score and inputs-only examples retain Ppara and Tbx2 with those same two DGs. Relaxed coverage examples add Ahr, Hnf4a, Nfe2l2 and Tbx15; pre-score RNA C03 also adds Smad4. All three-lipid examples add DG(20:0_18:1). Small inferred TF scores at the relaxed cutoffs remain weak model support even when their footprint activities are large. These networks do not restore Esr1, Nfic, or Dbp, which are absent from the supplied sign-filtered score tables.

An additional matched-cutoff audit (`matched_threshold_sensitivity.csv`) compares the same 1.85/0.55 setting across all states. The inputs-only protection depth4 result has exactly the same 388 signed edges as baseline C03 (edge Jaccard = 1), while the pre-score-RNA depth4 result has 368 edges and the same four-TF/three-DG endpoint set (edge Jaccard = 0.913924). Radius4/depth8 gives 429 edges with that same endpoint set. Thus the balanced baseline endpoint hypothesis persists under these alternative preprocessing choices; the increased global compression in the inputs-only variant does not improve that particular extracted graph.

`selected_validation.csv` passes every independent export check for all seven candidates. `selected_endpoints_and_margins.csv` preserves exact inferred versus input activities and secondary margins; `selected_shortest_routes.csv` gives representative paths. Each candidate directory includes the original-ID SIF/ATT, GraphML, exact solution RDS, and extraction audit. Main-review provenance enrichment may add further exports without changing the scoring state.

Reproduce from the repository root:

```sh
Rscript scripts/moon/search_compressed_threshold_review.R --sensitivity
```

The script discovers supplied run-state files under the main review's `topology_sensitivity`, `prescore_sensitivity`, and `inputs_only_sensitivity` directories and reuses per-state extraction caches. Source-state and script hashes are recorded in `search_manifest.json`.
