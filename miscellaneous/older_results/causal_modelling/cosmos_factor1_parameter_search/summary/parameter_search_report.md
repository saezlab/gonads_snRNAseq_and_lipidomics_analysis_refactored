# Factor 1 MOON parameter-search comparison

The grid contains **3,792 unique configurations/cutoff rows** (3,792 input rows). **152 rows** meet the declared validity, convergence, endpoint and ≤500-edge criteria; **132 rows** additionally retain at least two TFs and two lipid endpoints. These are candidate networks for comparison, not a single selected biological model.

## Declared ranking and validity

Threshold priority is lexicographic: higher secondary cutoff first, then higher primary cutoff, then more lipid endpoints, more upstream TFs and fewer edges. The interpretable ranking applies the same rule after requiring ≥2 TFs and ≥2 lipids. A complete, converged, valid, nonempty network with an upstream TF and downstream lipid is required. Missing or false validity certification fails eligibility. Materialization is reported separately because the scoring workflow may write candidate solutions after grid selection.

Depth, compression and topology radius change the scoring problem. Priority ranks across those configurations reflect the requested threshold preference; they do not make the scores or mechanisms biologically equivalent. Pareto frontiers and first-ranked shortlists are therefore computed separately within each configuration.

## Configuration-specific shortlist

| Candidate | Compression | Depth / radius | Secondary / primary | TFs / lipids | Edges |
|---|---|---:|---:|---:|---:|
| radius6_depth5_uncompressed__ecb3bfb8d1 | uncompressed | 5 / 6 | 1.65 / 1.95 | 1 / 2 | 282 |
| radius7_depth5_uncompressed__71d50aa9b6 | uncompressed | 5 / 7 | 1.65 / 1.95 | 1 / 2 | 282 |
| radius6_depth5_uncompressed__ea052cc19b | uncompressed | 5 / 6 | 1.6 / 2.05 | 2 / 2 | 318 |
| radius7_depth5_uncompressed__a27233dc7c | uncompressed | 5 / 7 | 1.6 / 2.05 | 2 / 2 | 318 |
| radius6_depth4_uncompressed__0ab7db0691 | uncompressed | 4 / 6 | 1.6 / 1.95 | 1 / 2 | 190 |
| radius7_depth4_uncompressed__e063f2f2ec | uncompressed | 4 / 7 | 1.6 / 1.95 | 1 / 2 | 190 |
| radius6_depth4_uncompressed__00147faa29 | uncompressed | 4 / 6 | 1.6 / 1.9 | 2 / 2 | 208 |
| radius7_depth4_uncompressed__f587487279 | uncompressed | 4 / 7 | 1.6 / 1.9 | 2 / 2 | 208 |
| radius6_depth4_compressed__1687fd4b77 | compressed | 4 / 6 | 0.85 / 2 | 2 / 2 | 255 |

## Pareto comparison and sensitivity

There are **48 nondominated rows within their own configurations**. The objectives retain cutoff stringency, lipid/TF coverage and network size explicitly, without collapsing them into an arbitrary composite score. `eligibility_stability_counts.csv` counts how many tested depth/compression/radius configurations pass at each matched cutoff pair. Its denominators can differ when coarse and refined searches cover different cells. Those counts measure reproducibility of eligibility, not agreement in TF identities, lipid identities, edges, signs or mechanisms.

## Candidate delivery and provenance

`candidate_export_plan.csv` lists materialization status and the exact exporter command for each configuration-specific shortlist row. Commands are plans only; this summarizer never fits models or writes into solution directories. Once the corresponding SIF/ATT files exist, the exporter independently checks exact PKN membership and evidence, score-sign consistency, upstream controllability and downstream observability. It writes Cytoscape-ready GraphML with distinct Factor 1 weight, inferred MOON score and TF-activity fields, route tables and figures. Do not call an unmaterialized or unaudited candidate a validated network.

Measured-species reaction rules remain substrate-specificity hypotheses. Positive reaction traversal is not a flux or stoichiometric production prediction. Scoring depth does not cap every TF-to-lipid route in a merged graph; compare bounded catalyst-only routes as well as unrestricted connectivity.

## Files

- `grid_assessment.csv`: every deduplicated row and explicit exclusion reasons.
- `ranked_eligible.csv` and `ranked_interpretable.csv`: threshold-priority rankings, with configuration-specific ranks.
- `pareto_within_configuration.csv` and `configuration_shortlist.csv`: within-configuration tradeoffs and shortlisted candidates.
- `eligibility_stability_counts.csv`: matched-cutoff eligibility counts and tested denominators.
- `candidate_export_plan.csv`: candidate artifact locations and auditable export commands.
- `summary_manifest.json`: exact source/script hashes, criteria and row counts.
