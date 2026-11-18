# Refined Factor 1 candidate comparison

All three candidates pass the independent exact-PKN provenance, edge-sign and upstream/downstream connectivity audits. They retain the same two observed DG endpoints but different upstream TF endpoints. This is stability of a downstream reaction hypothesis, not stable identification of its upstream TF control.

| Candidate | Scoring | Primary / secondary | Nodes / edges | Upstream TFs | Shortest routes |
|---|---|---:|---:|---|---:|
| U4 | uncompressed, depth 4 | 1.90 / 1.60 | 83 / 208 | Nfe2l2, Nr4a2 | 7-10 edges |
| U5 | uncompressed, depth 5 | 2.05 / 1.60 | 123 / 318 | Klf10, Nfe2l2 | 9-11 edges |
| C4 | compressed, depth 4 | 2.00 / 0.85 | 110 / 255 | Ppara, Tbx2 | 9-9 edges |

## Shared terminal hypothesis

The retained observations are DG(24:0_18:2) and DG(24:1_18:2). Both terminate in Smpd1-associated, rule-instantiated PC-to-DG reactions supported by the generic RHEA:10605 template. The exact PC/DG acyl-chain specificity remains assumed. Each retained enzyme-instance score has only one directly linked measured observation in a background of 74 measured nodes; the two separate enzyme instances do not provide independent validation of the reaction rules.

All three candidate graphs contain the Pdpk1 → Prkcd → Smpd1 → reaction-instance → DG tail. The complete selected shortest TF-to-DG routes are 7–11 edges long. None has a route of at most six edges, even though the score propagation depths are four or five; depth and route length are distinct quantities.

## What changes between candidates

- U4 versus U5: 78 shared nodes (Jaccard 0.609), 188 shared signed edges (Jaccard 0.556); TF endpoint Jaccard 0.333.
- U4 versus C4: 40 shared nodes (Jaccard 0.261), 70 shared signed edges (Jaccard 0.178); TF endpoint Jaccard 0.000.
- U5 versus C4: 42 shared nodes (Jaccard 0.220), 71 shared signed edges (Jaccard 0.141); TF endpoint Jaccard 0.000.

The three-way intersection contains 40 nodes and 70 signed edges, but no common upstream TF endpoint. Score-distribution and shared-node score comparisons are exported separately; they do not imply numerical equivalence between differently compressed or propagated runs.

All 2,161 inferred nodes shared by the saved U4 and U5 score tables have identical scores and levels. U5 adds scored nodes and uses a different primary cutoff; the TF endpoint change therefore cannot be described as changed scores on the shared nodes or attributed to depth alone.

U5 meets the higher primary-cutoff preference, whereas U4 is smaller and retains the same two lipid endpoints. Klf10 in U5 lies only 0.000186 above its secondary cutoff in absolute score. This threshold margin is a numerical sensitivity warning, not a confidence interval. U4's Nfe2l2 score also lies close to the 1.60 secondary cutoff. The compressed candidate changes the upstream TF endpoints completely and requires a much lower secondary cutoff.

## Recommended interpretation

Keep U4 and U5 as distinct uncompressed candidate hypotheses and C4 as a compression-sensitivity comparison. U4 is a useful compact network for inspection; U5 satisfies the stronger primary-cutoff preference. Neither is validated as the upstream mechanism. The common Smpd1-associated PC-to-DG tail is the reproducible part of these three extracted graphs, with measured-species reaction specificity still hypothetical. Treat the TF assignments and long regulatory routes as candidates requiring stronger context and perturbation evidence.

## Auditable outputs

- `candidate_summary.csv`, `endpoint_identities.csv`, `pairwise_overlap.csv`: exact identities, sizes, path bounds and signed overlaps.
- `candidate_paths.csv`, `reaction_observation_support.csv`: auditable paths and immediate measured support.
- `score_distributions.csv`, `candidate_node_scores.csv`, `paired_candidate_scores.csv`, `node_membership_and_scores.csv`: score scales, margins, memberships and comparisons.
- `shared_signed_edges_all_three.csv`: three-way overlap; it is an intersection for comparison, not a newly fitted network.
- `source_artifact_hashes.csv`: SHA-256 hashes of the saved source artifacts.

GraphML and complete edge/reaction provenance are alongside each original SIF/ATT, using its exact solution stem. Original scores, grids and networks were not changed.
