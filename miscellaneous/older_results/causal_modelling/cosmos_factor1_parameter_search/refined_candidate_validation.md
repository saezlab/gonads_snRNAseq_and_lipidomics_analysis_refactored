# Independent validation of the refined radius-6 candidates

Audit date: 2026-09-09. All six saved candidate networks passed the independent
checks below. The audit uses saved run-state RDS files and SIF/ATT exports only;
it does not rerun scoring or change the search. Reproduce with
`Rscript scripts/moon/audit_refined_candidates.R`.

| Scoring depth / representation | Primary / secondary | Nodes | Edges | TFs | Lipids | All checks |
|---|---|---:|---:|---:|---:|---|
| 4 / compressed | 2.00 / 0.75 | 114 | 261 | 3 | 2 | Pass |
| 4 / compressed | 2.00 / 0.85 | 110 | 255 | 2 | 2 | Pass |
| 4 / uncompressed | 1.70 / 1.50 | 156 | 461 | 4 | 3 | Pass |
| 4 / uncompressed | 1.90 / 1.60 | 83 | 208 | 2 | 2 | Pass |
| 5 / uncompressed | 2.00 / 1.50 | 171 | 499 | 3 | 2 | Pass |
| 5 / uncompressed | 2.05 / 1.60 | 123 | 318 | 2 | 2 | Pass |

Every graph has strictly `0 < edges < 500`; 499 is eligible. Exported nodes
exactly equal the edge endpoints, IDs and edges are unique, and every signed
edge occurs in its run's final PKN. Scores and levels match the exact saved
state. Every inferred node exceeds its secondary threshold, every edge has a
primary inferred endpoint, observed weights match the input vector, and
observations have no outgoing edges. All retained TF signs agree with the
intended activity inputs. Every node is reachable from an intended TF and can
reach a measured lipid. The declared extraction was independently reconstructed
without calling `moon_extract_paths()` and matched each saved edge set exactly.
Saved grid counts and final zero-removal coherence iterations also agree.

These are graph/implementation checks. They do not validate the biological
mechanisms or provide independent evidence for the reaction hypotheses.

## Depth 4 versus depth 5

The two refined uncompressed representatives use the same radius-6 topology,
49 TF inputs and 74 measured lipid inputs:

- **Depth 4, 1.90 / 1.60:** Nfe2l2 and Nr4a2; 208 edges.
- **Depth 5, 2.05 / 1.60:** Nfe2l2 and Klf10; 318 edges.

Both reach exactly the same two measured lipids, `DG(24:0_18:2)` and
`DG(24:1_18:2)`. Their selected networks share 188 edges, with an edge Jaccard of
0.5562. Comparing different primary cutoffs confounds the depth effect with
threshold selection. At the matched cutoff 1.90 / 1.60, the depth-4 network's
208 edges are entirely contained in the depth-5 network's 363 edges; depth 5
adds Klf10 and Nr2f1 but no measured lipid. At matched 2.05 / 1.60, depth 4 has
185 edges and one TF, Nfe2l2; all 185 edges occur in the 318-edge depth-5 result.

**All 2,235 depth-4 node scores are present and exactly identical at depth 5.**
Depth 5 adds 3,454 newly scored level-5 nodes. Klf10 was already scored at
level 4; the extra level supplies additional surviving routes rather than a
stronger Klf10 score. Its absolute score is 1.600186347, only **0.000186347 above
the secondary threshold of 1.60**. Its retention as a second TF is therefore
particularly threshold-sensitive. Nfe2l2 exceeds this cutoff by 0.028999776;
Nr4a2 in the depth-4 representative exceeds it by 0.394785614.

Selecting depth 5 is defensible under an explicitly declared preference for the
higher primary threshold while retaining two TFs and two lipids. It should not
be described as stronger biological support or improved lipid coverage. Depth 4
is the more compact representative and has less marginal retention of its
second TF. Reporting the pair makes that tradeoff visible. Neither result
establishes an optimal mechanistic model.

## Scoring depth is not retained route length

The refined depth-4 graph's shortest TF-to-lipid routes contain seven edges for
Nfe2l2 and ten for Nr4a2. In the refined depth-5 graph they contain nine edges
for Nfe2l2 and eleven for Klf10. Every TF can reach both displayed lipids.

This is consistent with the extraction contract: `n_layers` limits score
propagation, while the graph extractor retains signed, reachable edges without
requiring score levels to decrease along a displayed route. Do not call these
four-step or five-step mechanisms. For example, a shortest depth-5 route is:

`Nfe2l2 → Tle5 → Ctnnb1 → Foxo1 → Igf1r → Pdpk1 → Prkcd → Smpd1 → reaction → DG(24:1_18:2)`

The Klf10 route reaches the shared Ctnnb1 branch through Cdkn1a, Cdk1 and
Csnk2b. Such multi-step routes remain PKN hypotheses, and graph reachability
alone does not test their tissue specificity or biochemical feasibility.

## Output tables

- `refined_candidate_validation.csv`: all per-candidate checks, counts,
  threshold margins and shortest-route summaries.
- `refined_candidate_routes.csv`: every retained TF/lipid pair, its shortest
  directed distance, TF score, input activity and secondary-threshold margin.
- `refined_candidate_endpoint_edges.csv`: observed lipid endpoint edges and
  their inferred enzyme scores.
- `depth4_depth5_refined_comparison.csv`: selected and matched-threshold
  comparisons, including edge overlap and TF/lipid identities.
- `depth4_depth5_score_stability.csv`: exact score equality and newly scored
  node counts across these two depths.
