# Lipid-to-TF solution interpretation

Factor 1 lipid-to-TF solution. Level-zero score is an observed RNA-derived TF activity, not a lipid weight. moon_score, upstream_lipid_weight, downstream_tf_activity and RNA_input are separate quantities. Substrate participation is a qualitative metabolic projection and does not establish signaling or flux causality. Species-specific reaction rules are hypotheses; chemical-to-protein regulatory edges are explicitly identified. Evidence endpoints were projected using the explicit supplied alias table; original endpoints and mapping bases remain in exported provenance. All solution edges were verified against the exact supplied prior: results/causal_modelling/cosmos_mouse_pkn/pkn_lipid_curated_ph73.csv.

No exact-mapped measured lipid has an eligible direct regulatory edge in the full PKN. Any selected routes with a later chemical-to-protein transition nominate indirect metabolic availability hypotheses, not direct signaling lipids or confirmed biochemical precursors.

One shortest path per reachable input pair; greatest minimum absolute inferred MOON score; deterministic signed-edge tie breaking in the shortest-path DAG. Level-zero observed TF activities excluded from score bottleneck. Not all alternative routes are enumerated.

- A directed signed route is a prior-knowledge hypothesis, not evidence of causality, flux or biochemical feasibility in hepatocytes.
- Reaction instances encode substrate participation, product formation and catalyst availability; they do not implement mass balance or enzyme kinetics.
- The actual regulatory chemical may be an unmeasured metabolite or co-product; its abundance or availability is not independently established.
- Shortest-route classifications do not summarize every alternative walk through the solution.
- Lipid weights and RNA-derived TF activities are observational measurements on different scales; MOON scores are inferred statistics.
