# Expanded Factor 1 TF→lipid MOON sensitivity

Start with [the full report](../../../docs/cosmos_factor1_expanded_moon_rerun.md) and `selected_network.json`. The preferred conservative output is `currency_v2_expanded_layers4_lower_thresholds/`: 138 edges, 94 nodes, eight TFs and 13 measured features at 0.75/0.75, pruning 6 and scoring 4. `currency_v2_expanded_layers6_lower_thresholds/` has 177 edges and the same 13 features, adding Hoxa10. All displayed edges follow consecutive score layers, all catalysts are retained, and no currency bridges remain. Scores are modest and not calibrated significance.

The `expanded_network_audit/` within either selected extraction contains labeled networks, GraphML, exact edge evidence, complete reaction masters, participant/currency audits and route tables. Illustrative PNG/PDF routes are in `figures/`.

`expanded_inputs/` and `expanded_coverage_1.85_1.0/` retain initial unrestricted 101/206-edge candidates as diagnostics. Their Ar/Arnt→Cer displays lack recursive scoring witnesses; the 101-edge candidate depends entirely on a water bridge. Do not treat these as the preferred result. `fixed_thresholds/` contains each initial broad topological comparator; it is separate from the consecutive-layer selected extraction.

`frozen_inputs_control/`, currency-v1 diagnostics, and `expanded_layer_consistent/` preserve controlled comparisons and threshold searches. V1 and V2 give identical prepared graphs/scores in these data. `inputs/` preserves the 155 frozen mappings and 268 new one-to-one assignments. Previous manuscript artifacts remain separate.

`independent_validation.json` passed seven fits and seven extractions. Complete commands, settings, mapping distinctions and limitations are in the report; individual fit/grid/export directories contain source snapshots and hash manifests.
