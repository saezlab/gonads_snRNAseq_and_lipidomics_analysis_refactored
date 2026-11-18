# Relaxed lipid-to-PKN mapping audit

This audit quantifies the noise introduced if a measured lipid score were copied to multiple existing PKN metabolite nodes. It reads the frozen `pkn_lipid_rules` graph; it does not change the PKN or run MOON.

## Main files

- `per_feature_policy_summary.csv`: one row per measured feature and explicit mapping policy; use this for candidate metabolite nodes, product-side candidate metabolite nodes, reaction instances, product-side reaction instances, and substrate-side reaction instances.
- `candidate_chemical_identity_inventory.csv`: every source candidate identity, including identities that have no active PKN node.
- `candidate_pkn_node_inventory.csv`: each supported feature-to-metabolite-node assignment.
- `candidate_reaction_instance_inventory.csv`: each feature-to-candidate-reaction link and whether the lipid node is a reaction product or substrate.
- `pkn_node_reuse_by_policy.csv`: identifies generic nodes to which many measured features would collapse.
- `population_duplication_summary.csv` and `mapping_tier_summary.csv`: panel-level burden.

## Mapping tiers

The positional/structural tier (`compatible_more_specific`) is kept apart from aggregate composition, generic diacyl template, and hierarchy mappings. It can include more than sn position, such as double-bond or stereochemical detail. A generic Rhea class term is reported only in `generic_diacyl_template_stress_test`; it is not treated as a species-specific or sn-relaxed assignment.

`native_position_or_structure_relaxed` reports only direct ChEBI-addressable candidates and deliberately excludes the existing measured-species rule nodes. The `underscore_molecular_species_features` population is the narrowest relevant subset for the question of unreported sn position: it includes only molecular-species measurements with unordered acyl chains.

The audit is exhaustive over the candidate chemical identities addressable in the frozen SwissLipids-to-ChEBI mapping and the active `pkn_lipid_rules` nodes. It is not an exhaustive enumeration of all chemically possible positional isomers: detailed SwissLipids records without a ChEBI ID cannot be mapped to this PKN without adding new, uncurated nodes.

## MOON scoring implication

Candidate mappings must not be passed to MOON as independent copies of the same lipid score: the ULM score treats duplicated rows as separate observations and would over-weight an ambiguous feature. A relaxed sensitivity should instead retain the feature-to-candidate mapping table and evaluate one supported node per feature in repeated mapping draws, or use an explicitly weighted aggregation whose interpretation is validated separately.
