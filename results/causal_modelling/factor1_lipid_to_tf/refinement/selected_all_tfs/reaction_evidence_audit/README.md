# Complete reaction and regulatory evidence for the selected solution

Generated without fitting or changing the selected network. `selected_reaction_evidence.csv` includes every selected enzyme-reaction node, grouped by accession, complete frozen Rhea equations and all annotated mouse catalyst alternatives. JSON fields preserve the original participant stoichiometries, source annotations and selected incident edges. Graph presence does not establish measured availability or full multi-substrate reaction feasibility.

`selected_regulatory_evidence.csv` recovers the original frozen GtoPdb records and their binding affinities. The record enumeration is not an official persistent interaction accession. Ki is not EC50 or agonist efficacy. The interpretive primary-source audit is in `docs/cosmos_factor1_lipid_tf_selected_evidence.md`.

Reproduce with `python3 scripts/moon/ancillary/audit_lipid_tf_selected_reactions.py --results=results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs --output=NEW_EMPTY_DIRECTORY`. Occupied output directories are refused.
