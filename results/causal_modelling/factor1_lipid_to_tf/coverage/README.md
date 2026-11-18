# Factor 1 lipid-to-TF coverage audit

This is a read-only audit of frozen resource identities and topology. No MOON fit is run. The original resources and accepted measurement mapping are unchanged.

`topology_summaries.csv` counts endpoint reachability by shortest directed path with the indicated TF observations made terminal. These distances differ from iterative forward/backward MOON pruning; they are not activity scores or significance results. `endpoint_reachability.csv` preserves every pair.

`pkn_lipid_precursor.csv` and `pkn_ph73_lipid_precursor.csv` contain the filtered original and existing acid/base-normalized curated graph, respectively. The latter is a sensitivity, includes no measured-species reaction additions, and uses `normalized_edge_provenance.csv` to recover original edge IDs. Both omit nonlipid metabolites and protein-to-reaction catalysis so traversal cannot loop from signaling back into metabolism. Frozen classification admits explicit lipid descendants and their pH7.3 forms.

`exact_regulatory_edges.csv` contains accepted measured nodes with direct curated regulatory edges; an empty table is a valid result. `curated_regulatory_edges.csv` gives the complete current regulatory layer. `any_mapping_opportunities.csv` reports exact, compatible, class and alias candidates separately, including cached evidence tiers that are not used in the exported graphs.

The normalized LPC16-to-DHA route is acyl-remodeling participation: LPC plus unmeasured DHA-CoA makes PC, and phospholipase hydrolysis releases DHA while regenerating LPC. It does not prove net DHA production or a direct signaling effect of measured LPC. `example_paths.csv` also records precursor dead ends and unsafe cofactor shortcuts for inspection. Full reaction stoichiometry and catalyst provenance remain in the frozen PKN `rhea_reactions.jsonl` and `nodes.jsonl`.

`dependency_hashes.csv` includes the script and every consumed dependency.
