# Factor 1 MOON pilot outputs

Main report and reproduction: `../../docs/cosmos_factor1_moon_pilot.md`.
Input provenance and schemas: `inputs/README.md`.

The compressed tutorial-style run is in this directory. The independently
converged uncompressed comparison is in `uncompressed_sensitivity/`.
The primary compressed network is empty; it was not replaced by a relaxed
network under the same filename. The uncompressed primary network has 473
nodes, 1,676 edges, eight TFs and four observed lipid endpoints.

Start with `uncompressed_sensitivity/solution_primary.graphml` in Cytoscape,
and the annotated figures and path tables in the same subdirectory. The compact
`solution_catalyst_bounded6.graphml` is instead a selected route union from the
compressed **relaxed** graph. Neither is a validated causal mechanism.

`score_support_README.md` and `fullcomparison.json` explain why compression,
sparse reaction support and threshold choices materially affect the networks.
All measurement weights are kept separate from inferred ULM score magnitudes.
