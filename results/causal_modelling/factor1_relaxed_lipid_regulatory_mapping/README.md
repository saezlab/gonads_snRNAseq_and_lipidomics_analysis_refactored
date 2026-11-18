# Relaxed Factor 1 lipid-to-regulatory-PKN mapping
This pre-scoring ledger assigns each Factor 1 lipidomic measurement to all
supported candidate **chemical source nodes** in the directional exploratory
PKN. It does not change the PKN, choose a mapping draw, aggregate measurements,
or run MOON.
`default_relaxed_measurement_assignments.csv` contains only exact identity,
direct acid/base, same-formula structural-descendant, and compatible-more-
specific mappings. Every listed source has at least one signed direct
chemical-to-protein edge surviving the final PKN conflict filter.
`all_candidate_regulatory_assignments.csv` also contains broad class/ancestor,
lower-resolution, and headgroup-isomer annotations. Those rows have
`assignment_role=annotation_only` and must not be added to an upstream MOON
vector unless a separate policy is explicitly selected.
`mapping_paths.csv` preserves every SwissLipids/ChEBI route collapsed into an
assignment. `candidate_assignment_regulatory_edges.csv` provides direct
source-to-protein evidence and resource provenance. `measurement_mapping_summary.csv`
covers all 1,254 Factor 1 measurements. `candidate_node_collisions.csv` lists
the rare cases where different measurements reach the same PKN source node.
`default_relaxed_unique_topology_seeds.csv` has one row per default source node;
its value is blank for a cross-feature collision rather than inventing an
aggregation.
`candidate_measurement_value` is the unchanged oriented Factor 1 lipid weight.
Copying it to several candidate source nodes creates candidate signalling
hypotheses; it does not establish ligand potency, receptor occupancy, flux, or
causality. Feature-to-many-node copies are intentional here. A source node that
receives values from multiple *different* measurements should be scored once as
a unique topology seed; retain its collision rows for assignment-specific
sign-consistency and path-extraction review rather than counting it twice.
MRCLinksDB and the nonpeptidic CellPhoneDB adapter are reported in
`resource_policy_summary.csv` as excluded annotations. STITCH v5 is directionally
auditable but historical, and remains exploratory.
