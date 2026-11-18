# Measured lipid → directional STITCH mapping audit

This audit compares all **1,254** lipidomic features and the
frozen **155** Factor 1 PKN inputs with the raw,
directionality-audited STITCH v5 chemical-to-mouse-protein overlay. It does
not alter a PKN or run MOON.

The overlay contains **7,481** retained source records across
**1,986** chemical sources. Its own validation report
passed all checks. The existing Factor 1 inputs contain
**149** measured-species rule nodes and
**6** exact ChEBI nodes; none is
a directional STITCH chemical-source node.

## Interpretation contract

`exact_CHEBI_identity` is the only tier labelled direct measured chemical
identity. `explicit_conjugate_acid_base` is a charge-state counterpart and
`same_formula_is_a_descendant` resolves additional structural detail through
the ChEBI ontology; neither is a direct measured chemical claim. The latter is
limited to an exact reported-resolution lipid and same-formula/same-charge
descendants, but can still assume double-bond position, geometry, or
stereochemistry that the measurement does not resolve.

`feature_stitch_chemical_matches.csv` records feature-to-chemical mappings;
`feature_stitch_protein_sign_matches.csv` collapses duplicate source records
into feature/source/protein/sign triples and marks triples retained after the
assembled graph's conflict filter. `mapping_tier_summary.csv` reports all
SwissLipids mapping tiers separately for all measurements and for the existing
Factor 1 inputs. `feature_inventory.csv` retains the zero-match majority.

The STITCH source files are v5.0 and were last modified in 2016. Their raw
directionality is auditable in this overlay, but it does not establish current
evidence, direct mouse experimental support, hepatocyte context, or causal
regulation in this dataset.
