# Layer-constrained relaxed lipid-to-TF MOON extractions

This is a derived extraction from the completed score-first relaxed lipid-to-TF fits. It does not refit the PKN, pruning, compression, RNA coherence, or ULM scores. Before extracting routes, it checks the original dependency manifest.

Only sign-concordant feature-to-chemical assignments are eligible as source seeds. Each retained edge must descend exactly one recursive ULM scoring level and agree with the signed PKN and node scores. This distinguishes score-supporting routes from the broader sign-consistent connectivity screen in the parent analysis.

`score_layer_threshold_grid.csv` retains every threshold result. A compact result is listed in `eligible_score_layer_solutions.csv` only if it has fewer than 500 edges, at least one unique direct-edge fingerprint group, at least two observed TFs, and a nonzero secondary threshold. It is not a statistical model selection or proof of mechanism. Each materialized solution retains raw candidate-node count and fingerprint-group count so alternative mappings with identical direct PKN profiles are not counted as independent regulators.

The source measurement may be copied to several chemically plausible nodes. Mapping tier, Factor 1 loading/rank, and direct-edge provenance are saved with each materialized network. Bulk abundance, ligand exposure, receptor engagement, temporal ordering and causal regulation remain unmeasured.
