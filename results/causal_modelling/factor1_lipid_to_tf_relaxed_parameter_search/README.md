# Relaxed Factor 1 lipid-to-TF parameter search

This ancillary workflow uses the same score-first candidate semantics and strict consecutive score-layer extraction as the selected LPS refit.
Candidate measurements are copied only as post-score sign-consistency evidence; they are not numerical ULM inputs.

`fit_results.csv` records every requested compressed fit. Each successful `runs/<run_id>/` directory contains its fit state, compression map, coherence ledger, scores and candidate-status tables.
`grid_results.csv` records every threshold pair, including invalid, failed and empty extractions. Ordinary sign-consistent connectivity is intentionally not generated as a solution output.
`eligible_under_500.csv` requires a nonempty strict solution with fewer than 500 edges, a nonzero secondary cutoff, at least one independent direct-edge fingerprint group and at least two observed TFs.
`score_layer_stability.csv` groups identical signed SIFs across parameter settings. Only the first occurrence of each topology is materialized; later rows point to the same `solution_ref`.
Threshold and coverage frontiers describe tradeoffs only and do not define an optimal or high-confidence network.
