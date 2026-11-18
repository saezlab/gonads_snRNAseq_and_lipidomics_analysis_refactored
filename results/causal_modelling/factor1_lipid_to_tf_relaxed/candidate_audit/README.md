# Relaxed lipid-to-TF candidate audit

This ancillary audit reads the completed selected refit and does not repeat PKN preparation, pruning, compression, RNA coherence or ULM scoring.
`candidate_hypothesis_summary.csv` retains every feature-to-chemical mapping, its copied Factor 1 loading, candidate score/status and score-explanation flags.
`candidate_edge_fingerprint_summary.csv` prevents chemically unresolved candidates with identical signed direct PKN neighborhoods from being counted as independent evidence.
`candidate_score_explanations/` contains the compressed, layerwise `get_moon_scoring_network()` trace and direct mapping evidence for each scored candidate.
These traces support conditional hypotheses only; they do not resolve lipid identity or demonstrate signaling in hepatocytes.
