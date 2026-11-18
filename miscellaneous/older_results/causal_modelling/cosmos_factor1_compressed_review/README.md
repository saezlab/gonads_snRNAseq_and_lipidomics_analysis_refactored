# Compressed Factor 1 MOON review

The integrated report is [docs/cosmos_factor1_compressed_review.md](../../docs/cosmos_factor1_compressed_review.md).

- **Recommended broader exploratory network:** `threshold_search/C03/network.graphml` (compressed radius 6/depth 4, cutoffs 1.85/0.55, 153 nodes/388 edges, Ahr/Hnf4a/Ppara/Tbx2 and three DGs). Stronger compression reproduces the exact signed graph at the same cutoffs.
- **Stricter baseline:** `threshold_search/C01/network.graphml` (2.0/0.85, 255 edges, Ppara/Tbx2 and two DGs).
- **Wider-topology comparison:** `sensitivity_threshold_search/C01/network.graphml` (compressed radius 12/depth 8, 2.15/1.075, 468 edges, Klf10/Tbx2 and two DGs). Changed regression context prevents interpreting higher cutoffs as greater biological confidence.
- **Why Esr1, Nfic and Dbp are missing:** `tf_audit/README.md`. They survive preprocessing, are scored, then omitted because weak positive scores oppose their negative supplied TF activities. All 66 audited TF/run comparisons preserve this exclusion reason.
- **Complete score visibility:** `score_retention/radius6_depth4_compressed_scores_all.csv` distinguishes raw diagnostic scores, returned scores and final-candidate membership; diagnostic rows do not bypass sign filtering.
- **Compression implementation:** `compression_method_audit.md`. Previous uncompressed runs really skipped same-children compression; preserve them as diagnostics.
- **Search inventories and validation:** `threshold_search/README.md`, `sensitivity_threshold_search/README.md`, each directory's `grid_complete.csv` and `selected_validation.csv`.
- **Figures:** `figures/compressed_review_summary.pdf` and each main candidate's `figures/candidate_paths_network_network.pdf`.

7,563 additional state/cutoff pairs were tested. Candidate IDs are local to their search directory. All ten exported candidates pass exact score, PKN, sign, cutoff and connectivity checks. Biological specificity remains limited by the measured-species reaction hypotheses and the narrow set of retained DGs. No cutoffs are calibrated significance thresholds.
