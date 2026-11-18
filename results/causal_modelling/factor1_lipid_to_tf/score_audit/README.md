# Saved lipid-to-TF score audit

Scores were recomputed only on each saved final compressed scoring graph, with upstream_input=NULL to expose the package's final sign-filter exclusions. Every originally returned score and level was compared with the recomputation. No coherence loop, input mapping, parameter fitting or network extraction was rerun.

`lipid_input_audit.csv` diagnoses all 155 frozen accepted lipid inputs in every saved scored run. `two_lpc_score_summary.csv` isolates the two LPCs that support normalized precursor topology. Absent policy nodes, topology loss, final-graph loss, missing layer support, and upstream sign exclusion are distinct statuses.

`focus_tf_audit.csv` distinguishes Esr1, Nfic, Dbp and Rxra supplied as level-0 observations from intermediates with inferred scores. An observed activity is not an independently predicted TF activity.

`focal_node_scores.csv` and `scoring_support_edges.csv` report the relevant LPC/PC/DHA intermediates, reaction instances, and TFs, using final compressed score-support edges. Counts refer to signed children in the previous scoring layer and its compressed universe. They are not independent evidence counts or statistical confidence intervals.

Low or conflicting inferred lipid scores cannot be repaired by lowering extraction cutoffs alone. Precursor paths can depend on unmeasured cosubstrates; the LPC16/PC/DHA acyl-remodeling route regenerates LPC and does not establish net DHA production.
