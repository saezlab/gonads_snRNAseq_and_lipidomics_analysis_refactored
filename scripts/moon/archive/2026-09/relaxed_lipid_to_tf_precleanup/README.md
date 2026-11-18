# Pre-cleanup relaxed lipid-to-TF scripts

This snapshot preserves the four scripts that generated the original relaxed
lipid-to-TF discovery output in
`miscellaneous/older_results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory/`.  They were
copied here before the active workflow was reorganized on 2026-09-16.

The historical discovery output remains available for provenance, but these
scripts are not supported reproduction entry points.  In particular, the old
grid could materialize ordinary sign-consistent connectivity, whereas the
active parameter search materializes only strict consecutive-score-layer
solutions.  The active entry points are:

* `scripts/moon/run_factor1_lipid_to_tf_relaxed.R` for the fixed selected LPS
  sensitivity;
* `scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R` for
  the independent parameter search;
* `scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_candidates.R` and
  `scripts/moon/ancillary/audit_factor1_lipid_to_tf_relaxed_score_layer_evidence.py`
  for saved-result audits.

Do not run the archived scripts as current analysis commands.  They retain
their former paths and semantics solely to make the legacy discovery result
auditable.
