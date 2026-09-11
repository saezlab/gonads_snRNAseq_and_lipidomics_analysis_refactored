# Factor 1 MOON search summary inputs

Required columns: `run_id`, `run_dir`, `compression`, `n_layers`, `topology_radius`, `primary_threshold`, `secondary_threshold`, `status`, `converged`, `valid`, `nodes`, `edges`, `upstream_tfs`, `downstream_nodes`, `solution_stem`.

Optional columns: `phase`, `realized_max_level`, `candidate_edges`, `consistent_edges`, `short_route_pairs`, `route_pairs`, `run_manifest`.

A row describes one score configuration and threshold pair. `valid` certifies the scoring/extraction pipeline's provenance, sign and connectivity checks; `converged` certifies its own coherence loop. Only `status=complete` passes. Missing certification does not pass. `run_dir` and `solution_stem` locate eventual `STEM_SIF.csv` and `STEM_ATT.csv` files. Solutions may be materialized after selection; the summary records that availability separately.

Rank secondary cutoff first, primary cutoff second, then lipid and TF coverage, then fewer edges. Require nonempty valid/converged solutions with at least one TF and one lipid and <=500 edges. A separate interpretable subset requires at least two TFs and two lipids. Pareto frontiers are computed within each compression/depth/topology-radius configuration. Cross-configuration counts describe eligibility reproducibility, not biological network stability.
