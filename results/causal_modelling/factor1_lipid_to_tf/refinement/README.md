# Lipid-to-TF threshold refinement

This targeted refinement reads the frozen fits listed in dependency_manifest.csv; no MOON scores, observations, mappings or PKN edges are changed.
Only runs with nonempty coarse solutions are explored. Bounds and first-edge bottlenecks are explicit in each runs/<id>/ directory.
All extractions use moon_extract_paths and moon_validate_solution. Eligible means 0 < edges < 500, >=1 measured upstream lipid and >=2 observed downstream TFs.
all_tfs is the primary observation panel. strong_tfs is an observation-set sensitivity and is not pooled with it for biological ranking.
A compact output is descriptive; it does not establish direct lipid signaling, metabolic flux or a causal relationship.

The selected exploratory all_tfs network is in selected_all_tfs/: primary 1.65, secondary 1.10, 455 edges, 233 nodes, one LPC input and 64 observed TF endpoints. The neighborhood table documents threshold sensitivity. Reproduce the refinement and selection with `Rscript scripts/moon/ancillary/refine_lipid_tf_thresholds.R --output=NEW_EMPTY_DIRECTORY`; use --select-only only when selected_all_tfs does not yet exist in an existing refinement. Existing selected outputs are refused.

The fine grid contains 1,431 extractions: 26 eligible all_tfs combinations and no eligible strong_tfs combinations. See provenance/ for the exact initial grid-execution code snapshot and selected_all_tfs/selection_code_manifest.csv for final code hashes.
