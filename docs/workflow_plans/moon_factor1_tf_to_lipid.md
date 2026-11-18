# Optional MOON workflow: selected Factor 1 TF-to-lipid network

## Purpose and scope

`scripts/moon/run_factor1_tf_to_lipid.R` regenerates the selected manuscript network from the existing, audited Factor 1 inputs and lipid-expanded mouse PKN. It is separate from workflows 00–10 and does not retrain MOFA or rebuild the PKN. The lipid-to-TF analysis has its own entry point and plan, `moon_factor1_lipid_to_tf.md`.

Read this plan and `scripts/moon/README.md` before changing the workflow.
Preserve the selected configuration and reference comparison while adding
ancillary analyses. The base prior remains under
`results/causal_modelling/cosmos_mouse_pkn/`; historical pilot/search/selection
references now live under `miscellaneous/older_results/causal_modelling/`.
Those reference folders are inputs, not destinations for a new selected fit.

## Entry points and shared code

| Component | Path | Role |
|---|---|---|
| Main workflow | `scripts/moon/run_factor1_tf_to_lipid.R` | Executes the selected run, diagnostics, reference comparison and export. |
| Selected settings | `scripts/moon/config/factor1_tf_to_lipid.json` | Records input/output locations, the existing direction and selected numerical settings. |
| Shared operations | `scripts/moon/lib/workflow.R` | Configuration/input validation, topology preparation, compression, scoring/coherence, extraction and audits shared by main and search runners. |
| Compatibility helpers | `scripts/moon/moon_safe_helpers.R` | Existing safeguarded compression, complete score decompression and guarded dual-threshold path extraction; preserved unchanged during cleanup. |
| Ancillary parameter grid | `scripts/moon/ancillary/search_factor1_parameters.R` | Reuses the same compressed scoring core for explicit pruning/layer/threshold grids and records every combination. |
| Export/plots | `scripts/moon/ancillary/export_factor1_network.py`, `scripts/moon/ancillary/plot_factor1_network.R` | Export readable metadata, exact PKN evidence, GraphML, route illustrations and figures. |
| Historical implementations | `scripts/moon/archive/2026-09/` | Byte-preserved pilot, searches, sensitivities, audits, templates and migration manifest. Not active reproduction commands. |

## Upstream dependencies

| Input | Required content |
|---|---|
| `results/causal_modelling/cosmos_mouse_pkn/pkn_lipid_rules.rds` | Signed mouse PKN with measured-species lipid reaction hypotheses. |
| `results/causal_modelling/cosmos_mouse_pkn/` | Node catalogue and exact edge/reaction provenance for annotation and export. |
| `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds` | Jointly oriented TF activities, exact-mapped lipid loadings, hepatocyte RNA loadings, TF regulons and mapping audit. |
| `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_lipid_node_mapping.csv` | Accepted exact observed-feature/node mapping for exports. |
| `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth4_compressed/` | Frozen selected scoring state and final graphs for comparison. |
| `miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03/` | Frozen manuscript solution, not a candidate with the same name in another search directory. |

The accepted input bundle was prepared from workflow 04 outputs, the saved MOFA HDF5, sample metadata, CollecTRI/hormone regulons and the PKN crosswalk. The main workflow consumes this frozen bundle; it does not repeat the upstream data preparation or claim to reconstruct missing source data. Original input-construction scripts and manifests remain available in the archive and frozen input folder.

R dependencies include cosmosR, decoupleR, igraph and jsonlite; export uses Python and the R plotting dependencies declared by the active scripts. Consult the run's recorded package versions and setup checks. The archived baseline used cosmosR 1.21.2, decoupleR 2.17.0 and igraph 2.3.3; a version change requires examining the replay comparison, not assuming numerical equivalence.

## Main workflow stages

1. Validate settings, files and signed PKN/input contracts. TF footprints, lipid loadings and RNA coherence evidence retain separate semantics. `--validate-only` stops after validation without writing output.
2. Restrict upstream inputs to PKN-present TFs with `abs(activity) > 2`; retain mapped, nonzero lipid observations. Remove outgoing observation edges, then alternate directed upstream/downstream reachability pruning until stable. The selected pruning distance is six steps per search.
3. Compress nodes sharing signed child sets after protecting upstream/downstream input IDs and TF-regulon IDs, including the entire duplicate group of a protected parent. Save and verify the complete original-to-compressed mapping.
4. Score the compressed PKN using ULM with four MOON layers. Apply the original RNA-coherence filtering after scoring and repeat for at most ten passes; preserve the fixed initial compression map through the selected baseline loop. Record edge removals and convergence.
5. Decompress scores, project the final scoring graph back to original IDs, and extract the network at strict primary/secondary absolute-score thresholds 1.85/0.55. Preserve observed weights at level 0, require sign consistency and directed TF-to-observation connectivity.
6. Diagnose expected TF retention using final-network scores before the upstream sign filter. Keep diagnostic scores separate from the production solution. Esr1/Nfic/Dbp should be classified by the stage at which they are lost, not assumed pre-pruned.
7. Compare the refit with the frozen selected scores/graph, write provenance and export node/edge tables, GraphML and route figures. Report execution success only when the generated checks pass.

## Output contract and replay checks

Default output root: `results/causal_modelling/factor1_tf_to_lipid/`. Existing nonempty output directories are refused; use `--output=NEW_DIRECTORY` for a replay. The workflow keeps settings and software/input provenance, pruning and coherence ledgers, compression mapping, final scoring and original-ID graphs, score tables, extracted SIF/ATT, expected-TF diagnostics, reference comparison, evidence exports and figures together.

Key entry files are `results_summary.json`, `reference_comparison.json`, `run_manifest.json`, `dependency_manifest.csv`, `output_manifest.csv`, `network_SIF.csv`, `network_ATT.csv`, `network.graphml`, `moon_scores.csv`, `moon_scores_all.csv` and `focus_tf_retention.csv`. Exact graph and mapping RDS objects, coherence/pruning ledgers, reaction/path evidence and `figures/` provide the detailed audit trail.

The complete main workflow was executed successfully on 2026-09-15. The refit matched the frozen C03 edge/node sets and compression mapping, all score IDs and levels, and numerical scores within 1e-10. Three coherence iterations, 49 TF inputs, 74 lipid inputs, 2,238 pre-sign-filter scores and 2,207 returned scores were reproduced. Expected-TF exclusions, evidence export and plotting passed. This verifies the selected refit from frozen mapped inputs/PKN, not reconstruction of raw inputs or the PKN resource.

The frozen reference has 388 signed edges and 153 original-ID nodes connecting Ahr, Hnf4a, Ppara and Tbx2 to three DG observations. Its baseline converged in three passes. The initial 44,127-edge topology compresses to 36,569 edges; the final scoring graph has 27,746 edges and projects to 35,304 original-ID edges. There are 1,680 compressed returned score rows and 2,207 decompressed returned rows. These are replay expectations from saved results, not assertions that every future environment or parameter choice will reproduce them.

Required checks include finite unique named inputs/scores, signed unique PKN edges, unchanged protected identities, exact graph projection, final-PKN membership of solution edges, edge-wise score/sign agreement, terminal observations and directed TF-to-observation connectivity. Compare signed edge and node identities, aligned score values and scoring levels with the frozen reference; counts alone are insufficient. Neither six pruning steps nor four scoring layers is a hard bound on final displayed route length.

## Ancillary search

The reusable compressed-only search defaults to pruning distance 6 and layers 4,5,6,8. An explicit small grid can be run from the repository root as:

```sh
Rscript scripts/moon/ancillary/search_factor1_parameters.R \
  --pruning-steps=6,12 --layers=4,8 \
  --primary=1.85,2.15 --secondary=0.55,1.075 \
  --output=results/causal_modelling/factor1_parameter_exploration
```

Every valid threshold combination is retained, including empty solutions; scoring is reused within a pruning/layer configuration. The historical 7,563 additional threshold evaluations include distinct sensitivity policies and adaptive refinements. This generic runner does not claim to reconstruct that historical search or to replace its selection evidence. Historical grids remain under `miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/`.

## Maintenance and known limitations

- Keep this entry point's direction TF-to-lipid. The separate lipid-to-TF entry point has its own measurement/PKN decision, settings and output identity; it is not a source/target-column swap.
- Preserve frozen historical RDS/CSV/JSON path strings and source-hash keys. The new workflow uses current explicit input paths; the archive's old `results/cosmos_*` paths are not active commands.
- The archive preserves the final exploratory source files, not a claim of byte-identical source snapshots for every historical execution. `path_migration_manifest.csv` verifies the cleanup relocation independently.
- The selected PKN contains explicitly hypothetical measured-species reactions. Keep their evidence labels through export; these hypotheses and MOON scores are not experimental validation.
- Run any changed R script immediately, and update this plan, `scripts/moon/README.md` and the static result-query map if its input/output contract changes.

## 16 November 2026 repository cleanup

Historical causal result folders now live under
`miscellaneous/older_results/causal_modelling/`. Active configurations and
consumer defaults explicitly use their frozen inputs/reference states there.
The per-file old/new path and SHA256 ledger is
`docs/provenance/causal_cleanup_20261116.csv`; the active/frozen/archive dependency
map is `docs/CAUSAL_ANALYSIS_MAP.md`. No numerical settings, graph policy, mapping,
score or fitted result changed. Saved-artifact extraction/reference replay passed
without fitting or scoring. Full-fit version gates remain unchanged; the installed
decoupleR 2.9.7 differs from the pinned 2.17.0 and blocks full-fit validation.
