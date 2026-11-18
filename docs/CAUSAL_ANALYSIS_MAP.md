# Causal analysis: active workflows and frozen references

This optional branch uses saved MOFA Factor 1 inputs and mouse prior-knowledge
networks. It is independent of the expanded-marker zonation release. Repository
cleanup on 16 November 2026 moved historical exploration; it did not refit MOFA
or MOON, compute new scores, change thresholds, remap observations or change a
selected network.

## Active entry points and results

Run commands from the repository root. Read the corresponding `moon_*.md` plan
in `docs/workflow_plans/` and [the MOON guide](../scripts/moon/README.md) first.

| Analysis | Entry point | Current result directory under `results/causal_modelling/` |
|---|---|---|
| Selected TF activities → lipid weights | `scripts/moon/run_factor1_tf_to_lipid.R` | `factor1_tf_to_lipid/` |
| Exact-mapped lipid → TF hypotheses | `scripts/moon/run_factor1_lipid_to_tf.R` | `factor1_lipid_to_tf/` |
| Relaxed regulatory lipid → TF sensitivity | `scripts/moon/run_factor1_lipid_to_tf_relaxed.R` | `factor1_lipid_to_tf_relaxed/` |
| Relaxed threshold/depth exploration | `scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R` | `factor1_lipid_to_tf_relaxed_parameter_search/` |
| Chain-compatible TF → lipid sensitivity | `scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R` | `factor1_tf_to_lipid_expanded_species/` |
| Reviewed experimental target proposal | `scripts/moon/ancillary/propose_factor1_ko_validation.py` | `factor1_ko_validation_targets/alternating_prioritization/` |

Active implementation remains in `scripts/moon/{run*,config,lib,ancillary}` and
`scripts/pkn/`. `scripts/moon/moon_safe_helpers.R` retains its frozen bytes and
hash. The base PKN is `results/causal_modelling/cosmos_mouse_pkn/`; expanded
species and directional-STITCH overlays remain separate sibling resources.
Mapping/coverage/evidence audits and their current result directories are
retained. The optional KO expression audit reads its four historical culture
controls from `miscellaneous/older_inputs/legacy_rna/merged_counts_matrix_250508_VH00211_433_AACVFM3HV_with_gene_symbols_and_groups.csv`; this archived input remains a live dependency. Optional sensitivity outputs are not automatic replacements for the
selected TF → lipid analysis.

## Archived results that remain live inputs

The following complete folders moved from `results/causal_modelling/` to
[`miscellaneous/older_results/causal_modelling/`](../miscellaneous/older_results/causal_modelling/).
Active configurations and consumer defaults now name the archive explicitly.

| Folder | Why retained | Live dependency within it |
|---|---|---|
| `cosmos_factor1_pilot/` | Original pilot, input preparation and evidence | `inputs/factor1_inputs.rds`, mapping CSV and other audited input resources |
| `cosmos_factor1_parameter_search/` | Original parameter search and representation diagnostics | `runs/radius6_depth4_compressed/` and `topology_radius_6/compressed/compression.rds` |
| `cosmos_factor1_compressed_review/` | Historical selection and sensitivity evidence | `threshold_search/C03/` and `score_retention/radius6_depth4_compressed_scores_all.csv` |
| `factor1_lipid_to_tf_relaxed_regulatory/` | Superseded discovery topology/search | `runs/direct_strong_tfs_pruning4_layers4/run_state.rds` and its strict `score_layer_constrained/solutions/...__p0.5__s0.5/` reference |

These are frozen resources, not supported output destinations. Historical
`run_dir`, absolute paths, source hashes and manifests inside them are
unchanged. They describe original execution and must not be edited merely to
match the new directory names. Current code does not follow the obsolete
embedded paths to locate these inputs.

[The migration ledger](provenance/causal_cleanup_20261116.csv) records every old
and new relative file path, byte count, SHA256, role and dependency. Its
[directory summary](provenance/causal_cleanup_20261116.json) describes the move
units. Four result folders contain 1,283 files totalling 226,359,066 bytes.

[The minimal public dependency inventory](provenance/active_frozen_dependencies.json)
identifies 13 compact archived files (4,074,361 bytes total): 11 for active main
entry points/reference checks and two further files for expanded-species input
preparation. These retain full-precision inputs and frozen reference states
that current outputs and zonation inputs do not replace. Other discovery
payloads can remain outside the public checkout under the documented artifact
policy.

Historical implementation remains in `scripts/moon/archive/2026-09/` with its
earlier 29-file SHA256 migration manifest. Those scripts retain original path
assumptions and are not current commands. The separate initial PKN feasibility
probe moved unchanged to `miscellaneous/quick_tests/code/audit_cosmos_pkn.py`;
it is also historical and is included in the cleanup ledger.

The archived KO proposal PDF is now at
[`miscellaneous/quick_tests/output/pdf/cosmos_factor1_ko_validation_targets.pdf`](../miscellaneous/quick_tests/output/pdf/cosmos_factor1_ko_validation_targets.pdf).
Frozen KO result README/source snapshots can retain the former `output/pdf/`
link as historical provenance; this map gives its current physical location.

## Verification without refitting

```sh
Rscript scripts/moon/ancillary/validate_frozen_archives.R
Rscript scripts/pkn/validate_export.R --validate-only
```

The first command checks frozen input hashes, saved score/level/compression
identities and re-extracts the selected networks from saved scores. The selected
388-edge/153-node network and relaxed 5-edge/6-node network match their archived
references. It performs no MOON fitting or scoring. The second checks all five
base PKN CSVs against the installed COSMOS cleaner and writes nothing.

[Saved-artifact replay](provenance/causal_cleanup_20261116_saved_replay.json)
records the exact checks and environment. Full-fit `--validate-only` contracts
were deliberately preserved: the current environment has **decoupleR 2.9.7**,
whereas the frozen full replay requires **2.17.0**. The three main MOON
`--validate-only` commands fail at this pre-existing software gate. This cleanup
does not claim that a new fit is validated under the installed version; no
package was upgraded and no gate was relaxed.

## PKN rebuilding and provenance

Base PKN builders, exporters and their audit-cache references consistently use
`results/causal_modelling/cosmos_mouse_pkn/` and
`results/causal_modelling/cosmos_pkn_audit/`. This repairs the former mixture of
`results/cosmos_*` and relocated paths. Raw sources in
`support/cosmos_mouse/raw/` and `cosmos_pkn_audit/raw/` remain available.

`scripts/pkn/run_build.py` is a resource rebuild, not a cleanup or smoke test;
it writes its declared resource directory. The cleanup ran compilation,
small tests and saved-export validation, without rerunning the full builder or
changing frozen PKNs. Existing output/dependency manifests retain original
execution paths and source hashes; new code changes do not retroactively alter
those historical records.

Scores and routes remain exploratory prior-guided hypotheses. Supplied TF
observations, conditional chemical mappings, measured bulk-liver abundance and
inferred node scores retain their distinct meanings and original limitations.
