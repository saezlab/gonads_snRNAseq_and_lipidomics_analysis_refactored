# Repository organization and provenance

The cleanup on 16 November 2026 promotes the approved expanded-marker analysis without changing the approved scientific result files. Moves preserve bytes; copies are hash-verified. Nothing was deleted, committed or pushed.

## Maintained locations

| Location | Role |
| --- | --- |
| `scripts/workflows/` | Numbered entry points; 17 is the current zonation/diet chain |
| `scripts/core/`, `config/`, `visualization/` | Shared implementations, release registry and plotting code |
| `scripts/core/engines/` | Frozen original 05/07 source bodies used by validated adapters |
| `scripts/moon/`, `scripts/pkn/` | Active optional causal workflows and resources |
| `scripts/maintenance/` | Organization, path resolution and GitHub checks |
| `tests/` | Maintained tests and independent numerical QA; generated test logs are not source |
| `data/RNA/pseudobulk_csv/` | Public pseudobulk CSV parts under 50 MB, with verified source values and original-object provenance |
| `data/RNA/zonation/current/` | Ignored local expanded-marker source objects, scoring scripts and small compatibility CSVs |
| `results/zonation_marker_update_20261116/` | Immutable current zonation, diet panels, source tables and completed revision audit |
| Other main `results/` roots | Unchanged cell-type, lipid, MOFA, validation and optional active causal outputs |
| `miscellaneous/` | Historical results, comparison branches, source, refactoring, scratch outputs and local state |

The dated current output root is intentional: renaming its embedded paths would obscure the original production provenance. Current consumers resolve it through [analysis_release.json](../scripts/config/analysis_release.json), rather than inferring a default from folder age or old filenames.

## Archive map

| Original location | New location | Dependency/status |
| --- | --- | --- |
| `results/zonation/` | `miscellaneous/older_results/zonation_two_marker/` | Original FCG model/networks; frozen old/new and 12/13 references |
| `results/cf453_{standalone,diet_impact,diet_rewiring}_figure6/` | Same basenames under `miscellaneous/older_results/` | Superseded Figure6 versions; old/new comparison only |
| `results/wdsd_network_comparison/`, `wd_crosscohort_robustness/` | Same basenames under `miscellaneous/comparisons/` | Original-compartment cross-cohort comparisons, not current Figure6 |
| Old `data/RNA/*zonation*` and WDSD zonation RDS | `miscellaneous/older_inputs/two_marker_zonation/RNA/` | Historical inputs kept separate from new defaults |
| `data/RNA/legacy/` | `miscellaneous/older_inputs/legacy_rna/` | Includes culture counts still used by optional KO expression audit |
| `results/mofa/RNA/old_*/` | `miscellaneous/older_results/mofa/` | Earlier trained-model generations; current unversioned model set stays in place |
| `results/multi-omic/legacy/` | `miscellaneous/older_results/multi_omic/` | Earlier structural analysis |
| Old zonation-specific validation exports | `miscellaneous/older_results/validation_zonation_two_marker/` | Current zonal validation lives in approved release `fcg07/`; experimental KD tables stay active |
| Four historical causal runs | `miscellaneous/older_results/causal_modelling/` | Some contain mandatory frozen inputs/reference states; see causal map |
| `scripts/legacy/`, `scripts/old/`, `unused/` | `miscellaneous/older_code/` | Original source retained for interpretation, not current entry points |
| `refactoring_docs/` | `miscellaneous/refactoring/` | Historical notes/debug code/logs |
| Ad hoc `quick_*.R`, output/tmp exports | `miscellaneous/quick_tests/` | Scratch analyses and generated outputs; useful `quick_lookup.R` remains maintained |
| Workflow/package/test logs | `miscellaneous/logs/`, `quick_tests/logs/` | Execution history, not numerical source |
| Machine settings/history | `miscellaneous/local/` | Local only; excluded from Git |
| Earlier README/PIPELINE/agent instructions | `miscellaneous/documentation/pre_cleanup_20261116/` | Documentation snapshots before promotion |

Exact file paths override this compact summary. Use the machine-readable ledgers below.

## Dependency and checksum ledgers

- [Cleanup verification summary](provenance/repository_cleanup_20261116_qa.json): preservation, fresh branch equivalence and public-tree setup checks.
- [Repository migration](provenance/repository_migration_20261116.csv): every moved/snapshotted file, old/new repository-relative path, SHA256, size and role.
- [Causal migration](provenance/causal_cleanup_20261116.csv): historical causal artifacts and the original exploratory PKN probe.
- [Public current inputs](provenance/current_input_manifest.csv): distributed CSV/JSON payloads, SHA256 and size. [The original local-object ledger](provenance/current_input_manifest_20261116.csv) retains the previous byte-identical source-copy identities. See [the CSV contract](PSEUDOBULK_CSV_INPUTS.md).
- [Approved release inventory](provenance/current_release_artifacts.csv): approved result files protected against accidental modification during cleanup.
- [Git-history recovery map](provenance/archive_recovery.csv): exact Git object references where pre-cleanup tracked bytes match; otherwise an explicit local-archive requirement.
- [Causal dependency map](CAUSAL_ANALYSIS_MAP.md): active workflows versus frozen pilot/selection evidence and optional sensitivities.
- [Analysis catalog](analysis_catalog.csv): current/optional/historical workflow roles and output roots.

Run `python3 scripts/maintenance/verify_repository_organization.py` to verify the moves, stable inputs, historical reference resolution and approved release. Run `python3 scripts/maintenance/audit_github_readiness.py --strict` for the prospective Git tree. These checks do not stage or commit files.

## Resolving original paths

Frozen manifests keep the paths and hashes recorded at production time. They are evidence, not an editable directory index. The registry's `relocations` maps old prefixes to current archive paths. R uses `analysis_resolve_recorded_path()` from `scripts/config/analysis_release.R`; Python uses `resolve_recorded_path()` from `scripts/maintenance/repository_paths.py`.

A relocated file must still match its original SHA256. The resolver does not reinterpret the original two-marker reference as an expanded-marker result. Machine-specific historical records are retained; new default input paths are repository-relative.

## GitHub and restoration

Archive payloads remain in the local checkout; selected code, documentation and manifests are included in Git. Large generated payloads, cached models, contact sheets, repeated bundles, local settings and private manuscript source are excluded. The archive is not discarded just because `.gitignore` hides a payload.

To restore an excluded file in another checkout, locate it by its `old_path`/`new_path` in the migration ledger, copy the original bytes from the retained archive or a separately distributed artifact, then verify the SHA256 before running a historical branch. No public download URL is invented: an external data deposit has not been configured by this cleanup. Current default scientific runs use the manifest-selected expanded-marker CSV exports; manuscript auditing and older comparative branches may need additional archived/private inputs.

The optional MOON full-fit environment still enforces its original decoupleR version pin. Saved-artifact relocation/extraction checks passed; the local installed version does not meet that separate full-fit contract. See the causal map for exact versions and checks. The gate was not weakened and no package upgrade was made.
