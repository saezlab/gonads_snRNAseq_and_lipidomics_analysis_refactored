# Repository organization and provenance

The cleanup on 16 November 2026 promotes the approved expanded-marker analysis without changing the approved scientific result files. Moves preserve bytes; copies are hash-verified. The cleanup procedure itself did not delete files or perform Git publication; later user-authorized commits and pushes are recorded in Git history.

## Maintained locations

| Location | Role |
| --- | --- |
| `scripts/workflows/` | Numbered entry points; 17 is the current zonation/diet chain |
| `scripts/core/`, `scripts/config/`, `scripts/visualization/` | Shared implementations, release registry and plotting code |
| `scripts/core/engines/` | Frozen original 05/07 source bodies used by validated adapters |
| `scripts/moon/`, `scripts/pkn/` | Active optional causal workflows and resources |
| `scripts/maintenance/` | Organization, path resolution and GitHub checks |
| `tests/` | Maintained tests and independent numerical QA; generated test logs are not source |
| `data/RNA/pseudobulk_csv/` | Public pseudobulk CSV parts under 50 MB, with verified source values and original-object provenance |
| `data/RNA/zonation/current/` | Local expanded-marker source objects (excluded from Git), alongside retained scoring scripts and small compatibility CSVs |
| `results/current/` | Published results by Results section |
| `results/causal_modelling/` | Active optional causal outputs and priors |
| `results/legacy/` | Git-ignored complete prior result roots and access to the historical archives |
| `miscellaneous/` | Historical results, comparison branches, source, refactoring, scratch outputs and local state |

## Manuscript layout promotion

On 20 November 2026, the verified current results were published under
[the six Results sections](../results/README.md). The 16 November release remains
the scientific baseline. The [promotion ledger](provenance/result_layout_20261120.csv)
records each published path, original path, size and SHA-256. It also identifies
the verified replay additions and the two Figure 5C tables exported with only
the current rows. Preserved scientific files were copied byte for byte.

The previous flat results, including the complete 16 November bundle, were
moved into the ignored `results/legacy/pre_reorganization_baseline/` archive.
The local relocation inventory verifies 807 archived files and links. Output
versions, reports and comparisons omitted from the new tree remain available
there. Former tracked paths can also be recovered from Git commit
`8883c449f1a34675c2059e4d950264b19dffc994`. Private reviewer materials and meeting exports remain local.

Current workflow paths are documented in the [registry](../scripts/config/analysis_release.json).
Recorded paths resolve to a byte-identical published artifact when present
and to the legacy copy when a file was kept locally. Archived records and
their embedded paths retain their original bytes.

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
- [Complete workflow reference](WORKFLOW_REFERENCE.md) and [dependency inventory](workflow_dependencies.csv): exact inputs, producers, execution modes and availability, maintained through [the JSON contract](workflow_dependencies.json).
- [Documentation audit](provenance/workflow_documentation_audit_20261118.json): source review, preserved documentation snapshots, guard-only code changes and verification.

Run `python3 scripts/maintenance/verify_repository_organization.py` to verify the moves, stable inputs, historical reference resolution and approved release. Run `python3 scripts/maintenance/audit_github_readiness.py --strict` for the prospective Git tree. These checks do not stage or commit files.

## Resolving original paths

Frozen manifests keep the paths and hashes recorded at production time. They are evidence, not an editable directory index. The registry's `relocations` maps old prefixes to current archive paths. R uses `analysis_resolve_recorded_path()` from `scripts/config/analysis_release.R`; Python uses `resolve_recorded_path()` from `scripts/maintenance/repository_paths.py`.

A relocated file must still match its original SHA256. The resolver does not reinterpret the original two-marker reference as an expanded-marker result. Machine-specific historical records are retained; new default input paths are repository-relative.

## GitHub and restoration

Archive payloads remain in the local checkout; selected code, documentation and manifests are included in Git. Large generated payloads, cached models, contact sheets, repeated bundles, local settings and private manuscript source are excluded. The archive is not discarded just because `.gitignore` hides a payload.

To restore an excluded file in another checkout, locate it by its `old_path`/`new_path` in the migration ledger, copy the original bytes from the retained archive or a separately distributed artifact, then verify the SHA256 before running a historical branch. No public download URL is invented: an external data deposit has not been configured by this cleanup. Current default scientific runs use the manifest-selected expanded-marker CSV exports; manuscript auditing and older comparative branches may need additional archived/private inputs.

The optional MOON full-fit environment still enforces its original decoupleR version pin. Saved-artifact relocation/extraction checks passed; the local installed version does not meet that separate full-fit contract. See the causal map for exact versions and checks. The gate was not weakened and no package upgrade was made.
