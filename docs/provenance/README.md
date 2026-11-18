# Provenance records

All CSV paths in the cleanup manifests are repository-relative. Original production manifests can contain absolute paths; retain those identities and use the release resolver to find relocated files.

| Record | Meaning |
| --- | --- |
| `repository_cleanup_20261116_qa.json` | Consolidated cleanup, preservation, scientific equivalence and public-tree setup evidence |
| `repository_migration_20261116.csv` | Applied repository moves and documentation snapshots, with original/new path, role, size and SHA256 |
| `causal_cleanup_20261116.csv` / `.json` | Applied causal archive moves and whole-directory mappings |
| `archive_recovery.csv` | Exact historical Git object recovery versus local archive required; working-file hashes are never replaced with a nonidentical Git blob |
| `pseudobulk_csv_distribution_20261118.json` | CSV-only distribution gates, exact source round-trip and source-RDS-free public-tree smoke |
| `pseudobulk_csv_roundtrip.json`, `pseudobulk_csv_transport_qa.json` | Source-object equivalence and split/precision/metadata/corruption tests |
| `current_input_manifest.csv` | Current public pseudobulk CSV/JSON payload hashes and small preprocessing/compatibility interfaces |
| `current_input_manifest_20261116.csv` | Original local source-object/copy identities before CSV-only distribution; RDS objects are excluded from Git |
| `current_release_artifacts.csv` | Original approved release files protected by cleanup verification; subsequently added navigation is not an analysis modification. Finder .DS_Store is explicitly classified local_metadata, retaining its original hash but excluded from scientific verification |
| `current_panel_manifest.csv` | Portable public panel/source index derived from production `panel_trace.csv` |
| `causal_cleanup_20261116_qa.json`, `_saved_replay.json` | Causal path/byte/extraction replay checks and retained environment gate |
| `active_frozen_dependencies.json` | Minimal frozen reference files retained for active causal workflows, when required by their configs |

`repository_migration_plan_20261116.csv` is the initial preflight inventory; the applied migration ledger is authoritative and includes clarified dependency roles. Large payloads remain local when excluded from Git. SHA256 proves identity, not public availability; no external deposit was invented during this cleanup.

Regenerate local verification with `python3 scripts/maintenance/verify_repository_organization.py --require-full-archive`. Definitive current cleanup/readiness summaries are linked from the organization map. No data or model fitting is performed by these checks.
