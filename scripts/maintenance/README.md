# Repository maintenance

Run from the repository root. None of these tools stages, commits or pushes.

| Tool | Purpose |
| --- | --- |
| `audit_github_readiness.py --strict` | Prospective Git tree, file-size/local-state gates and nonshipping dependency/link ledgers |
| `verify_repository_organization.py --require-full-archive` | Verify every retained local archive/result/input record and approved scientific bytes |
| `verify_repository_organization.py` | Verify available records and mandatory current inputs; report absent external payloads separately |
| `build_analysis_index.py` | Rebuild portable current panel navigation and workflow catalog from approved manifests |
| `export_pseudobulk_csv.R` | Export or verify seven declared local pseudobulks as CSV parts under 50 MB; source RDS stays local |
| `build_pseudobulk_input_index.py` | Rebuild the public CSV/JSON input hash index; never includes RDS objects |
| `repository_paths.py` | Shared Python resolver for original paths in immutable provenance |
| `organize_repository.py` | Explicit dated migration implementation; do not rerun `--apply` on an already migrated checkout |

The release selection and prefix relocations are in `scripts/config/analysis_release.json`. The R counterpart is `scripts/config/analysis_release.R`. Exact moves, approved input copies and result inventories are under `docs/provenance/`.

The cleanup changes paths, defaults and documentation, not scientific results. To adopt a future analysis release, first validate it and record its provenance, then intentionally update the release registry and indices. Do not overwrite the approved production files or rewrite frozen manifests to conceal changes.
