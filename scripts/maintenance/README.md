# Repository maintenance

Run from the repository root. None of these tools stages, commits or pushes.

| Tool | Purpose |
| --- | --- |
| `audit_github_readiness.py --strict` | Prospective Git tree, file-size/local-state gates and nonshipping dependency/link ledgers |
| `verify_repository_organization.py --require-full-archive` | Verify every retained local archive/result/input record and approved scientific bytes |
| `verify_repository_organization.py` | Verify available records and mandatory current inputs; report absent external payloads separately |
| `build_analysis_index.py` | Rebuild current panel navigation from approved manifests, then delegate the complete workflow catalog/reference to the documentation builder |
| `build_workflow_documentation.py` / `--check` | Render / check all workflow contracts, source coverage, registry keys, input availability and generated documentation; never runs scientific entries |
| `export_pseudobulk_csv.R` | Export or verify seven declared local pseudobulks as CSV parts under 50 MB; source RDS stays local |
| `build_pseudobulk_input_index.py` | Rebuild the public CSV/JSON input hash index; never includes RDS objects |
| `repository_paths.py` | Shared Python resolver for original paths in immutable provenance |
| `organize_repository.py` | Explicit dated migration implementation; do not rerun `--apply` on an already migrated checkout |

The release selection and prefix relocations are in `scripts/config/analysis_release.json`. The R counterpart is `scripts/config/analysis_release.R`. Exact moves, approved input copies and result inventories are under `docs/provenance/`.

The cleanup changes paths, defaults and documentation, not scientific results. To adopt a future analysis release, first validate it and record its provenance, then intentionally update the release registry and indices. Do not overwrite the approved production files or rewrite frozen manifests to conceal changes.

## Workflow documentation maintenance

`docs/workflow_dependencies.json` is the maintained, code-audited contract. Update it together with the matching workflow plan whenever inputs, producers, execution modes or outputs change. Then run:

```sh
python3 scripts/maintenance/build_workflow_documentation.py
python3 scripts/maintenance/build_workflow_documentation.py --check
```

This generates `docs/WORKFLOW_REFERENCE.md`, `docs/analysis_catalog.csv` and `docs/workflow_dependencies.csv`. Availability reflects tracked files and the local checkout, including the payloads of CSV transport manifests. The check detects missing entry-point documentation, broken source/plan paths, unknown registry keys or producer/call IDs, and stale generated files. It does not infer undocumented reads from source code, install packages or guarantee a full run; source review and per-workflow validation are still required. Git inclusion changes can change the availability inventory and require regeneration.
