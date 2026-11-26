# Testing the maintained analyses

Run checks from the repository root. Test names alone do not establish that a command is read-only: some legacy runners execute full workflows and write their normal result paths. Read the matching [workflow plan](../docs/workflow_plans/) and [execution contract](../docs/WORKFLOW_REFERENCE.md) before running a workflow or numerical QA script.

## Documentation and release checks

These commands do not fit scientific models or regenerate panels. Repository audits write inventories under ignored `miscellaneous/local/`; the guard fixture creates and removes a temporary miniature repository, with model/render/subprocess calls stubbed.

```sh
python3 scripts/maintenance/build_workflow_documentation.py --check
python3 tests/test_workflow_documentation.py
python3 tests/test_github_readiness_links.py
Rscript tests/test_approved_output_guard.R
python3 scripts/maintenance/verify_repository_organization.py
python3 scripts/maintenance/audit_github_readiness.py --strict
git diff --check
```

The documentation check covers entry-point documentation, declared source/plan paths, registry keys and generated-file consistency. It does not infer every read from code or certify full scientific reproduction. Update `docs/workflow_dependencies.json` and the matching plan together, then run the documentation builder without `--check` to regenerate its three outputs.

Organization verification compares inputs and approved scientific artifacts with preserved hashes. Add `--require-full-archive` when checking the complete local archive; a public clone intentionally lacks some external historical payloads. Strict GitHub readiness checks prospective files, RNA object exclusions, CSV transport completeness/values and size policy. Neither audit stages, commits or pushes.

## Targeted execution and saved-result checks

| Check | Scope and prerequisites |
| --- | --- |
| `Rscript scripts/workflows/05_zonation_analysis.R --validate-only` and `--cohort=gdx --validate-only` | Check the selected expanded-marker CSV input. FCG also checks the temporary engine inputs and runtime syntax. |
| `Rscript scripts/workflows/17_expanded_marker_zonation_update.R --validate-only` | Reads/verifies the three current zone CSV transports, engine files and saved KD files; no fitting or writes. It does not check every full-run03/06/resource dependency. |
| `Rscript tests/test_analysis_release.R` | Checks registry, frozen engines, relocation and raw-input identity for the saved current diet model cache; requires ignored `cf453_models/models.rds` and adjacent provenance. Optional `--with-historical` also needs restored historical caches. |
| `Rscript tests/test_validation_global_entry.R` | With no extra argument, checks the07 global entry's prerequisites and generated runtime text; no fits. Its optional fresh-output argument compares an already-generated independent validation rerun; it does not create that rerun. |
| `tests/qa_*.R` | Workflow-specific numerical/export checks. Some replay calculations, render files or compare exact historical code/plan hashes. Inspect the script and plan before execution. |
| `tests/test_*.R` | Mixed unit/integration tests with different input/runtime requirements. Run the tests appropriate to the changed code, after reading them. |

Workflow11 `--validate-only` fits its observed weighted limma model in memory. Workflow16 `--models-only` performs its full numerical/resampling stage. Workflows14/15/16 have no `--validate-only` flag. The older01–04/06/08 entries execute on source and have fixed output roots; `--help` does not make them safe to inspect by execution.

## Scientific reproductions

Use the fresh-output commands and dependency order in [PIPELINE.md](../PIPELINE.md).
Run FCG 05 before zonal 07, or choose the paired GDX mode in 05. Diet-only
runs pass the new 14 and 15 cache paths explicitly to 16. Workflow 17 offers
batch reproduction. Use isolated checkouts for older fixed-root workflows.

Historical source-snapshot QA can flag intentional code/plan changes even when preserved numerical files and normal parent-cache checks pass. Keep the original recorded hashes and interpret the stated mismatch; do not refresh frozen provenance to hide drift.

## Legacy runners and documentation

`run_tests.R`, `run_workflow_tests.R`, `test_all_workflows.R`, `test_config.csv` and older examples were created for earlier repository structures. They are not the default release-validation commands and can execute workflows or rely on stale paths/reference counts. Review and adapt them in an isolated reproduction before use; do not run them indiscriminately in the approved checkout.

The previous README and workflow testing guide are byte-preserved in [the documentation audit snapshots](../miscellaneous/documentation/dependency_audit_20261118/README.md). Historical generated test logs are under `miscellaneous/quick_tests/logs/`; new local audit logs remain under `miscellaneous/local/`.
