# Workflow testing guide

The maintained testing entry is [tests/README.md](README.md). It lists checks that do not refit scientific models, cache-consuming checks, workflow-specific execution limits and the status of legacy runners.

For a change to a workflow:

1. Read its [workflow plan](../docs/workflow_plans/) and the [complete execution contract](../docs/WORKFLOW_REFERENCE.md). Identify inputs, outputs, parent caches and the actual behavior of validation flags.
2. Exercise the changed branch promptly with a meaningful targeted test. Stub fitting/rendering when testing dispatch or overwrite guards. For scientific logic changes, use an appropriate fresh run and compare estimates with frozen references.
3. Preserve approved inputs/results. Use explicit fresh output paths where supported and an isolated checkout/copy for older fixed-root entries. Do not source a script merely to inspect it unless it has a source guard.
4. Update the matching plan and `docs/workflow_dependencies.json`, regenerate documentation and run the checks in the maintained README.
5. Review input joins, matrix orientation, paired mouse identities, contrast order, denominators and multiple-testing families whenever these are affected. Successful file creation alone does not validate an analysis.

The former guide contained obsolete paths, counts and full-workflow execution examples. Its exact previous bytes are retained as [a historical snapshot](../miscellaneous/documentation/dependency_audit_20261118/tests_TESTING_GUIDE.md); they are not current execution instructions.
