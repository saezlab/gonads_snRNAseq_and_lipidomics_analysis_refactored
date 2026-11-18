# Historical and exploratory material

This archive preserves earlier work while the main README and release registry identify the current pipeline. Archived does not mean unimportant: some earlier files are frozen dependencies of explicitly historical comparisons or causal-model replay.

| Folder | Contents |
| --- | --- |
| `older_results/` | Original two-marker zonation, older Figure6 presentations, old MOFA generations, earlier structural/validation and causal runs |
| `older_inputs/` | Original zonation inputs and legacy RNA/culture data; explicit optional dependencies are documented |
| `comparisons/` | Original WDSD/CF459 cross-cohort and earlier reviewer-response analyses |
| `older_code/` | Original unnumbered scripts and previously unused code |
| `refactoring/` | Refactor/debug history |
| `quick_tests/` | Scratch code, temporary exports and old test reports |
| `logs/` | Historical workflow and package logs |
| `documentation/` | Previous repository instructions and pipeline documentation |
| `local/` | Machine settings and generated cleanup/readiness checks; ignored by Git |

Find exact dependencies and original paths in [the organization map](../docs/REPOSITORY_ORGANIZATION.md), [the release registry](../scripts/config/analysis_release.json), [the migration ledger](../docs/provenance/repository_migration_20261116.csv) and [the causal map](../docs/CAUSAL_ANALYSIS_MAP.md).

Large archive payloads stay local and are excluded from ordinary commits. Their location, size and SHA256 remain versioned. Do not replace a missing archived reference with the current expanded-marker data; restore the original bytes instead. Run archived notebook-style code only after reviewing its historical input/output assumptions.
