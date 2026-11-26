# Repository guidance for analysis agents

## Read first

- `README.md`: maintained overview and current default.
- `PIPELINE.md`: workflow order and dependency graph.
- `docs/AGENT_GUIDE.md`: static result lookup map.
- `docs/REPOSITORY_ORGANIZATION.md`: archive, provenance and restoration.
- `docs/workflow_plans/_meta_cross_workflow_overview.md`: cross-workflow details.
- `docs/workflow_plans/NN_*.md`: mandatory before changing the matching workflow.

The previous long instruction file is preserved at `miscellaneous/documentation/pre_cleanup_20261116/AGENTS.md`. Its old two-marker zonation defaults, unnumbered script paths and old result counts are historical; use generated results for quantitative claims.

## Current scientific defaults

This repository analyses gonadal versus chromosome effects in Four Core Genotypes mouse liver, with RNA/lipid integration and independent validation. The current zonation and standalone diet release is **expanded_markers_20261116**, selected in `scripts/config/analysis_release.json`.

- Main zonation entry: `scripts/workflows/05_zonation_analysis.R` for FCG and `--cohort=gdx` for paired gonadectomy. Zonal Workflow 07 uses the FCG checkpoint and its adjacent networks. Workflow 17 is an optional batch and audit entry.
- Public pseudobulk inputs: `data/RNA/pseudobulk_csv/`, using the seven manifest paths in the release registry. Each RNA CSV part must be strictly smaller than 50,000,000 bytes.
- Original RNA/Seurat RDS objects remain local and must not be allowed into the commit; readers prefer verified CSV manifests. Preserve source-object SHA256 separately from CSV payload hashes.
- Current result root: `results/current/`. Consult `analysis_current(component)` for release components and the relocation ledger for frozen original paths.
- Standalone diet analysis: Workflow 14 reads the expanded-marker diet transport directly and feeds 15 and 16.
- Current manuscript Figure 6B and 6E: `scripts/reports/cf453_manuscript_panels.R` uses saved Workflow15 lipid contrasts and Workflow16 gene–lipid pairs. Workflow17 renders both panels. The anchor layouts retain their native Workflow16 filenames.
- Old 05/07 bodies are implementation engines under `scripts/core/engines/`; use the maintained public entries/adapters, not an unreviewed legacy source call.
- Historical two-marker inputs/results and 12/13 comparisons resolve through explicit historical registry keys. Never point those keys at current data.
- Main whole-cell RNA, lipid, MOFA, structural and independent validation results retain their documented roots. Workflow11 is the independent Cyp2c38 in-vivo branch.
- COSMOS/MOON is an active optional extension; read `docs/CAUSAL_ANALYSIS_MAP.md`. Historical pilot folders can be live frozen dependencies.

Default scientific runs do not require the private manuscript or a Codex-specific runtime. Manuscript-audit/report replay is explicit opt-in. Full reruns require fresh output directories. Preserve the approved release and historical artifacts.

## Mandatory workflow change protocol

1. Read the matching `docs/workflow_plans/NN_*.md` before modifying a workflow. Read the meta overview when shared data/dependencies change.
2. Run each modified R script immediately to catch errors. Use a meaningful validate-only/source-guard or targeted execution when appropriate; exercise changed execution branches and compare estimates to frozen references when logic changes. Do not defer all testing until many R changes accumulate.
3. Use explicit namespaces for tidyverse functions: `dplyr::select`, `dplyr::filter`, `dplyr::mutate`, `dplyr::arrange`, `tidyr::pivot_longer`, `tidyr::pivot_wider`. Avoid MASS/stats conflicts.
4. Check matrix/data-frame format, feature-by-sample orientation, row names, sample identifiers/suffixes, paired mouse identities and contrast order. Print classes, dimensions and a short ID preview while debugging.
5. Update the corresponding workflow plan after changes, including outputs, parameters, dependencies and vulnerabilities.
6. Use `scripts/support_functions.R` for established shared utilities and preserve RNA split-file loading where applicable. Free large temporary matrices.

Common pitfalls: MOFA helpers can return data frames rather than lists; decoupleR expects features × samples; sample names can contain cell-type/gonad suffixes; limma numeric contrasts depend on actual level order. Remove interactive `View()` calls from batch execution. Check matching and dimensions first when errors occur.

## Scientific and provenance requirements

- Preserve the distinctions between effect size, t statistic, FDR family, animal count, lipid count and edge count. Correlated edges are not independent animals.
- Zones are paired RNA profiles per mouse; lipid measurements are whole-liver.
- The standard/Western-diet cohort is a standalone analysis; age differences prevent simple causal bridging to the original cohort.
- Current diet association/network panels are exploratory. Do not convert threshold-defined gains/losses or intervals including zero into confirmed differential correlations.
- Figure4A/S4B–C and associated nucleus counts/composition tests remain James's tasks. Pseudobulk expression cannot reconstruct per-nucleus composition QC.
- Supplied marker candidates/scoring scripts do not establish the exact cohort-specific executed retained marker lists. Keep that Methods qualification until James confirms.
- Never fabricate counts, citations, coefficients or mechanistic interpretations. Use exact saved tables and primary evidence.
- Preserve gonad/zone palette conventions from the plotting config: ovaries red, testes blue; darker PC, lighter PP.

## Archive and Git hygiene

`miscellaneous/` holds older results, comparisons, original code, quick tests and logs. Some payloads are excluded from Git but remain on disk. Their per-file old/new paths and SHA256 are in `docs/provenance/`.

Do not rewrite a frozen embedded provenance path/hash to hide a move. Resolve recorded paths with `analysis_resolve_recorded_path()` (R) or `resolve_recorded_path()` (Python). Add a relocation record for any subsequent move, and verify bytes.

Before handing off changes, run `python3 scripts/maintenance/verify_repository_organization.py` and `python3 scripts/maintenance/audit_github_readiness.py --strict`, plus tests appropriate to code changes. No automatic staging, commit, push, publication or destructive cleanup is implied. Keep private manuscript source, local histories/settings, caches, oversized downloads and generated bundles out of ordinary commits.
