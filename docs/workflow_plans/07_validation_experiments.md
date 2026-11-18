# Workflow 07: siRNA Validation Experiments

## Authoritative execution contract (audited 18 November 2026)

Public entry: `scripts/workflows/07_validation_experiments.R` (23-line source-safe
wrapper). `source()` defines functions; `Rscript` invokes
`zonation_validation_run()`. The original full body is immutable at
`scripts/core/engines/07_validation_experiments_legacy.R` (`engine_validation07`),
not a maintained direct execution entry. There are two intentionally separate modes.

| Mode | What runs | Destination and protection |
| --- | --- | --- |
| Default/current zonal | `core/zonation_update_fcg.R::zmu_run07()`: frozen engine8/9/11 with global genes empty; reuse unchanged KD effects, adaptive zonal contexts plus prespecified Mtmr11 male-PP. | `--output=ROOT` defaults to `current_root`; writes `ROOT/fcg07/` and `ROOT/runtime/`. A completed `fcg07_checkpoint.rds` blocks the wrapper. |
| `--global-only` | `core/validation_global.R`: frozen engine1-7 and10-11; refit original KD data, global MOFA validation and global structural consistency; no zonation network read. | Actual execution requires explicit `--output=FRESH_DIR`; any nonempty directory is refused. Files are written directly under this directory. |

Supported flags are `--output=...`, `--global-only`, and `--validate-only`.
`--validate-only` alone checks package/checkpoint/KD-file prerequisites without
writing; it does not verify the structural reference later opened by the fit.
`--global-only --validate-only` verifies its raw KD/global/structural inputs and
prepared runtime syntax without writing. The checkpoint guard on zonal runs does
not promise to preserve a partial destination lacking its completion checkpoint.
Full zonal execution also calls `analysis_assert_not_approved_output()` before
checking/fitting, rejecting the approved root and descendants even when model
caches are absent. Validation is allowed there; actual fits need an explicit
fresh root outside it. Cache guards additionally protect completed fresh runs.

Workflow17 coordinates current05 then zonal07. The independent in-vivo Cyp2c38
Workflow11 is a different experiment and is not called by this entry; former07b/07c
are phase compatibility wrappers for11.

## Actual inputs and ordering by mode

| Input | Zonal mode | Global-only mode | Producer |
| --- | --- | --- | --- |
| `<ROOT>/fcg05/fcg05_checkpoint.rds` and four `fcg05/networks/cormat_*_full_SIF.csv` | Required current reference; the exact mixed names are below. | Not read. | Current05 in the same run root. |
| `results/validation/knockdown_t_table.csv` | Reused directly for scoring and plots. | Not an input; a new table is fitted/written in the fresh output. | Retained original siRNA result. |
| `results/validation/knockdown_ttop_list.RData` | Entry-level prerequisite; the score calculations read the CSV. | Not an input; newly fitted. | Retained original siRNA fit. |
| `data/metabolomics/lipidomics_20250423transfections.csv` | Not read. | Required original KD measurements. | Curated experimental input. |
| `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | Not read. | Required permissive global gene-lipid reference. | Workflow04. |
| `results/lipidomic/structural_model/lipids_with_structure.csv` | Required for structural reference axes. | Required for structural reference axes. | Workflow06 after03. |

The curated KD/global/structural resources are public-tree artifacts with the
necessary R packages; neither mode requires a source RNA/Seurat object. Zonal
validation also requires its generated same-root05 checkpoint, which is an ignored
model cache. A public clone must first fit05/17 into a fresh root (or explicitly
restore that cache); plain07 is not a way to recreate it. Recompute04/06 first only when
rebuilding their global network/structural inputs. A new global-only run writes
its KD fit to its fresh output; it does not silently replace the curated canonical
`results/validation/` inputs reused by current zonal validation.

Exact current full-network basenames are
`cormat_testes_portal_full_SIF.csv`, `cormat_testes_central_full_SIF.csv`,
`cormat_Ovaries_portal_full_SIF.csv`, `cormat_ova_central_full_SIF.csv`.
All are below the selected root's `fcg05/networks/`. Old `results/zonation/`
paths occur only inside the frozen engine and are redirected by the zonal adapter.

## Implementation and outputs

Both adapters source the registry and shared FCG adapter, whose `zmu_require()`
checks the full declared package set. Generated runtimes additionally source
`analysis_config.R`, `datasets.R`, `support_functions.R`, and structural helpers.
The global adapter removes workspace clearing, binds an empty zonation candidate
set to avoid inherited interactive state, redirects only validation output paths,
and supplies the t-based mean-CI plotting helper. No statistical model changes
are introduced by those transport/execution adapters.

Zonal outputs under `<ROOT>/fcg07/` include:

- `zonation_network_validation_scores.csv`, on-target/adaptive context results in
  `fcg07_checkpoint.rds`, `warnings.csv`, and zonal validation scatter PDFs.
- `Mtmr11_male_PP_prespecified_ULM.csv`,
  `Mtmr11_male_PP_prespecified_source.csv`,
  `Mtmr11_male_PP_ULM_background_source.csv`,
  `Fig5F_Mtmr11_male_PP_prespecified.pdf`.
- `structural_consistency/[gene]_structural_consistency.pdf` and
  `structural_consistency/consistency_structural_regression_summary.csv` for zonal
  candidates. Runtime: `runtime/07_zonation_validation_versioned.R` at root level.

Global-only outputs directly under its explicit fresh directory include:

- `knockdown_PCA.pdf`, `knockdown_ttop_list.RData`, `knockdown_t_table.csv`,
  `knockdown_heatmap.pdf`, `global_network_validation_scores.csv`.
- Global scatter PDFs for available genes among Cyp2d9/Cyp2b9/Cyp2c38;
  `validation_summary.csv`, global `structural_consistency/` plots/regressions,
  `warnings.csv`, `runtime/07_global_validation_versioned.R` and its section ledger.

Additional old/new comparison tables and unchanged-KD-fit audits belong to the
retrospective FCG report driver/Workflow17 audit mode, not to the plain07 wrapper.

## Method and interpretation contracts

Global candidates are Cyp2d9/Cyp2b9/Cyp2c38; zonal adaptive candidates are
Sult3a1/Mtmr11 when available. Adaptive best-context scores remain distinct from
prespecified Mtmr11 male-PP validation. Structural consistency is correlation times
KD t-statistic times minus1, with SM/BMP/PG/TG focus classes. The ULM background
and measured-scatter feature count must remain distinct; zones are paired RNA
profiles and lipids are measured at whole-liver level. Preserve lowercase gene
matching and the exact network filenames. The original parenthesis-based target
filter is historical scoring behavior, not a general lipid-name identity rule.

The base-R mean-CI helper is a plotting-dependency replacement for Hmisc's normal-
theory t interval; it does not change ULM or regression calculations. Actual
independent global reproduction previously matched all KD fits, global ULM values
and global structural regressions at1e-12 after matching tied lipid rows by ID.

## Historical engine map and production record

The frozen engine's original1113-line map was:1 configuration,2 KD import,3 PCA,
4 limma KD fitting,5 heatmap,6 global ULM,7 global scatters,8 zonal ULM,9 zonal
scatters,10 combined summary,11 structural consistency. Its original fixed
`results/validation/` outputs and automatic combined execution are historical;
the maintained mode table above defines current behavior.

## Versioned zonation-reference update, 16 November 2026

`scripts/reports/run_zonation_update_fcg.R` and its isolated
`scripts/core/zonation_update_fcg.R` adapter rerun only the zonation-dependent
validation and structural-consistency sections against the revised FCG WF05
networks. Outputs are under `results/zonation_marker_update_20261116/fcg07/`;
historical validation outputs are preserved. This records the production
comparison; current defaults are the zonal wrapper mode described above. The global
MOFA validation and knockdown fit are unchanged. The saved knockdown t-table
is reused directly and independently checked against all five saved limma
toptables to numerical precision.

The historical adaptive best-context plots/scores are preserved as a separate
versioned branch. An additional prespecified **Mtmr11 male-periportal** result
and panel prevent an adaptive context choice from changing the manuscript's
stated validation target. All four historical Figure 5D genes continue to meet
the revised gene-selection rule, including Mtmr11.

Native `decoupleR::run_ulm()` uses the full 661-feature knockdown matrix and
zero prediction weights for features absent from the reference network. The
scatter contains the 345 matched lipids with measured in-vivo correlations;
it must not imply that all 661 features have measured correlation zero. A
separate background ledger distinguishes unmatched zero ULM weights from
measured correlations. Preserving this background exactly reproduces the
archived Mtmr11 score -7.211103 (`p=1.527607e-12`); the revised reference gives
-6.791136 (`p=2.487504e-11`). Both references use 22 male animals, not the
pooled FCG sample size of 43. Matched-only ULM is a different statistic and is
not substituted for the historical method.

Comparison tables include native scores/p-values, explicitly labelled
on-target BH summaries, matched correlation vectors, selected-context
structural regression summaries and the unchanged-fit audit. The original
parenthesis-based WF07 target filter is retained for faithful replay; the
three measured FCG features without parentheses do not overlap the knockdown
data, so their omission does not change these validation scores.

The legacy `stat_summary(mean_cl_normal)` plotting layer has a local base-R
fallback because Hmisc was unavailable from the active R package repositories.
It uses the documented normal-theory t interval (mean plus/minus the
t-quantile times its standard error) and was checked against `stats::t.test`.
This is a rendering dependency repair; ULM and regression calculations are
unchanged. Final warning logs contain no failed statistical plotting layers.
