# Workflow 07: Validation Experiments

## Current and independent execution entries - 16 November 2026

The source-safe Workflow07 entry defaults to current expanded-marker zonation
validation. It requires the selected output root's `fcg05/fcg05_checkpoint.rds`
and the unchanged siRNA KD table/fit, writes `fcg07/`, and refuses to overwrite
completed validation. `--validate-only` checks that zonation entry's requirements.

Independent KD/global/structural reproduction remains available explicitly:
`Rscript scripts/workflows/07_validation_experiments.R --global-only --output=results/runs/global_validation_reproduction`.
This mode requires a fresh directory, refits the original KD data, and executes
unchanged engine sections1-7 and10-11 without zonation networks. Its own
`--global-only --validate-only` checks the independent input scope. The adapter
is `scripts/core/validation_global.R`; output redirection and a t-based plot CI
helper are the only transformations beyond section selection. Actual execution
reproduced every KD fit, global ULM value and global structural regression at
1e-12 tolerance after matching tied lipid rows by ID. Overwrite guards were tested.

The original full body is frozen byte-for-byte in
`scripts/core/engines/07_validation_experiments_legacy.R`. The historical section
map below refers to that engine. Current zonation execution separately preserves
Mtmr11 male-periportal prespecified validation and the adaptive-context result;
neither changes the independent KD measurements or global MOFA networks.

---

The October 2026 meeting's heterozygous **in vivo Cyp2c38 KO** WP1 is an
independent extension, documented in `11_cyp2c38_in_vivo_validation.md` (with the original 07b/07c
phase plans retained) and implemented in
`scripts/workflows/11_cyp2c38_in_vivo_validation.R`. This
workflow retains the original primary-hepatocyte siRNA analyses.

**Script:** `scripts/workflows/07_validation_experiments.R`  
**Lines:** ~1113  
**Purpose:** Validate MOFA and zonation gene-lipid correlation predictions against siRNA knockdown lipidomics, then test whether knockdown-response consistency varies with lipid structure.

## Source Dependencies

```r
source("scripts/config/analysis_config.R")
source("scripts/config/datasets.R")
source("scripts/support_functions.R")

# Section 11 only:
source("scripts/core/structural_analysis.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 51-69 | Load config/datasets/support functions and create `results/validation/`. |
| 2 | Load Knockdown Data | 72-103 | Read `lipidomics_20250423transfections.csv`, clean lipid data, infer condition from sample names. |
| 3 | QC PCA | 106-122 | PCA of knockdown lipidomics. |
| 4 | Differential Analysis | 125-190 | Build every knockdown-vs-control contrast, run limma, write ttop list and t table. |
| 5 | Heatmap | 193-222 | Heatmap of lipid t-statistics across knockdowns. |
| 6 | Global Network Validation | 227-274 | Load MOFA full gene-lipid network and score prediction/knockdown consistency with decoupleR ULM. |
| 7 | Global Scatter Plots | 279-398 | Generate per-gene predicted-vs-observed scatter plots for available global genes. |
| 8 | Zonation Network Validation | 404-481 | Load four workflow 05 full zonation networks, combine them, and score on-target validation. |
| 9 | Zonation Scatter Plots | 483-620 | Generate zonation-specific predicted-vs-observed scatter plots for available zonation genes. |
| 10 | Summary Statistics | 628-660 | Write validation summary table. |
| 11 | Structural Consistency | 665-1079 | Parse lipid structures, regress KO-response consistency on `C`, `DB`, and `C:DB`, and write plots/table. |
| - | Complete | 1081-1113 | Print output summary. |

## Input Files

**Knockdown data**
- `data/metabolomics/lipidomics_20250423transfections.csv`

**Prediction networks**
- `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` - Global MOFA-derived gene-lipid network from workflow 04.
- `results/zonation/networks/cormat_testes_portal_full_SIF.csv`
- `results/zonation/networks/cormat_testes_central_full_SIF.csv`
- `results/zonation/networks/cormat_Ovaries_portal_full_SIF.csv`
- `results/zonation/networks/cormat_ova_central_full_SIF.csv`

**Structural reference**
- `results/lipidomic/structural_model/lipids_with_structure.csv` - Workflow 06 output used for shared axis limits in section 11.

## Output Files

- `results/validation/knockdown_PCA.pdf`
- `results/validation/knockdown_ttop_list.RData`
- `results/validation/knockdown_t_table.csv`
- `results/validation/knockdown_heatmap.pdf`
- `results/validation/global_network_validation_scores.csv`
- `results/validation/[gene]_validation_scatter.pdf` - Global network scatter plots for available genes among `cyp2d9`, `cyp2b9`, `cyp2c38`.
- `results/validation/zonation_network_validation_scores.csv`
- `results/validation/[gene]_zonation_validation_scatter.pdf` - Zonation scatter plots for available genes among `sult3a1`, `mtmr11`.
- `results/validation/validation_summary.csv`
- `results/validation/structural_consistency/[gene]_structural_consistency.pdf`
- `results/validation/structural_consistency/consistency_structural_regression_summary.csv`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `lipidomic_raw` | data.frame | Raw knockdown lipidomics file |
| `lipidomic` | data.frame | Cleaned lipid abundance matrix, lipids x samples |
| `targets` | data.frame | Sample and inferred condition metadata |
| `comparisons` | named list | All knockdown conditions compared against detected control (`neg`, `control`, or `scramble`) |
| `ttop_list` | list | Limma top tables by knockdown |
| `t_table` | data.frame | Lipid t-statistic table used for validation |
| `lipid_network` | data.frame | Global MOFA gene-lipid network with `source`, `target`, `mor` |
| `sszb_network` | data.frame | Combined zonation-specific network |
| `focus_classes` | character | `c("SM", "BMP", "PG", "TG")` for structural consistency |

## Genes Validated

| Gene set | Hardcoded candidates | Notes |
|----------|----------------------|-------|
| Global MOFA network | `cyp2d9`, `cyp2b9`, `cyp2c38` | Only genes present in `names(comparisons)` are plotted/scored as available. |
| Zonation network | `sult3a1`, `mtmr11` | Best gonad-zone context is chosen by highest absolute decoupleR score. |
| Heatmap/limma table | All knockdowns in the 2025-04-23 file | Additional knockdowns can be present even if they do not get dedicated scatter plots. |

## Coherence Assessment

**Status: GOOD.** The workflow is linear. It conditionally skips zonation validation if workflow 05 full networks are missing, but global validation expects the workflow 04 full network to exist.

## Vulnerabilities to Keep in Mind When Modifying

1. **Network filename coupling:** WF07 expects the exact mixed WF05 filenames listed above, including `Ovaries` capitalization and `ova` abbreviation.
2. **Control detection:** Section 4 searches for `neg`, `control`, or `scramble`; otherwise the first condition becomes the reference.
3. **Structural dependency:** Section 11 requires WF06's `lipids_with_structure.csv`; a fresh run of WF07 after deleting structural outputs will fail there.
4. **Case normalization:** Validation filters compare lowercased gene names. Preserve this when adding genes with mixed-case mouse symbols.

## Notes for Context Restoration

1. Validation consistency is computed as predicted correlation times knockdown t-statistic times `-1`; positive consistency means the knockdown direction agrees with the FCGL prediction.
2. Dedicated scatter plots are not automatically made for every knockdown; the hardcoded global and zonation candidate lists control those plots.
3. Structural consistency plots use the same focus classes and axis-range logic as workflow 06.

## Versioned zonation-reference update, 16 November 2026

`scripts/reports/run_zonation_update_fcg.R` and its isolated
`scripts/core/zonation_update_fcg.R` adapter rerun only the zonation-dependent
validation and structural-consistency sections against the revised FCG WF05
networks. Outputs are under `results/zonation_marker_update_20261116/fcg07/`;
historical validation outputs and workflow defaults are preserved. The global
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
