# Workflow 08: Variance Analysis

**Script:** `scripts/workflows/08_variance_analysis.R`
**Lines:** ~782
**Purpose:** Compare lipid variance (not mean) between ovaries and testes groups. Tests whether gonadal hormones affect lipid variability in addition to abundance.

## Source Dependencies

```r
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/config/analysis_config.R")
source("scripts/config/plotting_themes.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 1-40 | Load config, set output dir |
| 2 | Data Loading | 41-100 | Load lipidomic data, log2 transform |
| 3 | Sample Separation | 101-150 | Split into ovaries vs testes sample groups |
| 4 | CV Calculation | 151-250 | Coefficient of variation per lipid per group |
| 5 | Levene's Test | 251-350 | Test equality of variances (lipid by lipid) |
| 6 | Global Paired Test | 351-420 | Paired Wilcoxon: CVs ovaries vs testes overall |
| 7 | Identify Differential Variance | 421-500 | Filter significant Levene's test results |
| 8 | Visualizations | 501-600 | CV scatter, density plots, volcano-style plot |
| 9 | Class-Level Analysis | 601-700 | Aggregate variance differences by lipid class |
| 10 | Save Results | 701-750 | Write CSVs and RData |
| 11 | Summary | 751-782 | Print summary statistics |

## Input Files

- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` - Lipidomic data

## Output Files

- `results/lipidomic/variance_analysis/cv_by_group.csv` - CV per lipid per group
- `results/lipidomic/variance_analysis/levene_test_results.csv` - Levene's test results
- `results/lipidomic/variance_analysis/differential_variance_lipids.csv` - Significant lipids
- `results/lipidomic/variance_analysis/class_level_variance.csv` - Class-level summary
- `results/lipidomic/variance_analysis/*.pdf` - Visualization plots

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `lipidomic` | data.frame | Log2-transformed lipid abundances |
| `cv_df` | data.frame | CV values per lipid, split by gonad group |
| `levene_results` | data.frame | Levene's test p-values per lipid |
| `diff_var_lipids` | data.frame | Lipids with significantly different variance |

## Coherence Assessment

**Status: GOOD** - Self-contained, no cross-workflow dependencies.

- Clean linear flow: load -> compute -> test -> visualize -> save
- No dependencies on other workflow outputs
- Only needs the raw lipidomic data file
- Levene's test is appropriate for testing variance equality

## Statistical Methods

1. **Coefficient of Variation (CV):** `CV = SD / mean` per lipid per gonad group
2. **Levene's Test:** Non-parametric test for equality of variances (robust to non-normality)
3. **Paired Wilcoxon:** Global test comparing CVs across all lipids (ovaries vs testes)
4. **Multiple testing correction:** BH (Benjamini-Hochberg) on Levene's p-values

## Notes for Context Restoration

When returning to this workflow:
1. This is a NEW workflow (untracked in git) - independent of the main analysis chain
2. Tests a different hypothesis: are lipids more VARIABLE in one gonad group vs the other?
3. This is complementary to workflow 03 which tests mean differences
4. Self-contained: only needs raw lipidomic CSV, no dependencies on other workflow outputs
5. Results go to `results/lipidomic/variance_analysis/` subdirectory
