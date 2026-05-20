# Workflow 08: Variance Analysis

**Script:** `scripts/workflows/08_variance_analysis.R`  
**Lines:** ~781  
**Purpose:** Compare lipid variability between ovaries and testes groups, separately from mean differential abundance tested in workflow 03.

## Source Dependencies

```r
source("scripts/config/analysis_config.R")
source("scripts/config/plotting_themes.R")
source("scripts/core/data_loading.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 41-58 | Load config, themes, data loaders, and create `results/lipidomic/variance_analysis/`. |
| 2 | Load Data | 61-94 | Load log2 lipidomics and sample metadata, subset metadata to lipidomic samples. |
| 3 | Separate Samples | 98-119 | Split columns into ovaries and testes sample sets. |
| 4 | CV Calculation | 123-178 | Compute per-lipid mean, SD, variance, CV, CV difference, and CV ratio. |
| 5 | Levene's Test | 181-241 | Brown-Forsythe Levene test per lipid; BH adjustment. |
| 6 | Global Paired Test | 246-289 | Paired Wilcoxon and paired t-test on lipid CVs. |
| 7 | Differential Variance Labels | 293-342 | Mark significant lipids and report top more-variable lipids by group. |
| 8 | Lipid-Level Plot | 349-414 | Volcano-style plot of CV difference versus Levene p-value. |
| 9 | Class-Level Analysis | 417-634 | Aggregate lipids to class means, repeat CV/Levene/global tests, and plot class volcano. |
| 10 | Save Results | 638-716 | Write result CSVs. |
| 11 | Final Summary | 720-781 | Print interpretation and output paths. |

## Input Files

- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv`
- `data/Suppl_Table_1_OverviewData.xlsx` via `load_metadata()`

## Output Files

- `results/lipidomic/variance_analysis/variance_analysis_results.csv` - Per-lipid CVs, Levene statistics, adjusted p-values, and direction labels.
- `results/lipidomic/variance_analysis/significant_differential_variance_lipids.csv` - FDR-significant subset from the per-lipid table.
- `results/lipidomic/variance_analysis/variance_analysis_by_class.csv` - Class-level variance analysis.
- `results/lipidomic/variance_analysis/variance_analysis_summary.csv` - Text-like metric/value summary table.
- `results/lipidomic/variance_analysis/cv_difference_volcano.pdf`
- `results/lipidomic/variance_analysis/cv_difference_volcano_by_class.pdf`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `lipidomic` | data.frame | Log2-transformed lipid abundance matrix |
| `metadata_lipid` | data.frame | Metadata for samples present in lipidomics |
| `ovaries_samples`, `testes_samples` | character | Sample names for the two gonadal groups |
| `cv_results` | data.frame | Per-lipid mean/SD/variance/CV and Levene results |
| `sig_var_lipids` | data.frame | FDR-significant per-lipid differential variance subset |
| `cv_class_results` | data.frame | Class-level CV and Levene results |
| `summary_stats` | data.frame | Output metric/value summary |

## Coherence Assessment

**Status: GOOD.** The workflow is standalone and only depends on raw lipidomics plus metadata helpers. The plan has been updated to current `variance_analysis_*` output names.

## Statistical Methods

1. **CV:** `SD / abs(mean) * 100`; the absolute mean is intentional because data are log2-transformed and means can be negative.
2. **Levene's test:** `car::leveneTest(..., center = "median")`, i.e. Brown-Forsythe variant.
3. **Multiple testing:** BH adjustment on per-lipid and per-class Levene p-values.
4. **Global tests:** Paired Wilcoxon and paired t-test comparing CVs across lipids/classes.

## Notes for Context Restoration

1. This workflow complements workflow 03: it tests variability, not mean abundance.
2. It does not consume workflow 03 outputs; it reloads lipidomics directly.
3. Volcano plot colors use FDR-significant Levene results and direction of `cv_testes - cv_ovaries`.
