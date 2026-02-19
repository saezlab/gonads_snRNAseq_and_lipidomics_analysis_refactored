# Workflow 06: Structural Analysis

**Script:** `scripts/workflows/06_structural_analysis.R`
**Lines:** ~1685
**Purpose:** Analyze how lipid structural features (chain length, double bond position) modulate gene-lipid correlations. Regression models test whether gene-lipid associations depend on structural properties.

## Source Dependencies

```r
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/core/structural_analysis.R")
source("scripts/config/analysis_config.R")
source("scripts/config/plotting_themes.R")
source("scripts/visualization/heatmaps.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 1-50 | Load config, set focus genes/classes |
| 2 | MOFA Weights | 51-150 | Load weights, convert to long format |
| 3 | RNA + Lipidomic Data | 151-250 | Load expression and lipidomic matrices |
| 4 | Lipid Structural Features | 251-400 | Load chain length, DB position from `for_bubbles_with_names/` |
| 5 | Lipid-Only Analysis | 401-550 | Structural analysis WITHOUT gene correlations |
| 5B | Lipid-Only Visualization | 551-650 | Heatmap of lipid-only structural effects |
| 6 | Lipid Metabolism Gene Set | 651-720 | Filter to GO:CELLULAR_LIPID_METABOLIC_PROCESS genes |
| 7 | Gene-Lipid Correlations | 721-850 | Pearson correlations (hepatocyte genes x lipids) |
| 8 | Global Regression | 851-1000 | Per-gene: `cor ~ class + abundance + chain + DB + interactions` |
| 9 | Class-Specific Regression | 1001-1150 | Same model but within each lipid class (min 10 lipids) |
| 10 | Identify Key Genes | 1151-1250 | Genes with chain length, DB, or interaction effects |
| 11 | Visualization | 1251-1450 | 9-panel plots per gene per class (correlation vs structural features) |
| 12 | Summary Statistics | 1451-1520 | Print counts of significant associations |
| 13 | Biological Interpretation | 1521-1600 | Categorize genes by structural preference type |
| 14 | Summary Figure Panels | 1601-1685 | Panels A (overview heatmap), B (example gene), C (class summary) |

## Input Files

- `results/mofa/main_model_RNA_lipid_weight.RData` - MOFA weights (from workflow 04)
- `data/RNA/scaled_RNA_{1,2,3}.csv` - Expression data
- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` - Lipidomic data
- `data/metabolomics/for_bubbles_with_names/*.csv` - Structural data per lipid class
- `support/m5.all.v2023.2.Mm.symbols.gmt` - GO gene sets (for lipid metabolism filter)

## Output Files

**Data:**
- `results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv` - Full correlation matrix
- `results/multi-omic/structural_model/lipids_with_structure.csv` - Lipid structural features
- `results/multi-omic/structural_model/global_regression_results.csv` - Global model results
- `results/multi-omic/structural_model/class_specific_regression_results.csv` - Class-specific results
- `results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv` - All significant
- `results/multi-omic/structural_model/genes_associated_with_chain_length.csv`
- `results/multi-omic/structural_model/genes_associated_with_db_position.csv`
- `results/multi-omic/structural_model/genes_with_significant_interactions.csv`
- `results/multi-omic/structural_model/genes_multiple_structural_associations.csv`
- `results/multi-omic/structural_model/analysis_summary.csv`

**Lipid-Only Results:**
- `results/multi-omic/structural_model/lipid_only_class_summaries.csv`
- `results/multi-omic/structural_model/lipid_only_class_coefficients.csv`
- `results/multi-omic/structural_model/lipid_only_chain_length_effects.csv`
- `results/multi-omic/structural_model/lipid_only_db_position_effects.csv`
- `results/multi-omic/structural_model/lipid_only_interaction_effects.csv`
- `results/multi-omic/structural_model/lipid_only_summary_heatmap.pdf`
- `results/multi-omic/structural_model/lipid_only_plots/*.pdf`

**Figures:**
- `results/multi-omic/structural_model/gene_plots/{gene}_lipid_structure.pdf` - Per-gene plots (4 genes: Cyp2b13, Cyp2b9, Cyp2c38, Cyp2d9)
- `results/multi-omic/structural_model/summary_figure_panelA.pdf`
- `results/multi-omic/structural_model/summary_figure_panelB.pdf`
- `results/multi-omic/structural_model/summary_figure_panelC.pdf`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `weights` | data.frame | MOFA weights (long format after section 2 conversion) |
| `cor_mat` | matrix | Gene-lipid Pearson correlation matrix |
| `lipids_with_structure` | data.frame | Lipid name, class, chain_length, db_position, logFC |
| `global_results` | data.frame | Regression results across all lipid classes |
| `class_results` | data.frame | Class-specific regression results |
| `focus_genes` | character | `c("Cyp2d9", "Cyp2b13", "Cyp2c38", "Cyp2b9")` |
| `focus_classes` | character | `c("SM", "BMP", "PE", "TG")` |

## Coherence Assessment

**Status: GOOD** — No open bugs. Weights format conversion is intentional and correct.

## Vulnerabilities to Keep in Mind When Modifying

1. **Focus genes are hardcoded (section 1):** `focus_genes = c("Cyp2d9", "Cyp2b13", "Cyp2c38", "Cyp2b9")` and `focus_classes = c("SM", "BMP", "PE", "TG")`. These are the validated Cyp2 family members and are correct. If upstream MOFA or validation analysis identifies different key genes, update here manually.

2. **Structural data file matching (section 4):** Reads one CSV per lipid class from `data/metabolomics/for_bubbles_with_names/` (e.g., `PC.csv`, `TG.csv`). If lipid class naming changes or new classes are added, the file list needs updating.

## Regression Model

```r
# Global model (section 8):
gene_correlation ~ lipid_class + lipid_abundance + chain_length + db_position +
                   lipid_abundance:chain_length + lipid_abundance:db_position

# Class-specific model (section 9, per lipid class with >= 10 lipids):
gene_correlation ~ lipid_abundance + chain_length + db_position +
                   lipid_abundance:chain_length + lipid_abundance:db_position
```

## Notes for Context Restoration

When returning to this workflow:
1. Depends on MOFA weights from workflow 04 (`main_model_RNA_lipid_weight.RData`)
2. The lipid structural features come from `data/metabolomics/for_bubbles_with_names/` - one CSV per class
3. The regression asks: "For gene X, does its correlation with lipids depend on chain length or double bonds?"
4. Section 5 (lipid-only) is a NEW addition - analyzes structural effects WITHOUT gene involvement
5. Visualization: green regression line = significant (p < 0.05), black = non-significant
6. Panels A/B/C in section 14 are designed for the manuscript figure
7. The 4 focus genes are the Cyp2 family members validated by knockdown experiments
