# Workflow 07: Validation Experiments

**Script:** `scripts/workflows/07_validation_experiments.R`
**Lines:** ~672
**Purpose:** Validate gene-lipid correlations from MOFA and zonation analyses using siRNA knockdown experiments in primary hepatocytes.

## Source Dependencies

```r
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/core/differential_analysis.R")
source("scripts/config/analysis_config.R")
source("scripts/config/plotting_themes.R")
source("scripts/visualization/volcano_plots.R")
source("scripts/visualization/heatmaps.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 1-40 | Load config, set validation output dir |
| 2 | Load Knockdown Data | 41-120 | Load siRNA knockdown lipidomic data |
| 3 | QC PCA | 121-200 | PCA of knockdown samples, check clustering |
| 4 | DE Analysis | 201-300 | limma: knockdown vs control for each gene |
| 5 | Heatmap | 301-380 | Heatmap of knockdown effects across lipids |
| 6 | Global Network Validation | 381-470 | Compare knockdown FC with MOFA network correlations |
| 7 | Scatter Plots | 471-530 | Scatter: predicted (MOFA cor) vs observed (knockdown FC) |
| 8 | Zonation Network Validation | 531-590 | Same as 6-7 but using zonation-specific networks |
| 9 | Zonation Scatter Plots | 591-640 | Scatter plots for Sult3a1, Mtmr11 |
| 10 | Summary | 641-672 | Print validation statistics |

## Input Files

**Knockdown Data:**
- Knockdown lipidomic data (loaded via `data_loading.R` functions)
- Knockdown experiment metadata

**Network Files (from previous workflows):**
- `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` - Global MOFA network (from workflow 04)
- `results/zonation/networks/cormat_testes_central_full_SIF.csv` - Testes pericentral network (from workflow 05)
- `results/zonation/networks/cormat_testes_portal_full_SIF.csv` - Testes periportal network (from workflow 05)
- `results/zonation/networks/cormat_ova_central_full_SIF.csv` - Ovaries pericentral network (from workflow 05)
- `results/zonation/networks/cormat_Ovaries_portal_full_SIF.csv` - Ovaries periportal network (from workflow 05)

## Output Files

- `results/validation/knockdown_PCA.pdf` - PCA of knockdown samples
- `results/validation/knockdown_heatmap.pdf` - Knockdown effects heatmap
- `results/validation/cyp2d9_validation_scatter.pdf` - Cyp2d9 validation
- `results/validation/cyp2b9_validation_scatter.pdf` - Cyp2b9 validation
- `results/validation/cyp2c38_validation_scatter.pdf` - Cyp2c38 validation
- `results/validation/sult3a1_zonation_validation_scatter.pdf` - Sult3a1 zonation validation
- `results/validation/mtmr11_zonation_validation_scatter.pdf` - Mtmr11 zonation validation

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `knockdown_data` | data.frame | Knockdown lipidomic measurements |
| `knockdown_meta` | data.frame | Knockdown experiment metadata |
| `kd_ttop_list` | list | DE results per knockdown gene |
| `global_network` | data.frame | MOFA gene-lipid network edges |
| `zonation_networks` | list | 4 zone x gonad network data.frames |

## Genes Validated

| Gene | Network | Zone | Validation Type |
|------|---------|------|-----------------|
| Cyp2d9 | Global (MOFA) | - | Global correlation |
| Cyp2b9 | Global (MOFA) | - | Global correlation |
| Cyp2c38 | Global (MOFA) | - | Global correlation |
| Sult3a1 | Zonation | OPC | Zone-specific correlation |
| Mtmr11 | Zonation | TPP | Zone-specific correlation |

## Coherence Assessment

**Status: GOOD** — Clean validation pipeline, no open issues.

## Vulnerabilities to Keep in Mind When Modifying

1. **Network file naming is inconsistent but matched (section 6, 8):** The zonation network filenames use mixed conventions (`testes_portal`, `Ovaries_portal` with capital O, `ova_central` abbreviated). These match exactly what workflow 05 writes, so it works. But if you rename any network output in workflow 05, you must update the paths in this workflow too (lines 405-408).

## Notes for Context Restoration

When returning to this workflow:
1. This is a validation workflow - it TESTS predictions from workflows 04 and 05
2. Validation logic: if gene X correlates with lipid Y in the network, then knockdown of gene X should change lipid Y abundance
3. Expected: negative correlation between MOFA network weight and knockdown fold-change (gene removed -> lipid changes in opposite direction)
4. Global validation uses MOFA networks; zonation validation uses zone x gonad networks
5. Sult3a1 and Mtmr11 are the zonation-specific validation genes (not in MOFA top genes)
6. The knockdown data is from in vitro primary hepatocyte experiments
7. Also validates Cyp2b13 (not shown in scatter but included in heatmap)
