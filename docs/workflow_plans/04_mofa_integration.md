# Workflow 04: MOFA Integration

**Script:** `scripts/workflows/04_mofa_integration.R`  
**Lines:** ~1473  
**Purpose:** Integrate RNA pseudobulk views and lipidomics with MOFA, enrich factors for gonadal type, extract factor weights, infer TF activity, and build gene-lipid and TF-gene-lipid networks.

## Current execution contract (audited 18 November 2026)

The full body runs under `source()` and `Rscript`, with no CLI parser, source guard,
validate-only or fresh-output option. It always prepares/saves `MOFA_data.RData`
and can overwrite weights, enrichment, network and figure outputs under
`results/mofa/`. Use an isolated checkout/copy for reproduction without changing
curated results. The entry's literal/default paths control loading; neither
`DATASETS` nor `MOFA_SETTINGS` redirects its hardcoded model loop/options.

It skips training only when at least11 filenames match
`results/mofa/RNA/mofa_res_.*factor.hdf5`. All eleven5-15 factor models are retained
in the public tree, so the normal supplied-cache run reuses them. If this count
gate fails, the script automatically writes `data/mofa/MOFA_RNA_data.csv` and trains
missing factors5-15 via `python scripts/mofa/mofa2.py options_MOFA.csv` (GPU attempt,
CPU retry). A cached weights RData does not bypass this earlier training step.
The fallback weights branch is reached later if no model is selected/available.

`core/multiomics.R` unconditionally loads `MOFA2`; therefore MOFA2 is an actual
package prerequisite despite the entry's later unavailable-package/fallback
message. Python/mofapy2 is needed only for the automatic training contingency.
Public main CSVs, lipid data, workbook, support resources and retained models are
sufficient; no Seurat object, Workflow02 DE or current zonation output is consumed.
Current05 does not depend on04. Global07 consumes the permissive gene-lipid network;
11/12 and optional causal workflows use their separately declared frozen branches.

## Source Dependencies

```r
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/core/multiomics.R")
source("scripts/core/enrichment_analysis.R")
source("scripts/config/plotting_themes.R")
source("scripts/visualization/heatmaps.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Prepare MOFA Input | 78-117 | Load RNA/lipidomics, harmonize samples in `prepare_mofa_input()`, and save `MOFA_data.RData`. |
| 2 | MOFA Model Training | 120-197 | Use existing 5-15 factor models or run Python `scripts/mofa/mofa2.py`. |
| 3 | Model Selection | 200-246 | Prefer `mofa_res_12factor.hdf5`. |
| 4 | Extract Weights and TF Enrichment | 249-369 | Load MOFA model or existing weights, compute TF scores, and save weights RData. |
| 5 | Factor Enrichment | 372-672 | Enrich factor scores for gonadal type and plot variance/enrichment figures. |
| 6 | Analyze TF Enrichment | 675-753 | Summarize TF enrichment in hepatocyte weights. |
| 7 | Gene-Lipid Networks | 756-1390 | Select Factor 1 features, compute gene-lipid correlations, save full/sparse/TF/compressed networks. |
| - | Complete | 1438-1473 | Print output summary. |

## Input Files

- `data/RNA/scaled_RNA_{1,2,3}.csv` - Workflow 01 main RNA matrix.
- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv`
- `results/mofa/RNA/mofa_res_*factor.hdf5` - Retained pre-trained models; a count below11 triggers automatic training of missing5-15 factors.
- `results/mofa/main_model_RNA_lipid_weight.RData` - Later weight fallback if no model is available; it does not suppress the earlier training gate or the MOFA2 package requirement.
- `data/Suppl_Table_1_OverviewData.xlsx` via `load_metadata()`.
- `support/m5.all.v2023.2.Mm.symbols.gmt` - Lipid-metabolism gene filter for networks.
- `support/TF_regulons.RData`
- `support/Gonadal_Hormone_Genes.csv`

## Output Files

- `results/mofa/MOFA_data.RData`
- `data/mofa/MOFA_RNA_data.csv` - Written only if model training is needed.
- `results/mofa/RNA/mofa_res_12factor.hdf5` - Expected optimal model location.
- `results/mofa/main_model_RNA_lipid_weight.RData`
- `results/mofa/factor_analysis/factor_enrichment_gonadal.RData`
- `results/mofa/factor_analysis/variance_explained_pheatmap.pdf`
- `results/mofa/factor_analysis/gonadal_enrichment_pheatmap.pdf`
- `results/mofa/factor_analysis/mofa_variance_explained_with_gonadal_enrichment.pdf`
- `results/mofa/networks/gene_lipid_correlations_factor1.RData`
- `results/mofa/networks/gene_lipid_correlations_factor1.csv`
- `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_full_SIF.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_full_ATT.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv`
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_ATT.csv`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `MOFA_data` | data.frame | Long-format MOFA input with six views |
| `model_dir` | character | `results/mofa/RNA` |
| `optimal_model_file` | character | `results/mofa/RNA/mofa_res_12factor.hdf5` |
| `weights` | list | MOFA weights by view, plus TF score view when available |
| `factor_enrichment` | data.frame | Factor association with gonadal type |
| `factor1_features` | list | Top positive/negative Factor 1 features by view |
| `hepatocyte_cor_mat_RNA_lip` | matrix | Factor 1 gene-lipid correlations for sparse network |
| `hepatocyte_cor_mat_RNA_lip_permissive` | matrix | Top 33% feature correlation matrix for full network |

## Coherence Assessment

**Status: GOOD.** Workflow 04 is independent of workflow 02 DE outputs. It depends on workflow 01 RNA data and raw lipidomics, then feeds workflow 07 through the full gene-lipid network.

## Network Types Generated

| Type | Feature filter | Network resolution | Output |
|------|----------------|--------------------|--------|
| Sparse gene-lipid | Top 2.5% positive and negative Factor 1 features | Individual lipid species; top 3 genes per lipid | Internal input for TF-gene-lipid network |
| Full gene-lipid | Top 33% positive and negative Factor 1 features | Individual lipid species; all gene-lipid pairs | `hepatocyte_gene_lipid_network_full_SIF.csv` |
| TF-gene-lipid full | Sparse gene-lipid plus TF-gene edges | TFs, genes, lipid species | `hepatocyte_TF_gene_lipid_network_full_*` |
| TF-gene-lipid compressed | Aggregates lipid species to classes | TFs, genes, lipid classes | `hepatocyte_TF_gene_lipid_network_compressed_*` |

## Notes for Context Restoration

1. Factor 1 is the main gonad-associated factor in the current analysis narrative.
2. `prepare_mofa_input()` now filters RNA and lipidomics to their sample intersection.
3. `results/mofa/main_model_RNA_lipid_weight.RData` is important for lookup, but current workflow 06 no longer consumes it.
4. Workflow 07 consumes `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`.
5. Workflow 12 pooled phase (former 05c) consumes the sparse species-level Figure 3 branch instead:
   192 lipid-gene pairs from `hepatocyte_TF_gene_lipid_network_full_SIF.csv`
   and its 139-by-64 `gene_lipid_correlations_factor1.csv` precursor.
   The revised figure's "FULL SPECIES-LEVEL NETWORK" label does not refer
   to the permissive 33%-tail network. Workflow 12 replays frozen values,
   preserves top-three membership (minimum |r| approximately 0.5847), and
   never retrains MOFA or imposes a new 0.65 threshold.
