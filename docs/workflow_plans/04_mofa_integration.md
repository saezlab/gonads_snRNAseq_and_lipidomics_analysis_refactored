# Workflow 04: MOFA Integration

**Script:** `scripts/workflows/04_mofa_integration.R`
**Lines:** ~1438
**Purpose:** Multi-Omics Factor Analysis integrating RNA expression (by cell type) with lipidomics, factor enrichment, TF analysis, gene-lipid correlation networks.

## Source Dependencies

```r
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/core/enrichment_analysis.R")
source("scripts/core/multiomics.R")
source("scripts/config/analysis_config.R")
source("scripts/config/plotting_themes.R")
source("scripts/visualization/heatmaps.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | MOFA Input Prep | 1-150 | Combine RNA (per cell type) + lipidomics into MOFA format |
| 2 | Model Training | 151-300 | Train MOFA models (5-15 factors), optional skip if pre-trained |
| 3 | Model Selection | 301-400 | Compare models, select model12 (12 factors) |
| 4 | Weight Extraction | 401-600 | Extract W matrix, enrich for TFs using CollecTRI |
| 5 | Factor Enrichment | 601-800 | Visualize factor-gonad associations, publication figures |
| 6 | TF Analysis | 801-950 | Nfic, Dbp, Esr1 enrichment in factor weights |
| 7 | Gene-Lipid Networks | 951-1350 | Sparse (r>0.65), full (r>0), compressed (by class) networks |
| 8 | Summary | 1351-1438 | Print summary, list output files |

## Input Files

- `data/RNA/scaled_RNA_{1,2,3}.csv` - Expression data (from workflow 01)
- `data/RNA/meta_data_df.csv` - Metadata
- `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` - Lipidomic data
- `data/mofa/` - Pre-computed MOFA models (optional)
- `support/TF_regulons.RData` - TF regulons
- `support/Gonadal_Hormone_Genes.csv` - ESR1/AR targets
- `support/m5.all.v2023.2.Mm.symbols.gmt` - GO gene sets (for lipid metabolism filter)

## Output Files

- `results/mofa/MOFA_data.RData` - MOFA input data
- `results/mofa/main_model_RNA_lipid_weight.RData` - **KEY OUTPUT**: MOFA weights with TF scores
- `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` - Full gene-lipid network
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_*_SIF.csv` - TF-gene-lipid networks
- `results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata` - Correlation matrix for structural analysis

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `mofa_model` | MOFA | Trained MOFA model object |
| `weights` | list | MOFA W matrix (list of views -> data.frames) |
| `TF_scores_df` | data.frame | TF activity scores (TFs x samples) |
| `cor_mat` | matrix | Gene-lipid Pearson correlations |
| `threshold_cor` | numeric | Correlation threshold (0.65 for sparse, 0 for full) |

## Sample Harmonization (2026-03-13)

`prepare_mofa_input()` in `scripts/core/multiomics.R` now filters **both** RNA and lipidomic
samples to their intersection (43 common samples). Previously, only RNA was filtered to match
lipidomic, leaving 48 lipidomic samples (5 extra) in the MOFA input. Old models archived in
`results/mofa/RNA/old_20260313/`.

## Coherence Assessment

**Status: GOOD** — No real issues found. Notes on design patterns below.

### Model loading has proper fallback logic
- Three paths: (1) train new model, (2) load pre-trained `.hdf5`, (3) load pre-computed `weights` from `.RData`
- Both paths 2 and 3 produce the same `weights` list object (list of views -> data.frames)
- If neither model nor weights file exists, the script calls `stop()` — no silent failure
- Downstream code uses `if (!is.null(model))` to branch correctly

### Network generation is sequential by design
- Three versions (sparse, full, compressed) are built in order within section 7
- The compressed version intentionally overwrites `cormat_*` variables after the uncompressed versions have been saved to CSV
- This is safe because the script runs sequentially top-to-bottom

## Network Types Generated

| Type | Threshold | Node Resolution | Purpose |
|------|-----------|-----------------|---------|
| Sparse | r > 0.65 | Individual lipids | Cytoscape visualization |
| Full | r > 0 | Individual lipids | Validation (workflow 07) |
| Compressed | r > 0.65 (pre-agg) | Lipid classes | Summary figures |

## Notes for Context Restoration

When returning to this workflow:
1. `main_model_RNA_lipid_weight.RData` is the single most important output - used by workflows 05, 06, 07
2. Factor 1 shows strongest gonad association - this is the key biological factor
3. Top 50% filtering: features are selected by absolute weight magnitude from Factor 1
4. Gene filter: Only lipid metabolism genes (GO:CELLULAR_LIPID_METABOLIC_PROCESS) are included in networks
5. `make_cor_network()` and `expand_attributes()` are defined in `scripts/support_functions.R`
6. Model12 (12 factors) was selected after cross-validation
7. Key TFs identified: Nfic, Dbp, Esr1 (correlated with Cyp2b9, Cyp2d9, Cyp2b13, Cyp2c38)
