# Plan: Switch to New FCG Pseudobulk and Re-run Full Pipeline

## Context

A collaborator regenerated all four pseudobulk datasets. The new `fcg_pseudobulks.rds` is identical to the 20241114 version but differs from the 20241023 version previously used by workflow 01:
- Same dimensions (20,123 genes x 215 samples), same metadata
- Expression values differ: mean abs diff = 0.17, max = 12.03, median gene correlation r = 0.94

Switching FCG pseudobulk input changes all downstream outputs: RNA DE results, MOFA model, correlation networks, zonation analysis, structural analysis, and validation. This requires a careful, sequential re-run of the entire pipeline (workflows 01-08).

**Branch:** `update-fcg-pseudobulk`
**Strategy:** Step-by-step execution with a commit after each workflow, documenting changes as we go.

## Dependency Chain

```
01 (new FCG input) -> scaled_RNA_{1,2,3}.csv, meta_data_df.csv
  |-> 02 (RNA DE) -> ttop_full_list.RData, tval_full_df.csv
  |-> 03 (lipid DE) -> unchanged FCG lipidomic input, re-run for regenerated plots
  |-> 04 (MOFA) -> RE-TRAIN models -> new weights, networks
  |    |-> 06 (structural) -> new gene-lipid correlations
  |    |-> 07 (validation) -> new network files
  |-> 05 (zonation) -> also update zonation pseudobulk input
  |    |-> 07 (validation) -> new zonation networks
08 (variance) -> self-contained, re-run for completeness
```

---

## Steps 1-4: COMPLETED

### Step 1: Update FCG Main Pseudobulk + Run Workflow 01
- Switched `seurat_path` to `data/RNA/fcg_pseudobulks.rds`
- Formalized all Seurat-to-CSV conversions (FCG, zonation from `fcg_zonation_pseudobulks.rds`, GDX from `gdx_pseudobulks.rds`)
- Removed unused legacy data files from `data/RNA/` (moved to `data/RNA/legacy/`)
- **Commit:** `d34ce7d` "Switch FCG pseudobulk to new version and regenerate CSVs"
- **Commit:** `7478d22` "Formalize Seurat-to-CSV conversions in workflow 01 for all datasets"

### Step 2: Re-run Workflow 02 (RNA DE)
- No code changes needed
- **Commit:** `ddbac2d` "Regenerate RNA DE results with new FCG pseudobulk"

### Step 3: Re-run Workflow 03 (Lipidomic DE)
- No code changes needed. FCG lipidomic input unchanged.
- **Commit:** `746e499` "Re-run lipidomic analysis (FCG lipidomics unchanged, confirms clean state)"

### Step 4: Re-train MOFA Models + Run Workflow 04
- MOFA models re-trained with new FCG data (11 models, factors 5-15)
- Factor 1 still shows strongest gonadal association
- **Fix applied:** Replaced `load("lipidomic.RData")` with `load_lipidomic_data()` in WF04 and WF05 (eliminates hidden dependency on legacy file)
- **Commit:** `fd28b19` "Re-train MOFA models with new FCG data and regenerate weights/networks"
- **Commit:** `65278bd` "Fix lipidomic loading vulnerability and re-run WF01-04 with harmonized inputs"

### Additional fixes during Steps 1-4
- Added `source("scripts/config/plotting_themes.R")` to WF05 for `ZONATION_GONAD_COLORS`
- Added `source("scripts/core/data_loading.R")` to WF05 for `load_lipidomic_data()`
- **Commit:** `3b9ea49` "Regenerate zonation CSVs from 20241114 RDS and add conversion to workflow 01"
- **Commit:** `ee7f14f` "Re-run zonation analysis (workflow 05) with updated pseudobulk data"

---

## Steps 5-8: RE-RUN WITH RESULT COMPARISON

WF05-08 outputs currently on disk are stale (generated before lipidomic fix and before WF01-04 clean re-run). We need to re-run them with the fully harmonized pipeline and compare results before/after.

### Baseline Metrics (Before re-run)

| Metric | Current Value | Source |
|--------|---------------|--------|
| Zonation sex-specific genes | 138 | `sex_specific_zonation_genes_summary.csv` |
| Zonation gonad DEGs (testes_vs_ova) | 5,744 | `summary_statistics.csv` |
| Zonation interaction genes | 297 | `summary_statistics.csv` |
| Structural associations (total rows) | 22,757 | `significant_gene_lipid_class_associations.csv` |
| Chain length genes | 1,162 | `genes_associated_with_chain_length.csv` |
| Validation: genes validated (global) | 3 | `validation_summary.csv` |
| Validation: gene-lipid edges (global) | 643,406 | `validation_summary.csv` |
| Variance: sig diff variance lipids (FDR<0.05) | 2 | `variance_analysis_summary.csv` |
| Variance: mean CV Ovaries | 3.724% | `variance_analysis_summary.csv` |
| Variance: mean CV Testes | 3.567% | `variance_analysis_summary.csv` |

### Step 5: Run WF05 (Zonation Analysis)

**Script:** `scripts/workflows/05_zonation_analysis.R`

**Key inputs (from WF01-04):**
- `results/pseudobulk/ttop_full_list.RData` (WF02)
- `results/pseudobulk/tval_full_df.csv` (WF02)
- Lipidomic data via `load_lipidomic_data()` (fixed)
- Zonation CSVs from `data/RNA/scaled_RNA_zonation.csv`

**Comparison after run:**
1. Count sex-specific zonation genes (was 138)
2. Count testes_vs_ova DEGs in zonation (was 5,744)
3. Count interaction genes (was 297)
4. Verify all `results/zonation/*.csv` and `*.RData` updated

**Acceptance:** Counts should be in same ballpark (within ~10-15%).

### Step 6: Run WF06 (Structural Analysis)

**Script:** `scripts/workflows/06_structural_analysis.R`

**Key inputs (from WF04):**
- `results/mofa/main_model_RNA_lipid_weight.RData`
- Lipidomic data and structural files from `data/metabolomics/`
- RNA data from `data/RNA/scaled_RNA_{1,2,3}.csv`

**Comparison after run:**
1. Total significant associations (was 22,757 rows)
2. Chain length genes count (was 1,162)
3. Whether key genes (Cyp2b9, Cyp2d9, Cyp2c38, etc.) still appear
4. Check all `results/multi-omic/structural_model/*.csv` updated

**Acceptance:** Key genes should persist. Total counts may shift but same order of magnitude.

### Step 7: Run WF07 (Validation Experiments)

**Script:** `scripts/workflows/07_validation_experiments.R`

**Key inputs:**
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv` (WF04)
- `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_ATT.csv` (WF04)
- Zonation network files from `results/zonation/` (from WF05 re-run)
- siRNA knockdown data from `data/metabolomics/sirna_*.csv`

**Comparison after run:**
1. Genes validated (global network) (was 3)
2. Gene-lipid edges (global) (was 643,406)
3. Validation scatter plots regenerated

**Acceptance:** Validation gene count should be 3 (Cyp2b9, Cyp2c38, Cyp2d9). Edge counts may vary.

### Step 8: Run WF08 (Variance Analysis)

**Script:** `scripts/workflows/08_variance_analysis.R`

**Key inputs:**
- `results/lipidomic/ttop_full_list.RData` (WF03)
- Lipidomic data via `load_lipidomic_data()`

**Comparison after run:**
1. Sig differential variance lipids FDR<0.05 (was 2)
2. Mean CV values (Ovaries: 3.724%, Testes: 3.567%)

**Acceptance:** Since lipidomic data is unchanged, expect identical or near-identical results.

---

## Step 9: Final Verification and Documentation

### Metrics Comparison Table

| Metric | Before (20241023) | After (new) | Acceptable range | Status |
|--------|-------------------|-------------|------------------|--------|
| Hepatocyte gonad DEGs (FDR<0.05) | ~4000 | 4,103 | >3000 | PASS |
| GDX RNA validation r | -0.534 | -0.5351 | < -0.3 | PASS |
| GDX lipid validation r | -0.721 | -0.7206 | < -0.5 | PASS |
| MOFA Factor 1 gonadal R² | ~0.8 | ~0.8 (confirmed) | > 0.6 | PASS |
| Gonad-dependent lipids | 773 | 773 | Should be identical | PASS |
| Zonation sex-specific genes | ~137 | 138 | > 100 | PASS |
| Structural associations | 22,757 | 21,989 (-3.4%) | Same order of magnitude | PASS |
| Validation genes (global) | 3 | 3 | 3 | PASS |
| Variance sig lipids (FDR<0.05) | 2 | 2 | Identical | PASS |

### Final actions
1. Commit all WF05-08 regenerated results - DONE (commit `a849904`)
2. Produce final comparison table with filled-in values - DONE (see above)
3. Review git log for complete branch history - DONE
4. Merge `update-fcg-pseudobulk` branch into `main` - PENDING (user decision)

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| MOFA factor structure shifts | Medium | High | Compare variance explained, check Factor 1 still gonadal |
| Zonation threshold changes | Medium | Medium | May change gene selection; document |
| Network file naming breaks | Low | High | Same code, just new data values |
| Python mofapy2 install fails | Low | High | pip install; fallback: use old weights |
| Expression r=0.94 causes major result shifts | Low | Medium | Most results should be qualitatively similar |

## Execution Order

All steps are sequential (each depends on the previous):
Step 1 → Step 2 → Step 3 → Step 4 → Step 5 → Step 6 → Step 7 → Step 8 → Step 9
