# GDX Data Update Plan

**Date:** 2026-02-25
**Context:** Collaborator regenerated all four pseudobulk datasets. GDX pseudobulks now contain only 4 groups (sham + GDX for males/females) instead of previous 8 groups (which also included estrogen/testosterone replacement). GDX lipidomics must be filtered to match the 22 samples in the new pseudobulks.

## Decision Required: FCG Pseudobulk Version

The new `fcg_pseudobulks.rds` is **identical** to `fcgl_pseudobulks_gonad_20241114.rds` but **differs** from `fcgl_pseudobulk_gonad_20241023.rds` (the version currently used by workflow 01).

| Metric | New vs 20241023 | New vs 20241114 |
|--------|-----------------|-----------------|
| Genes | Same (20,123) | Same (20,123) |
| Samples | Same (215) | Same (215) |
| Metadata | Identical | Identical |
| Expression values | Mean diff = 0.17, max = 12.03 | Identical (diff < 2.2e-15) |
| Gene correlations | Median r = 0.94, min r = -0.20 | r = 1.0 |

---

## Option A: Keep Current FCG (20241023) — GDX-Only Update

**Pros:** MOFA model untouched, minimal changes, faster
**Cons:** FCG pseudobulk in repo doesn't match what collaborator considers "current"

### Steps

#### A1. Convert new GDX pseudobulk to CSVs

- Input: `original_datasets/snRNAseq_seurat_objects/gdx_pseudobulks.rds`
- Adapt `scripts/legacy/gdxhr_seurat_to_df.R`:
  - The new data has ~110 pseudobulk samples (22 mice × 5 cell types) instead of ~190
  - Hardcoded column index on line 21 (`scaled_RNA[,c(191,1:190)]`) must be updated dynamically
  - Split into 3 CSVs proportionally
- Output: `data/RNA/scaled_RNA_gdxhr_[1,2,3].csv` + `data/RNA/meta_data_df_gdxhr.csv`
- Backup old files first

#### A2. Filter GDX lipidomics to 22 matching samples

- Load new GDX pseudobulk metadata → extract unique `sample_name` values (should be 22)
- Load `original_datasets/metabolomics/lipidomics_gdx.csv` (44 samples)
- Match sample names: RNA metadata names (e.g., `csb177_01`) to lipidomics column names (e.g., `csb177-01_m_s_v`)
- Keep only the 22 matching sample columns + all metadata columns (1-34)
- Save as `data/metabolomics/lipidomics_gdxhr.csv` (replacing old 38-sample file)

#### A3. Update GDX contrasts in workflow 02 (Section 3B, lines ~381-390)

Old contrasts (8 groups → indices 1-8):
```r
gdx_contrasts <- list(
    "testes_gonadosectomy" = c(5, -1),
    "testes_oestro" = c(8, -5),
    "testes_testo" = c(7, -5),
    "ova_gonadosectomy" = c(6, -3),
    "ova_oestro" = c(4, -6),
    "ova_testo" = c(2, -6)
)
```

New contrasts (4 groups, alphabetical: f-g-v=1, f-s-v=2, m-g-v=3, m-s-v=4):
```r
gdx_contrasts <- list(
    "testes_gonadosectomy" = c(3, -4),  # m-g-v vs m-s-v
    "ova_gonadosectomy" = c(1, -2)      # f-g-v vs f-s-v
)
```
**Note:** Exact indices must be verified after loading the metadata — print `unique(targets$condition)` to confirm order.

Also update surrounding text/comments that reference 6 contrasts and hormone replacement groups.

#### A4. Update GDX contrasts in workflow 03 (Section 3.5, lines ~244-284)

Two changes:
1. **Line 245** — Column selection `lipidomic_gdx[, c(4, 35:72)]` must change. Column 4 is "Metabolite name", columns 35:72 were the old 38 sample columns. With 22 samples and potentially different metadata column count, indices need recalculation. Better approach: select column 4 + all columns matching sample name patterns.
2. **Lines 277-284** — Reduce from 6 contrasts to 2, with indices matching the 4 lipidomics groups (alphabetical: f_g_v=1, f_s_v=2, m_g_v=3, m_s_v=4):
```r
comparisons_gdx <- list(
    "testes_gonadosectomy" = c(3, -4),  # m_g_v vs m_s_v
    "ova_gonadosectomy" = c(1, -2)      # f_g_v vs f_s_v
)
```
**Note:** Verify indices after loading — the condition extraction (`sub("^[^_]*_", "", ...)`) determines the group names.

#### A5. Re-run workflows 02 and 03

- `Rscript scripts/workflows/02_fcg_pseudobulk_analysis.R` — FCG sections produce identical output (same input CSVs), GDX sections regenerated
- `Rscript scripts/workflows/03_fcg_lipidomic_analysis.R` — Same: FCG unchanged, GDX regenerated
- Do **NOT** re-run workflow 01 (FCG CSVs must stay the same)

#### A6. Verify MOFA inputs unchanged

- Confirm `data/RNA/scaled_RNA_[1,2,3].csv` are byte-identical (not touched)
- Confirm `data/RNA/meta_data_df.csv` is identical
- Confirm `data/metabolomics/lipidomic.RData` was never modified

#### A7 (Optional). Convert new GDX zonation pseudobulk

- Input: `original_datasets/snRNAseq_seurat_objects/gdx_zonation_pseudobulks.rds`
- Adapt `scripts/legacy/gdxhr_seurat_to_df_zonation.R` similarly
- Output: `data/RNA/scaled_RNA_gdxhr_zonation_[1,2,3].csv` + metadata
- Not currently used by any workflow (05-08 have no GDX sections), but good to have ready

### Expected GDX Output Changes

| File | Change |
|------|--------|
| `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/gdxrh/gdhxr_groups_ttop_list.RData` | Regenerated: 2 contrasts instead of 6 |
| `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv` | Regenerated: 2 columns instead of 6 |
| `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/gdx_validation_comparison.pdf` | Regenerated: may show different correlation |
| `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gdxhr_ttop_list.RData` | Regenerated: 2 contrasts instead of 6 |
| `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gdxhr_t_table.csv` | Regenerated: 2 columns instead of 6 |
| `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/fcg_gdx_correlation.pdf` | Regenerated: may show different correlation |

---

## Option B: Switch to New FCG (matches 20241114) — Full Re-run

**Pros:** Everything consistent with collaborator's current data, clean slate
**Cons:** MOFA must be re-trained (or at minimum re-verified), all downstream results change

### Steps

#### B1. Update FCG pseudobulk reference in workflow 01

- Copy `original_datasets/snRNAseq_seurat_objects/fcg_pseudobulks.rds` to `data/RNA/`
- Update `scripts/workflows/01_data_preparation_qc.R` line 91:
  ```r
  seurat_path <- "data/RNA/fcg_pseudobulks.rds"  # was fcgl_pseudobulk_gonad_20241023.rds
  ```

#### B2. Re-run workflow 01

- Generates new `scaled_RNA_[1,2,3].csv` and `meta_data_df.csv`
- These will differ from current (matches 20241114 scaling, not 20241023)

#### B3-B7. Same as Option A steps A1-A5

(Convert GDX data, filter lipidomics, update contrasts, re-run workflows 02+03)

#### B8. Re-run workflow 04 (MOFA)

- MOFA input data will change (new FCG expression values)
- Two sub-options:
  - **B8a.** Re-train MOFA models (slow, requires MOFA2 + Python backend)
  - **B8b.** Load existing pre-trained model but regenerate weights/networks from new data
    - Risk: Factor structure may not match if expression differences are large enough
- Verify Factor 1 still shows strong gonadal association (R² ~ 0.8)
- Regenerate all network files in `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/`

#### B9. Re-run workflows 05-08

Since MOFA outputs changed, all downstream workflows need re-running:
- **05** (Zonation) — Uses MOFA weights for TF enrichment
- **06** (Structural analysis) — Uses MOFA-derived gene-lipid correlations
- **07** (Validation) — Uses MOFA network files
- **08** (Variance analysis) — Uses pseudobulk results

#### B10. Also update FCG zonation pseudobulk

- Copy `original_datasets/snRNAseq_seurat_objects/fcg_zonation_pseudobulks.rds`
- Convert to CSVs for workflow 05
- The new FCG zonation file (23.3 MB) matches `fcgl_heps_zonationsubsetted_pseudobulks_gonad_20241114.rds` (23.3 MB) — verify identity like we did for the main FCG

### Expected Impact

| Component | Affected? |
|-----------|-----------|
| FCG RNA CSVs | Yes — expression values change |
| FCG DE results (workflow 02) | Yes — t-values/p-values change |
| FCG lipidomic results (workflow 03) | No — FCG lipidomics input is independent |
| MOFA model (workflow 04) | Yes — input changed, needs re-run |
| MOFA networks | Yes — regenerated |
| Zonation results (workflow 05) | Yes — depends on MOFA + RNA |
| Structural analysis (workflow 06) | Yes — depends on MOFA |
| Validation (workflow 07) | Yes — depends on networks |
| Variance analysis (workflow 08) | Yes — depends on DE results |
| All GDX results | Yes — updated data + new FCG comparison baseline |

---

## Summary

| | Option A (Keep 20241023) | Option B (Switch to new FCG) |
|---|---|---|
| Workflows to re-run | 02, 03 (GDX sections only) | 01-08 (everything) |
| MOFA re-training | No | Yes |
| FCG results change | No | Yes |
| GDX results change | Yes | Yes |
| Estimated effort | Small (contrast updates + re-run 2 workflows) | Large (full pipeline re-run) |
| Risk | FCG pseudobulk mismatch with collaborator | MOFA factor structure might shift |
