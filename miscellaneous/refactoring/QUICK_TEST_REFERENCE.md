# Quick Test Reference Card

## Run All Workflow Tests

From R console in project root:

```r
source("tests/run_workflow_tests.R")
```

This will run all 6 workflows and report:
- ✓ Which workflows passed
- ✗ Which workflows failed
- File outputs created
- Key result validations

**Expected runtime:** 70-135 minutes total

---

## Test Individual Workflows

### Workflow 2: Pseudobulk Analysis (Most Important)

```r
source("workflows/02_fcg_pseudobulk_analysis.R")

# Quick validation
load("results/pseudobulk/ttop_full_list.RData")
names(ttop_full_list)  # Should show 6 contrasts

# Check PRIMARY finding
hepa <- ttop_full_list$gonadal_type$Hepatocytes
sum(hepa$adj.P.Val <= 0.05)  # Should be ~4,000 DEGs
```

**Expected:** 3,000-5,000 DEGs (manuscript reports >4,000)

### Workflow 3: Lipidomic Analysis

```r
source("workflows/03_fcg_lipidomic_analysis.R")

# Quick validation
load("results/lipidomic/ttop_full_list.RData")

# Check PRIMARY finding
gonad <- ttop_full_list$gonadal_type
sum(gonad$adj.P.Val <= 0.05)  # Should be ~733 lipids
```

**Expected:** 650-850 lipids (manuscript reports 733, 53.9%)

---

## What Should Work Now

✓ **Fixed Issues:**
1. `group_variable` correctly extracted from contrast config
2. All 6 contrasts use the right metadata column
3. FCG group comparisons (XXT_vs_XXO, etc.) now work

✓ **All Contrasts:**
1. `gonadal_type` → uses column `gonadal_type`
2. `allosome_type` → uses column `allosome_type`
3. `XXT_vs_XXO` → uses column `group`
4. `XYT_vs_XYO` → uses column `group`
5. `XYO_vs_XXO` → uses column `group`
6. `XYT_vs_XXT` → uses column `group`

---

## Key Validation Checks

### ✓ Pseudobulk Analysis (Workflow 2)
- [ ] All 6 contrasts analyzed
- [ ] All 5 cell types analyzed
- [ ] ~4,000 hepatocyte gonadal DEGs
- [ ] <100 hepatocyte allosome DEGs

### ✓ Lipidomic Analysis (Workflow 3)
- [ ] ~733 gonad-dependent lipids (53.9%)
- [ ] ~0 chromosome-dependent lipids
- [ ] Male-biased: TG, DG, phospholipids
- [ ] Female-biased: Sphingolipids

### ✓ MOFA Integration (Workflow 4)
- [ ] Factor 1 associates with gonad
- [ ] Cyp2b9, Cyp2d9, Cyp2c38, Cyp2b13 identified
- [ ] Nfic, Dbp, Esr1 TFs identified

---

## Common Errors (SOLVED)

### ❌ OLD ERROR:
```
Error in get_contrast(contrast_type) :
  Contrast 'group' not found
```

### ✓ FIXED:
Now extracts `group_variable` from contrast config:
```r
contrast_config <- get_contrast("XXT_vs_XXO")
group_variable <- contrast_config$group_variable  # Returns "group"
```

---

## Files Created

### Test Scripts:
- `tests/run_workflow_tests.R` - Interactive testing (recommended)
- `tests/test_all_workflows.R` - Comprehensive testthat suite
- `tests/TESTING_GUIDE.md` - Detailed testing documentation (this file)

### Reference Documents:
- `ANALYSIS_SCOPE_SUMMARY.md` - Complete analysis scope from manuscript
- `QUICK_TEST_REFERENCE.md` - This quick reference card

---

## Quick Troubleshooting

**Workflow fails immediately:**
```r
# Check data files exist
file.exists("data/RNA/scaled_RNA_1.csv")
file.exists("data/RNA/meta_data_df.csv")

# If missing, regenerate:
source("scripts/seurat_to_df.R")
```

**Wrong number of DEGs:**
```r
# Check metadata
source("scripts/core/data_loading.R")
meta <- load_metadata(rna_metadata_file = "data/RNA/meta_data_df.csv")
table(meta$gonadal_type)  # Should show Ovaries and Testes
table(meta$group)         # Should show XXO, XXT, XYO, XYT
```

---

**For detailed troubleshooting, see:** `tests/TESTING_GUIDE.md`

**For analysis scope details, see:** `ANALYSIS_SCOPE_SUMMARY.md`
