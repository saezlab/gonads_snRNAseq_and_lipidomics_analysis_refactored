# Refactoring Guide

This guide provides step-by-step instructions for safely refactoring this codebase while ensuring scientific reproducibility.

## Quick Start

```bash
# 1. Save current outputs as reference
Rscript tests/save_reference_outputs.R

# 2. Make your changes to the code
# ... edit scripts ...

# 3. Re-run the analysis
Rscript run_pipeline.R

# 4. Verify outputs haven't changed
Rscript tests/test_output_consistency.R
```

If step 4 passes, your refactoring is safe! If it fails, review the differences to ensure they're expected.

## Prerequisites

Before starting refactoring, ensure you have:

1. **R and required packages installed**
   ```r
   install.packages(c("readr", "readxl", "ggplot2", "ggrepel", "pheatmap",
                      "limma", "reshape2", "cowplot", "gridExtra", "decoupleR",
                      "dplyr", "Seurat", "MOFA2", "ComplexHeatmap", "circlize",
                      "digest", "jsonlite", "GSEABase"))
   ```

2. **All required data files in place**
   - Check that `data/` directory contains all input files
   - See `PIPELINE.md` for complete list of required files

3. **Baseline reference outputs**
   ```bash
   Rscript tests/save_reference_outputs.R
   ```

4. **Git repository initialized** (recommended)
   ```bash
   git init
   git add .
   git commit -m "Initial commit before refactoring"
   ```

## Refactoring Strategy

### Recommended Order

Refactor scripts in reverse dependency order (outputs to inputs):

1. **Visualization scripts** (safest - don't affect data)
   - Figure generation scripts
   - Plotting code in analysis scripts

2. **Downstream analysis**
   - `scripts/pseudobulk_pathways.R` - Pathway analysis
   - `scripts/MOFA_lipidcell_downstream.R` - MOFA downstream

3. **Core analysis**
   - `scripts/pseudobulk_DA.R` - RNA differential analysis
   - `scripts/lipidomic_DA.R` - Lipidomics differential analysis

4. **Data preparation** (most risky - affects everything)
   - `scripts/seurat_to_df.R` - RNA data preparation

5. **Support functions**
   - `scripts/support_functions.R` - Only after testing downstream scripts

6. **MOFA** (special handling required)
   - `scripts/MOFA_lipidcell_RNA.R` - See MOFA section below

### Common Refactoring Tasks

#### 1. Extract Repeated Code into Functions

**Before:**
```r
# In pseudobulk_DA.R
volcano_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1,
                            pValIndex = 5, label = T, nlabels = 40)
name_of_file <- paste("results/pseudobulk/volcanos/",
                      paste(names(comparisons)[1],".pdf",sep = ""),sep = "")
pdf_size <- 9
ggsave(filename = name_of_file, plot = volcano_plt, device = "pdf",
       width = pdf_size, height = pdf_size)

# Repeated 5 more times...
```

**After:**
```r
# Add to support_functions.R
save_volcano_plot <- function(ttop, comparison_name, output_dir = "results/pseudobulk/volcanos/",
                              FCIndex = 2, IDIndex = 1, pValIndex = 5,
                              label = TRUE, nlabels = 40, pdf_size = 9) {
  volcano_plt <- volcano_nice(ttop, FCIndex = FCIndex, IDIndex = IDIndex,
                              pValIndex = pValIndex, label = label, nlabels = nlabels)
  filename <- file.path(output_dir, paste0(comparison_name, ".pdf"))
  ggsave(filename = filename, plot = volcano_plt, device = "pdf",
         width = pdf_size, height = pdf_size)
  return(volcano_plt)
}

# In pseudobulk_DA.R
volcano_plt <- save_volcano_plot(ttop_list[[1]], names(comparisons)[1])
```

**Test:**
```bash
Rscript tests/test_output_consistency.R
```

#### 2. Standardize File Paths

**Before:**
```r
# Scattered throughout scripts
read_csv("data/RNA/scaled_RNA_1.csv")
read_csv("results/pseudobulk/ttop_full_list.RData")
```

**After:**
```r
# Add to top of each script or support_functions.R
DATA_DIR <- "data"
RESULTS_DIR <- "results"

RNA_DATA_DIR <- file.path(DATA_DIR, "RNA")
LIPID_DATA_DIR <- file.path(DATA_DIR, "metabolomics")
PSEUDOBULK_RESULTS <- file.path(RESULTS_DIR, "pseudobulk")
LIPID_RESULTS <- file.path(RESULTS_DIR, "lipidomic")

# In scripts
scaled_RNA_1 <- read_csv(file.path(RNA_DATA_DIR, "scaled_RNA_1.csv"))
```

#### 3. Add Progress Indicators

**Before:**
```r
for (name in names(csv_files)) {
  df <- read_csv(csv_files[[name]])
  # process...
}
```

**After:**
```r
cat(sprintf("Processing %d files...\n", length(csv_files)))
for (i in seq_along(csv_files)) {
  name <- names(csv_files)[i]
  cat(sprintf("[%d/%d] Processing %s...\n", i, length(csv_files), name))
  df <- read_csv(csv_files[[name]])
  # process...
}
cat("Done!\n")
```

#### 4. Add Input Validation

**Before:**
```r
runLimma <- function(measurements, targets, comparisons) {
  # ... limma code ...
}
```

**After:**
```r
runLimma <- function(measurements, targets, comparisons) {
  # Validate inputs
  if (!is.data.frame(measurements) && !is.matrix(measurements)) {
    stop("measurements must be a data frame or matrix")
  }
  if (!is.data.frame(targets)) {
    stop("targets must be a data frame")
  }
  if (!is.list(comparisons)) {
    stop("comparisons must be a list")
  }
  if (ncol(measurements) != nrow(targets)) {
    stop(sprintf("Dimensions mismatch: %d columns in measurements, %d rows in targets",
                 ncol(measurements), nrow(targets)))
  }

  # ... limma code ...
}
```

#### 5. Modularize Large Scripts

For scripts over 200 lines, consider breaking into modules:

**Before:** `pseudobulk_DA.R` (400+ lines)

**After:** Create multiple files:
- `scripts/pseudobulk/01_load_data.R`
- `scripts/pseudobulk/02_run_limma.R`
- `scripts/pseudobulk/03_visualize_results.R`
- `scripts/pseudobulk/04_save_outputs.R`

Then in `pseudobulk_DA.R`:
```r
source("scripts/pseudobulk/01_load_data.R")
source("scripts/pseudobulk/02_run_limma.R")
source("scripts/pseudobulk/03_visualize_results.R")
source("scripts/pseudobulk/04_save_outputs.R")
```

## Special Cases

### MOFA (Multi-Omics Factor Analysis)

MOFA requires special handling due to computational cost:

**What you CAN refactor:**
- Code that prepares MOFA input (BEFORE line 60 in MOFA_lipidcell_RNA.R)
- Code that analyzes MOFA outputs (AFTER line 79)
- Visualization code

**What you CANNOT easily refactor:**
- MOFA training code (lines 63-67, commented out)
- MOFA input data structure

**Critical constraint:**
- `data/mofa/MOFA_RNA_data.csv` must remain IDENTICAL
- Tests will FAIL if this file changes
- If it must change, you need to retrain all MOFA models (very time-consuming)

**Testing MOFA refactoring:**
```r
# After refactoring MOFA input preparation
source("scripts/MOFA_lipidcell_RNA.R")  # Only runs up to line 60

# Check MOFA input hasn't changed
Rscript tests/test_output_consistency.R  # Will fail if MOFA input changed
```

### Support Functions

`scripts/support_functions.R` is used by many scripts. Refactoring it affects the entire pipeline:

**Strategy:**
1. Refactor one function at a time
2. Test ALL scripts that use it
3. Consider versioning functions if major changes needed

**Functions and their users:**
- `runLimma()` - Used by: pseudobulk_DA.R, lipidomic_DA.R, gdxhr_pseudobulk_DA.R, etc.
- `volcano_nice()` - Used by: pseudobulk_DA.R, lipidomic_DA.R, etc.
- `nicePCA()` - Used by: lipidomic_DA.R, pseudobulk_DA.R
- `createLinearColors()` - Used by: many scripts for heatmaps

**Testing support functions:**
```bash
# After modifying support_functions.R, test entire pipeline
Rscript run_pipeline.R
Rscript tests/test_output_consistency.R
```

## Testing Workflow

### 1. Before Each Refactoring Session

```bash
# Ensure reference exists
ls tests/reference_outputs/reference_checksums.rds

# If not, create it
Rscript tests/save_reference_outputs.R

# Create a git branch
git checkout -b refactor-descriptive-name
```

### 2. During Refactoring

```bash
# Make small, incremental changes
# Test frequently

# After each change:
Rscript scripts/your_modified_script.R  # Or run_pipeline.R for full test
Rscript tests/test_output_consistency.R
```

### 3. If Tests Fail

```bash
# Review what changed
cat tests/reference_outputs/reference_checksums.json | grep -A 5 "failed_file"

# Options:
# A) Bug in refactoring - fix it
# B) Expected change (e.g., bug fix) - document and update reference
# C) Numerical precision issue - investigate tolerance

# If expected change, update reference:
Rscript tests/save_reference_outputs.R
git add tests/reference_outputs/
git commit -m "Update reference after fixing bug in X"
```

### 4. After Refactoring

```bash
# Run full pipeline
Rscript run_pipeline.R

# Verify all tests pass
Rscript tests/test_output_consistency.R

# Commit changes
git add .
git commit -m "Refactor: descriptive message"
git push origin refactor-descriptive-name
```

## Common Issues

### Issue: "Reference checksums not found"

**Cause:** Haven't created reference yet

**Solution:**
```bash
Rscript tests/save_reference_outputs.R
```

### Issue: "MOFA input has changed"

**Cause:** Refactoring affected data fed into MOFA

**Solutions:**
1. Revert changes that affected MOFA input
2. Accept that MOFA models need retraining (time-consuming)
3. Investigate if change is due to bug fix (document carefully)

### Issue: Tests pass but plots look different

**Cause:** Tests don't check visual appearance

**Solution:** Visual inspection of plots is still required. Tests only verify numerical consistency.

### Issue: Small numerical differences in floating point

**Cause:** Order of operations or rounding differences

**Example:**
```r
# These may give slightly different results
sum(x) / length(x)  # vs
mean(x)
```

**Solution:**
- Review code for order-dependent operations
- Consider using `all.equal()` with tolerance instead of exact equality
- Document if acceptable precision loss

### Issue: File exists but checksum differs

**Cause:**
- Line ending changes (Windows/Unix)
- Encoding changes
- White space changes

**Solution:** Check content checksums (Test 2) - they're more robust to metadata changes.

## Best Practices

1. **Make small, incremental changes**
   - Refactor one function or section at a time
   - Test after each change
   - Commit frequently

2. **Document changes**
   - Update comments when refactoring
   - Add docstrings to new functions
   - Update CLAUDE.md if major changes

3. **Preserve scientific reproducibility**
   - Don't change random seeds
   - Don't change statistical methods without documentation
   - Keep track of any differences

4. **Use version control**
   - Branch for each refactoring task
   - Commit before and after tests
   - Tag major milestones

5. **Test on real data**
   - Don't just test on toy examples
   - Use the full pipeline data
   - Check edge cases

## Getting Help

If you encounter issues:

1. Check `PIPELINE.md` for pipeline details
2. Check `tests/README.md` for testing details
3. Review `CLAUDE.md` for project context
4. Check git history for recent changes
5. Run individual scripts with verbose output to debug

## Checklist

Before submitting refactored code:

- [ ] All tests pass (`test_output_consistency.R` returns 0)
- [ ] Full pipeline runs successfully (`run_pipeline.R` completes)
- [ ] Code is more readable/maintainable than before
- [ ] Comments and documentation are updated
- [ ] No unnecessary files added
- [ ] Git commits are logical and well-described
- [ ] MOFA input remains unchanged (unless intentionally retrained)
- [ ] Visual inspection of key plots completed
- [ ] Performance hasn't significantly degraded

## Additional Resources

- `PIPELINE.md` - Detailed pipeline documentation
- `tests/README.md` - Testing framework documentation
- `CLAUDE.md` - Project overview and context
- `scripts/support_functions.R` - Core utility functions
- Original paper - See README.md for citation
