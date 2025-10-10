# Testing Framework Setup Summary

## What Was Created

A comprehensive testing framework has been created to support safe refactoring of this codebase. The framework ensures that code changes don't inadvertently alter scientific results.

### Files Created

1. **`tests/save_reference_outputs.R`**
   - Saves checksums and statistics from current analysis outputs
   - Creates baseline for comparison
   - **Status:** ✅ Executed successfully

2. **`tests/test_output_consistency.R`**
   - Compares current outputs against saved reference
   - Verifies numerical consistency
   - Reports any differences

3. **`tests/reference_outputs/reference_checksums.rds`**
   - Binary R object containing all reference data
   - **Status:** ✅ Created (16 files tracked)

4. **`tests/reference_outputs/reference_checksums.json`**
   - Human-readable JSON summary of reference data
   - Useful for reviewing what's tracked

5. **`tests/reference_outputs/mofa_input_sample.csv`**
   - Sample of MOFA input data for inspection
   - Helps verify MOFA input format

6. **`tests/README.md`**
   - Documentation for the testing framework
   - Usage instructions and troubleshooting

7. **`PIPELINE.md`**
   - Comprehensive pipeline documentation
   - Execution order for all analysis scripts
   - Dependencies and outputs for each step

8. **`REFACTORING_GUIDE.md`**
   - Step-by-step guide for safe refactoring
   - Best practices and common patterns
   - Examples of refactoring tasks

9. **`run_pipeline.R`**
   - Script to run the complete analysis pipeline
   - Supports running individual steps or ranges
   - Provides progress reporting

## Files Tracked by Testing Framework

The testing framework monitors 16 key output files:

### RNA Analysis (7 files)
- `data/RNA/scaled_RNA_1.csv` - Scaled RNA expression (samples 1-79)
- `data/RNA/scaled_RNA_2.csv` - Scaled RNA expression (samples 80-159)
- `data/RNA/scaled_RNA_3.csv` - Scaled RNA expression (samples 160-215)
- `data/RNA/meta_data_df.csv` - Sample metadata
- `results/pseudobulk/ttop_full_list.RData` - Differential expression results
- `results/pseudobulk/tval_full_df.csv` - T-values for all genes
- `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` - PCA coordinates

### Lipidomics Analysis (6 files)
- `data/metabolomics/lipidomic.RData` - Processed lipidomics matrix
- `results/lipidomic/gonade_type_ttop.csv` - Full differential abundance results
- `results/lipidomic/gonade_type_ttop_20.csv` - Top 20 significant lipids
- `results/lipidomic/PCA_sample_coord_lipid.csv` - PCA coordinates
- `results/lipidomic/factor_enrichment_long.csv` - PC enrichment
- `results/lipidomic/ttop_full_list.RData` - All contrasts

### Multi-Omics Integration (2 files)
- `data/mofa/MOFA_RNA_data.csv` - **CRITICAL: MOFA input (must not change!)**
- `results/mofa/main_model_RNA_lipid_weight.RData` - MOFA weights

### Pathway Analysis (1 file)
- `results/pseudobulk/pathways/pathway_scores_df_all.RData` - Pathway scores

## How to Use the Testing Framework

### Before Refactoring

```bash
# Save current outputs as reference (already done!)
Rscript tests/save_reference_outputs.R
```

✅ **This has already been executed successfully.**

### During Refactoring

```bash
# 1. Make your code changes
vim scripts/pseudobulk_DA.R  # or any other script

# 2. Re-run the analysis
Rscript run_pipeline.R

# 3. Test consistency
Rscript tests/test_output_consistency.R
```

### If Tests Fail

Review the output to see what changed:
- **File checksum mismatch** - File content changed
- **Content checksum mismatch** - Data values changed
- **Dimension mismatch** - Number of rows/columns changed
- **MOFA input changed** - Critical issue, MOFA models need retraining

Options:
1. Fix the bug in your refactored code
2. If change is expected (e.g., bug fix), document and update reference
3. Investigate numerical precision issues

## Pipeline Execution Order

The main analysis pipeline consists of 5 steps:

1. **seurat_to_df.R** - RNA data preparation (Seurat → CSV)
2. **pseudobulk_DA.R** - RNA differential analysis (limma)
3. **lipidomic_DA.R** - Lipidomics differential analysis
4. **pseudobulk_pathways.R** - Pathway enrichment
5. **MOFA_lipidcell_RNA.R** - Multi-omic integration (loads pre-trained models)

Run the complete pipeline:
```bash
Rscript run_pipeline.R
```

Run specific steps:
```bash
Rscript run_pipeline.R --step 2      # Only RNA differential analysis
Rscript run_pipeline.R --from 3      # Start from lipidomics
Rscript run_pipeline.R --to 4        # Stop before MOFA
Rscript run_pipeline.R --skip-mofa   # Skip MOFA (faster for testing)
```

## Special Considerations

### MOFA (Multi-Omics Factor Analysis)

⚠️ **CRITICAL:** The MOFA input data (`data/mofa/MOFA_RNA_data.csv`) must remain unchanged.

- MOFA model training is computationally expensive (~hours)
- Pre-trained models are stored in `results/mofa/RNA/mofa_res_*factor.hdf5`
- `MOFA_lipidcell_RNA.R` loads these pre-trained models (no retraining)
- If MOFA input changes, all models need to be retrained

**What you can refactor:**
- Code that prepares MOFA input (verify input stays identical)
- Code that analyzes MOFA outputs (after model loading)
- Visualization code

**What you cannot easily refactor:**
- MOFA training code (commented out in script)
- MOFA input data structure

### Support Functions

`scripts/support_functions.R` is used by many scripts. Changes affect the entire pipeline.

**Strategy:**
- Refactor one function at a time
- Test the complete pipeline after changes
- Consider versioning functions if major changes are needed

**Key functions:**
- `runLimma()` - Differential analysis wrapper
- `volcano_nice()` - Volcano plot generation
- `nicePCA()` - PCA visualization
- `createLinearColors()` - Heatmap color palettes

## What Tests Do NOT Cover

The testing framework only verifies numerical outputs. It does NOT test:
- Plots and figures (only the data behind them)
- PDF outputs
- Console messages
- Performance/timing
- Memory usage

Visual inspection is still required for plots and figures.

## Reference Statistics

From the initial reference capture:

- **Total files tracked:** 16
- **CSV files with content checksums:** 7
- **RData files tracked:** 5
- **MOFA input:**
  - Rows: 297,980
  - Samples: 215
  - Views: 6 (Hepatocytes, Endothelial, Kupffer, LSECs, Stellate, Lipidomic)
  - Features: 1,386 (lipids) + variable (RNA genes by cell type)

## Next Steps for Refactoring

Now that the testing framework is in place, you can proceed with refactoring:

### Suggested Refactoring Order

1. **Start with visualization scripts** (safest, don't affect data)
   - Figure generation
   - Plotting code

2. **Downstream analysis scripts**
   - `pseudobulk_pathways.R`
   - `MOFA_lipidcell_downstream.R`

3. **Core analysis scripts**
   - `pseudobulk_DA.R`
   - `lipidomic_DA.R`

4. **Data preparation scripts** (most risky)
   - `seurat_to_df.R`

5. **Support functions** (after testing downstream)
   - `support_functions.R`

### Common Refactoring Tasks

See `REFACTORING_GUIDE.md` for detailed examples of:
- Extracting repeated code into functions
- Standardizing file paths
- Adding progress indicators
- Adding input validation
- Modularizing large scripts

## Documentation

All documentation is now in place:

- **[README.md](README.md)** - Project overview
- **[CLAUDE.md](CLAUDE.md)** - Comprehensive project guide for AI assistants
- **[PIPELINE.md](PIPELINE.md)** - Pipeline execution order and details
- **[REFACTORING_GUIDE.md](REFACTORING_GUIDE.md)** - Refactoring best practices
- **[tests/README.md](tests/README.md)** - Testing framework documentation

## Quick Reference Commands

```bash
# Run complete pipeline
Rscript run_pipeline.R

# Run specific step
Rscript run_pipeline.R --step 2

# Save reference (before refactoring)
Rscript tests/save_reference_outputs.R

# Test consistency (after refactoring)
Rscript tests/test_output_consistency.R

# Get help
Rscript run_pipeline.R --help
```

## Status

✅ **Testing framework is complete and operational**

- Reference outputs saved successfully
- 16 files tracked
- All documentation created
- Ready for refactoring to begin

You can now safely refactor the codebase with confidence that any changes to scientific results will be detected and reported.
