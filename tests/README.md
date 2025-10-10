# Testing Framework for FCG Analysis Scripts

This directory contains the testing infrastructure for validating analysis scripts in this project.

## Overview

The testing framework is designed to:
1. Verify that analysis scripts run without errors
2. Capture outputs (dataframes, results) for validation
3. Track script execution logs for debugging
4. Provide a lightweight way to ensure reproducibility

## Directory Structure

```
tests/
├── README.md              # This file
├── outputs/               # Test outputs (CSV files from analysis)
├── logs/                  # Execution logs for each test run
├── expected_outputs/      # Reference outputs for comparison (optional)
├── run_tests.R            # Main test runner script
└── test_config.csv        # Configuration for which scripts to test
```

## Testing Philosophy

### Core Analysis Scripts
These scripts generate results files (RData, CSV) and should be tested by:
- Running the script
- Checking that expected output files exist
- Optionally comparing output dimensions/structure to expected values

Examples: `pseudobulk_DA.R`, `lipidomic_DA.R`, `MOFA_lipidcell_RNA.R`

### Visualization Scripts
These scripts generate plots (PDF files). To enable testing, they have been modified to:
- Save underlying dataframes to CSV files in `tests/outputs/`
- Continue generating plots as before
- Only save test outputs when `TESTING_MODE` environment variable is set

The modification is minimal and looks like:
```r
# At the end of the plotting code, before or after ggsave()
if (Sys.getenv("TESTING_MODE") == "TRUE") {
  write_csv(plot_data, "tests/outputs/scriptname_plotdata.csv")
}
```

Examples: `lipid_structure_visualisation.R`, `figure_1_UMAP.R`

## Running Tests

### Run All Tests
```r
source("tests/run_tests.R")
```

### Run Specific Script
```r
source("tests/run_tests.R")
test_single_script("scripts/pseudobulk_DA.R")
```

### Set Testing Mode
Before running visualization scripts individually:
```r
Sys.setenv(TESTING_MODE = "TRUE")
source("scripts/lipid_structure_visualisation.R")
```

## Test Configuration

The `test_config.csv` file specifies which scripts to test and their expected outputs:

| script_path | test_type | expected_output | timeout_min |
|-------------|-----------|-----------------|-------------|
| scripts/pseudobulk_DA.R | analysis | results/pseudobulk/ttop_full_list.RData | 10 |
| scripts/lipid_structure_visualisation.R | visualization | tests/outputs/lipid_structure_vis_*.csv | 5 |

## Interpreting Results

### Successful Test
- Script runs without errors
- Expected outputs are generated
- Log file shows "PASS"

### Failed Test
- Script throws an error (captured in log)
- Expected outputs missing
- Log file shows "FAIL" with error message

## Adding New Tests

To add a test for a new script:

1. **For analysis scripts:** Add entry to `test_config.csv` with expected output path
2. **For visualization scripts:**
   - Add dataframe saving code at the end of the script
   - Add entry to `test_config.csv`
   - Run once to generate reference output

## Test Output Files

### Analysis Scripts
Outputs are saved to their normal locations:
- `results/pseudobulk/`
- `results/lipidomic/`
- `results/mofa/`
- etc.

### Visualization Scripts
Dataframe outputs are saved to:
- `tests/outputs/[script_name]_[data_description].csv`

Examples:
- `tests/outputs/lipid_structure_vis_cyp2c38_data.csv`
- `tests/outputs/figure1_umap_coords.csv`
- `tests/outputs/pseudobulk_da_volcano_data.csv`

## Modified Scripts

The following visualization scripts have been modified to save dataframes when `TESTING_MODE=TRUE`:

1. **`scripts/lipid_structure_visualisation.R`**
   - Saves plot data for each gene's structural association plots
   - Output: `tests/outputs/lipid_structure_vis_[gene]_data.csv`
   - Contains: lipid, lipid_class, chain_length, db_position, fold_change, gene_correlation

2. **`scripts/lipidomic_bubbles.R`**
   - Saves bubble plot data for each lipid class
   - Output: `tests/outputs/lipidomic_bubbles_[lipid_class]_data.csv`
   - Contains: lipid structural features (C, DB, LogFC) for each lipid class

(Additional scripts can be modified following the same pattern)

## Notes

- **Non-invasive:** Test modifications only activate when `TESTING_MODE=TRUE`
- **Lightweight:** No heavy testing dependencies (testthat, etc.)
- **Practical:** Focuses on "does it run" rather than statistical correctness
- **Flexible:** Can be expanded to include output validation, dimension checks, etc.

## Future Enhancements

Potential improvements:
- Compare outputs to reference files
- Check dataframe dimensions and column names
- Validate statistical results (e.g., number of significant genes)
- Add performance benchmarking
- Integrate with CI/CD pipeline
