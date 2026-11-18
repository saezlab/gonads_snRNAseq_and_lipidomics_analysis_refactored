# Current release and organization checks

The expanded-marker default is selected by `scripts/config/analysis_release.json`. Current saved-result checks use the approved release; older WP2/WP3 tests explicitly retain their original two-marker references. Do not enforce the former two-marker nucleus gate on new multi-gene assignments.

Run `python3 scripts/maintenance/verify_repository_organization.py --require-full-archive` for local migration/input/result hashes and `python3 scripts/maintenance/audit_github_readiness.py --strict` for the proposed commit. Cleanup smoke outputs are under ignored `miscellaneous/local/repository_cleanup/`. Old generated test logs and historical TEST_RESULTS reports are under `miscellaneous/quick_tests/logs/`; test implementations remain here.

---

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

## Cyp2c38 In Vivo Extension

The independent Cyp2c38 in vivo extension has focused statistical and identity
checks: `Rscript tests/test_cyp2c38_in_vivo.R`. It verifies known simulated
effects with unequal male group sizes, equal-sex pooled contrasts, interaction
and perturbation signs, precision-weight invariance of effects, mouse-order
alignment and molecular-annotation distinctions. Run
`Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --validate-only`
for actual input/reference replay checks. The full workflow additionally checks
the observed contrast identities and one-to-one lipid overlap.

The ULM follow-up has two additional focused checks:
`Rscript tests/test_cyp2c38_ulm_bootstrap.R` verifies native decoupleR ULM
equivalence, constrained weighted null means/leverage, nuisance dimensions,
shared whole-mouse multipliers and reproducible paired reference resampling.
`Rscript tests/test_cyp2c38_sex_networks.R` independently checks the full
feature universe, direct within-gonad Pearson coefficients, pairing/order,
historical replay and six network exports. Workflow 11 also verifies
unchanged KO t-statistics, all observed ULM scores and strict 856/552 joins.

## WDSD Network Comparison (WP2)

Run `Rscript tests/test_wdsd_intake.R` for updated CF453 labels, Seurat layer
alignment and aggregate marker/scaling checks; `Rscript tests/test_wdsd_reference.R`
for independent original Pearson replay and branch selection;
`Rscript tests/test_wdsd_comparison_setup.R` for strict molecular joins, frozen
eligibility and intentional mispairing failures; and
`Rscript tests/test_wdsd_network_statistics.R` for native ULM/decoupleR equivalence,
all 5,040 exact pairings, shared whole-mouse null/stability draws, original
genotype strata and all 49 SD/WD mouse omissions. Workflow 12's
`--validate-only` checks the actual inputs without writing analysis outputs.

## Pooled Figure 3 and Molecular Overlap (WP2 Follow-up)

`Rscript tests/test_wdsd_lipid_overlap_audit.R` verifies annotation resolution,
literal alias enumeration, adduct-corrected masses and global collision guards.
`Rscript tests/test_wdsd_pooled_reference.R` independently replays all archived
Figure 3 species and precursor coefficients, top-three membership and class
compression. `Rscript tests/test_wdsd_pooled_statistics.R` checks native ULM,
complete mouse strata/draws, identity audit and conditional null behavior with
a pure between-sex control. Workflow 12's `--validate-only` checks
actual inputs. `tests/qa_wdsd_pooled_actual.R` independently checks completed
output coefficients, ULM, all null tails and saved mouse-draw alignment.

## Work Package Consolidation and Local Inputs

Run `Rscript tests/test_wp1_wp2_local_inputs.R` for imported input hashes,
`Rscript tests/test_work_package_workflows.R` for phase order, CLI routing,
fresh-input validation and overwrite guards, and the focused statistics tests
above for model invariants. Both new entry points passed actual `--validate-only`.
A complete Workflow 11 run with reduced resampling counts tested every phase
in a fresh temporary directory; its observed differential and native ULM
outputs reproduce the existing results. Workflow 12 is also tested end to end
with the production resampling counts. These runs do not replace saved primary
results. The 12 WP2 inference/coefficient/sensitivity tables reproduce within
1e-10, including empirical p-values, and 16 primary resampling/index/LOCO
banks match exactly. Repeat these actual-data checks after a fresh run:

```sh
Rscript tests/test_work_package_workflows.R --wp1-output=/path/to/wp1_run --wp2-output=/path/to/wp2_run
```

WP1 replay checks compare observed results, so reduced-draw smoke outputs are
suitable for that check; WP2 inference replay requires the original default
settings. Migration/refactor evidence is under `tests/outputs/wp1_wp2_migration/`;
PDF preservation/layout QA is under `results/wdsd_network_comparison/QA/`.

## Western-Diet Cross-Cohort Robustness (WP3)

`Rscript tests/test_wd_intake.R` checks verified source hashes, decimal-comma/BOM
import, actual Seurat layers, mouse/sex/zone alignment and deliberate failures.
`Rscript tests/test_wd_crosscohort_setup.R` checks three-way strict identity,
ambiguous-candidate exclusion, independent direct Pearson values, frozen
membership and independence from the unanalysed CF453 SD arm.
`Rscript tests/test_wd_crosscohort_statistics.R` checks native ULM against
OLS/decoupleR, exact24/5040 index banks/tails, n=4 degenerate draws, shared indices,
paired stability, all28+11 omissions and seeded replay.
`Rscript tests/test_wd_crosscohort_workflow.R` checks CLI/source guards and output
protection. Actual-input validation is the new Workflow 13 `--validate-only`.

Independent completed-output verification recomputes all coefficients and exact
permutation banks with base cor/lm, plus selected bootstrap and all omission draws:

```sh
Rscript tests/qa_wd_crosscohort_actual.R results/wd_crosscohort_robustness/analysis.rds results/wd_crosscohort_robustness/contexts.rds
```

`Rscript tests/test_wd_fcg_comparison.R` tests the original FCG anchor,
identical-edge delta, all 24 CF459 pairings, independent genotype-stratified
FCG seed stream, shared CF459 paired draw, 61 omissions and exact contingency
certification of all-zero null vectors. Weak varying vectors remain eligible.
The default CLI extension replay is identical to production; the complete
19-draw workflow smoke also reproduces observed and exact comparison outputs.

The extension's independent verifier reconstructs FCG directly from original
CSVs and checks every exact/null/paired draw, mouse index and omission, including
the corrected historical CF453-anchor null bank:

```sh
Rscript tests/qa_wd_fcg_comparison_actual.R results/wd_crosscohort_robustness/fcg_comparison/analysis.rds results/wd_crosscohort_robustness/contexts.rds results/wd_crosscohort_robustness/analysis.rds
```

Production verification passes 1,642 checks (maximum discrepancy 1.73e-12).
The original 500-check audit did not detect the exact-zero marginal-null
roundoff issue; historical files remain archived and separately corrected
extension tables supersede their bootstrap p-values for interpretation.

## Future Enhancements

Potential improvements:
- Compare outputs to reference files
- Check dataframe dimensions and column names
- Validate statistical results (e.g., number of significant genes)
- Add performance benchmarking
- Integrate with CI/CD pipeline
