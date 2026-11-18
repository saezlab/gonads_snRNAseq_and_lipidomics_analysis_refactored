################################################################################
# Workflow 13: Western-Diet Cross-Cohort Robustness (WP3)
################################################################################
# Run from repository root. Read docs/workflow_plans/13_wd_crosscohort_robustness.md.
# Rscript scripts/workflows/13_wd_crosscohort_robustness.R --validate-only
# Defaults: seed20261014; 4999 marginal null draws; 999 paired stability draws;
# exhaustive 24 CF459 and 5040 CF453 mouse-pairing permutations.
# Full runs include the original male FCG SD anchor and direct three-cohort
# closeness stability. Full runs require a fresh --output. --report-only refreshes
# report Markdown/extension figures from saved data. Add the extension to a saved
# base run with --fcg-comparison-only [--comparison-output=FRESH_PATH].
################################################################################

################################################################################
# 1. CONFIGURATION AND SHARED FUNCTIONS
################################################################################
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/wd_crosscohort_workflow.R")

################################################################################
# 2. INTAKE, FIXED REFERENCES AND STRICT MATCHING
# 3. CROSS-COHORT NETWORKS, EXACT NULLS AND STABILITY
# 4. PRESPECIFIED QC AND RNA-LAYER SENSITIVITIES
# 5. ORIGINAL MALE FCG SD ANCHOR AND DIRECT CLOSENESS
# 6. SCIENTIFIC FIGURES, REPORT AND PROVENANCE
################################################################################
# The runner executes these sections in order; source() defines functions only.
if (sys.nframe() == 0L) wp3_run()
