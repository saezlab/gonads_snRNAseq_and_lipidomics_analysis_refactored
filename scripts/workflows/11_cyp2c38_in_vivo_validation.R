################################################################################
# Workflow 11: Cyp2c38 In Vivo Validation (WP1)
################################################################################
# Run from the repository root. Read docs/workflow_plans/11_cyp2c38_in_vivo_validation.md.
# Independent raw-input/reference checks: Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --validate-only
# Reproduction settings: --permutations=1999 --bootstraps=4999 --base-bootstraps=999 --workers=3
# Completed outputs are preserved; select a fresh --output for analysis reruns.
# Shared helper functions retain the original statistical implementations.
################################################################################

################################################################################
# 1. CONFIGURATION AND SHARED FUNCTIONS
################################################################################
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) {
  .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
}
source("scripts/core/wp1_workflow.R")

################################################################################
# 2. INTAKE AND REFERENCE REPLAY
# 3. WORK PACKAGE ANALYSES
# 4. REPORTS AND REPRODUCIBILITY EXPORTS
################################################################################
# The package runner executes these sections in the documented order and prints
# phase boundaries. CLI settings and phase provenance accompany the result tables.
if (sys.nframe() == 0L) wp1_package_run()
