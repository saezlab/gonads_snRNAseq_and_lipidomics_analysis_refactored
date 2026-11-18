################################################################################
# Workflow 12: WDSD Pooled and Zonation Network Comparison (WP2)
################################################################################
# Run from the repository root. Read docs/workflow_plans/12_wdsd_network_comparison.md.
# Independent raw-input/reference checks: Rscript scripts/workflows/12_wdsd_network_comparison.R --validate-only
# Reproduction settings: --permutations=19999 --bootstraps=4999 --paired-bootstraps=999
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
source("scripts/core/wp2_workflow.R")

################################################################################
# 2. INTAKE AND REFERENCE REPLAY
# 3. WORK PACKAGE ANALYSES
# 4. REPORTS AND REPRODUCIBILITY EXPORTS
################################################################################
# The package runner executes these sections in the documented order and prints
# phase boundaries. CLI settings and phase provenance accompany the result tables.
if (sys.nframe() == 0L) wp2_package_run()
