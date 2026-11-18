# Compatibility entry for one historical phase of the work package.
# Use scripts/workflows/12_wdsd_network_comparison.R for the complete package.
# The matching legacy plan documents phase-specific defaults and outputs.
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) {
  .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
}
source("scripts/core/wp2_workflow.R")
if (sys.nframe() == 0L) wp2_run()
