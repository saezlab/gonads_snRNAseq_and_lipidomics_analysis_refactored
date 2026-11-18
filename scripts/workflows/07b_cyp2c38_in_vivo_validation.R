# Compatibility entry for one historical phase of the work package.
# Use scripts/workflows/11_cyp2c38_in_vivo_validation.R for the complete package.
# The matching legacy plan documents phase-specific defaults and outputs.
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) {
  .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
}
source("scripts/core/wp1_workflow.R")
if (sys.nframe() == 0L) wp1_run()
