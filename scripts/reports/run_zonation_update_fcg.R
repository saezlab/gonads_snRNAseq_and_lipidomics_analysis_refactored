#!/usr/bin/env Rscript
if (dir.exists(".local/R-library")) .libPaths(c(".local/R-library", .libPaths()))
source("scripts/core/zonation_update_fcg.R")
zonation_update_fcg_run <- function(args=commandArgs(trailingOnly=TRUE)) {
value <- function(key, default) {
  x <- args[startsWith(args, paste0(key, "="))]
  if (length(x)) sub(paste0("^", key, "="), "", x[1]) else default
}
output <- value("--output", analysis_current())
input <- value("--input", if(file.exists(file.path(output,"inputs/fcg/scaled_RNA_zonation.csv"))) file.path(output,"inputs/fcg") else dirname(analysis_path("current_fcg_zonation_csv")))
if ("--validate-only" %in% args) {
  zmu_require(); zmu_validate(input)
  invisible(parse(text = zmu_runtime05(input, file.path(output, "fcg05"))))
  cat("Versioned FCG input/runtime validation passed; no plots or models generated.\n")
} else if ("--postprocess-only" %in% args) {
  zmu_finish(output,rerun07="--rerun-zonal-validation" %in% args)
} else {
  analysis_assert_not_approved_output(output)
  zmu_run05(input, output)
  zmu_run07(output)
  if("--science-only"%in%args)zmu_science_finish(output) else zmu_finish(output)
}

}
if(sys.nframe()==0L)zonation_update_fcg_run()
