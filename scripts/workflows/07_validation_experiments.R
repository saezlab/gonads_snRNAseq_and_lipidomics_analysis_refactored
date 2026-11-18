# Current expanded-marker zonation validation; global/KD results remain unchanged.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
source("scripts/core/zonation_update_fcg.R")
zonation_validation_run <- function(args=commandArgs(TRUE)) {
  if(any(!grepl("^--output=.+|^--validate-only$|^--global-only$",args)))stop("Use --output=, --global-only or --validate-only")
  x<-args[startsWith(args,"--output=")]
  output<-if(length(x))sub("^--output=","",x[1]) else analysis_current()
  if("--global-only"%in%args){
    source("scripts/core/validation_global.R")
    if(!length(x)&&!"--validate-only"%in%args)stop("Use --global-only with a fresh --output=<directory>")
    return(validation_global_run(output,validate_only="--validate-only"%in%args))
  }
  if(!"--validate-only"%in%args)analysis_assert_not_approved_output(output)
  required<-c(file.path(output,"fcg05/fcg05_checkpoint.rds"),
    "results/validation/knockdown_t_table.csv","results/validation/knockdown_ttop_list.RData")
  zmu_require();stopifnot(all(file.exists(required)))
  if("--validate-only"%in%args){cat("Current zonation reference and unchanged KD fits validated; no files changed.\n");return(invisible(TRUE))}
  if(file.exists(file.path(output,"fcg07/fcg07_checkpoint.rds")))
    stop("Preserve approved validation: use a fresh --output with an FCG05 checkpoint, or the explicit revision-audit postprocess command")
  zmu_run07(output)
}
if(sys.nframe()==0L)zonation_validation_run()
