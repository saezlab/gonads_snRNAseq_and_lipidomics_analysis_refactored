# Current expanded-marker zonation validation; global/KD results remain unchanged.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
source("scripts/core/zonation_update_fcg.R")
zonation_validation_run <- function(args=commandArgs(TRUE)) {
  if(any(!grepl("^(--(output|parent)=.+|--validate-only|--global-only)$",args)) ||
      anyDuplicated(sub("=.*$", "", args)))
    stop("Use --output=, --parent= for zonal validation, --global-only or --validate-only")
  x<-args[startsWith(args,"--output=")]
  parent <- args[startsWith(args, "--parent=")]
  output<-if(length(x))sub("^--output=","",x[1]) else analysis_current()
  validate_only<-"--validate-only"%in%args
  if("--global-only"%in%args){
    if(length(parent))stop("--parent is used only for zonal validation")
    if(!length(x)&&!validate_only)stop("Use --global-only with a fresh --output=<directory>")
    source("scripts/core/validation_global.R")
    return(validation_global_run(output,validate_only=validate_only))
  }
  if(!validate_only)analysis_assert_not_approved_output(output)
  checkpoint <- if (length(parent)) sub("^--parent=", "", parent) else
    if (length(x)) file.path(output, "fcg05/fcg05_checkpoint.rds") else analysis_current("fcg05/fcg05_checkpoint.rds")
  reference <- c(checkpoint, zmu_zonal_network_paths(checkpoint))
  if (!all(file.exists(reference)))
    stop("Missing Workflow 05 FCG reference: ", paste(reference[!file.exists(reference)], collapse = ", "),
      ". Run Workflow 05 into this output root or pass --parent=PATH to a completed FCG05 checkpoint.")
  required<-c(
    "results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv","results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_ttop_list.RData")
  zmu_require();stopifnot(all(file.exists(required)))
  if(validate_only){cat("Current zonation reference and unchanged KD fits validated; no files changed.\n");return(invisible(TRUE))}
  if(file.exists(file.path(output,"fcg07/fcg07_checkpoint.rds")))
    stop("Preserve approved validation: use a fresh --output with an FCG05 checkpoint, or the explicit revision-audit postprocess command")
  cat("Using Workflow 05 reference:", normalizePath(checkpoint, mustWork = TRUE), "\n")
  zmu_run07(output, checkpoint_path = checkpoint)
}
if(sys.nframe()==0L)zonation_validation_run()
