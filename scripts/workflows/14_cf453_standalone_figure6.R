source("scripts/config/analysis_release.R")
# Standalone CF453 Figure 6. See docs/workflow_plans/14_cf453_standalone_figure6.md.
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/wdsd_intake.R")
source("scripts/core/cf453_standalone.R")

cf6_run <- function(args = commandArgs(TRUE)) {
  allowed<-c("--models-only","--finish","--plot-only")
  values<-args[startsWith(args,"--output=")|startsWith(args,"--bootstraps=")|startsWith(args,"--seed=")|startsWith(args,"--input-zonation=")]
  if(any(!args %in% c(allowed,values)) || anyDuplicated(sub("=.*$","",args))) stop("Unknown or duplicate CLI argument")
  mode<-args[args %in% allowed];if(length(mode)>1L) stop("Select one mode")
  if(any(startsWith(args,"--input-zonation=")) && any(mode %in% c("--finish","--plot-only")))
    stop("--input-zonation applies to a fresh model fit, not a saved-result mode")
  value<-function(key,default) {x<-args[startsWith(args,paste0(key,"="))];if(length(x)) sub("^[^=]*=","",x) else default}
  output<-value("--output",analysis_current("cf453_models"))
  b<-as.integer(value("--bootstraps","999"));seed<-as.integer(value("--seed","20261030"))
  if(is.na(b)||b<19L||is.na(seed))stop("Invalid resampling settings")
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  if("--plot-only" %in% mode) {
    a<-readRDS(file.path(output,"analysis.rds"))
  } else {
    if("--finish" %in% mode) {
      if(file.exists(file.path(output,"analysis.rds")))stop("Preserve completed analysis; --finish requires models without completed analysis")
      a<-readRDS(file.path(output,"models.rds"))
    } else {
      if(file.exists(file.path(output,"models.rds"))||file.exists(file.path(output,"analysis.rds")))stop("Preserve completed results: choose a fresh --output")
      a<-cf6_models(seed=seed,bootstraps=b,input_zonation=value("--input-zonation",NULL))
      cf6_save_models(a,output)
    }
    if("--models-only" %in% mode)return(invisible(a))
    a<-cf6_networks(a);a<-cf6_density_summaries(a);a<-cf6_representative_draws(a);cf6_save(a,output)
  }
  source("scripts/core/cf453_figure6_plots.R")
  manifest<-cf6_make_panels(a,output)
  cf6_write(manifest,file.path(output,"panel_manifest.csv"))
  if(file.exists("scripts/core/cf453_figure6_report.R")) {
    source("scripts/core/cf453_figure6_report.R")
    cf6_report(a,manifest,output)
  }
  invisible(a)
}
if(sys.nframe()==0L) cf6_run()
