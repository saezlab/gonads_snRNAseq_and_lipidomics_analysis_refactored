source("scripts/config/analysis_release.R")
# Three independent diet-rewiring panels. Read Workflow 16 plan before editing.
source("scripts/support_functions.R")
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/core/cf453_standalone.R")
source("scripts/core/cf453_diet_rewiring.R")
cfdr_run <- function(args=commandArgs(TRUE)) {
  allowed<-c("--models-only","--plot-only")
  values<-args[startsWith(args,"--output=")|startsWith(args,"--parent=")|startsWith(args,"--direct=")|startsWith(args,"--parent-provenance=")]
  if(any(!args %in% c(allowed,values))||anyDuplicated(sub("=.*$","",args)))stop("Unknown/duplicate CLI options")
  if(all(allowed %in% args))stop("Choose one execution mode")
  value<-function(key,default){x<-args[startsWith(args,paste0(key,"="))];if(length(x))sub("^[^=]*=","",x) else default}
  output<-value("--output",analysis_current("cf453_rewiring"))
  if("--plot-only" %in% args) a<-readRDS(file.path(output,"rewiring_analysis.rds")) else {
    if(dir.exists(output)&&length(setdiff(list.files(output,all.files=FALSE),"QA")))stop("Choose a fresh output directory")
    a<-cfdr_analyse(value("--parent",analysis_current("cf453_models/models.rds")),
      value("--direct",analysis_current("cf453_direct_diet/diet_analysis.rds")),
      provenance_path=value("--parent-provenance",NULL))
    cfdr_save(a,output)
  }
  if("--models-only" %in% args)return(invisible(a))
  source("scripts/core/cf453_diet_rewiring_plots.R")
  manifest<-cfdr_make_panels(a,output)
  cf6_write(manifest,file.path(output,"panel_manifest.csv"))
  source("scripts/core/cf453_diet_rewiring_report.R")
  cfdr_report(a,manifest,output);cfdr_provenance(a,output)
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  invisible(a)
}
if(sys.nframe()==0L) cfdr_run()
