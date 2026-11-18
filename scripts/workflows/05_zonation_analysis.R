# Current expanded-marker FCG zonation; frozen historical engine lives in core/engines.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
source("scripts/core/zonation_update_fcg.R")
fcg_zonation_run <- function(args=commandArgs(TRUE)) {
  if(any(!grepl("^--(output|input)=.+|^--validate-only$",args)))stop("Use --output=, --input= or --validate-only")
  value<-function(key,default){x<-args[startsWith(args,paste0(key,"="))];if(length(x))sub("^[^=]*=","",x[1]) else default}
  input<-value("--input",dirname(analysis_path("current_fcg_zonation_csv")))
  output<-value("--output",analysis_current())
  if("--validate-only"%in%args){
    zmu_require();zmu_validate(input)
    invisible(parse(text=zmu_runtime05(input,file.path(output,"fcg05"))))
    cat("Current expanded-marker FCG inputs and isolated engine validated; no files changed.\n")
    return(invisible(TRUE))
  }
  zmu_run05(input,output)
}
if(sys.nframe()==0L)fcg_zonation_run()
