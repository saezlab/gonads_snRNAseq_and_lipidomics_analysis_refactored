source("scripts/config/analysis_release.R")
source("scripts/support_functions.R")
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/core/zonation_update_gdx.R")
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(any(!startsWith(args,"--output=")&args!="--science-only"))stop("Use --output=PATH and optional --science-only")
  path<-args[startsWith(args,"--output=")]
  zu_gdx_run(if(length(path))sub("^--output=","",path[1]) else analysis_current(),science_only="--science-only"%in%args)
}
