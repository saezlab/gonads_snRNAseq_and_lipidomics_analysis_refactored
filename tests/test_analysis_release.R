# Registry, archived identities and both current/historical CF453 cache linkage.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
registry<-analysis_release()
stopifnot(identical(registry$release_id,"expanded_markers_20261116"))
public_keys<-c("current_root","engine_fcg05","engine_validation07",grep("_transport$",names(registry$paths),value=TRUE))
for(key in public_keys)stopifnot(file.exists(analysis_path(key)))
stopifnot(identical(analysis_resolve_recorded_path("results/zonation_marker_update_20261116"),analysis_current()),
  identical(analysis_resolve_recorded_path("/prior/checkout/results/zonation/networks"),file.path(analysis_path("historical_fcg05"),"networks")))
stopifnot(digest::digest(analysis_path("engine_fcg05"),file=TRUE,algo="sha256")=="c560572332047ce3ce0cb065662d32a112241bcc2d46600824b424232ae25410",
  digest::digest(analysis_path("engine_validation07"),file=TRUE,algo="sha256")=="730420ace073e6bd94ce2e15ca7b90d063b082cb9ebf73c58fa0d68f1e26b66b")
args<-commandArgs(TRUE)
if(any(args!="--with-historical"))stop("Use optional --with-historical for restored archived model caches")
roots<-analysis_current("cf453_models")
if("--with-historical"%in%args)roots<-c(roots,analysis_path("historical_cf453_models"))
for(root in roots){
  a<-readRDS(file.path(root,"models.rds"));p<-utils::read.csv(file.path(root,"provenance_sha256.csv"))
  recorded<-unname(a$input_paths);expected<-p$sha256[match(recorded,p$path)]
  actual<-analysis_verify_inputs(recorded,expected)
  stopifnot(!anyNA(expected),identical(recorded,actual$path),all(actual$verified))
  cat("Verified raw-input bytes and recorded provenance identity:",root,"\n")
}
cat("Current release registry, frozen engines, relocations and CF453 cache linkage passed.\n")
