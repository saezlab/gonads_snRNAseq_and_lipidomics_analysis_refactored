# Changed explicit inputs must not inherit a same-basename export identity.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
test_analysis_provenance_identity <- function(){
  manifest<-analysis_path("wd_zonation_transport",mustWork=TRUE)
  m<-jsonlite::read_json(manifest,simplifyVector=FALSE)
  directory<-tempfile("changed_pseudobulk_identity_");dir.create(directory)
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  fixture<-file.path(directory,basename(m$source$path))
  saveRDS(list(changed=TRUE,fixture="different explicit source bytes"),fixture)
  actual_sha<-digest::digest(fixture,file=TRUE,algo="sha256")
  stopifnot(!identical(actual_sha,m$source$sha256))
  trace("readRDS",where=baseenv(),print=FALSE,tracer=quote(stop("Identity test must not open an RDS object")))
  on.exit(untrace("readRDS",where=baseenv()),add=TRUE)
  stopifnot(is.null(analysis_pseudobulk_bridge(fixture)))
  provenance<-analysis_provenance(fixture)
  stopifnot(nrow(provenance)==1L,provenance$path==fixture,
    provenance$resolved_path==normalizePath(fixture),provenance$sha256==actual_sha,
    provenance$sha256_scope=="file_bytes",provenance$verification=="direct_file_bytes")
  # An explicitly expected frozen identity still selects its verified CSV,
  # independent of the changed local file at the recorded input name.
  bound<-analysis_verify_inputs(fixture,m$source$sha256)
  stopifnot(bound$verified,bound$source_RDS_sha256==m$source$sha256,
    bound$verification=="verified_lossless_pseudobulk_CSV_export",
    bound$opened_path==normalizePath(manifest),bound$actual_file_sha256!=actual_sha,
    digest::digest(fixture,file=TRUE,algo="sha256")==actual_sha)
  cat("Changed same-basename input retains actual-byte provenance; expected-SHA CSV verification remains source-RDS-free.\n")
  invisible(TRUE)
}
test_analysis_provenance_identity()
