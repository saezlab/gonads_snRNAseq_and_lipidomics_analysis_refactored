# Verify public CSV preference, source-RDS bridges and unchanged intake values.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
source("scripts/core/pseudobulk_csv.R")
source("scripts/core/wdsd_intake.R")
source("scripts/core/wd_intake.R")

qa_pseudobulk_csv_consumers <- function(baseline_dir=NULL,output=NULL) {
  rows<-list()
  check<-function(label,value){rows[[length(rows)+1L]]<<-data.frame(check=label,passed=isTRUE(value));if(!isTRUE(value))stop(label)}
  # Source objects must never be opened, even while local originals are present.
  trace("readRDS",where=baseenv(),print=FALSE,tracer=quote({
    if(is.character(file)&&length(file)==1L&&
      grepl("(/data/RNA/|/inputs/(fcg|gdx|wdsd)/)",normalizePath(file,mustWork=FALSE)))
      stop("CSV smoke forbids opening a source pseudobulk RDS: ",file)
  }))
  on.exit(untrace("readRDS",where=baseenv()),add=TRUE)
  for(role in c("wdsd_current","wdsd_historical","wd")){
    x<-switch(role,wdsd_current=wp2_load_inputs(),
      wdsd_historical=wp2_load_inputs(input_zonation=analysis_path("historical_wdsd_zonation_transport")),wd=wd_load_inputs())
    check(paste(role,"all provenance inputs are real files"),all(file.exists(x$input_paths)))
    check(paste(role,"CSV preferred while local RDS exists"),!any(grepl("[.]rds$",x$input_paths,ignore.case=TRUE)))
    if(!is.null(baseline_dir)){
      old<-readRDS(file.path(baseline_dir,paste0("baseline_",role,".rds")))
      for(n in setdiff(names(old),"input_paths")){
        a<-old[[n]];b<-x[[n]]
        if(is.data.frame(a))b<-b[,names(a),drop=FALSE]
        check(paste(role,n,"identical preserved values"),isTRUE(all.equal(a,b,tolerance=0,check.attributes=TRUE)))
      }
    }
    rm(x);invisible(gc())
  }
  five <- new.env(parent = globalenv())
  sys.source("scripts/workflows/05_zonation_analysis.R", envir = five)
  selected_transports <- character()
  five$analysis_path <- function(key, ...) {
    selected_transports <<- c(selected_transports, key)
    analysis_path(key, ...)
  }
  check("Workflow 05 FCG CSV-only validation", isTRUE(five$fcg_zonation_run("--validate-only")))
  check("Workflow 05 GDX CSV-only validation", isTRUE(five$fcg_zonation_run(c("--cohort=gdx", "--validate-only"))))
  check("Workflow 05 directly selects both expanded-marker manifests",
    identical(selected_transports, c("current_fcg_zonation_transport", "current_gdx_zonation_transport")))
  source("scripts/workflows/17_expanded_marker_zonation_update.R",local=TRUE)
  check("Workflow17 CSV-only validation",isTRUE(zu17_run("--validate-only")))
  a<-readRDS(analysis_current("cf453_models/models.rds"))
  prior<-utils::read.csv(analysis_current("cf453_models/provenance_sha256.csv"))
  p<-unname(a$input_paths);expected<-prior$sha256[match(p,prior$path)]
  verification<-analysis_verify_inputs(p,expected)
  check("Frozen approved model sources verified without opening RDS",all(verification$verified))
  check("Frozen raw RDS identities use verified CSV bridges",all(verification$verification[grepl("[.]rds$",p,ignore.case=TRUE)]=="verified_lossless_pseudobulk_CSV_export"))
  bad<-tryCatch(analysis_verify_inputs(p[1L],paste(rep("0",64L),collapse="")),error=identity)
  check("Changed expected input identity is rejected",inherits(bad,"error"))
  ledger<-analysis_provenance(p)
  check("Original source RDS and opened CSV hashes remain distinct",all(c("recorded_source_RDS_identity","file_bytes")%in%ledger$sha256_scope))
  source("tests/qa_current_zonation_release.R",local=TRUE)
  check("Approved current scientific invariants accept CSV bridges",all(qa_current_zonation_release()$passed))
  result<-do.call(rbind,rows)
  if(!is.null(output)){dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE);utils::write.csv(result,output,row.names=FALSE)}
  cat(nrow(result),"CSV consumer, provenance and no-source-RDS checks passed\n")
  invisible(result)
}
if(sys.nframe()==0L){
  args<-commandArgs(TRUE)
  if(any(!grepl("^--(baseline|output)=.+",args)))stop("Use optional --baseline= and --output=")
  value<-function(key){z<-args[startsWith(args,paste0(key,"="))];if(length(z))sub("^[^=]*=","",z[1L]) else NULL}
  qa_pseudobulk_csv_consumers(value("--baseline"),value("--output"))
}
