# Current expanded-marker scientific pipeline; optional retrospective manuscript audit.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
zu17_run <- function(args=commandArgs(TRUE)) {
  flags<-c("--validate-only","--report-only","--with-manuscript-audit")
  values<-args[grepl("^--(output|source-root|python)=",args)]
  if(any(!args%in%c(flags,values))||anyDuplicated(sub("=.*$","",args))||all(c("--validate-only","--report-only")%in%args))stop("Unknown, duplicate or mutually exclusive CLI options")
  value<-function(key,default){v<-args[startsWith(args,paste0(key,"="))];if(length(v))sub("^[^=]*=","",v[1]) else default}
  output<-value("--output",analysis_current())
  source_root<-value("--source-root",analysis_path("pseudobulk_csv_root"))
  staging<-new.env(parent=globalenv());sys.source("scripts/reports/stage_zonation_marker_update.R",envir=staging)
  py<-value("--python",Sys.getenv("ZONATION_PYTHON",unset=Sys.which("python3")))
  audit<-any(c("--with-manuscript-audit","--report-only")%in%args)
  if("--validate-only"%in%args){
    source("scripts/core/pseudobulk_csv.R",local=TRUE)
    for(cohort in c("fcg","gdx","wdsd")){
      input<-staging$zu_source_input(cohort,source_root);x<-pb_components(input)
      stopifnot(all(c("counts","data","scale.data")%in%names(x$layers)),all(table(x$meta.data$zonation_type)==length(x$cells)/2L))
      rm(x)
    }
    stopifnot(file.exists(analysis_path("engine_fcg05")),file.exists(analysis_path("engine_validation07")),
      file.exists("results/validation/knockdown_t_table.csv"),file.exists("results/validation/knockdown_ttop_list.RData"))
    if(audit)stopifnot(nzchar(py),file.exists(py),file.exists(file.path(analysis_path("historical_fcg05"),"zonation_ttop_list.RData")),
      file.exists(file.path(analysis_path("historical_cf453_models"),"models.rds")))
    cat("Current scientific inputs, frozen engines and unchanged KD fits validated; no outputs written.\n")
    if(audit)cat("Optional retrospective-audit historical references and Python also validated.\n")
    return(invisible(TRUE))
  }
  run<-function(command,arguments,name){
    dir.create(file.path(output,"logs"),recursive=TRUE,showWarnings=FALSE)
    log<-file.path(output,"logs",paste0(name,".log"));cat("Running",name,"->",log,"\n")
    status<-system2(command,vapply(arguments,shQuote,character(1L)),stdout=log,stderr=log)
    if(status!=0L)stop("Stage failed: ",name,"; see ",log)
  }
  if(!"--report-only"%in%args){
    if(file.exists(file.path(output,"cf453_models/models.rds"))||file.exists(file.path(output,"fcg05/zonation_ttop_list.RData")))stop("Approved results are immutable: choose a fresh --output=results/runs/NAME or explicit --report-only")
    dir.create(file.path(output,"QA"),recursive=TRUE,showWarnings=FALSE)
    if(audit){
      files<-c(analysis_path("historical_fcg_zonation_csv"),analysis_path("historical_fcg_zonation_metadata"),
        analysis_path("historical_fcg_zonation_rds"),analysis_path("historical_gdx_zonation_rds"),analysis_path("historical_wdsd_zonation_rds"),
        list.files(analysis_path("historical_fcg05"),recursive=TRUE,full.names=TRUE),
        list.files("results/validation",recursive=TRUE,full.names=TRUE))
      files<-sort(unique(files[file.exists(files)&grepl("[.](csv|rds|RData|pdf|docx)$",files)&!grepl("cyp2c38_in_vivo",files)]))
      original<-data.frame(path=normalizePath(files),sha256=unname(vapply(files,digest::digest,character(1L),file=TRUE,algo="sha256")),bytes=file.info(files)$size)
      utils::write.csv(original,file.path(output,"QA/historical_sha256.csv"),row.names=FALSE)
    }
    staged<-staging$zu_stage_inputs(output,source_root)
    mode<-if(audit)character() else "--science-only"
    run("Rscript",c("scripts/reports/run_zonation_update_fcg.R",paste0("--output=",output),mode),"workflow17_fcg")
    run("Rscript",c("scripts/reports/run_zonation_update_gdx.R",paste0("--output=",output),mode),"workflow17_gdx")
    run("Rscript",c("scripts/workflows/14_cf453_standalone_figure6.R","--models-only",paste0("--input-zonation=",staged$sources[["wdsd"]]),paste0("--output=",output,"/cf453_models")),"workflow17_cf453_14")
    run("Rscript",c("scripts/workflows/15_cf453_diet_impact_figure6.R","--models-only",paste0("--parent=",output,"/cf453_models/models.rds"),paste0("--output=",output,"/cf453_direct_diet")),"workflow17_cf453_15")
    run("Rscript",c("scripts/workflows/16_cf453_diet_rewiring_figure6.R",paste0("--parent=",output,"/cf453_models/models.rds"),paste0("--direct=",output,"/cf453_direct_diet/diet_analysis.rds"),paste0("--output=",output,"/cf453_rewiring")),"workflow17_cf453_16")
    run("Rscript",c("scripts/reports/cf453_figure6a_layouts.R",paste0("--parent=",output,"/cf453_rewiring/rewiring_analysis.rds"),paste0("--output=",output,"/cf453_layouts")),"workflow17_cf453_layouts")
    if(audit)run("Rscript",c("scripts/reports/compare_cf453_zonation_update.R",output),"workflow17_cf453_comparisons")
  }
  run("Rscript",c("tests/qa_current_zonation_release.R",paste0("--root=",output)),"workflow17_scientific_QA")
  if(nzchar(py)&&file.exists(py))run(py,c("scripts/reports/prepare_zonation_update_report_sources.py","--output",output,"--tables-only"),"workflow17_scientific_tables")
  if(audit){
    if(!nzchar(py)||!file.exists(py))stop("Optional manuscript audit requires Python3; set --python= or ZONATION_PYTHON")
    run(py,c("scripts/reports/prepare_zonation_update_report_sources.py","--output",output,"--source-root",source_root),"workflow17_report_sources")
    run("Rscript",c("tests/qa_zonation_marker_update.R",output),"workflow17_independent_audit_QA")
    run(py,c("scripts/reports/zonation_update_manuscript_audit.py","--output",file.path(output,"report")),"workflow17_manuscript_inventory")
    run(py,c("scripts/reports/build_zonation_update_report.py","--output",output),"workflow17_report")
    run(py,c("scripts/reports/package_zonation_marker_update.py","--output",output),"workflow17_bundle")
  }
  cat("Current expanded-marker scientific pipeline complete in",output,"\n")
  invisible(output)
}
if(sys.nframe()==0L)zu17_run()
