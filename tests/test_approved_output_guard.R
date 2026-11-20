# Cache-absent public-clone fixtures; fitting/rendering functions are stubbed.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
test_approved_output_guard <- function(){
  root<-tempfile("approved_guard_repo_");dir.create(file.path(root,"scripts/config"),recursive=TRUE)
  root<-normalizePath(root,mustWork=TRUE)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  jsonlite::write_json(list(schema_version=1L,paths=list(current_root="results/approved",legacy_results_root="results/legacy")),
    file.path(root,"scripts/config/analysis_release.json"),auto_unbox=TRUE)
  approved<-file.path(root,"results/approved")
  dir.create(file.path(approved,"fcg05"),recursive=TRUE)
  artifact<-file.path(approved,"fcg05/t_table_zonation.csv");writeLines("curated,bytes",artifact)
  original<-digest::digest(artifact,file=TRUE,algo="sha256")
  expect_guard<-function(expr){
    error<-tryCatch(force(expr),error=identity)
    stopifnot(inherits(error,"error"),grepl("Preserve the approved release",conditionMessage(error),fixed=TRUE))
  }
  for(path in c(approved,file.path(approved,"missing_child"),"results/approved", "results/approved/new/../fcg05"))
    expect_guard(analysis_assert_not_approved_output(path,root))
  expect_guard(analysis_assert_not_approved_output("results/legacy/old_run",root))
  alias<-file.path(root,"approved_alias")
  if(file.symlink(approved,alias)){
    expect_guard(analysis_assert_not_approved_output(file.path(alias,"new_child"),root))
    expect_guard(analysis_assert_not_approved_output(file.path(root,"missing","..","approved_alias","new_child"),root))
    outside<-analysis_assert_not_approved_output(file.path(root,"missing","..","approved_alias","..","fresh_child"),root)
    stopifnot(identical(outside,file.path(dirname(approved),"fresh_child")))
  }
  dangling<-file.path(root,"dangling_alias")
  if(file.symlink(file.path(approved,"missing_target"),dangling)){
    for(path in c(file.path(dangling,"new_child"),
      file.path(approved,"missing_target","..","..","..","dangling_alias","new_child"))){
      error<-tryCatch(analysis_assert_not_approved_output(path,root),error=identity)
      stopifnot(inherits(error,"error"),grepl("Cannot safely resolve dangling output symlink",conditionMessage(error),fixed=TRUE),
        !dir.exists(file.path(approved,"missing_target")))
    }
  }
  if(.Platform$OS.type=="unix"){
    home_parts<-strsplit(normalizePath(path.expand("~"),mustWork=TRUE),"/",fixed=TRUE)[[1L]]
    # Reach only the temporary fixture from~; never create files in the home directory.
    tilde_root<-paste(c("~",rep("..",sum(nzchar(home_parts))),sub("^/","",root)),collapse="/")
    tilde_approved<-file.path(tilde_root,"results/approved")
    stopifnot(identical(normalizePath(tilde_approved,mustWork=TRUE),approved))
    expect_guard(analysis_assert_not_approved_output(tilde_approved,root))
    expect_guard(analysis_assert_not_approved_output(file.path(tilde_approved,"new_child"),root))
  }
  stopifnot(identical(analysis_assert_not_approved_output("results/approved_extra",root),file.path(root,"results/approved_extra")))
  stopifnot(!file.exists(file.path(approved,"fcg05/fcg05_checkpoint.rds")),
    !file.exists(file.path(approved,"fcg07/fcg07_checkpoint.rds")))

  five<-new.env(parent=globalenv());sys.source("scripts/workflows/05_zonation_analysis.R",envir=five)
  five$analysis_current<-function(...)approved
  five$analysis_assert_not_approved_output<-function(output)analysis_assert_not_approved_output(output,root)
  five$zmu_require<-function()TRUE;five$zmu_validate<-function(...)TRUE
  five$zmu_runtime05<-function(...)"TRUE"
  fits05<-0L;five$zmu_run05<-function(...){fits05<<-fits05+1L;invisible(TRUE)}
  expect_guard(five$fcg_zonation_run())
  expect_guard(five$fcg_zonation_run(paste0("--output=",approved)))
  stopifnot(fits05==0L,isTRUE(five$fcg_zonation_run("--validate-only")))
  fresh05<-file.path(root,"results/fresh05")
  five$fcg_zonation_run(paste0("--output=",fresh05))
  stopifnot(fits05==1L,!dir.exists(fresh05))

  driver<-new.env(parent=globalenv());sys.source("scripts/reports/run_zonation_update_fcg.R",envir=driver)
  driver$analysis_current<-function(...)approved
  driver$analysis_assert_not_approved_output<-function(output)analysis_assert_not_approved_output(output,root)
  driver$zmu_run05<-function(...)stop("No report-driver fit may run in this guard test")
  driver$zmu_require<-function()TRUE;driver$zmu_validate<-function(...)TRUE
  driver$zmu_runtime05<-function(...)"TRUE"
  driver$zmu_finish<-function(...)invisible(TRUE)
  expect_guard(driver$zonation_update_fcg_run())
  expect_guard(driver$zonation_update_fcg_run("--science-only"))
  driver$zonation_update_fcg_run("--validate-only")
  stopifnot(isTRUE(driver$zonation_update_fcg_run("--postprocess-only")))

  seven<-new.env(parent=globalenv());sys.source("scripts/workflows/07_validation_experiments.R",envir=seven)
  seven$analysis_current<-function(component="",...)if(nzchar(component))file.path(approved,component) else approved
  seven$analysis_assert_not_approved_output<-function(output)analysis_assert_not_approved_output(output,root)
  seven$zmu_require<-function()TRUE
  fits07<-0L;seven$zmu_run07<-function(...){fits07<<-fits07+1L;invisible(TRUE)}
  expect_guard(seven$zonation_validation_run())
  expect_guard(seven$zonation_validation_run(paste0("--output=",approved)))
  stopifnot(fits07==0L)
  # Validation remains allowed when its declared prerequisites exist.
  writeLines("fixture checkpoint, never read",file.path(approved,"fcg05/fcg05_checkpoint.rds"))
  stopifnot(isTRUE(seven$zonation_validation_run("--validate-only")))
  fresh07<-file.path(root,"results/fresh07")
  dir.create(file.path(fresh07,"fcg05"),recursive=TRUE)
  writeLines("fixture checkpoint, never read",file.path(fresh07,"fcg05/fcg05_checkpoint.rds"))
  seven$zonation_validation_run(paste0("--output=",fresh07))
  stopifnot(fits07==1L,!dir.exists(file.path(fresh07,"fcg07")))

  current<-function(component="",...)if(nzchar(component))file.path(approved,component) else approved
  seventeen<-new.env(parent=globalenv());sys.source("scripts/workflows/17_expanded_marker_zonation_update.R",envir=seventeen)
  seventeen$analysis_current<-current
  seventeen$analysis_assert_not_approved_output<-function(output)analysis_assert_not_approved_output(output,root)
  stage_calls<-0L;commands<-character()
  seventeen$sys.source<-function(file,envir,...){
    stage_calls<<-stage_calls+1L
    envir$zu_source_input<-function(...)"fixture manifest"
    envir$zu_stage_inputs<-function(...)stop("No scientific staging may run in this guard test")
  }
  seventeen$source<-function(file,local=FALSE,...){
    stopifnot(basename(file)=="pseudobulk_csv.R")
    assign("pb_components",function(...)list(layers=list(counts=1,data=1,scale.data=1),
      cells=c("PP","PC"),meta.data=data.frame(zonation_type=c("PP","PC"))),envir=parent.frame())
    invisible(NULL)
  }
  seventeen$dir.create<-function(...)TRUE
  seventeen$system2<-function(command,args,...){commands<<-c(commands,command);0L}
  expect_guard(seventeen$zu17_run())
  expect_guard(seventeen$zu17_run("--with-manuscript-audit"))
  stopifnot(stage_calls==0L,isTRUE(seventeen$zu17_run("--validate-only")))
  # Full-run report refresh needs a saved work directory outside published/legacy paths.
  error<-tryCatch(seventeen$zu17_run("--report-only"),error=identity)
  stopifnot(inherits(error,"error"),grepl("saved Workflow 17 run",conditionMessage(error),fixed=TRUE))
  expect_guard(seventeen$zu17_run(c("--report-only",paste0("--output=",approved))))
  seventeen$zu17_run(c("--report-only",paste0("--output=",root,"/results/saved17")))
  stopifnot(length(commands)>0L)

  fourteen<-new.env(parent=globalenv());sys.source("scripts/workflows/14_cf453_standalone_figure6.R",envir=fourteen)
  fourteen$analysis_current<-current
  fourteen$analysis_assert_not_approved_output<-function(output)analysis_assert_not_approved_output(output,root)
  fourteen$dir.create<-function(...)TRUE
  fourteen$readRDS<-function(...)stop("Explicit saved-result mode reached read stub")
  fourteen$cf6_models<-function(...)stop("No model may be fitted in this guard test")
  expect_guard(fourteen$cf6_run())
  expect_guard(fourteen$cf6_run("--models-only"))
  for(mode in c("--finish","--plot-only")){
    error<-tryCatch(fourteen$cf6_run(mode),error=identity)
    stopifnot(inherits(error,"error"),conditionMessage(error)=="Explicit saved-result mode reached read stub")
  }

  layouts<-new.env(parent=globalenv());sys.source("scripts/reports/cf453_figure6a_layouts.R",envir=layouts)
  layouts$analysis_current<-current
  layouts$analysis_assert_not_approved_output<-function(output)analysis_assert_not_approved_output(output,root)
  layouts$dir.create<-function(...)TRUE
  layouts$cf6a_original_files<-function(...)character()
  layouts$cf6a_prepare<-function(...)stop("Fresh layout reached prepare stub")
  expect_guard(layouts$cf6a_run())
  expect_guard(layouts$cf6a_run(paste0("--output=",approved,"/another_layout")))
  error<-tryCatch(layouts$cf6a_run(paste0("--output=",root,"/results/fresh_layout")),error=identity)
  stopifnot(inherits(error,"error"),conditionMessage(error)=="Fresh layout reached prepare stub",
    identical(original,digest::digest(artifact,file=TRUE,algo="sha256")))
  cat("Approved-root guards passed for05/07/14/17/layouts with missing caches and stubbed models/renders; explicit validation/refresh modes stay allowed.\n")
  invisible(TRUE)
}
if(sys.nframe()==0L)test_approved_output_guard()
