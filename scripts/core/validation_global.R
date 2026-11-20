# Independent siRNA/global-network branch of the frozen Workflow 07 engine.
# Only section selection, output redirection and the plot CI helper are adapted.
source("scripts/core/zonation_update_fcg.R")

validation_global_runtime <- function(output) {
  text <- readLines(analysis_path("engine_validation07",mustWork=TRUE),warn=FALSE)
  at8 <- grep('^# 8\\. VALIDATION:',text)[1L]
  at10 <- grep('^# 10\\. SUMMARY',text)[1L]
  end11 <- grep('^cat\\("Output directory: results/validation/structural',text)[1L]
  stopifnot(all(is.finite(c(at8,at10,end11))),at8<at10,at10<end11)
  selected <- c(seq_len(at8-1L),seq.int(at10,end11))
  branch <- text[selected]
  branch <- branch[!grepl('^rm\\(list = ls\\(\\)\\)',branch)]
  # Bind after support_functions loads, before any plotting. This implements
  # Hmisc's mean +/- t * SE without making Hmisc a required plotting package.
  at <- grep('^source\\("scripts/support_functions.R"\\)',branch)[1L]
  branch <- append(branch,c('mean_cl_normal <- zmu_mean_cl_normal',
    'available_zonation <- character(0) # Prevent inherited interactive zonation state.'),after=at)
  branch <- gsub("results/validation",output,branch,fixed=TRUE)
  branch <- gsub("results/mofa/", "results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/", branch, fixed=TRUE)
  branch <- gsub("results/lipidomic/", "results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/", branch, fixed=TRUE)
  stopifnot(!any(grepl('^# [89]\\.',branch)),!any(grepl('"results/zonation/',branch,fixed=TRUE)))
  attr(branch,"engine_lines") <- selected
  branch
}

validation_global_run <- function(output,validate_only=FALSE) {
  zmu_require()
  required <- c(analysis_path("engine_validation07",mustWork=TRUE),
    "data/metabolomics/lipidomics_20250423transfections.csv",
    "results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv",
    "results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipids_with_structure.csv")
  if(!all(file.exists(required)))stop("Missing independent validation input: ",paste(required[!file.exists(required)],collapse=", "))
  output <- file.path(normalizePath(dirname(output),mustWork=FALSE),basename(output))
  runtime <- validation_global_runtime(output)
  invisible(parse(text=runtime))
  if(validate_only){cat("Independent KD/global/structural inputs and frozen-engine extraction validated; no files changed.\n");return(invisible(TRUE))}
  if(dir.exists(output)&&length(list.files(output,all.files=TRUE,no..=TRUE)))
    stop("Preserve existing validation: --global-only requires a fresh --output")
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  output <- normalizePath(output)
  runtime <- validation_global_runtime(output)
  runtime_path <- file.path(output,"runtime/07_global_validation_versioned.R")
  zmu_write_runtime(runtime,runtime_path)
  utils::write.csv(data.frame(engine_line=attr(runtime,"engine_lines")),
    file.path(output,"runtime/engine_section_selection.csv"),row.names=FALSE)
  e <- new.env(parent=globalenv());warnings <- list()
  on.exit({
    w <- if(length(warnings))do.call(rbind,warnings) else data.frame(message=character(),call=character())
    utils::write.csv(w,file.path(output,"warnings.csv"),row.names=FALSE)
  },add=TRUE)
  withCallingHandlers(sys.source(runtime_path,envir=e),warning=function(w){
    warnings[[length(warnings)+1L]] <<- data.frame(message=conditionMessage(w),call=paste(deparse(conditionCall(w)),collapse=" "))
    invokeRestart("muffleWarning")
  })
  invisible(output)
}
