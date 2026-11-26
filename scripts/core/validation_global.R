# Independent siRNA/global-network branch of the frozen Workflow 07 engine.
# Reuse the frozen global fit and render Figure 3F–G from matched lipids.
source("scripts/core/zonation_update_fcg.R")
source("scripts/core/validation_fig3fg.R")

validation_global_network_path <- paste0(
  "results/current/03_gonadal_sex_biased_regulatory_circuitry/",
  "mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv"
)

validation_fig3fg_run <- function(output, validate_only = FALSE) {
  knockdown_path <- paste0("results/current/03_gonadal_sex_biased_regulatory_circuitry/",
                           "in_vitro_validation/knockdown_t_table.csv")
  inputs <- c(knockdown_path, validation_global_network_path)
  if (!all(file.exists(inputs))) stop("Missing Figure 3F–G input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  t_table <- utils::read.csv(knockdown_path, check.names = FALSE)
  sif <- utils::read.csv(validation_global_network_path)
  network <- stats::setNames(sif[c("node1", "node2", "mor")], c("source", "target", "mor"))
  panels <- fig3fg_prepare(t_table, network)
  print(do.call(rbind, lapply(panels, "[[", "summary")), row.names = FALSE)
  if (validate_only) return(invisible(TRUE))

  if (file.exists(output) && (!dir.exists(output) || length(list.files(output, all.files = TRUE, no.. = TRUE)))) {
    stop("Preserve existing panels: --fig3fg-only requires a fresh --output")
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  fig3fg_write(panels, output)
}

validation_global_runtime <- function(output) {
  text <- readLines(analysis_path("engine_validation07",mustWork=TRUE),warn=FALSE)
  at8 <- grep('^# 8\\. VALIDATION:',text)[1L]
  at10 <- grep('^# 10\\. SUMMARY',text)[1L]
  end11 <- grep('^cat\\("Output directory: results/validation/structural',text)[1L]
  stopifnot(all(is.finite(c(at8,at10,end11))),at8<at10,at10<end11)
  selected <- c(seq_len(at8-1L),seq.int(at10,end11))
  branch <- text[selected]
  branch <- branch[!grepl('^rm\\(list = ls\\(\\)\\)',branch)]
  # Keep all global candidates for scoring and structural checks. Only the
  # first loop is the historical scatter, now needed solely for Cyp2b9.
  plot_loop <- match("for (gene_name in available_global) {", branch)
  stopifnot(!is.na(plot_loop))
  branch[plot_loop] <- "for (gene_name in setdiff(available_global, tolower(fig3fg_genes))) {"
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
    validation_global_network_path,
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
  withCallingHandlers({
    sys.source(runtime_path,envir=e)
    fig3fg_write(fig3fg_prepare(e$t_table, e$lipid_network), output)
  },warning=function(w){
    warnings[[length(warnings)+1L]] <<- data.frame(message=conditionMessage(w),call=paste(deparse(conditionCall(w)),collapse=" "))
    invokeRestart("muffleWarning")
  })
  invisible(output)
}
