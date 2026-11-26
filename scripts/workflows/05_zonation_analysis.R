# Current FCG and paired gonadectomy zonation from the released pseudobulks.
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/zonation_update_fcg.R")

fcg_zonation_run <- function(args = commandArgs(TRUE)) {
  if (any(!grepl("^(--(output|input)=.+|--cohort=(fcg|gdx)|--validate-only)$", args)) ||
      anyDuplicated(sub("=.*$", "", args)))
    stop("Use --cohort=fcg or --cohort=gdx, --input=, --output= and optionally --validate-only")
  value <- function(key, default) {
    x <- args[startsWith(args, paste0(key, "="))]
    if (length(x)) sub("^[^=]*=", "", x) else default
  }
  cohort <- value("--cohort", "fcg")
  input <- value("--input", analysis_path(paste0("current_", cohort, "_zonation_transport")))
  output <- value("--output", analysis_current())
  validate_only <- "--validate-only" %in% args
  if (!validate_only) analysis_assert_not_approved_output(output)

  if (cohort == "gdx") {
    source("scripts/core/zonation_update_gdx.R", local = TRUE)
    if (!validate_only) return(zu_gdx_run(output, science_only = TRUE, input = input))
    if (!requireNamespace("limma", quietly = TRUE)) stop("Missing R package: limma")
    gdx <- zu_gdx_read_input(input)
    cat("GDX input validated:", length(gdx$features), "genes x",
      length(gdx$cells), "profiles from", length(unique(gdx$meta.data$sample_name)), "mice.\n")
    return(invisible(TRUE))
  }

  if (!validate_only && file.exists(file.path(output, "fcg05/fcg05_checkpoint.rds")))
    stop("A completed FCG05 destination cannot be overwritten: ", file.path(output, "fcg05"))
  input_dir <- input
  if (!dir.exists(input) || (file.exists(file.path(input, "manifest.json")) &&
      !file.exists(file.path(input, "scaled_RNA_zonation.csv")))) {
    source("scripts/core/pseudobulk_csv.R", local = TRUE)
    source("scripts/core/zonation_inputs.R", local = TRUE)
    pseudobulk <- pb_components(input)
    if (validate_only) {
      input_dir <- tempfile("fcg_zonation_")
      on.exit(unlink(input_dir, recursive = TRUE), add = TRUE)
    } else {
      input_dir <- analysis_assert_not_approved_output(file.path(output, "inputs/fcg"))
    }
    engine_inputs <- zu_write_scaled_inputs(pseudobulk, input_dir, "fcg")
    if (!validate_only) {
      utils::write.csv(analysis_provenance(c(pseudobulk$source_paths, engine_inputs)),
        file.path(input_dir, "input_provenance.csv"), row.names = FALSE)
    }
    rm(pseudobulk)
  }
  if (validate_only) {
    zmu_require()
    zmu_validate(input_dir)
    invisible(parse(text = zmu_runtime05(input_dir, file.path(output, "fcg05"))))
    cat("FCG inputs and isolated engine validated.\n")
    return(invisible(TRUE))
  }
  zmu_run05(input_dir, output, current_panels = TRUE)
}
if (sys.nframe() == 0L) fcg_zonation_run()
