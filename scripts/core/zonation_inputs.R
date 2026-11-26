# Adapt supplied scaled pseudobulk profiles to the zonation engine's CSV format.
zu_write_scaled_inputs <- function(pseudobulk, destination, cohort) {
  stopifnot(cohort %in% c("fcg", "gdx"))
  x <- as.matrix(pseudobulk$layers$scale.data)
  metadata <- pseudobulk$meta.data
  stopifnot(identical(colnames(x), rownames(metadata)),
    all(c("orig.ident", "sample_name") %in% names(metadata)))
  metadata$sample <- gsub("-", "_", as.character(metadata$orig.ident), fixed = TRUE)
  metadata$sample_name <- gsub("-", "_", metadata$sample_name, fixed = TRUE)
  stopifnot(!anyNA(metadata$sample), !anyDuplicated(metadata$sample))
  metadata <- metadata[, c("sample", setdiff(names(metadata), c("sample", "orig.ident")))]
  colnames(x) <- metadata$sample
  tag <- if (cohort == "fcg") "" else "gdx_"
  paths <- file.path(destination, paste0(c("scaled_RNA_", "meta_data_df_"), tag, "zonation.csv"))
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(gene = rownames(x), x, check.names = FALSE), paths[1L], row.names = FALSE)
  utils::write.csv(metadata, paths[2L], row.names = FALSE)
  invisible(paths)
}
