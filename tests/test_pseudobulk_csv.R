#!/usr/bin/env Rscript
# Targeted transport tests; temporary aggregate profiles only, no model fitting.
source("scripts/support_functions.R")
source("scripts/core/pseudobulk_csv.R")
stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
root <- tempfile("pseudobulk_csv_test_"); dir.create(root)
ids <- paste0("mouse", 1:6, "_Hepatocytes_group")
counts <- matrix(as.double(seq_len(1200) %% 27), 200, 6,
  dimnames = list(paste0("Gene", seq_len(200)), ids))
md <- data.frame(orig.ident = ids, sample_name = paste0("mouse", 1:6), group = "group",
  literal = c("NA", NA, "", "comma,quote\"", "two\nlines", "literal"),
  factor_col = factor(c("NA", "", NA, "high", "NA", "high"),
    levels = c("", "NA", "high", "unused"), ordered = TRUE),
  integer_col = c(1L, NA_integer_, 3L, 4L, 5L, 6L),
  logical_col = c(TRUE, FALSE, NA, TRUE, FALSE, TRUE),
  row.names = ids, check.names = FALSE)
o <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = md)
values <- matrix(rep(jsonlite::fromJSON(paste0("[0.24913119230963265,",
  "0.076555803664553529,1.0000000000000002,-1.0000000000000002,1e-100,1e100]")),
  length.out = length(counts)), nrow(counts), ncol(counts), dimnames = dimnames(counts))
SeuratObject::LayerData(o, assay = "RNA", layer = "data") <- values
SeuratObject::LayerData(o, assay = "RNA", layer = "scale.data") <- values / 3
source <- file.path(root, "aggregate_profiles.rds"); saveRDS(o, source)
output <- file.path(root, "transport")
pb_csv_write(source, output, max_bytes = 10000L)
a <- pb_components(source); b <- pb_read(file.path(output, "manifest.csv"))
stopifnot(identical(a$layers, b$layers), identical(a$meta.data, b$meta.data),
  identical(a$identities, b$identities), identical(a$identity_levels, b$identity_levels),
  any(vapply(b$manifest$layers, function(z) length(z$chunks) > 1L, logical(1))),
  all(file.info(pb_source_files(output))$size[grepl("[.]csv$", pb_source_files(output))] < 10000L),
  identical(b$meta.data$literal, md$literal), identical(b$meta.data$factor_col, md$factor_col))
bad <- file.path(output, b$manifest$layers$data$chunks[[1]]$file)
cat("\n", file = bad, append = TRUE)
rejected <- try(pb_read(output), silent = TRUE)
stopifnot(inherits(rejected, "try-error"))
# Nucleus-style metadata repeats a mouse ID instead of one aggregate-profile ID.
o$orig.ident <- rep("mouse1", ncol(o))
nuclei_style <- file.path(root, "nucleus_style_metadata.rds"); saveRDS(o, nuclei_style)
rejected <- try(pb_csv_write(nuclei_style, file.path(root, "must_not_export")), silent = TRUE)
stopifnot(inherits(rejected, "try-error"), !dir.exists(file.path(root, "must_not_export")))
unlink(root, recursive = TRUE)
cat("Pseudobulk transport tests passed: exact binary64, forced split, metadata types/NA, tamper and nucleus-style rejection.\n")
