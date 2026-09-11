#!/usr/bin/env Rscript
# Read-only bridge for the Python candidate auditor; no scoring or extraction.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript read_saved_candidate_pkn.R FINAL_PKN_RDS CANDIDATE_ATT_CSV")
pkn <- readRDS(args[[1L]])
attributes <- utils::read.csv(args[[2L]], check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("source", "interaction", "target") %in% names(pkn)), "source" %in% names(attributes))
focus <- unique(c(attributes$source, "Gene__Nfe2l2", "Gene__Nr4a2", "Gene__Klf10", "Gene__Ppara", "Gene__Tbx2"))
pkn <- pkn[pkn$source %in% focus, c("source", "interaction", "target"), drop = FALSE]
utils::write.csv(pkn, stdout(), row.names = FALSE, na = "")
