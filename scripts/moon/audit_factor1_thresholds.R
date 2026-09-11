#!/usr/bin/env Rscript
# Extraction sensitivity only: reuse the frozen final scores and scoring graph.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
out <- "results/cosmos_factor1_pilot"
scores <- utils::read.csv(file.path(out, "moon_scores.csv"))
annotated <- utils::read.csv(file.path(out, "moon_scores_annotated.csv"), check.names = FALSE)
pkn <- readRDS(file.path(out, "pkn_final_uncompressed.rds"))
up <- utils::read.csv(file.path(out, "upstream_input.csv"))
upstream <- stats::setNames(up$score, up$source)
rows <- list()
for (secondary in c(0, 0.25, 0.5, 0.75, 1)) {
  sol <- moon_extract_paths(scores, pkn, upstream, primary_thresh = 1.7, secondary_thresh = secondary)
  rows[[length(rows) + 1L]] <- data.frame(primary = 1.7, secondary = secondary,
    nodes = nrow(sol$ATT), edges = nrow(sol$SIF),
    upstream_tfs = sum(sol$ATT$type == "upstream_input"),
    downstream_nodes = sum(sol$ATT$level == 0))
  if (secondary == 0.5) {
    utils::write.csv(sol$SIF, file.path(out, "solution_intermediate_SIF.csv"), row.names = FALSE)
    att <- annotated[match(sol$ATT$source, annotated$source), , drop = FALSE]
    utils::write.csv(att, file.path(out, "solution_intermediate_ATT.csv"), row.names = FALSE, na = "")
  }
}
grid <- do.call(rbind, rows)
utils::write.csv(grid, file.path(out, "secondary_threshold_sensitivity.csv"), row.names = FALSE)
print(grid, row.names = FALSE)
