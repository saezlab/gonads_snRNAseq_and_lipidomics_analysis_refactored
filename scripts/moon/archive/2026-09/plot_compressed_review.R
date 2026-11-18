#!/usr/bin/env Rscript
# Summary of the finalized compressed candidates and the three focused TFs.
# Read-only with respect to scoring graphs, score tables and parameter grids.
source("scripts/support_functions.R")
base <- "results/cosmos_factor1_compressed_review"
selected <- utils::read.csv(file.path(base, "threshold_search/selected_candidates.csv"),
                            stringsAsFactors = FALSE)
audit <- utils::read.csv(file.path(base, "tf_audit/focus_tf_stage_audit.csv"),
                         stringsAsFactors = FALSE)
selected <- selected[match(c("C01", "C03", "C02"), selected$candidate_id), ]
focus <- audit[audit$run_id == "radius6_depth4_compressed", ]
focus <- focus[match(paste0("Gene__", c("Esr1", "Nfic", "Dbp")), focus$source), ]
stopifnot(nrow(selected) == 3L, !anyNA(selected$candidate_id), all(selected$edges < 500),
          all(selected$depth == 4L), all(selected$valid), all(selected$converged),
          all(focus$in_pruned_pkn), all(focus$compressed_id_unchanged),
          all(focus$in_final_scoring_pkn), all(is.finite(focus$raw_final_score)),
          all(focus$activity < 0), all(focus$raw_final_score > 0),
          all(!focus$in_saved_scores),
          all(focus$exclusion_reason == "post_scoring_upstream_sign_conflict"))
figures <- file.path(base, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

draw_summary <- function() {
  grid::grid.newpage()
  txt <- function(label, x, y, size = 11, color = "#29343D", face = "plain", just = "left") {
    grid::grid.text(label, x, y, just = just,
      gp = grid::gpar(fontfamily = "sans", fontsize = size, col = color, fontface = face))
  }
  ln <- function(x0, y0, x1, y1, color = "#DFE4E8", width = 1, type = 1) {
    grid::grid.lines(c(x0, x1), c(y0, y1), gp = grid::gpar(col = color, lwd = width, lty = type))
  }
  bar <- function(x0, x1, y, color, height = 0.025) {
    grid::grid.rect(x = (x0 + x1) / 2, y = y, width = abs(x1 - x0), height = height,
                    gp = grid::gpar(col = NA, fill = color))
  }
  txt("Compressed MOON solutions and TF exclusion audit", 0.055, 0.956, 20, face = "bold")
  txt("Factor 1: hepatocyte TF activities -> liver lipid weights | Compressed network, depth 4, topology radius 6",
      0.055, 0.920, 11, color = "#596773")
  txt("A  Lower extraction thresholds retain compact compressed solutions", 0.055, 0.873, 13, face = "bold")
  txt("Primary / secondary cutoffs are shown for each candidate; bars count original-ID edges after decompression.",
      0.055, 0.842, 10, color = "#596773")
  bar_start <- 0.326
  edge_x <- function(v) bar_start + 0.468 * v / 550
  ys <- c(0.774, 0.694, 0.614)
  labels <- c("C01  Higher thresholds", "C03  Balanced candidate", "C02  Broader coverage")
  fills <- c("#A9BBC7", "#245D79", "#A9BBC7")
  for (tick in seq(0, 500, 100)) {
    ln(edge_x(tick), 0.589, edge_x(tick), 0.801, color = "#EDF0F2")
    txt(as.character(tick), edge_x(tick), 0.565, 9, color = "#596773", just = "center")
  }
  ln(edge_x(500), 0.589, edge_x(500), 0.801, color = "#A3603D", type = 2, width = 1.3)
  txt("500-edge cap", edge_x(500), 0.814, 9, color = "#925532", just = "center")
  for (i in seq_len(nrow(selected))) {
    x <- selected[i, ]
    txt(labels[i], 0.055, ys[i] + 0.015, 11, face = if (i == 2L) "bold" else "plain")
    txt(paste0("Thresholds ", format(x$primary, trim = TRUE), " / ", format(x$secondary, trim = TRUE)),
        0.055, ys[i] - 0.015, 10, color = "#596773")
    bar(bar_start, edge_x(x$edges), ys[i], fills[i], height = 0.032)
    txt(paste(x$edges, "edges"), 0.945, ys[i] + 0.013, 11, just = "right", face = "bold")
    txt(paste0(x$upstream_tfs, " TFs; ", x$downstream_lipids, " lipids"),
        0.945, ys[i] - 0.015, 10, color = "#596773", just = "right")
  }
  ln(0.055, 0.534, 0.945, 0.534)
  txt("B  RNA-derived TF activities", 0.055, 0.495, 13, face = "bold")
  txt("C  Raw MOON scores before sign filtering", 0.545, 0.495, 13, face = "bold")
  txt("Negative along the oriented factor", 0.055, 0.464, 10, color = "#596773")
  txt("Positive, on an independent score scale", 0.545, 0.464, 10, color = "#596773")
  ys <- c(0.397, 0.317, 0.237)
  act_x <- function(v) 0.456 + 0.27 * v / 9
  moon_x <- function(v) 0.648 + 0.245 * v / 0.9
  for (tick in c(-8, -4, 0)) {
    ln(act_x(tick), 0.211, act_x(tick), 0.426, color = if(tick == 0) "#93A0AA" else "#EDF0F2")
    txt(as.character(tick), act_x(tick), 0.187, 10, color = "#596773", just = "center")
  }
  for (tick in c(0, 0.4, 0.8)) {
    ln(moon_x(tick), 0.211, moon_x(tick), 0.426, color = if(tick == 0) "#93A0AA" else "#EDF0F2")
    txt(format(tick, nsmall = 1), moon_x(tick), 0.187, 10, color = "#596773", just = "center")
  }
  for (i in seq_len(nrow(focus))) {
    gene <- sub("Gene__", "", focus$source[i], fixed = TRUE)
    txt(gene, 0.055, ys[i], 12, face = "italic")
    txt(gene, 0.545, ys[i], 12, face = "italic")
    bar(act_x(focus$activity[i]), act_x(0), ys[i], "#E54141")
    bar(moon_x(0), moon_x(focus$raw_final_score[i]), ys[i], "#4173D8")
    txt(sprintf("%.2f", focus$activity[i]), act_x(focus$activity[i]) - 0.01,
        ys[i], 10, just = "right")
    txt(sprintf("+%.3f", focus$raw_final_score[i]), moon_x(focus$raw_final_score[i]) + 0.01,
        ys[i], 10)
  }
  txt("RNA activity (ULM t-statistic)", 0.326, 0.151, 10, just = "center")
  txt("Inferred MOON score (ULM t-statistic)", 0.780, 0.151, 10, just = "center")
  ln(0.055, 0.114, 0.945, 0.114)
  txt("All three TFs survive pruning and compression, but their MOON signs oppose their RNA activity signs.",
      0.5, 0.080, 11, face = "bold", just = "center")
  txt("They are removed after scoring. Lowering extraction thresholds cannot restore them while retaining the sign rule.",
      0.5, 0.047, 10, color = "#596773", just = "center")
}

grDevices::pdf(file.path(figures, "compressed_review_summary.pdf"), width = 12, height = 8,
               useDingbats = FALSE)
draw_summary()
invisible(grDevices::dev.off())
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(file.path(figures, "compressed_review_summary.png"), width = 12, height = 8,
                units = "in", res = 220)
} else {
  grDevices::png(file.path(figures, "compressed_review_summary.png"), width = 2640, height = 1760,
                 res = 220)
}
draw_summary()
invisible(grDevices::dev.off())
cat("Saved compressed review PDF and PNG. Input assertions passed for three candidates and three TFs.\n")
