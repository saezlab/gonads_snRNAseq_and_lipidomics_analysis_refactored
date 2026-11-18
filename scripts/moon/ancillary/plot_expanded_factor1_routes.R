#!/usr/bin/env Rscript
# Illustrate three audited adjacent-layer routes; no fitting or selection occurs here.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/plot_expanded_factor1_routes.R AUDIT_DIR NEW_FIGURE_DIR\n")
  quit(save = "no")
}
if (length(args) != 2L) stop("Use --help")
audit <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
if (dir.exists(out) && length(list.files(out))) stop("Figure directory must be empty")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
nodes <- utils::read.csv(file.path(audit, "network_labeled_ATT.csv"), check.names = FALSE)
paths <- utils::read.csv(file.path(audit, "candidate_paths_score_layer_aligned.csv"), check.names = FALSE)
requests <- data.frame(tf = c("Gene__Ppara", "Gene__Nr2f1", "Gene__Gata6"),
  lipid = c("DG(16:0_22:6)", "PA(18:1_18:2)", "LPC(20:5)"),
  chemistry = c("PA dephosphorylation", "PC headgroup hydrolysis", "PC acyl hydrolysis"))
chosen <- lapply(seq_len(nrow(requests)), function(i) {
  row <- paths[paths$upstream_tf == requests$tf[i] & paths$downstream_label == requests$lipid[i], , drop = FALSE]
  stopifnot(nrow(row) == 1L, tolower(as.character(row$score_layer_aligned)) == "true",
    tolower(as.character(row$contains_currency_bridge)) == "false")
  row
})
draw <- function() {
  graphics::par(mar = rep(0.5, 4), family = "sans")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  graphics::text(.035, .975, "Expanded Factor 1 MOON: audited example routes", adj = c(0, 1), cex = 1.3, font = 2)
  graphics::text(.035, .925, "Compressed ULM; pruning 6, scoring 4; strict cutoffs 0.75 / 0.75; generic cofactor traversal blocked",
    adj = c(0, 1), cex = .84)
  ys <- c(.76, .53, .30)
  for (i in seq_along(chosen)) {
    path <- chosen[[i]]
    ids <- strsplit(path$node_ids, "|", fixed = TRUE)[[1]]
    att <- nodes[match(ids, nodes$node_id), , drop = FALSE]
    stopifnot(!anyNA(att$node_id), all(diff(att$level) == -1))
    x <- seq(.08, .92, length.out = length(ids)); y <- ys[i]
    edges <- jsonlite::fromJSON(path$signed_edges_json)
    graphics::text(.035, y + .11, paste0(i, ". ", path$upstream_label, " to ", path$downstream_label,
      "   (minimum |inferred score| ", formatC(path$minimum_absolute_inferred_score, digits = 2, format = "f"), ")"),
      adj = c(0, .5), cex = .90, font = 2)
    for (j in seq_len(length(ids) - 1L)) {
      graphics::arrows(x[j] + .035, y, x[j + 1L] - .035, y, length = .065, col = "#555555", lwd = 1.5)
      graphics::text(mean(x[c(j, j + 1L)]), y + .023,
        if (as.numeric(edges[j, 2]) > 0) "+" else "-", cex = 1)
    }
    for (j in seq_along(ids)) {
      row <- att[j, ]
      color <- if (row$score > 0) "#4173d8" else "#e54141"
      enzyme <- row$node_type == "enzyme_reaction_instance"
      if (enzyme) graphics::polygon(x[j] + c(-.030, 0, .030, 0), y + c(0, .025, 0, -.025),
        border = color, col = "white", lwd = 2, lty = 2)
      else if (row$level == 0) graphics::points(x[j], y, pch = 23, bg = color, col = color, cex = 3)
      else graphics::points(x[j], y, pch = 21, bg = color, col = color, cex = 3)
      label <- if (enzyme) paste0(gsub("Gene__", "", row$catalysts, fixed = TRUE), " reaction\n", requests$chemistry[i]) else row$label
      graphics::text(x[j], y - .047, paste(strwrap(label, width = 22), collapse = "\n"), adj = c(.5, 1), cex = .78)
      graphics::text(x[j], y + .047, paste0(if (row$level == 0) "W = " else "score = ",
        formatC(row$score, digits = 2, format = "f"), "; level ", row$level), cex = .70, col = color)
    }
  }
  graphics::text(.035, .105, "Each arrow follows a consecutive scoring layer; all catalysts are present in the selected graph.", adj = c(0, 1), cex = .80)
  graphics::text(.035, .070, "Blue: positive Factor 1 sign; red: negative. Terminal W values are measured lipid loadings, not inferred scores.", adj = c(0, 1), cex = .78)
  graphics::text(.035, .035, "Hypotheses only: assumed species specificity and unmeasured co-substrates; no flux or joint substrate availability is inferred.",
    adj = c(0, 1), cex = .76)
}
grDevices::pdf(file.path(out, "audited_example_routes.pdf"), width = 14, height = 9, useDingbats = FALSE)
draw(); grDevices::dev.off()
grDevices::png(file.path(out, "audited_example_routes.png"), width = 2100, height = 1350, res = 150)
draw(); grDevices::dev.off()
utils::write.csv(do.call(rbind, chosen), file.path(out, "illustrated_routes.csv"), row.names = FALSE)
utils::write.csv(data.frame(path = c(file.path(audit, "network_labeled_ATT.csv"), file.path(audit, "candidate_paths_score_layer_aligned.csv")),
  md5 = unname(tools::md5sum(c(file.path(audit, "network_labeled_ATT.csv"), file.path(audit, "candidate_paths_score_layer_aligned.csv"))))),
  file.path(out, "input_manifest.csv"), row.names = FALSE)
cat("PASS: three audited routes rendered to", out, "\n")
