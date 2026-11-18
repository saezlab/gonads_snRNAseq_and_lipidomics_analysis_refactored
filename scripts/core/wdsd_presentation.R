# WP2 scientific plots and output-derived report. Source support_functions.R
# first. Sourcing this file defines functions without creating output artifacts.

wp2_context_label <- function(context) {
  labels <- c(male_periportal = "Male periportal", male_pericentral = "Male pericentral",
              female_periportal = "Female periportal", female_pericentral = "Female pericentral")
  result <- unname(labels[as.character(context)])
  result[is.na(result)] <- as.character(context)[is.na(result)]
  result
}

wp2_format_number <- function(x, digits = 3L) {
  vapply(x, function(value) if (!is.finite(value)) "NA" else
    formatC(value, format = "f", digits = digits), character(1L))
}

wp2_format_p <- function(x) {
  vapply(x, function(value) {
    if (!is.finite(value)) return("NA")
    if (value == 0) return("0 (numerical underflow)")
    if (value < 0.001) formatC(value, format = "e", digits = 2L) else
      formatC(value, format = "f", digits = 4L)
  }, character(1L))
}

wp2_inference_row <- function(inference, context, universe, gene = NULL) {
  selected <- inference$context == context & inference$universe == universe
  selected <- selected & if (is.null(gene)) inference$scope == "network" else
    inference$scope == "gene" & !is.na(inference$gene) & inference$gene == gene
  result <- inference[which(selected), , drop = FALSE]
  if (nrow(result) > 1L) stop("Duplicate presentation inference row.")
  result
}

wp2_plot_coefficient_panel <- function(edges, inference, context, universe,
                                      stage = c("baseline", "diet"), gene = NULL,
                                      compact = FALSE, tag = NULL) {
  stage <- match.arg(stage)
  e <- edges[edges$context == context & edges$eligible &
               (universe == "full" | edges$primary_edge), , drop = FALSE]
  if (!is.null(gene)) e <- e[e$gene == gene, , drop = FALSE]
  label <- wp2_context_label(context)
  title <- if (is.null(gene)) label else paste(gene, "-", label)
  if (!is.null(tag)) title <- paste0(tag, "  ", title)
  if (!nrow(e)) {
    graphics::plot.new(); graphics::title(main = title)
    graphics::text(0.5, 0.5, "No eligible fixed pairs", cex = 0.95)
    return(invisible(NULL))
  }
  x <- if (stage == "baseline") e$r_reference else e$r_SD
  y <- if (stage == "baseline") e$r_SD else e$r_WD
  colour <- if (grepl("^male_", context)) "#4173d8" else "#e54141"
  alpha <- if (nrow(e) > 1000L) 0.14 else 0.55
  pointsize <- if (nrow(e) > 1000L) 0.33 else 0.60
  xlab <- if (stage == "baseline") "Original FCG gene-lipid Pearson r" else
    "WDSD standard-diet Pearson r"
  ylab <- if (stage == "baseline") "WDSD standard-diet Pearson r" else
    "WDSD Western-diet Pearson r"
  graphics::plot(x, y, xlim = c(-1, 1), ylim = c(-1, 1),
                 pch = 16, cex = pointsize, col = grDevices::adjustcolor(colour, alpha.f = alpha),
                 xlab = xlab, ylab = ylab, main = title,
                 cex.main = if (compact) 0.96 else 1.02, cex.lab = 0.93,
                 axes = FALSE, xaxs = "i", yaxs = "i")
  graphics::axis(1, at = seq(-1, 1, by = 0.5)); graphics::axis(2, at = seq(-1, 1, by = 0.5), las = 1)
  graphics::box(col = "#555555")
  graphics::abline(h = 0, v = 0, col = "#dddddd", lty = 3)
  graphics::abline(0, 1, col = "#777777", lty = 2)
  if (stats::sd(x) > 0) graphics::abline(stats::lm(y ~ x), col = colour, lwd = 1.6)
  universe_label <- if (universe == "primary") "original displayed pairs" else "full eligible pairs"
  reference_n <- if (grepl("^male_", context)) 22L else 21L
  subtitle <- if (stage == "baseline")
    paste0(nrow(e), " ", universe_label, "; FCG n=", reference_n, ", SD n=7 mice") else
    paste0(nrow(e), " ", universe_label, "; SD n=7, WD n=7 mice")
  graphics::mtext(subtitle, side = 3, line = 0.18, cex = if (compact) 0.67 else 0.78)
  if (stage == "baseline") {
    row <- wp2_inference_row(inference, context, universe, gene)
    if (nrow(row)) {
      annotation <- c(paste0("Coefficient-vector r = ", wp2_format_number(row$coefficient_vector_r)),
                      paste0("ULM = ", wp2_format_number(row$ULM_score),
                             "; nominal p = ", wp2_format_p(row$ULM_p_value)),
                      paste0("Exact pairing p = ", wp2_format_p(row$pairing_permutation_p)),
                      paste0("Null-bootstrap p = ", wp2_format_p(row$null_bootstrap_p)))
    } else annotation <- "ULM inference unavailable for this fixed subset"
  } else {
    rho <- if (nrow(e) >= 3L && stats::sd(x) > 0 && stats::sd(y) > 0) stats::cor(x, y) else NA_real_
    diet_stats <- wp2_vector_stats(x, y)
    annotation <- c(paste0("SD-WD coefficient-vector r = ", wp2_format_number(rho)),
                    paste0("ULM = ", wp2_format_number(diet_stats[["ULM_score"]]),
                           "; nominal p = ", wp2_format_p(diet_stats[["ULM_p_value"]])),
                    paste0("Median |change in r| = ", wp2_format_number(stats::median(e$abs_delta_r))),
                    paste0("Sign reversals = ", wp2_format_number(100 * mean(e$sign_change), 1), "%"),
                    "Descriptive diet comparison; no change p-value")
  }
  graphics::text(-0.94, 0.95, paste(annotation, collapse = "\n"), adj = c(0, 1),
                 cex = if (compact) 0.70 else 0.83,
                 col = "#222222", xpd = FALSE)
  invisible(e)
}

wp2_plot_actual_pair <- function(edge, contexts, tag = NULL) {
  context <- as.character(edge$context[1L]); x <- contexts[[context]]
  gene <- as.character(edge$gene[1L]); feature <- as.character(edge$feature_id[1L])
  lipid_label <- as.character(edge$reference_lipid[1L])
  xs <- as.numeric(x$rna_SD[gene, ]); xw <- as.numeric(x$rna_WD[gene, ])
  ys <- as.numeric(x$lipid_SD[feature, ]); yw <- as.numeric(x$lipid_WD[feature, ])
  if (length(xs) != 7L || length(xw) != 7L || length(ys) != 7L || length(yw) != 7L) {
    stop("Selected-pair figures require seven exactly paired mice in each diet group.")
  }
  xr <- range(c(xs, xw)); yr <- range(c(ys, yw))
  pad <- function(z) z + c(-0.13, 0.13) * max(diff(z), 0.3)
  title <- paste(gene, "-", lipid_label)
  if (!is.null(tag)) title <- paste(tag, title)
  graphics::plot(c(xs, xw), c(ys, yw), type = "n", xlim = pad(xr), ylim = pad(yr),
                 xlab = "",
                 ylab = "Supplied lipid intensity (log2)", main = title,
                 cex.main = 0.90, cex.lab = 0.89)
  graphics::points(xs, ys, pch = 16, col = "#365b8c", cex = 0.9)
  graphics::points(xw, yw, pch = 17, col = "#d1812d", cex = 0.9)
  if (stats::sd(xs) > 0) graphics::abline(stats::lm(ys ~ xs), col = "#365b8c", lwd = 1.5)
  if (stats::sd(xw) > 0) graphics::abline(stats::lm(yw ~ xw), col = "#d1812d", lwd = 1.5)
  ids <- function(z) sub("^cf453-", "", colnames(z))
  graphics::text(xs, ys, labels = ids(x$rna_SD), pos = 3, offset = 0.25,
                 cex = 0.57, col = "#365b8c")
  graphics::text(xw, yw, labels = ids(x$rna_WD), pos = 1, offset = 0.25,
                 cex = 0.57, col = "#a05c12")
  graphics::legend("topleft", c("Standard diet", "Western diet"), pch = c(16, 17),
                   col = c("#365b8c", "#d1812d"), bty = "n", cex = 0.72)
  subtitle <- paste0(wp2_context_label(context), "; 1 fixed pair; n=7 SD + 7 WD mice")
  graphics::mtext(subtitle, side = 3, line = 0.12, cex = 0.73)
  graphics::mtext(paste(gene, "supplied scaled RNA expression"), side = 1, line = 2.65, cex = 0.89)
  if ("loco_same_delta_direction_fraction" %in% names(edge)) {
    note <- paste0("r_SD=", wp2_format_number(edge$r_SD), "; r_WD=", wp2_format_number(edge$r_WD),
                   "; change=", wp2_format_number(edge$delta_r),
                   "; direction stable ",
                   wp2_format_number(100 * edge$loco_same_delta_direction_fraction, 0), "%")
    graphics::mtext(note, side = 1, line = 4.05, cex = 0.68)
  }
  if (all(c("supplied_fill_fraction", "qc_cv_percent") %in% names(edge))) {
    low_fill <- is.finite(edge$supplied_fill_fraction) && edge$supplied_fill_fraction < 0.5
    note <- paste0("Supplied Fill=", wp2_format_number(edge$supplied_fill_fraction),
                    "; QC CV=", wp2_format_number(edge$qc_cv_percent, 1), "%",
                    if (low_fill) "; low-Fill sensitivity excludes this species" else "")
    graphics::mtext(note, side = 1, line = 5.05, cex = 0.66,
                    col = if (low_fill) "#a05c12" else "#555555")
  }
  invisible(NULL)
}

wp2_plot_PCA <- function(matrix, metadata, title) {
  stopifnot(identical(colnames(matrix), metadata$mouse_id))
  variance <- apply(matrix, 1L, stats::var)
  matrix <- matrix[is.finite(variance) & variance > 0, , drop = FALSE]
  pca <- stats::prcomp(t(matrix), center = TRUE, scale. = FALSE)
  fraction <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  colour <- ifelse(metadata$sex == "male", "#4173d8", "#e54141")
  symbol <- ifelse(metadata$diet == "SD", 16L, 17L)
  xr <- range(pca$x[, 1L]); yr <- range(pca$x[, 2L])
  xlim <- xr + c(-0.08, 0.08) * max(diff(xr), 1)
  # Reserve a lower strip for the legend and upper room for mouse-ID labels.
  ylim <- yr + c(-0.25, 0.12) * max(diff(yr), 1)
  graphics::plot(pca$x[, 1L], pca$x[, 2L], pch = symbol, col = colour, cex = 0.9,
                 xlim = xlim, ylim = ylim,
                 xlab = paste0("PC1 (", wp2_format_number(fraction[1L], 1), "% variance)"),
                 ylab = paste0("PC2 (", wp2_format_number(fraction[2L], 1), "% variance)"),
                 main = title, cex.main = 1.00)
  graphics::text(pca$x[, 1L], pca$x[, 2L], labels = sub("^cf453-", "", metadata$mouse_id),
                 pos = 3, offset = 0.28, cex = 0.58, col = colour)
  graphics::legend("bottomleft", c("Female SD", "Female WD", "Male SD", "Male WD"),
                   col = c("#e54141", "#e54141", "#4173d8", "#4173d8"),
                   pch = c(16, 17, 16, 17), bty = "n", cex = 0.67)
  graphics::mtext(paste0("28 mice; n=7 per sex/diet; ", nrow(matrix), " varying features"),
                  side = 3, line = 0.13, cex = 0.76)
}

#' Generate scientific figures from the completed fixed-universe analysis.
wp2_make_plots <- function(out, data, contexts, analysis, sensitivities = NULL) {
  required <- c("context", "gene", "feature_id", "reference_lipid", "r_reference",
                "primary_edge", "r_SD", "r_WD", "delta_r", "abs_delta_r", "sign_change", "eligible")
  if (!all(required %in% names(analysis$edges))) stop("Plotting requires the completed edge-analysis schema.")
  if (!all(c("inference", "comparison_summary", "loco_aggregate_summary") %in% names(analysis))) {
    stop("Plotting requires inference, diet and LOCO output tables.")
  }
  plot_dir <- file.path(out, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  context_names <- c("male_periportal", "male_pericentral", "female_periportal", "female_pericentral")
  if (!identical(names(contexts), context_names)) stop("Plot contexts must have the expected fixed order.")
  edges <- analysis$edges; inference <- analysis$inference
  if (all(c("feature_id", "Fill %", "qc_cv_percent") %in% names(data$lipid_ledger))) {
    ledger_index <- match(edges$feature_id, data$lipid_ledger$feature_id)
    if (anyNA(ledger_index)) stop("Every plotted species requires its exact annotation/QC ledger row.")
    edges$supplied_fill_fraction <- data$lipid_ledger[["Fill %"]][ledger_index]
    edges$qc_cv_percent <- data$lipid_ledger$qc_cv_percent[ledger_index]
  }
  manifest <- list()
  pdf_file <- function(name, width, height, draw, description) {
    path <- file.path(plot_dir, paste0(name, ".pdf"))
    grDevices::pdf(path, width = width, height = height, family = "Helvetica",
                    title = description, useDingbats = FALSE, onefile = TRUE)
    tryCatch(draw(), finally = grDevices::dev.off())
    manifest[[length(manifest) + 1L]] <<- data.frame(
      file = paste0("plots/", name, ".pdf"), description = description,
      stringsAsFactors = FALSE)
  }
  pdf_file("Figure6_draft", 17, 9.2, function() {
    graphics::par(mfrow = c(2, 4), mar = c(4.6, 4.0, 3.4, 0.8), oma = c(1.7, 0, 2.0, 0), cex = 0.88)
    for (i in seq_along(context_names)) wp2_plot_coefficient_panel(
      edges, inference, context_names[i], "primary", "baseline", compact = TRUE,
      tag = LETTERS[i])
    for (i in seq_along(context_names)) wp2_plot_coefficient_panel(
      edges, inference, context_names[i], "primary", "diet", compact = TRUE,
      tag = LETTERS[i + 4L])
    graphics::mtext("Figure 6 draft - original displayed gene-lipid pairs", outer = TRUE, side = 3,
                    line = 0.5, cex = 1.20, font = 2)
    graphics::mtext("Top: original FCG versus WDSD SD. Bottom: WDSD SD versus WD. Dashed line: equal coefficients.",
                    outer = TRUE, side = 1, line = 0.4, cex = 0.86)
  }, "Figure 6 draft: frozen displayed pairs, baseline conservation and diet comparison")
  for (stage in c("baseline", "diet")) {
    pdf_file(if (stage == "baseline") "baseline_conservation" else "diet_network_comparison",
             12, 10, function() {
      for (universe in c("primary", "full")) {
        graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.4, 3.8, 1.0), oma = c(1.3, 0, 2.0, 0), cex = 0.97)
        for (context in context_names) wp2_plot_coefficient_panel(
          edges, inference, context, universe, stage)
        graphics::mtext(paste(if (stage == "baseline") "Standard-diet baseline conservation" else
                               "Western versus standard diet - descriptive association changes",
                              if (universe == "primary") "| Original displayed pairs" else "| Full eligible pairs"),
                        outer = TRUE, side = 3, line = 0.3, cex = 1.05, font = 2)
        graphics::mtext(if (stage == "baseline")
                         "Positive ULM indicates conservation. Nominal p-values treat dependent edges as independent; animal-level null tests are shown beside them." else
                         "Large coefficient changes alone do not establish rewiring. Seven independent mice in each diet group; omission sensitivity is reported separately.",
                        outer = TRUE, side = 1, line = 0.1, cex = 0.76)
      }
    }, if (stage == "baseline") "Baseline conservation: original displayed and full eligible pairs" else
      "Diet comparison: original displayed and full eligible pairs")
  }
  mtmr11 <- edges[edges$context == "male_periportal" & edges$gene == "Mtmr11" & edges$eligible, , drop = FALSE]
  selected <- list()
  if (nrow(mtmr11)) {
    pdf_file("Mtmr11_coefficients", 12, 5.5, function() {
      graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.4, 3.8, 1.0), oma = c(1.3, 0, 1.5, 0))
      for (universe in c("full", "primary")) wp2_plot_coefficient_panel(
        edges, inference, "male_periportal", universe, "diet", gene = "Mtmr11")
      graphics::mtext("Mtmr11 - predefined male periportal illustration", outer = TRUE, side = 3,
                      line = 0.2, font = 2, cex = 1.10)
      graphics::mtext("Species remain fixed by the original export. The two diets contain distinct mice; points are molecular pairs, not animals.",
                      outer = TRUE, side = 1, line = 0.1, cex = 0.80)
    }, "Mtmr11 diet coefficient comparison: full and original displayed pairs")
    priority <- order(!mtmr11$primary_edge, -mtmr11$loco_same_delta_direction_fraction,
                      -mtmr11$abs_delta_r, mtmr11$reference_lipid, na.last = TRUE)
    mt_selected <- mtmr11[head(priority, 6L), , drop = FALSE]
    mt_selected$figure <- "Mtmr11_selected_pairs"
    mt_selected$selection_reason <- "Predefined gene; displayed membership first, then omission direction stability and absolute change (descriptive)."
    selected[[length(selected) + 1L]] <- mt_selected
    pdf_file("Mtmr11_selected_pairs", 13, 12.4, function() {
      graphics::par(mfrow = c(3, 2), mar = c(6.5, 4.4, 3.5, 1.2), oma = c(1.1, 0, 1.5, 0), cex = 0.93)
      for (i in seq_len(nrow(mt_selected))) wp2_plot_actual_pair(mt_selected[i, , drop = FALSE], contexts, LETTERS[i])
      graphics::mtext("Mtmr11: individual paired mouse measurements", outer = TRUE, side = 3, line = 0.3,
                      cex = 1.15, font = 2)
      graphics::mtext("Labels are CF453 mouse numbers. Lipids are whole-liver measurements. Species selection is descriptive and recorded in the figure ledger.",
                      outer = TRUE, side = 1, line = 0.1, cex = 0.76)
    }, "Mtmr11 selected molecular pairs: seven SD and seven WD mice per panel")
  }
  changed_selected <- list()
  for (context in context_names) {
    e <- edges[edges$context == context & edges$eligible, , drop = FALSE]
    robust <- e$loco_valid_combinations == 49L & e$loco_same_delta_direction_fraction == 1
    e <- e[order(!robust, -e$loco_min_abs_delta, -e$abs_delta_r, e$gene, e$reference_lipid,
                  na.last = TRUE), , drop = FALSE]
    e <- head(e, 2L)
    e$figure <- "changing_pairs"
    e$selection_reason <- "Descriptive ranking: all 49 omissions defined with stable delta direction first, then minimum absolute omitted delta and observed absolute delta."
    changed_selected[[context]] <- e
  }
  changed_selected <- do.call(rbind, changed_selected)
  selected[[length(selected) + 1L]] <- changed_selected
  pdf_file("changing_pairs", 12, 10.5, function() {
    pages <- split(seq_len(nrow(changed_selected)), ceiling(seq_len(nrow(changed_selected)) / 4L))
    for (page in pages) {
      graphics::par(mfrow = c(2, 2), mar = c(6.4, 4.5, 3.5, 1.0), oma = c(1.0, 0, 1.5, 0), cex = 0.95)
      for (i in page) wp2_plot_actual_pair(changed_selected[i, , drop = FALSE], contexts)
      graphics::mtext("Descriptively selected changing associations - mouse-level inspection", outer = TRUE,
                      side = 3, line = 0.3, font = 2, cex = 1.10)
      graphics::mtext("Outcome-based illustrative selection; no formal differential-correlation significance is assigned.",
                      outer = TRUE, side = 1, line = 0.1, cex = 0.82)
    }
  }, "Descriptively selected changing pairs and their actual mouse measurements")
  selected <- do.call(rbind, selected)
  utils::write.csv(selected, file.path(plot_dir, "selected_edges_for_figures.csv"), row.names = FALSE)
  pdf_file("QC_diagnostics", 12, 10, function() {
    graphics::par(mfrow = c(2, 2), mar = c(4.7, 4.5, 3.5, 0.9), oma = c(1.0, 0, 1.5, 0), cex = 0.96)
    wp2_plot_PCA(data$lipid, data$metadata, "Lipidomics PCA - log2 supplied intensity")
    wp2_plot_PCA(data$rna$Periportal, data$metadata, "Periportal RNA PCA - supplied scaled layer")
    wp2_plot_PCA(data$rna$Pericentral, data$metadata, "Pericentral RNA PCA - supplied scaled layer")
    cv <- data$lipid_ledger$qc_cv_percent[data$lipid_ledger$keep_primary]
    graphics::hist(cv, breaks = 35L, col = "#8eabc5", border = "white", main = "Analytical QC variability",
                   xlab = "QC CV (%) across seven injections", ylab = "Eligible lipid features")
    graphics::abline(v = 30, col = "#ba5d36", lty = 2, lwd = 1.7)
    graphics::mtext(paste0(length(cv), " retained features; ", sum(cv <= 30), " with CV <=30%"),
                    side = 3, line = 0.1, cex = 0.79)
    graphics::mtext("Diagnostic overview; no outcome-driven mouse exclusion", outer = TRUE, side = 3,
                    line = 0.2, font = 2, cex = 1.15)
    graphics::mtext("PCAs center features without additional variance scaling. QC injections are technical controls, not mice.",
                    outer = TRUE, side = 1, line = 0.1, cex = 0.80)
    graphics::par(mfrow = c(2, 3), mar = c(6.2, 4.0, 3.0, 0.8), oma = c(1.0, 0, 1.8, 0), cex = 0.92)
    group_order <- c("female-sd", "female-wd", "male-sd", "male-wd")
    for (gene in c("Cyp2e1", "Cyp2f2", "Glul", "Gls2", "Hal", "Sds")) {
      m <- data$marker_by_profile[data$marker_by_profile$gene == gene, , drop = FALSE]
      m$label <- paste(m$group, ifelse(m$zonation_type == "Periportal", "PP", "PC"))
      levels <- as.vector(t(outer(group_order, c("PP", "PC"), paste)))
      m$label <- factor(m$label, levels = levels)
      graphics::boxplot(supplied_data_expression ~ label, data = m, las = 2,
                        col = rep(c("#c7dcee", "#819ab3"), 4L), outline = FALSE,
                        main = gene, xlab = "", ylab = "Supplied RNA data layer", cex.axis = 0.68)
      graphics::stripchart(supplied_data_expression ~ label, data = m, vertical = TRUE,
                           add = TRUE, method = "jitter", pch = 16, cex = 0.47, col = "#444444")
      graphics::mtext("28 mice; n=7 per sex/diet; paired PP and PC", side = 3, line = 0.1, cex = 0.71)
    }
    graphics::mtext("Zone marker consistency in processed pseudobulks", outer = TRUE, side = 3,
                    line = 0.2, font = 2, cex = 1.12)
    graphics::mtext("Aggregate marker distributions do not reconstruct or validate the original nucleus-level gate or nucleus QC.",
                    outer = TRUE, side = 1, line = 0.1, cex = 0.78)
  }, "Lipid/RNA PCAs, seven-injection QC CV and aggregate zonation marker diagnostics")
  manifest <- do.call(rbind, manifest)
  utils::write.csv(manifest, file.path(out, "plot_manifest.csv"), row.names = FALSE)
  cat("WP2 presentation: ", nrow(manifest), " figure PDFs and ", nrow(selected),
      " selected-pair audit rows saved.\n", sep = "")
  invisible(manifest)
}

wp2_markdown_table <- function(x, labels = names(x)) {
  if (!nrow(x)) return("No eligible rows are available for this table.")
  if (length(labels) != ncol(x)) stop("Markdown table labels do not match its columns.")
  cell <- function(z) {
    z <- as.character(z); z[is.na(z)] <- "NA"
    z <- gsub("[\r\n]+", " ", z)
    gsub("|", "/", z, fixed = TRUE)
  }
  values <- lapply(x, cell)
  rows <- vapply(seq_len(nrow(x)), function(i)
    paste0("| ", paste(vapply(values, function(z) z[i], character(1L)), collapse = " | "), " |"),
    character(1L))
  c(paste0("| ", paste(labels, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(labels)), collapse = " | "), " |"), rows)
}

#' Write a self-contained Markdown report from saved completed output tables.
#' No result is inferred from a desired biological narrative or a plot alone.
wp2_make_report <- function(out) {
  paths <- c(inference = "baseline_inference.csv", comparison = "diet_comparison_summary.csv",
             loco = "loco_aggregate_summary.csv", edges = "network_edges.csv",
             replay = "reference_replay_summary.csv", qc = "input_qc_summary.csv",
             coverage = "network_coverage.csv", identity = "lipid_identity_crosswalk.csv",
             eligibility = "edge_eligibility_ledger.csv", sensitivities = "sensitivities.csv",
             metadata = "metadata.csv", markers = "zonation_marker_summary.csv",
             lipid_ledger = "lipid_annotation_and_QC_ledger.csv")
  absent <- paths[!file.exists(file.path(out, paths))]
  if (length(absent)) stop("Report requires completed saved tables: ", paste(absent, collapse = ", "))
  tables <- lapply(paths, function(path) utils::read.csv(file.path(out, path),
                                                       check.names = FALSE, stringsAsFactors = FALSE))
  inference <- tables$inference; edges <- tables$edges
  required_inference <- c("context", "universe", "scope", "gene", "n_edges", "n_SD", "n_WD",
                          "coefficient_vector_r", "ULM_score", "ULM_p_value", "pairing_permutation_p",
                          "null_bootstrap_p", "exact_pairing_permutations", "permutation_valid_draws",
                          "bootstrap_valid_draws", "null_bootstrap_draws", "paired_valid_draws",
                          "paired_stability_draws", "paired_same_direction_fraction",
                          "paired_reference_valid_draws", "paired_reference_same_direction_fraction")
  if (!all(required_inference %in% names(inference))) stop("Report inference schema differs from the analysis contract.")
  contexts <- c("male_periportal", "male_pericentral", "female_periportal", "female_pericentral")
  primary <- inference[inference$scope == "network" & inference$universe == "primary", , drop = FALSE]
  full <- inference[inference$scope == "network" & inference$universe == "full", , drop = FALSE]
  primary <- primary[match(contexts, primary$context), , drop = FALSE]
  full <- full[match(contexts, full$context), , drop = FALSE]
  if (nrow(primary) != 4L || nrow(full) != 4L || anyNA(primary$context) || anyNA(full$context)) {
    stop("Four primary and four full network summaries must be complete before reporting.")
  }
  metric <- function(name) {
    index <- match(name, tables$qc$metric)
    if (is.na(index)) stop("Missing reported QC metric: ", name)
    as.character(tables$qc$value[index])
  }
  baseline_table <- function(x) {
    wp2_markdown_table(data.frame(
      Context = wp2_context_label(x$context), Edges = as.integer(x$n_edges),
      `Vector r` = wp2_format_number(x$coefficient_vector_r),
      ULM = wp2_format_number(x$ULM_score),
      `Nominal ULM p` = wp2_format_p(x$ULM_p_value),
      `Exact pairing p` = wp2_format_p(x$pairing_permutation_p),
      `Null-bootstrap p` = wp2_format_p(x$null_bootstrap_p), check.names = FALSE))
  }
  replay <- tables$replay
  replay <- replay[match(contexts, replay$context), , drop = FALSE]
  replay_table <- wp2_markdown_table(data.frame(
    Context = wp2_context_label(replay$context), `FCG mice` = replay$n_mice,
    Genes = replay$n_genes, `Full reference edges` = replay$n_full_edges,
    `Displayed reference edges` = replay$n_primary_edges,
    `Eligible full edges` = as.integer(full$n_edges),
    `Eligible displayed edges` = as.integer(primary$n_edges), check.names = FALSE))
  valid_table <- wp2_markdown_table(data.frame(
    Context = wp2_context_label(primary$context),
    `Pairings valid / total` = paste(primary$permutation_valid_draws, primary$exact_pairing_permutations, sep = " / "),
    `Null bootstrap valid / total` = paste(primary$bootstrap_valid_draws, primary$null_bootstrap_draws, sep = " / "),
    `Paired stability valid / total` = paste(primary$paired_valid_draws, primary$paired_stability_draws, sep = " / "),
    `Paired same-direction fraction` = wp2_format_number(primary$paired_same_direction_fraction),
    `Reference-paired same-direction fraction` = wp2_format_number(primary$paired_reference_same_direction_fraction),
    check.names = FALSE))
  diet_table <- function(universe) {
    rows <- lapply(contexts, function(context) {
      e <- edges[edges$context == context & edges$eligible &
                   (universe == "full" | edges$primary_edge), , drop = FALSE]
      c <- tables$comparison[tables$comparison$context == context &
                               tables$comparison$universe == universe & tables$comparison$scope == "network", , drop = FALSE]
      if (nrow(c) != 1L || !nrow(e)) stop("Missing diet summary for fixed eligible pairs.")
      stable_direction <- e$loco_valid_combinations == 49L & e$loco_same_delta_direction_fraction == 1
      stable_reversal <- e$loco_valid_combinations == 49L & e$loco_sign_reversal_fraction == 1
      data.frame(Context = wp2_context_label(context), Edges = nrow(e),
                 `Median absolute change` = wp2_format_number(stats::median(e$abs_delta_r)),
                 `Sign reversal %` = wp2_format_number(100 * mean(e$sign_change), 1),
                 `Direction stable in all 49 %` = wp2_format_number(100 * mean(stable_direction), 1),
                 `Reversal in all 49 %` = wp2_format_number(100 * mean(stable_reversal), 1), check.names = FALSE)
    })
    wp2_markdown_table(do.call(rbind, rows))
  }
  diet_association_table <- function(universe) {
    x <- tables$comparison[tables$comparison$universe == universe &
                             tables$comparison$scope == "network", , drop = FALSE]
    x <- x[match(contexts, x$context), , drop = FALSE]
    wp2_markdown_table(data.frame(Context = wp2_context_label(x$context), Edges = x$n_edges,
                                  `SD-WD vector r` = wp2_format_number(x$coefficient_vector_r_SD_WD),
                                  `Descriptive ULM` = wp2_format_number(x$ULM_score_SD_WD),
                                  `Nominal ULM p` = wp2_format_p(x$ULM_p_value_SD_WD), check.names = FALSE))
  }
  mt_inference <- inference[inference$scope == "gene" & !is.na(inference$gene) &
                              inference$gene == "Mtmr11", , drop = FALSE]
  mt_edges <- edges[edges$context == "male_periportal" & edges$gene == "Mtmr11" & edges$eligible, , drop = FALSE]
  mt_baseline_text <- if (nrow(mt_inference)) {
    unlist(lapply(c("primary", "full"), function(universe) {
      x <- mt_inference[mt_inference$universe == universe, , drop = FALSE]
      c("", paste0("### Mtmr11 baseline: ", if (universe == "primary") "original displayed pairs" else "full eligible pairs"),
        "", baseline_table(x))
    }), use.names = FALSE)
  } else "No Mtmr11 subset has sufficient eligible pairs for a native ULM summary."
  mt_comparison <- tables$comparison[tables$comparison$scope == "gene" &
                                      !is.na(tables$comparison$gene) & tables$comparison$gene == "Mtmr11", , drop = FALSE]
  mt_diet_association <- if (nrow(mt_comparison)) wp2_markdown_table(data.frame(
    Universe = ifelse(mt_comparison$universe == "primary", "Displayed", "Full"),
    Edges = mt_comparison$n_edges,
    `SD-WD vector r` = wp2_format_number(mt_comparison$coefficient_vector_r_SD_WD),
    `Descriptive ULM` = wp2_format_number(mt_comparison$ULM_score_SD_WD),
    `Nominal ULM p` = wp2_format_p(mt_comparison$ULM_p_value_SD_WD), check.names = FALSE)) else
      "No sufficient Mtmr11 diet coefficient-vector summary."
  mt_diet_changes <- if (nrow(mt_edges)) {
    rows <- lapply(c("primary", "full"), function(universe) {
      e <- mt_edges[universe == "full" | mt_edges$primary_edge, , drop = FALSE]
      stable <- e$loco_valid_combinations == 49L & e$loco_same_delta_direction_fraction == 1
      reversed <- e$loco_valid_combinations == 49L & e$loco_sign_reversal_fraction == 1
      data.frame(Universe = if (universe == "primary") "Displayed" else "Full",
                 Edges = nrow(e), `Median absolute change` = wp2_format_number(stats::median(e$abs_delta_r)),
                 `Sign reversals` = paste(sum(e$sign_change), nrow(e), sep = " / "),
                 `Direction stable in all 49` = paste(sum(stable), nrow(e), sep = " / "),
                 `Reversal in all 49` = paste(sum(reversed), nrow(e), sep = " / "), check.names = FALSE)
    })
    wp2_markdown_table(do.call(rbind, rows))
  } else "No eligible Mtmr11 pairs for diet-change summaries."
  mt_primary <- mt_edges[mt_edges$primary_edge, , drop = FALSE]
  mt_ledger_index <- match(mt_primary$feature_id, tables$lipid_ledger$feature_id)
  if (anyNA(mt_ledger_index)) stop("Mtmr11 species require complete annotation/QC identities.")
  mt_primary$supplied_fill_fraction <- tables$lipid_ledger[["Fill %"]][mt_ledger_index]
  mt_primary$qc_cv_percent <- tables$lipid_ledger$qc_cv_percent[mt_ledger_index]
  mt_species_table <- wp2_markdown_table(data.frame(
    Lipid = mt_primary$reference_lipid, `r SD` = wp2_format_number(mt_primary$r_SD),
    `r WD` = wp2_format_number(mt_primary$r_WD), `Change in r` = wp2_format_number(mt_primary$delta_r),
    `Supplied Fill` = wp2_format_number(mt_primary$supplied_fill_fraction),
    `QC CV %` = wp2_format_number(mt_primary$qc_cv_percent, 1),
    `Reversal all 49` = ifelse(mt_primary$loco_valid_combinations == 49L &
                                mt_primary$loco_sign_reversal_fraction == 1, "Yes", "No"),
    check.names = FALSE))
  mt_primary_inference <- mt_inference[mt_inference$universe == "primary", , drop = FALSE]
  mt_ordering_note <- if (nrow(mt_primary)) paste0(
    "The ", nrow(mt_primary), " original displayed Mtmr11 coefficients occupy a narrow range (",
    wp2_format_number(min(mt_primary$r_reference)), " to ", wp2_format_number(max(mt_primary$r_reference)),
    "); ", sum(mt_primary$r_reference > 0), " are positive originally and ", sum(mt_primary$r_SD > 0),
    " remain positive in WDSD SD. Its coefficient-vector statistic concerns relative ordering of these associations. ",
    "A negative vector association does not imply that all seven gene-lipid associations reverse sign, nor opposite biological regulation. ",
    if (nrow(mt_primary_inference)) paste0("The ordinary paired bootstrap preserves the observed vector direction in ",
      wp2_format_number(100 * mt_primary_inference$paired_same_direction_fraction, 1), "% of valid draws; including paired original-reference resampling gives ",
      wp2_format_number(100 * mt_primary_inference$paired_reference_same_direction_fraction, 1),
      "% (both stability descriptions, not null p-values).") else "") else
      "Mtmr11 has no eligible originally displayed pairs."
  mt_change_note <- if (nrow(mt_primary)) {
    stable_negative <- mt_primary$loco_valid_combinations == 49L &
      mt_primary$loco_same_delta_direction_fraction == 1 & mt_primary$delta_r < 0
    persistent_reversals <- mt_primary$reference_lipid[
      mt_primary$loco_valid_combinations == 49L & mt_primary$loco_sign_reversal_fraction == 1]
    paste0(sum(stable_negative), " of ", nrow(mt_primary),
           " displayed Mtmr11 pairs have negative SD-to-WD delta direction in all 49 omissions. ",
           length(persistent_reversals), " of ", nrow(mt_primary),
           " remain actual sign reversals in all 49 omissions",
           if (length(persistent_reversals)) paste0(": ", paste(persistent_reversals, collapse = "; ")) else "",
           ". Persistent direction changes and persistent sign reversals are different observations; neither is a formal rewiring test.")
  } else "No Mtmr11 displayed pairs are available for omission analysis."
  mt_low_fill_reversals <- mt_primary[mt_primary$loco_valid_combinations == 49L &
                                       mt_primary$loco_sign_reversal_fraction == 1 &
                                       mt_primary$supplied_fill_fraction < 0.5, , drop = FALSE]
  mt_high_fill_stable <- mt_primary[mt_primary$loco_valid_combinations == 49L &
                                     mt_primary$loco_same_delta_direction_fraction == 1 &
                                     mt_primary$delta_r < 0 & mt_primary$supplied_fill_fraction >= 0.5, , drop = FALSE]
  mt_fill_note <- if (nrow(mt_low_fill_reversals)) paste0(
    "QC qualification: ", paste(mt_low_fill_reversals$reference_lipid, collapse = "; "),
    " is the displayed Mtmr11 sign-reversal illustration with low supplied Fill (",
    paste(wp2_format_number(mt_low_fill_reversals$supplied_fill_fraction), collapse = "; "),
    ") despite low analytical QC CV (", paste(wp2_format_number(mt_low_fill_reversals$qc_cv_percent, 1), collapse = "; "),
    "%). It is excluded by the Fill >=0.5 sensitivity, so its apparent persistent reversal has an additional input-detection/provenance limitation. ",
    nrow(mt_high_fill_stable), " higher-Fill displayed PE pairs retain negative delta direction across all 49 omissions: ",
    paste(mt_high_fill_stable$reference_lipid, collapse = "; "),
    ". Their persistent delta direction does not mean all are persistent sign reversals. The low-Fill pair remains in the frozen primary analysis, with its quality information reported explicitly.") else
    "No displayed persistent Mtmr11 sign reversal falls below the supplied Fill >=0.5 sensitivity cutoff."
  selected_path <- file.path(out, "plots/selected_edges_for_figures.csv")
  selected <- if (file.exists(selected_path)) utils::read.csv(selected_path, check.names = FALSE) else NULL
  selected_table <- if (!is.null(selected) && any(selected$figure == "changing_pairs")) {
    s <- selected[selected$figure == "changing_pairs", , drop = FALSE]
    wp2_markdown_table(data.frame(
      Context = wp2_context_label(s$context), Gene = s$gene, Lipid = s$reference_lipid,
      `r SD` = wp2_format_number(s$r_SD), `r WD` = wp2_format_number(s$r_WD),
      `Change in r` = wp2_format_number(s$delta_r),
      `Omission direction fraction` = wp2_format_number(s$loco_same_delta_direction_fraction), check.names = FALSE))
  } else "Selected-pair figures have not been generated; no illustrative pairs are reported."
  limits_path <- file.path(out, "input_limitations.txt")
  if (!file.exists(limits_path)) stop("Report requires saved input limitations.")
  limitations <- readLines(limits_path, warn = FALSE)
  limitations <- limitations[nzchar(trimws(limitations))]
  n_positive <- sum(primary$coefficient_vector_r > 0)
  n_pairing <- sum(primary$coefficient_vector_r > 0 & primary$pairing_permutation_p <= 0.05)
  n_bootstrap <- sum(primary$coefficient_vector_r > 0 & primary$null_bootstrap_p <= 0.05)
  strongest <- primary[which.max(primary$coefficient_vector_r), , drop = FALSE]
  interpretation <- paste0(
    n_positive, " of four displayed-edge comparisons have positive coefficient-vector association. ",
    n_pairing, " of four have positive association and exact pairing p <=0.05; ",
    n_bootstrap, " of four have positive association and null-bootstrap p <=0.05. ",
    "These counts summarize four predefined contexts and do not define a universal conservation threshold. ",
    "The largest displayed-edge point estimate is ", wp2_context_label(strongest$context),
    " (vector r=", wp2_format_number(strongest$coefficient_vector_r), ", ULM=",
    wp2_format_number(strongest$ULM_score), ", nominal p=", wp2_format_p(strongest$ULM_p_value),
    ", exact pairing p=", wp2_format_p(strongest$pairing_permutation_p),
    ", null-bootstrap p=", wp2_format_p(strongest$null_bootstrap_p), ").")
  sensitivity_text <- if (nrow(tables$sensitivities)) {
    s <- tables$sensitivities[tables$sensitivities$scope == "network" &
                               tables$sensitivities$universe == "primary", , drop = FALSE]
    unlist(lapply(unique(s$sensitivity), function(name) {
      x <- s[s$sensitivity == name, , drop = FALSE]
      x <- x[match(contexts, x$context), , drop = FALSE]
      label <- switch(name, QC_CV_le30 = "QC CV <=30%", supplied_fill_ge50 = "Supplied Fill % >=0.5", name)
      c("", paste0("### ", label, ": original displayed pairs"), "", baseline_table(x))
    }), use.names = FALSE)
  } else "No eligible prespecified sensitivity rows were available."
  fill <- tables$lipid_ledger[["Fill %"]]
  if (is.null(fill) || !is.numeric(fill) || !any(is.finite(fill))) stop("Missing supplied Fill % audit.")
  fill_note <- paste0("The supplied `Fill %` field ranges from ", wp2_format_number(min(fill, na.rm = TRUE)),
                       " to ", wp2_format_number(max(fill, na.rm = TRUE)), "; ",
                       sum(is.finite(fill) & fill < 0.5), " supplied rows and ",
                       sum(tables$lipid_ledger$keep_primary & is.finite(fill) & fill < 0.5),
                       " retained rows have values below 0.5 despite positive exported intensities. ",
                       "Their underlying detection, missingness and possible imputation provenance is unavailable. ",
                       "The Fill % >=0.5 restriction was added after intake inspection and before association analysis as an additional QC sensitivity.")
  coverage <- tables$coverage[match(contexts, tables$coverage$context), , drop = FALSE]
  shared_lipids <- unique(coverage$n_eligible_full_lipids)
  if (length(shared_lipids) != 1L || any(coverage$n_missing_gene_edges != 0) ||
      any(coverage$n_original_full_genes != coverage$n_eligible_full_genes)) {
    stop("Report coverage statement requires the audited common lipid universe and complete reference-gene coverage.")
  }
  coverage_note <- paste0(shared_lipids, " of 1,254 original molecular annotations have strict eligible matches (",
                           1254L - shared_lipids, " lack a strict match), giving ", sum(full$n_edges),
                           " eligible full pairs and ", sum(primary$n_edges), " eligible displayed pairs. ",
                           "All originally selected context genes are measurable and varying in both diets. ",
                           "Here feature loss is driven by exact species matching, rather than absent genes or invariant observed measurements. ",
                           "Unmatched annotations are retained in the crosswalk; no lower-resolution species/chain collapse is used to force coverage.")
  stageB_primary <- tables$loco[tables$loco$universe == "primary", , drop = FALSE]
  stageB_primary <- stageB_primary[match(contexts, stageB_primary$context), , drop = FALSE]
  stageB_full <- tables$loco[tables$loco$universe == "full", , drop = FALSE]
  stageB_full <- stageB_full[match(contexts, stageB_full$context), , drop = FALSE]
  persistent_primary <- edges[edges$eligible & edges$primary_edge, , drop = FALSE]
  n_persistent_reversal <- sum(persistent_primary$loco_valid_combinations == 49L &
                                persistent_primary$loco_sign_reversal_fraction == 1)
  stageB_note <- paste0("Displayed-pair median absolute SD-to-WD changes are ",
                         wp2_format_number(stageB_primary$median_abs_delta_r[3L]), " (female periportal) and ",
                         wp2_format_number(stageB_primary$median_abs_delta_r[4L]), " (female pericentral), compared with ",
                         wp2_format_number(stageB_primary$median_abs_delta_r[1L]), " and ",
                         wp2_format_number(stageB_primary$median_abs_delta_r[2L]), " in the corresponding male contexts. ",
                         n_persistent_reversal, " of ", nrow(persistent_primary),
                         " displayed pairs remain sign reversals across all 49 joint mouse omissions. ",
                         "The corresponding full-universe median changes are ",
                         paste(wp2_format_number(stageB_full$median_abs_delta_r), collapse = ", "),
                         " in male PP, male PC, female PP, female PC order. ",
                         "The apparent female excess among displayed pairs therefore depends on the fixed feature selection. ",
                         "These are descriptive follow-up candidates, not evidence of sex-specific causal rewiring.")
  analysis_path <- file.path(out, "network_analysis.rds")
  if (!file.exists(analysis_path)) stop("Report requires saved analysis settings.")
  settings <- readRDS(analysis_path)$metadata
  report <- c(
    "# WP2: WDSD baseline conservation and Western-diet association changes",
    "", "Analysis date: 7 October 2026. This report implements work package 2 of the 5 October meeting plan.",
    "", "## Main findings", "", interpretation,
    "", stageB_note,
    "", "The analysis assesses whether the original Figure 5 gene-lipid coefficient patterns transfer to standard-diet CF453 mice, then describes SD-to-WD changes within CF453. Nominal ULM results are reported beside mouse-level null tests. Diet changes are descriptive and assessed for omission stability; they are not established differential correlations or causal rewiring.",
    "", "## Cohort, processing and QC", "",
    paste0("The updated mouse-ID crosswalk pairs all ", metric("biological_mice"),
           " CF453 mice between modalities, with seven female SD, seven female WD, seven male SD and seven male WD mice. ",
           "There are ", metric("hepatocyte_profiles"), " whole-hepatocyte profiles and ",
           metric("zonation_profiles"), " paired periportal/pericentral profiles. Each mouse supplies one whole-liver lipid measurement shared by its two RNA zones; 56 profiles do not create 56 independent animals."),
    "", "The meeting described approximately 10-month-old mice after nine months of dietary exposure, compared with an approximately three-month-old original cohort, and a different genetic background. These descriptions are meeting context; exact ages, strain/cross and exposure records are not present in the supplied analysis metadata.",
    "", paste0("The lipid export contains ", metric("supplied_lipid_feature_rows"),
               " annotation rows and ", metric("analytical_QC_injections"), " analytical QC injections. Exclusions comprise ",
               metric("internal_standard_exclusions"), " internal standards, ", metric("additional_isotope_label_exclusions"),
               " additional isotope-labeled features and ", metric("repository_class_exclusions"),
               " repository-excluded-class rows, leaving ", metric("primary_lipid_features"), " primary features. ",
               metric("QC_CV_le30_lipid_features"), " eligible features have QC CV <=30%; ",
               metric("primary_features_QC_CV_gt30"), " exceed 30%. Median QC CV is ",
               wp2_format_number(as.numeric(metric("primary_QC_CV_median")), 1), "% and maximum is ",
               wp2_format_number(as.numeric(metric("primary_QC_CV_maximum")), 1), "%. The QC-CV cutoff is a sensitivity, not a primary outcome-driven filter."),
    "", "Primary lipid values are log2 supplied positive intensities. No extra normalization or imputation is applied. Primary RNA values are the supplied `scale.data` layer. Row scaling the supplied `data` layer over its complete original profile roster reproduces `scale.data` numerically; this checks scaling, rather than validating the upstream normalization. PCA and aggregate marker distributions are diagnostics, not grounds for excluding animals because they weaken conservation.",
    "", fill_note,
    "", paste0("Expected Cyp2e1/Cyp2f2 marker ordering occurs in ", metric("profiles_expected_marker_order"),
               " of 56 aggregate zone profiles. ", metric("profiles_aggregate_marker_threshold"),
               " satisfy the conventional threshold when applied to the stored aggregate data layer. This aggregate check cannot reconstruct the original nucleus-level gate."),
    "", paste0("Per-mouse and per-zone nucleus counts are unavailable for all ",
               metric("mice_with_missing_nucleus_counts"), " mice. Supplied count-library totals and detected genes are reported, but are not nucleus contributions. The recorded command history contains only ScaleData; upstream nuclear QC, doublets, aggregation, missing-value treatment and normalization provenance remain unresolved."),
    "", "QC source tables: [sample pairing](metadata.csv), [input QC metrics](input_qc_summary.csv), [zone marker summary](zonation_marker_summary.csv), [input limitations](input_limitations.txt). [QC figures](plots/QC_diagnostics.pdf) show lipid/RNA PCAs, analytical variability and six marker distributions.",
    "", "## Frozen Figure 5 reference and molecular coverage", "", replay_table,
    "", coverage_note,
    "", paste0("Original Pearson coefficients are independently reconstructed from 22 testes-bearing or 21 ovary-bearing mice. All ",
               sum(replay$n_full_edges), " archived gene-lipid edges replay using log2 lipid intensities, with maximum absolute error ",
               format(max(replay$log2_max_abs_archive_difference), scientific = TRUE, digits = 3),
               ". Raw intensities do not replay (maximum context discrepancies ",
               paste(wp2_format_number(replay$raw_max_abs_archive_difference), collapse = ", "), ")."),
    "", paste0("The original full exports freeze context-specific genes and 1,254 lipidomics features; primary membership is taken directly from original displayed SIF files (documented absolute original r >0.65). The ",
               sum(replay$n_excluded_TF_gene_edges), " TF-to-gene edges are excluded. Pitx3 is also a regulon source and carries an ATT TF label, but its measured RNA-to-lipid branch is correctly retained."),
    "", "Matching preserves gene symbols and each original lipid's chain order, slash/underscore separator, ether prefix, oxygenation and structural resolution. Only the final pipe-delimited annotation and class formatting are harmonized. Formatting agreement does not independently establish molecular identity. Alignment ID is not treated as a unique key; every annotation row has a feature ID. Missing, excluded, ambiguous or invariant measurements remain in the complete audit ledgers.",
    "", "The original displayed-edge and full-edge universes are frozen before WDSD correlations are examined. Eligibility requires the same measured, varying RNA/lipid features in both SD and WD. No WDSD r cutoff is imposed. Original selection of strong associations makes displayed-edge conservation conditional on the original selection; the complementary full universe reduces this particular selection restriction, while retaining the originally selected genes.",
    "", "Coverage and identity audits: [reference replay](reference_replay_summary.csv), [network coverage](network_coverage.csv), [lipid crosswalk](lipid_identity_crosswalk.csv), [edge eligibility](edge_eligibility_ledger.csv).",
    "", "## Stage A: standard-diet conservation", "",
    "Each context compares original reference coefficients with the matching sex and zone in seven standard-diet mice. Positive coefficient-vector Pearson r and positive native ULM slope t indicate conserved direction. Negative values indicate inverse alignment. ULM p-values are nominal, because edges share genes, lipids and seven animals; the molecular edge count is not the biological sample size.",
    "", "### Primary: original displayed pairs", "", baseline_table(primary),
    "", "### Complementary: all eligible original full-network pairs", "", baseline_table(full),
    "", "The primary exact pairing test exhausts all 7!=5,040 RNA-lipid pairings within each sex's SD group, including the observed identity pairing. One pairing is shared across all genes, lipids and both zones. Its null is lack of mouse-pairing information, stronger than simply no alignment with this specific reference. The approximate empirical-marginal null-bootstrap sensitivity independently resamples complete RNA-mouse and lipid-mouse vectors and tests their independence conditional on the observed empirical margins. Both p-values use two-sided absolute coefficient-vector association; at fixed edge count its ordering is equivalent to absolute ULM. The bootstrap uses `(1 + exceedances)/(1 + valid draws)`; the exact permutation uses its exhaustive tail proportion.",
    "", "These animal-level tests retain dependence among molecular features. They remain conditional, small-cohort tests and do not establish regulation, substrate specificity or a mechanistic pathway. A small two-sided p-value with negative alignment would indicate inverse alignment, rather than positive conservation. Neither null allows arbitrary RNA-lipid associations while asking only whether they align with the reference. Exact pairings retain all seven observed mice per margin; empirical bootstrap samples contain about 4.62 distinct mice per modality on average. Its tails can be coarse, unstable or undefined at this group size. A higher or lower bootstrap p-value does not establish conservatism; type-I error and power have not been calibrated for this WDSD bootstrap procedure.",
    "", "Exact inference assumes exchangeability of RNA-lipid pairing under the independence null within each sex's SD group. Batch and litter covariates are unavailable, so pairing tests do not independently control those possible shared sources of variation. Four predefined network contexts and secondary focus-gene summaries are exploratory; nominal edge-based p-values and unadjusted mouse-level p-values do not establish multiplicity-adjusted discoveries.",
    "", "An ordinary paired-mouse bootstrap resamples matched RNA/lipid mice together and describes observed conservation stability. An additional paired bootstrap resamples original FCG mice within sex-chromosome/gonad groups and WDSD paired mice independently, sharing draws across zones. Neither supplies a null p-value. Original selected edges stay fixed in these draws, so uncertainty in their original selection is not propagated. The fractions below are descriptive stability proportions, not probabilities that the biological interpretation is true.",
    "", valid_table,
    "", "Undefined constant-feature draws invalidate a fixed-universe summary instead of silently changing its edges. Valid and invalid counts for every universe and focus gene are saved in [baseline inference](baseline_inference.csv). [Baseline figures](plots/baseline_conservation.pdf) retain original displayed and full eligible results on separate pages.",
    "", "## Stage B: Western versus standard diet", "",
    "The two diet groups contain different mice. Every predefined eligible pair has `delta_r = r_WD - r_SD`; absolute change and sign reversal are exported alongside both coefficients. All 49 combinations dropping one SD and one WD mouse are evaluated, with the same omission applied across both zones. Direction stability and reversal stability require all 49 combinations to be defined, rather than counting a fraction of one from a smaller valid subset.",
    "", "### Primary displayed-pair association", "", diet_association_table("primary"),
    "", "These native ULM statistics describe SD-WD coefficient-vector association. Their nominal p-values test the fitted edge-vector slope under independent residual assumptions; they are not differential-correlation or rewiring tests. No diet-bootstrap p-value is assigned.",
    "", "### Primary displayed-pair changes", "", diet_table("primary"),
    "", "### Full eligible-pair association", "", diet_association_table("full"),
    "", "### Full eligible-pair changes", "", diet_table("full"),
    "", "A large absolute change is a descriptive observation at n=7 per group. Stability under omissions makes that observation less dependent on one mouse; it does not make it statistically significant or prove rewiring. No differential-correlation p-value is assigned. A simple diet-label shuffle would also test diet differences in the joint expression/lipid distributions, including abundance, rather than specifically equality of correlations.",
    "", "The strongest illustrative pairs below are chosen descriptively: complete 49-omission delta-direction stability first, then minimum omitted absolute delta and observed absolute delta. This figure selection is recorded; it is not used to select the analyzed universe or to infer significance.",
    "", selected_table,
    "", "[Diet coefficient figures](plots/diet_network_comparison.pdf), [selected mouse measurements](plots/changing_pairs.pdf), [complete coefficient/change table](network_edges.csv) and [LOCO summaries](loco_aggregate_summary.csv) retain the species-level evidence.",
    "", "## Predefined Mtmr11 focus", "",
    paste0("Mtmr11 is present only in the frozen male periportal reference: 1,254 full original pairs and nine originally displayed pairs. ",
           nrow(mt_edges), " full pairs and ", sum(mt_edges$primary_edge), " displayed pairs are eligible in WDSD."),
    "", mt_baseline_text,
    "", mt_ordering_note,
    "", "### Mtmr11 diet coefficient-vector association", "", mt_diet_association,
    "", "### Mtmr11 diet changes and omission stability", "", mt_diet_changes,
    "", mt_change_note,
    "", "### Displayed Mtmr11 species and analytical quality", "", mt_species_table,
    "", mt_fill_note,
    "", "Mtmr11 was specified before examining WDSD results. Its coefficient plots compare full and displayed species without requiring their new correlations to exceed the original cutoff. Up to six species are illustrated using actual scaled RNA and log2 whole-liver lipid values from seven mice per diet, prioritizing original displayed membership and then descriptive omission stability/change. The source identities and selection rule are in [the figure ledger](plots/selected_edges_for_figures.csv).",
    "", if (nrow(mt_edges)) "[Mtmr11 coefficient comparison](plots/Mtmr11_coefficients.pdf) and [individual mouse measurements](plots/Mtmr11_selected_pairs.pdf) support visual inspection; they do not measure enzyme activity, lipid flux or zonal lipid abundance." else "Mtmr11 was unavailable for the paired measurement figures.",
    "", "## Prespecified sensitivities", "",
    "The QC-CV <=30% feature restriction was prespecified, and the supplied Fill % >=0.5 restriction was added after input inspection before association analysis. Both rerun native ULM, exact pairings and the independent-marginal null bootstrap on the explicitly restricted original pairs. Neither is independent biological validation. Different denominators can themselves change coefficient-vector statistics; outcomes do not choose the preferred filter. The supplied RNA data layer is an affine replay of the scaled layer, documented in [RNA layer sensitivity](RNA_data_layer_sensitivity.csv), rather than independent normalization evidence.",
    "", sensitivity_text,
    "", "The [complete sensitivity table](sensitivities.csv) records retained universes and statistics. No confidence interval, reversed concordance score or multiplicity-adjusted discovery claim is needed for this exploratory report.",
    "", "## Interpretation limits and manuscript use", "",
    "Baseline conservation and within-cohort diet changes answer different questions. Where baseline evidence is inconclusive, a difference from the original network cannot be attributed specifically to Western diet; age, background, processing and sampling differences remain possible. The within-CF453 SD/WD comparison addresses diet-associated coefficient changes more directly, but does not establish their causality or reproducibility beyond this cohort.",
    "", "The Figure 6 draft uses the frozen original displayed pairs in its baseline row and its SD-versus-WD row. Full-universe results, animal-level uncertainty and QC remain available beside the draft. Draft labels identify sample counts, feature coverage, ULM summaries and the explicit null tests, rather than treating edge-level nominal significance as biological replication.",
    "", "Recorded input limitations:", "", paste0("- ", limitations),
    "", "## Artifacts and reproduction", "",
    "- [Figure 6 draft](plots/Figure6_draft.pdf): original displayed baseline and diet comparisons across four contexts.",
    "- [Plot manifest](plot_manifest.csv): scientific figure descriptions and file locations.",
    "- [Baseline inference](baseline_inference.csv), [diet summaries](diet_comparison_summary.csv), [complete species edges](network_edges.csv), [LOCO summaries](loco_aggregate_summary.csv) and [sensitivity results](sensitivities.csv).",
    "- [Mouse metadata](metadata.csv), [lipid identities](lipid_identity_crosswalk.csv), [edge eligibility](edge_eligibility_ledger.csv), [network coverage](network_coverage.csv) and [reference replay](reference_replay_summary.csv).",
    "", paste0("Recorded seed: ", settings$seed, "; independent-marginal null draws: ", settings$B,
                 "; paired stability draws: ", settings$paired_B, "; exact pairings per sex: ",
                 settings$exact_permutations, ". Shared mouse indices and per-omission IDs are retained."),
    "", "Run from the repository root, choosing a fresh output directory:", "", "```sh",
    "Rscript scripts/workflows/05b_wdsd_network_comparison.R --validate-only",
    paste0("Rscript scripts/workflows/05b_wdsd_network_comparison.R --bootstraps=", settings$B,
           " --paired-bootstraps=", settings$paired_B, " --seed=", settings$seed,
           " --output=results/wdsd_network_comparison_reproduction"), "```", "",
    "The entry point records settings, seeds, shared mouse indices, input/code checksums and package/session versions. Original FCG results and supplied objects remain the fixed inputs.")
  # On the first workflow run these audit files are produced after the report;
  # on presentation-only regeneration they are available and linked directly.
  audit_files <- c("run_manifest.json", "input_and_code_SHA256.csv", "sessionInfo.txt")
  if (all(file.exists(file.path(out, audit_files)))) {
    report <- c(report, "", "Run audit: [settings and scope](run_manifest.json), [input/code SHA256](input_and_code_SHA256.csv), [session versions](sessionInfo.txt).")
  }
  # Resolve artifact links at authoring time, so report completion cannot mask
  # a missing plotted or tabulated deliverable.
  links <- regmatches(report, gregexpr("\\]\\(([^)]+)\\)", report, perl = TRUE))
  targets <- unlist(lapply(links, function(z) sub("^\\]\\((.*)\\)$", "\\1", z)), use.names = FALSE)
  local_targets <- targets[!grepl("^[a-z]+://", targets)]
  if (any(!file.exists(file.path(out, local_targets)))) {
    stop("Report references absent artifacts: ", paste(local_targets[!file.exists(file.path(out, local_targets))], collapse = ", "))
  }
  writeLines(report, file.path(out, "README.md"), useBytes = TRUE)
  cat("WP2 output-derived scientific report saved: ", file.path(out, "README.md"), "\n", sep = "")
  invisible(file.path(out, "README.md"))
}
