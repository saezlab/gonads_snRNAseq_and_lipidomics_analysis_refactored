# Output-derived plots/report for the lipid-overlap and Figure 3 follow-up.
# Source support_functions.R and wdsd_network_statistics.R first.
# Sourcing this helper creates no artifacts.

wp2p_number <- function(x, digits = 3L) vapply(x, function(z) {
  if (!is.finite(z)) "NA" else formatC(z, format = "f", digits = digits)
}, character(1L))

wp2p_pvalue <- function(x) vapply(x, function(z) {
  if (!is.finite(z)) return("NA")
  if (z == 0) return("<2.23e-308")
  if (z < 0.001) formatC(z, format = "e", digits = 2L) else formatC(z, format = "f", digits = 4L)
}, character(1L))

wp2p_inference <- function(x, universe = "primary", analysis = "pooled", gene = NULL) {
  ix <- x$universe == universe & x$analysis == analysis
  ix <- ix & if (is.null(gene)) x$scope == "network" else
    x$scope == "gene" & !is.na(x$gene) & x$gene == gene
  answer <- x[which(ix), , drop = FALSE]
  if (nrow(answer) > 1L) stop("Duplicate pooled presentation inference.")
  answer
}

wp2p_plot_coefficients <- function(e, row, label, centered = FALSE) {
  x <- e[[if (centered) "r_reference_withinsex" else "r_reference"]]
  y <- e[[if (centered) "r_SD_withinsex" else "r_SD"]]
  if (!nrow(e) || any(!is.finite(x)) || any(!is.finite(y))) {
    graphics::plot.new(); graphics::title(main = label)
    graphics::text(0.5, 0.5, "Fixed-vector summary undefined", cex = 0.85)
    return(invisible(NULL))
  }
  colour <- if (centered) "#7561a8" else "#4173d8"
  graphics::plot(x, y, xlim = c(-1, 1), ylim = c(-1, 1), pch = 16,
                 col = grDevices::adjustcolor(colour, if (nrow(e) > 1000L) 0.15 else 0.5),
                 cex = if (nrow(e) > 1000L) 0.4 else 0.65,
                 xlab = if (centered) "FCG within-sex-centered Pearson r" else "FCG pooled Pearson r",
                 ylab = if (centered) "SD within-sex-centered Pearson r" else "SD pooled Pearson r",
                 main = label, cex.main = 0.93, cex.lab = 0.85,
                 xaxs = "i", yaxs = "i", axes = FALSE)
  graphics::axis(1, seq(-1, 1, .5)); graphics::axis(2, seq(-1, 1, .5), las = 1)
  graphics::box(col = "#555555")
  graphics::abline(h = 0, v = 0, col = "#dddddd", lty = 3)
  graphics::abline(0, 1, col = "#888888", lty = 2)
  if (stats::sd(x) > 0) graphics::abline(stats::lm(y ~ x), col = colour, lwd = 1.5)
  z <- wp2_vector_stats(x, y)
  annotation <- c(paste0("Vector r = ", wp2p_number(z[["coefficient_vector_r"]])),
                   paste0("ULM = ", wp2p_number(z[["ULM_score"]]), "; nominal p = ", wp2p_pvalue(z[["ULM_p_value"]])))
  if (!centered && nrow(row)) annotation <- c(annotation,
    paste0("Overall pairing p = ", wp2p_pvalue(row$unrestricted_permutation_absolute_p)),
    paste0("Overall bootstrap p = ", wp2p_pvalue(row$unrestricted_bootstrap_absolute_p)),
    paste0("Within-sex pairing upper p = ", wp2p_pvalue(row$withinsex_permutation_upper_p)),
    paste0("Within-sex bootstrap upper p = ", wp2p_pvalue(row$withinsex_bootstrap_upper_p)))
  if (centered) annotation <- c(annotation, "Descriptive residual-network sensitivity")
  if (!centered && !nrow(row)) annotation <- c(annotation, "Descriptive selected-gene complement")
  graphics::text(-0.94, 0.95, paste(annotation, collapse = "\n"), adj = c(0, 1), cex = .64)
  graphics::mtext(paste0(nrow(e), " frozen pairs; 43 FCG / 14 SD mice"), side = 3, line = .10, cex = .64)
}

wp2p_plot_null <- function(draws, row, column, method, conditional = FALSE) {
  values <- draws[[column]]
  if ("draw" %in% names(draws)) values <- values[draws$draw != 0L]
  values <- values[is.finite(values)]
  if (!length(values)) { graphics::plot.new(); graphics::title(main = method); return(invisible(NULL)) }
  observed <- row$coefficient_vector_r
  window <- if (conditional) range(c(values, observed)) else c(-1, 1)
  if (conditional) window <- pmax(-1, pmin(1, window + c(-.12, .12) * max(diff(window), .001)))
  breaks <- if (conditional) pretty(window, n = 30L) else seq(-1, 1, length.out = 41L)
  graphics::hist(values, breaks = breaks, col = "#b8c8df", border = "white",
                 xlim = window, main = method, xlab = "Null coefficient-vector Pearson r", cex.main = 1,
                 cex.lab = .9, ylab = "Whole-mouse draws")
  graphics::abline(v = observed, col = "#c83737", lwd = 2)
  graphics::abline(v = stats::median(values), col = "#555555", lty = 2, lwd = 1.4)
  prefix <- if (conditional) "withinsex" else "unrestricted"
  bank <- if (grepl("bootstrap", method, fixed = TRUE)) "bootstrap" else "permutation"
  p <- row[[paste(prefix, bank, if (conditional) "upper_p" else "absolute_p", sep = "_")]]
  graphics::legend("topleft", c(paste0("Observed r = ", wp2p_number(observed)),
                                paste0("Null median = ", wp2p_number(stats::median(values))),
                                paste0(if (conditional) "Upper-tail p = " else "Absolute-tail p = ", wp2p_pvalue(p)),
                                paste0(length(values), " valid null draws")),
                   bty = "n", cex = .78)
  graphics::mtext(if (conditional) "Sex groups preserved; asks about excess alignment beyond sex" else
                   "Across all 14 mice; shared sex effects can cause rejection", side = 3, line = .12, cex = .66)
}

wp2p_plot_mouse_pair <- function(edge, context) {
  gene <- as.character(edge$gene); feature <- as.character(edge$feature_id)
  x <- as.numeric(context$rna_SD[gene, ]); y <- as.numeric(context$lipid_SD[feature, ])
  meta <- context$metadata_SD
  if (!identical(colnames(context$rna_SD), meta$mouse_id) ||
      !identical(colnames(context$lipid_SD), meta$mouse_id) || length(x) != 14L) {
    stop("Selected Figure 3 pair requires the same 14 ordered SD animals in both modalities.")
  }
  colours <- ifelse(meta$sex == "female", "#e54141", "#4173d8")
  symbols <- ifelse(meta$sex == "female", 16L, 17L)
  pad <- function(a) range(a) + c(-.13, .13) * max(diff(range(a)), .3)
  ylimits <- range(y) + c(-.30, .13) * max(diff(range(y)), .3)
  label <- paste(gene, "-", edge$reference_lipid)
  graphics::plot(x, y, col = colours, pch = symbols, cex = 1,
                 xlim = pad(x), ylim = ylimits, main = label,
                 xlab = paste(gene, "supplied scaled RNA"), ylab = "Whole-liver lipid intensity (log2)",
                 cex.main = .92, cex.lab = .88)
  # Give closely clustered animal labels separate vertical positions, with
  # leader lines retaining the actual measurement identity.
  for (sex in c("female", "male")) {
    ix <- which(meta$sex == sex)
    ix <- ix[order(y[ix], x[ix])]
    label_y <- y[ix] + .025 * diff(ylimits)
    for (j in seq_along(label_y)[-1L]) label_y[j] <- max(label_y[j], label_y[j - 1L] + .035 * diff(ylimits))
    label_y <- label_y - max(0, max(label_y) - (ylimits[2L] - .025 * diff(ylimits)))
    label_x <- x[ix] + .025 * diff(pad(x))
    graphics::segments(x[ix], y[ix], label_x, label_y, col = grDevices::adjustcolor(colours[ix][1L], .4))
    graphics::text(label_x, label_y, labels = sub("^cf453-", "", meta$mouse_id[ix]),
                   adj = c(0, .5), cex = .62, col = colours[ix])
  }
  graphics::abline(stats::lm(y ~ x), lwd = 1.5, col = "#444444", lty = 2)
  for (sex in c("female", "male")) {
    ix <- meta$sex == sex
    if (stats::sd(x[ix]) > 0) {
      model <- stats::lm(y[ix] ~ x[ix])
      xs <- range(x[ix])
      graphics::lines(xs, stats::coef(model)[1L] + stats::coef(model)[2L] * xs, lwd = 1.4,
                       col = if (sex == "female") "#e54141" else "#4173d8")
    }
  }
  graphics::legend("bottomleft", c("Female SD (7)", "Male SD (7)", "Pooled fit"),
                   col = c("#e54141", "#4173d8", "#444444"), pch = c(16, 17, NA),
                   lty = c(NA, NA, 2), bty = "n", cex = .71)
  graphics::mtext("14 SD mice; labels give CF453 mouse-ID suffixes", side = 3, line = .1, cex = .7)
  graphics::mtext(paste0("FCG pooled r=", wp2p_number(edge$r_reference), "; SD pooled r=", wp2p_number(edge$r_SD),
                        "; SD sex-centered r=", wp2p_number(edge$r_SD_withinsex)), side = 1, line = 3.7, cex = .71)
  graphics::mtext("Selected by original displayed membership and original |r| only", side = 1, line = 4.6, cex = .67)
}

#' Create figures from completed pooled and molecular-overlap analyses.
wp2_pooled_make_plots <- function(out, data, context, analysis, audit, sensitivities = NULL) {
  required <- c("gene", "feature_id", "reference_lipid", "r_reference", "r_SD", "primary_edge", "eligible",
                "r_reference_withinsex", "r_SD_withinsex")
  if (!all(required %in% names(analysis$edges))) stop("Pooled plots require completed fixed-edge coefficients.")
  dir.create(file.path(out, "plots"), recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(out, "plots")
  manifest <- list()
  add <- function(filename, description, pages) {
    manifest[[length(manifest) + 1L]] <<- data.frame(file = file.path("plots", filename),
                                                  description = description, pages = pages, stringsAsFactors = FALSE)
  }
  e <- analysis$edges[analysis$edges$eligible, , drop = FALSE]
  universes <- list(primary = e[e$primary_edge, ],
                    full = e)
  labels <- c(primary = "Figure 3 displayed pairs", full = "Full sparse-selection precursor")
  grDevices::pdf(file.path(plot_dir, "pooled_conservation.pdf"), width = 10.5, height = 8.9, useDingbats = FALSE)
  graphics::par(mfrow = c(2, 2), mar = c(4.0, 3.9, 3.3, 1), oma = c(1.8, .2, 2.5, .2))
  for (centered in c(FALSE, TRUE)) for (u in names(universes)) {
    row <- if (u %in% c("primary", "full")) wp2p_inference(analysis$inference, u,
                        if (centered) "withinsex_centered" else "pooled") else data.frame()
    wp2p_plot_coefficients(universes[[u]], row, labels[[u]], centered)
  }
  graphics::mtext("Pooled whole-hepatocyte Figure 3 replication - frozen molecular pairs", outer = TRUE,
                  side = 3, line = 1.1, cex = 1.1, font = 2)
  graphics::mtext("Top: pooled sexes. Bottom: each cohort's sex means removed; nominal ULM is descriptive under dependent edges.",
                  outer = TRUE, side = 1, line = 0, cex = .75)
  grDevices::dev.off()
  add("pooled_conservation.pdf", "Fixed Figure 3 displayed and full sparse universes; pooled and sex-centered coefficients", 1L)
  primary <- wp2p_inference(analysis$inference)
  grDevices::pdf(file.path(plot_dir, "pooled_null_checks.pdf"), width = 10.5, height = 8.7, useDingbats = FALSE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3.3, 1), oma = c(1.8, 0, 2.5, 0))
  key <- as.character(primary$specification)
  wp2p_plot_null(analysis$permutation_unrestricted_draws, primary, key, "Unrestricted pairing permutations")
  wp2p_plot_null(analysis$bootstrap_unrestricted_draws, primary, key, "Unrestricted null bootstrap")
  wp2p_plot_null(analysis$permutation_withinsex_draws, primary, key, "Within-sex pairing permutations", TRUE)
  wp2p_plot_null(analysis$bootstrap_withinsex_draws, primary, key, "Within-sex null bootstrap", TRUE)
  graphics::mtext(paste0("Original displayed pairs: ", primary$n_edges, " eligible pairs; 14 SD animals"),
                  outer = TRUE, side = 3, line = 1.1, cex = 1.1, font = 2)
  graphics::mtext("Red: observed. Grey: null median. Lower panels use a closer horizontal scale. Stability draws do not supply p-values.",
                  outer = TRUE, side = 1, line = 0, cex = .71)
  grDevices::dev.off()
  add("pooled_null_checks.pdf", "Whole-mouse null checks distinguish overall correspondence from alignment beyond shared sex structure", 1L)
  focus <- c("Cyp2b9", "Cyp2b13", "Cyp2c38", "Cyp2d9")
  grDevices::pdf(file.path(plot_dir, "Cyp_pooled_conservation.pdf"), width = 10.5, height = 8.7, useDingbats = FALSE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3.2, 1), oma = c(1.8, 0, 2.5, 0))
  for (gene in focus) {
    row <- wp2p_inference(analysis$inference, "primary", gene = gene)
    wp2p_plot_coefficients(e[e$primary_edge & e$gene == gene, , drop = FALSE], row,
                            paste(gene, "- original displayed pairs"))
  }
  graphics::mtext("Four predefined Cyp genes - original Figure 3 species membership", outer = TRUE,
                  side = 3, line = 1.1, cex = 1.1, font = 2)
  graphics::mtext("Native ULM and animal-level checks are exploratory; these are pooled association patterns, not evidence of causal regulation.",
                  outer = TRUE, side = 1, line = 0, cex = .73)
  grDevices::dev.off()
  add("Cyp_pooled_conservation.pdf", "Predefined four Cyp genes on original displayed species pairs", 1L)
  selected <- list()
  for (gene in focus) {
    candidate <- e[e$gene == gene & e$primary_edge, , drop = FALSE]
    if (!nrow(candidate)) next
    candidate <- candidate[order(-abs(candidate$r_reference), candidate$reference_lipid), , drop = FALSE]
    chosen <- candidate[1L, , drop = FALSE]
    chosen$selection_rule <- "Highest absolute ORIGINAL FCG correlation among matched originally displayed pairs for this predefined Cyp gene; alphabetical original lipid tie-break"
    selected[[gene]] <- chosen
  }
  selected <- do.call(rbind, selected)
  if (!is.null(selected) && nrow(selected)) {
    utils::write.csv(selected, file.path(plot_dir, "figure_pair_selection_ledger.csv"), row.names = FALSE)
    grDevices::pdf(file.path(plot_dir, "selected_pooled_pairs.pdf"), width = 10.5, height = 9.2, useDingbats = FALSE)
    graphics::par(mfrow = c(2, 2), mar = c(6.4, 4.1, 3.2, 1), oma = c(.3, 0, 2.3, 0))
    for (i in seq_len(nrow(selected))) wp2p_plot_mouse_pair(selected[i, , drop = FALSE], context)
    graphics::mtext("Individual SD animals for four predefined Cyp examples", outer = TRUE,
                    side = 3, line = 1, cex = 1.1, font = 2)
    grDevices::dev.off()
    add("selected_pooled_pairs.pdf", "Actual SD measurements selected from original Cyp pairs without using new correlation strength", 1L)
  }
  old <- audit$original_ledger
  old$plot_category <- ifelse(old$overlap_category == "exact_final_annotation", "Same final annotation",
                             ifelse(old$overlap_category == "no_exported_annotation_or_same_class_formula_candidate",
                                    "No recorded candidate", "Coarser or unresolved candidate"))
  categories <- c("Same final annotation", "Coarser or unresolved candidate", "No recorded candidate")
  cols <- c("#4173d8", "#d89b43", "#a9acb3")
  totals <- sort(table(old$lipid_class), decreasing = TRUE)
  top <- names(totals)[seq_len(min(14L, length(totals)))]
  old$plot_class <- ifelse(old$lipid_class %in% top, old$lipid_class, "Other classes")
  ct <- table(factor(old$plot_category, levels = categories),
               factor(old$plot_class, levels = c(top, "Other classes")))
  grDevices::pdf(file.path(plot_dir, "lipid_overlap_audit.pdf"), width = 10.5, height = 8.4, useDingbats = FALSE)
  graphics::par(mar = c(5.6, 6, 4, 1), oma = c(.8, 0, 0, 0))
  positions <- graphics::barplot(ct[, ncol(ct):1L, drop = FALSE], horiz = TRUE, las = 1, col = cols,
                                 border = NA, xlab = "Original measured lipid rows", cex.names = .82,
                                 main = "Original-to-WDSD overlap by lipid class",
                                 xlim = c(0, max(colSums(ct)) * 1.12))
  graphics::text(colSums(ct)[ncol(ct):1L], positions, labels = colSums(ct)[ncol(ct):1L], pos = 4, cex = .74, xpd = TRUE)
  graphics::legend("bottomright", categories, fill = cols, bty = "n", cex = .78)
  graphics::mtext("Candidate categories do not establish an additional one-to-one chemical identity match.", side = 1, line = 4.7, cex = .78)
  category <- audit$category_summary
  category <- category[order(category$n_original_rows), , drop = FALSE]
  short <- c(exact_final_annotation = "Same final annotation", shared_coarse_annotation_alias_only = "Shared coarse alias only",
              unordered_chain_CL_candidate_unresolved = "Unresolved CL chain assignment",
              position_resolution_candidate_only = "Position / resolution change",
              library_identifier_only_structure_mismatch = "Library ID; structure mismatch",
              same_class_formula_only_possible_isomers = "Same formula; possible isomers",
              no_exported_annotation_or_same_class_formula_candidate = "No exported candidate")
  graphics::par(mar = c(5.6, 13, 4, 1), oma = c(.8, 0, 0, 0))
  positions <- graphics::barplot(category$n_original_rows, names.arg = unname(short[category$overlap_category]),
                                 horiz = TRUE, las = 1, col = "#b8c8df", border = NA, cex.names = .8,
                                 xlab = "Original measured lipid rows", main = "Complete partition of the 1,254 original lipid labels",
                                 xlim = c(0, max(category$n_original_rows) * 1.13))
  graphics::text(category$n_original_rows, positions, category$n_original_rows, pos = 4, cex = .84)
  graphics::mtext(paste0(nrow(audit$safe_map), " strict same-annotation pairs; ",
                         nrow(audit$conditional_compatibility_map), " pairs in conditional coarse-alias sensitivity"),
                  side = 1, line = 4.7, cex = .8)
  grDevices::dev.off()
  add("lipid_overlap_audit.pdf", "Class coverage and complete molecular-annotation loss partition, separating safe joins from candidates", 2L)
  manifest <- do.call(rbind, manifest)
  utils::write.csv(manifest, file.path(out, "plot_manifest.csv"), row.names = FALSE)
  invisible(manifest)
}

wp2p_markdown_table <- function(x) {
  if (!nrow(x)) return("No eligible fixed summary.")
  x[] <- lapply(x, function(z) gsub("\\|", "/", as.character(z)))
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1L, function(z) paste0("| ", paste(z, collapse = " | "), " |")))
}

#' Write a self-contained scientific addendum from the completed output tables.
wp2_pooled_make_report <- function(out) {
  read <- function(file) utils::read.csv(file.path(out, file), check.names = FALSE, stringsAsFactors = FALSE)
  inference <- read("pooled_inference.csv")
  edges <- read("matched_pooled_edges.csv")
  replay <- read("pooled_reference_replay.csv")
  coverage <- read("pooled_coverage.csv")
  analysis <- readRDS(file.path(out, "pooled_analysis.rds"))
  audit_path <- file.path(out, "lipid_overlap_audit", "lipid_overlap_audit.rds")
  if (file.exists(audit_path)) audit <- readRDS(audit_path) else {
    read_audit <- function(nm) read(file.path("lipid_overlap_audit", paste0(nm, ".csv")))
    audit <- lapply(c("original_ledger", "new_ledger", "safe_map", "conditional_compatibility_map",
                      "category_summary", "stage_summary", "class_summary"), read_audit)
    names(audit) <- c("original_ledger", "new_ledger", "safe_map", "conditional_compatibility_map",
                      "category_summary", "stage_summary", "class_summary")
  }
  primary <- wp2p_inference(inference)
  full <- wp2p_inference(inference, "full")
  centered <- inference[inference$analysis == "withinsex_centered" & inference$scope == "network", ]
  old <- audit$original_ledger
  strict <- nrow(audit$safe_map); total <- nrow(old)
  aliases <- sum(old$recorded_alias_candidates > 0L)
  strict_coverage <- coverage[coverage$mapping == "strict" & coverage$universe == "primary", ]
  conditional_coverage <- coverage[coverage$mapping == "conditional_compatibility" & coverage$universe == "primary", ]
  coverage_note <- paste0("For the actual Figure 3 displayed network, the strict join retains ",
                          strict_coverage$matched_lipids, "/", strict_coverage$original_lipids,
                          " lipids, ", strict_coverage$matched_edges, "/", strict_coverage$original_edges,
                          " displayed pairs and ", strict_coverage$matched_genes, "/", strict_coverage$original_genes,
                          " displayed genes. The conditional compatibility sensitivity retains ",
                          conditional_coverage$matched_lipids, "/", conditional_coverage$original_lipids,
                          " lipids, ", conditional_coverage$matched_edges, "/", conditional_coverage$original_edges,
                          " pairs and ", conditional_coverage$matched_genes, "/", conditional_coverage$original_genes,
                          " genes. These Figure 3 denominators are separate from the 752/1,254 whole-lipid-survey overlap.")
  strict_qc <- if ("keep_qc30" %in% names(audit$new_ledger)) {
    sum(audit$safe_map$feature_id %in% audit$new_ledger$feature_id[audit$new_ledger$keep_qc30])
  } else NA_integer_
  qc_overlap_note <- paste0("Repository exclusions remove zero of the 752 strict matches. The later QC-CV <=30% sensitivity retains ",
                            strict_qc, " of those 752 survey-wide matches, excluding ", strict - strict_qc,
                            " already matched rows. All 36 strict Figure 3 lipids pass that cutoff, so its pooled network is unchanged. ",
                            "The 502 initially absent exact annotations therefore cannot be attributed to these analysis exclusions.")
  missing_classes <- audit$class_coverage[order(-audit$class_coverage$n_original_unmatched), ]
  missing_classes <- head(missing_classes, 5L)
  class_loss_note <- paste0("The largest missing-annotation contributions are ",
                            paste(paste0(missing_classes$lipid_class, " (", missing_classes$n_original_unmatched, ")"), collapse = ", "),
                            ": ", sum(missing_classes$n_original_unmatched), " of the 502 strict-unmatched rows. ",
                            "These counts include candidate and candidate-free categories; they are not experimentally confirmed detection losses.")
  denominator_note <- character()
  if (file.exists(file.path(out, "QA/raw_FCG_excluded_class_summary.csv"))) {
    partition <- read("QA/raw_FCG_excluded_class_summary.csv")
    excluded_classes <- sum(partition$n_raw_rows[partition$descriptive_partition == "class_matches_repository_exclusion_list"])
    other_rows <- sum(partition$n_raw_rows[partition$descriptive_partition == "other_raw_row_absent_from_cleaned_reference"])
    denominator_note <- paste0("The 1,254-species denominator is the pre-existing cleaned reference, rather than every row in the original 1,386-row export. ",
                              "The ", sum(partition$n_raw_rows), " raw rows absent from that reference include ", excluded_classes,
                              " whose classes match the repository exclusion list. The remaining ", other_rows,
                              " are 18 SM rows with ontology ASM, two Cer_EOS rows and one Cer_EODS row. ",
                              "This recorded inventory establishes which rows were removed before this follow-up, but does not establish every historical filtering reason. ",
                              "[Raw-to-cleaned ledger](QA/raw_FCG_1386_to_cleaned_1254_ledger.csv) and [excluded-class inventory](QA/raw_FCG_excluded_class_summary.csv) retain the full partition.")
  }
  sex_mean_note <- character()
  if (file.exists(file.path(out, "QA/sex_mean_only_concordance.csv"))) {
    means <- read("QA/sex_mean_only_concordance.csv")
    means <- means[match(c("primary", "full"), means$universe), ]
    sex_mean_note <- paste0("A [descriptive sex-mean diagnostic](QA/sex_mean_only_concordance.csv) keeps only the between-sex covariance in SD, normalized by observed pooled feature variation. ",
                            "Its alignment with the original FCG coefficient vector is ",
                            wp2p_number(means$sex_mean_only_expected_coefficient_vector_r[1L], 4),
                            " for displayed pairs and ", wp2p_number(means$sex_mean_only_expected_coefficient_vector_r[2L], 4),
                            " for the precursor, close to observed pooled alignment ",
                            wp2p_number(means$observed_vector_r[1L], 4), " and ", wp2p_number(means$observed_vector_r[2L], 4),
                            ". This shows why shared sex structure can generate high pooled correspondence. It supplies no new p-value and does not prove that within-sex associations are absent.")
  }
  category_labels <- c(exact_final_annotation = "Same final molecular annotation",
    unordered_chain_CL_candidate_unresolved = "Unresolved CL chain-order candidate",
    position_resolution_candidate_only = "Positional / structural-resolution candidate",
    shared_coarse_annotation_alias_only = "Shared coarse annotation alias only",
    library_identifier_only_structure_mismatch = "Library identifier only; structure mismatch",
    same_class_formula_only_possible_isomers = "Same class / formula only; possible isomers",
    no_exported_annotation_or_same_class_formula_candidate = "No exported annotation or same-class formula candidate")
  cat_table <- data.frame(Category = unname(category_labels[audit$category_summary$overlap_category]),
                         `Original rows` = audit$category_summary$n_original_rows, check.names = FALSE)
  table_stats <- function(rows) wp2p_markdown_table(data.frame(
    Universe = ifelse(rows$universe == "primary", "Figure 3 displayed", "Full sparse precursor"),
    Pairs = rows$n_edges, `Vector r` = wp2p_number(rows$coefficient_vector_r),
    ULM = wp2p_number(rows$ULM_score), `Nominal ULM p` = wp2p_pvalue(rows$ULM_p_value), check.names = FALSE))
  null_table <- function(rows) wp2p_markdown_table(data.frame(
    Universe = ifelse(rows$universe == "primary", "Figure 3 displayed", "Full sparse precursor"),
    `Overall pairing p` = wp2p_pvalue(rows$unrestricted_permutation_absolute_p),
    `Overall bootstrap p` = wp2p_pvalue(rows$unrestricted_bootstrap_absolute_p),
    `Within-sex pairing upper p` = wp2p_pvalue(rows$withinsex_permutation_upper_p),
    `Within-sex bootstrap upper p` = wp2p_pvalue(rows$withinsex_bootstrap_upper_p), check.names = FALSE))
  rows <- rbind(primary, full)
  core_note <- paste0("The original displayed Figure 3 coefficient pattern has pooled FCG-to-SD vector r=",
                      wp2p_number(primary$coefficient_vector_r), " and native ULM=", wp2p_number(primary$ULM_score),
                      " (nominal p=", wp2p_pvalue(primary$ULM_p_value), "). The overall mouse-pairing p is ",
                      wp2p_pvalue(primary$unrestricted_permutation_absolute_p), " and overall null-bootstrap p is ",
                      wp2p_pvalue(primary$unrestricted_bootstrap_absolute_p), ". Those overall tests allow shared sex differences to explain agreement. ",
                      "When sex margins are retained, the upper-tail pairing p is ", wp2p_pvalue(primary$withinsex_permutation_upper_p),
                      " and upper-tail bootstrap p is ", wp2p_pvalue(primary$withinsex_bootstrap_upper_p), ". ",
                      if (is.finite(primary$withinsex_permutation_upper_p) && primary$withinsex_permutation_upper_p < .05)
                        "There is exploratory evidence of positive alignment beyond pairing expected from sex margins in this fixed comparison; multiplicity and bootstrap calibration remain unresolved."
                      else "These conditional checks do not establish positive alignment beyond that expected from the shared sex structure.")
  strict_note <- paste0("The original 752/1,254 figure is reproducible: ", strict, " original species (",
                        wp2p_number(100 * strict / total, 1), "%) have the same final molecular annotation in the new export. ",
                        "This is a deliberately strict molecular-label overlap, not an estimate of how many lipids were physically detected. ",
                        aliases, " original labels share at least one recorded pipe-delimited alias with a retained new row; ",
                        "these aliases often describe only sum composition. The cumulative audit reaches ",
                        audit$stage_summary$original_rows_covered_retained_new[audit$stage_summary$rule == "plus_recorded_alias_candidates"],
                        " after also allowing chain-order and positional-resolution candidates. These broader counts are candidates, not validated additional molecular identities.")
  focus <- inference[inference$analysis == "pooled" & inference$scope == "gene" & inference$universe == "primary", ]
  focus <- focus[match(c("Cyp2b9", "Cyp2b13", "Cyp2c38", "Cyp2d9"), focus$gene), , drop = FALSE]
  focus <- focus[!is.na(focus$gene), , drop = FALSE]
  focus_table <- wp2p_markdown_table(data.frame(Gene = focus$gene, Pairs = focus$n_edges,
    `Vector r` = wp2p_number(focus$coefficient_vector_r), ULM = wp2p_number(focus$ULM_score),
    `Nominal ULM p` = wp2p_pvalue(focus$ULM_p_value),
    `Within-sex pairing upper p` = wp2p_pvalue(focus$withinsex_permutation_upper_p),
    `Within-sex bootstrap upper p` = wp2p_pvalue(focus$withinsex_bootstrap_upper_p), check.names = FALSE))
  sensitivity_text <- "The complete sensitivity outputs are linked below; none is used to replace the strict primary comparison."
  sensitivity_paths <- c("mapping_sensitivity_inference.csv", "pooled_QC_sensitivity.csv")
  sensitivity_summaries <- list()
  for (p in sensitivity_paths[file.exists(file.path(out, sensitivity_paths))]) {
    s <- read(p)
    if (all(c("scope", "analysis", "universe", "n_edges", "ULM_score", "ULM_p_value") %in% names(s))) {
      s <- s[s$scope == "network" & s$analysis == "pooled", , drop = FALSE]
      if (nrow(s)) {
        label_fields <- intersect(c("mapping", "sensitivity", "sensitivity_name", "restriction"), names(s))
        label <- if (length(label_fields)) as.character(s[[label_fields[1L]]]) else rep(sub(".csv", "", p, fixed = TRUE), nrow(s))
        human_labels <- c(strict = "Strict annotation", conditional_compatibility = "Conditional alias map",
                           QC_CV_le30 = "QC CV <=30%", Fill_ge0.5 = "Fill >=50%", fill_ge50 = "Fill >=50%",
                           Fill_ge50 = "Fill >=50%", supplied_fill_ge0.5 = "Fill >=50%", supplied_fill_ge50 = "Fill >=50%")
        replaced <- unname(human_labels[label])
        label[!is.na(replaced)] <- replaced[!is.na(replaced)]
        sensitivity_summaries[[p]] <- data.frame(Sensitivity = label,
          Universe = ifelse(s$universe == "primary", "Displayed", "Sparse precursor"), Pairs = s$n_edges,
          ULM = wp2p_number(s$ULM_score), `Nominal ULM p` = wp2p_pvalue(s$ULM_p_value),
          `Within-sex pairing upper p` = if ("withinsex_permutation_upper_p" %in% names(s)) wp2p_pvalue(s$withinsex_permutation_upper_p) else "Not rerun",
          `Within-sex bootstrap upper p` = if ("withinsex_bootstrap_upper_p" %in% names(s)) wp2p_pvalue(s$withinsex_bootstrap_upper_p) else "Not rerun",
          check.names = FALSE)
      }
    }
  }
  if (length(sensitivity_summaries)) sensitivity_text <- wp2p_markdown_table(do.call(rbind, sensitivity_summaries))
  limitations <- if (file.exists(file.path(out, "input_limitations.txt"))) readLines(file.path(out, "input_limitations.txt")) else character()
  full_intake_note <- grepl("Twenty-eight mice remain", limitations, fixed = TRUE)
  limitations[full_intake_note] <- paste("The complete WDSD intake contains 28 mice with two RNA zones and one whole-liver lipid profile each.",
                                        "This pooled follow-up uses 14 standard-diet mice, whole-hepatocyte RNA and no zone profiles.")
  crows <- coverage
  report <- c("# WP2 follow-up: lipid annotation overlap and pooled Figure 3 conservation", "",
    "Analysis date: 7 October 2026. This addendum addresses the two requested follow-ups to WP2. Original results and the preceding Figure 5 / diet report remain unchanged.",
    "", "## Main findings", "", strict_note, "", core_note,
    "", "## 1. Why does the lipid overlap appear limited?", "",
    "The two exports use differing combinations of molecular-species labels, sum-composition aliases, sn-position notation and library identifiers. The original exported table has 1,386 rows; its cleaned network reference retains 1,254 measured rows. WDSD supplies 1,240 rows; the established exclusion of 14 internal standards and 49 rows in repository-excluded classes leaves 1,177 primary features. Multiple measured/adduct rows are never collapsed or summed to improve overlap.",
    "", denominator_note,
    "", qc_overlap_note,
    "", wp2p_markdown_table(cat_table), "",
    "Categories are a priority partition, so each original row appears once. A row can have more than one candidate in the detailed ledger. In particular, the 202 coarse-alias-only rows are the largest potentially resolvable group after labels already assigned to the chain-order and positional categories. The 280 rows without a recorded annotation or same-class formula candidate represent a real exported-feature difference; their absence does not independently establish absence of a chromatographic peak or abundance below detection.",
    "", class_loss_note,
    "", "Three concrete cautions explain why a broad match is not automatically safe:", "",
    "- CL chain quartets can be written in a different underscore order, but several original or new CL rows share the same unordered quartet. CL pairing/group assignment and peak identity remain unresolved; chain sorting does not provide a globally unique rescue.",
    "- A shared alias can describe a sum such as a total carbon / double-bond composition while the final annotations resolve different acyl chains or ether/plasmalogen structure. It supports compatibility at that coarser resolution, rather than equality of molecular species.",
    "- Exported InChIKey or SMILES fields can agree despite conflicting or differently resolved final annotations. The CL example `CL(16:0_16:1_18:1_18:2)` versus `CL(16:0_18:1_18:1_16:2)` changes reported chain constituents. `PE(O-18:1_18:2)` versus `PE(P-18:0_18:2)` may describe the same plasmalogen composition at different structural resolution because P- notation counts the vinyl bond differently. Shared library identifiers alone do not settle these unresolved structural differences.",
    "", paste0("A separate compatibility sensitivity extends the join to ", nrow(audit$conditional_compatibility_map),
                " original/new pairs: the ", strict, " strict pairs plus ", nrow(audit$conditional_compatibility_map) - strict,
                " recorded coarse-alias candidates. Each additional candidate has the same lipid class and adduct, precursor agreement within 5 ppm, and bidirectional degree one across the complete candidate graph. Already assigned strict features are blocked. This protects against feature reuse; it does not validate molecular identity. Those pairs remain explicitly conditional matches at their reported annotation resolution."),
    "", "Mass, retention time, formula and library fields are exported diagnostics. Cross-run retention times have not been calibrated and raw spectra have not been reidentified. Resolving the remaining candidate groups requires the collaborator's original peak/annotation evidence, especially acyl-chain fragments, adduct provenance and calibrated retention matching. No experimental-detection conclusion is drawn from the annotation audit alone.",
    "", "The lipid notation distinguishes molecular species, positional resolution and more general composition. The matching rules follow [updated lipid shorthand nomenclature](https://pmc.ncbi.nlm.nih.gov/articles/PMC7707175/), the [LIPID MAPS structural hierarchy](https://www.lipidmaps.org/lipid_nomenclature/rules/hierarchy) and [glycerophospholipid notation rules](https://www.lipidmaps.org/lipid_nomenclature/glycerophospholipids/pl_lpl). The safe join changes formatting only; chemically broader matches stay separate.",
    "", "Audits: [complete original annotation ledger](lipid_overlap_audit/original_ledger.csv), [new annotation ledger](lipid_overlap_audit/new_ledger.csv), [every candidate pair](lipid_overlap_audit/candidate_pairs.csv), [class-level coverage](lipid_overlap_audit/class_summary.csv), [strict safe map](lipid_overlap_audit/safe_map.csv), [conditional compatibility map](lipid_overlap_audit/conditional_compatibility_map.csv) and [structure/identifier conflict examples](lipid_overlap_audit/identifier_structure_mismatch_examples.csv). [Overlap figures](plots/lipid_overlap_audit.pdf) show the loss partition and class distribution.",
    "", "## 2. Which Figure 3 network was compared?", "",
    "The current revised Figure 3 A states top three correlated genes per lipid, and panel B labels that species display FULL SPECIES-LEVEL NETWORK. Its Workflow 04 integrated export has 267 edges: 192 species-level gene-lipid edges involving 25 displayed genes and 64 lipids, plus 75 TF-gene edges. This analysis uses the 192 species pairs and excludes TF edges. The figure's use of FULL refers to molecular species resolution; it is distinct from the permissive 33%-tail archive containing 769 genes, 828 lipids and 636,732 pairs. That separate archive is inventoried, rather than substituted for Figure 3.",
    "", "Original selection uses MOFA Factor 1's top 2.5% in each tail, a lipid-metabolism gene filter, and the three largest absolute gene correlations per selected lipid. It does not impose a 0.65 correlation cutoff: the minimum absolute displayed coefficient is approximately 0.5847. The unthresholded precursor includes 139 selected genes by the same 64 lipids, or 8,896 pairs; it is a complementary universe preserving the original MOFA selection. Neither uses new SD correlations to choose features.",
    "", paste0("All 43 original FCG mouse pairings are retained: 21 ovaries and 22 testes, with XX/XY and group metadata preserved. Scaled whole-hepatocyte RNA and log2 original lipid intensities replay the displayed and precursor archives to maximum absolute error ",
                format(max(replay$log2_max_abs_archive_difference), scientific = TRUE, digits = 3),
                ". New SD uses 14 distinct CF453 animals, seven recorded females and seven males. There is one whole-hepatocyte RNA profile and one whole-liver lipid profile per mouse. No zone-specific RNA profiles enter this analysis."),
    "", "Matching is one-to-one through the strict final molecular annotations; missing or invariant features remain recorded in the coverage ledger. The same eligible, originally selected pairs enter every observed comparison and mouse draw. A new correlation is kept even if it is weak or reverses sign. [Reference replay](pooled_reference_replay.csv), [coverage](pooled_coverage.csv), [complete matched coefficients](matched_pooled_edges.csv) and [mouse metadata](metadata.csv) retain those identities.",
    "", coverage_note,
    "", "## Pooled coefficient conservation", "", table_stats(rows), "", null_table(rows),
    "", "Positive vector r / positive ULM indicate positive alignment between the original and new correlation patterns. Native ULM is the intercept-model slope t statistic. Its nominal p-value treats edge residuals as independent even though pairs share genes, lipids and mice. Edge count therefore does not represent the biological sample size. These are exploratory unadjusted summaries, not independent confirmation of every edge. A returned numerical p-value of zero is reported as <2.23e-308, indicating floating-point underflow rather than an exact zero probability.",
    "", "The unrestricted null permutes complete lipid mouse vectors relative to complete RNA vectors across all 14 SD animals. The independent-marginal null bootstrap resamples RNA mice and lipid mice separately across that pooled roster. Both test broad absence of pairing information; either may reject because females and males differ similarly in both modalities. They do not isolate regulation or correspondence within sex.",
    "", "The conditional null makes those rearrangements or independent resamples within the two seven-mouse sex groups. Within-sex permutations retain the observed sex means and empirical margins exactly. The independent within-sex bootstrap preserves sex membership and resamples from each sex's observed empirical distribution, while individual draw means vary. Both conditional reference-alignment distributions can consequently be centered above zero. The reported conditional p-value is the prespecified one-sided upper tail: does the positive observed alignment exceed what the shared sex structure can already generate? Lower and equal-tail diagnostics are saved without being relabeled as positive-conservation evidence.",
    "", sex_mean_note,
    "", paste0("The main run uses ", analysis$metadata$random_permutations,
                " random pairing draws per null and ", analysis$metadata$null_bootstrap_draws,
                " null-bootstrap draws. Each draw is shared across every gene and lipid. Identity pairing is saved separately as draw zero; Monte Carlo p-values use `(1 + exceedances)/(1 + valid draws)`. All tails and valid/invalid counts are in [pooled inference](pooled_inference.csv). Undefined features invalidate their fixed summary vector instead of shrinking its edge set within a draw. The bootstrap's type-I error and power at seven animals per sex have not been calibrated."),
    "", "Paired stability bootstraps resample matched SD animals within sex. Additional stability draws independently resample original FCG paired mice within their chromosome/gonad groups. They describe sensitivity to mouse sampling and reference uncertainty, while keeping historical feature selection fixed. They do not supply a null p-value and do not propagate uncertainty in the original feature-selection process.",
    "", "[Coefficient figures](plots/pooled_conservation.pdf) separate displayed and precursor universes. [Null figures](plots/pooled_null_checks.pdf) show why unrestricted and sex-conditioned checks can differ.",
    "", "## Removing the between-sex means", "", table_stats(centered),
    "", "This sensitivity subtracts each cohort's female or male mean from every RNA and lipid feature before recalculating Pearson correlations on unchanged pairs. It describes a different residual network, rather than replaying Figure 3. Reported residual ULM p-values remain nominal under dependent edges; no residual-network null p-value is assigned. Pooled conservation can coexist with weak residual conservation because between-sex effects and within-sex mouse variation answer different questions.",
    "", "## Four predefined Cyp genes", "", focus_table,
    "", "The Cyp2b9, Cyp2b13, Cyp2c38 and Cyp2d9 summaries were specified before examining new results and retain their original displayed species membership. Full precursor summaries for those genes are also saved. These secondary summaries are exploratory and their p-values are unadjusted across genes and universes. A gene's pooled coefficient pattern does not identify enzyme activity, substrate specificity or causal regulation.",
    "", "[Cyp coefficient figures](plots/Cyp_pooled_conservation.pdf) show all available originally displayed pairs. [Individual-mouse examples](plots/selected_pooled_pairs.pdf) show one pair per predefined gene: the matched originally displayed pair with the largest absolute original FCG coefficient, with original lipid-name ordering resolving ties. The new correlation does not enter selection. Seven females and seven males are colored separately, with all mouse IDs shown; pooled and within-sex fits help expose mean-driven associations. The [figure selection ledger](plots/figure_pair_selection_ledger.csv) records these choices.",
    "", "## Molecular mapping and analytical-quality sensitivities", "", sensitivity_text,
    "", "The 854-pair conditional molecular compatibility map is a sensitivity, rather than a revised identity claim. Its original and new coefficients correspond to the explicitly recorded, potentially different molecular annotations, so any agreement is interpreted at the supported coarse annotation resolution. The QC-CV <=30% and supplied Fill % >=0.5 restrictions test dependence on the exported quality fields. They alter the denominator and do not constitute independent biological replication. [Mapping sensitivity](mapping_sensitivity_inference.csv) and [QC sensitivities](pooled_QC_sensitivity.csv) preserve complete results and null definitions.",
    "", "Supplied positive lipid values are log2 transformed without additional normalization or imputation. Positivity does not establish that upstream missingness or imputation was absent; the supplied Fill % field records incomplete coverage for some rows. Supplied RNA scale.data is used without changing its upstream normalization. Row scaling the stored data layer replays scale.data, which validates the scaling relationship rather than independently validating normalization. Nucleus contributions, detailed upstream nuclear QC, batch and litter covariates remain unavailable.",
    "", "## Interpretation and reproduction", "",
    "This pooled test answers whether the sex-pooled, non-zonated Figure 3 pattern transfers to a new standard-diet cohort. It is complementary to WP2's within-sex zonation analysis and should not be described as rescuing a failed zonal result: the feature selection, biological profiles and between-sex contributions differ. Differences in age, background, diet duration and processing also limit direct cohort equivalence; meeting descriptions of those differences are not fully recorded in the supplied metadata.",
    "", "A defensible manuscript distinction is overall pooled pattern reproduction versus additional correspondence beyond sex margins. Neither establishes direct lipid regulation or a mechanistic pathway. Non-significant conditional tests do not prove lack of a biological relationship; they indicate that this cohort and this conditional test do not establish the specified extra correspondence. Exact-annotation gaps likewise do not prove that the new experiment failed to detect those lipids.",
    "", "Recorded input limitations:", "", paste0("- ", limitations),
    "", "Run from the repository root, using a fresh output directory:", "", "```sh",
    "Rscript scripts/workflows/05c_wdsd_lipid_overlap_and_pooled_network.R --validate-only",
    paste0("Rscript scripts/workflows/05c_wdsd_lipid_overlap_and_pooled_network.R --permutations=",
           analysis$metadata$random_permutations, " --bootstraps=", analysis$metadata$null_bootstrap_draws,
           " --paired-bootstraps=", analysis$metadata$paired_stability_draws,
           " --seed=", analysis$metadata$seed, " --output=results/wdsd_network_comparison_followup_reproduction"),
    "```", "", "The workflow retains all shared resampling indices and null draws, original feature membership, molecular candidate and exclusion ledgers, input/code SHA256, settings and package/session versions. [Plot manifest](plot_manifest.csv) lists the figure pages. The original analyses and preceding WP2 report are preserved.")
  report_links <- regmatches(report, gregexpr("\\]\\(([^)]+)\\)", report, perl = TRUE))
  targets <- unlist(lapply(report_links, function(z) sub("^\\]\\((.*)\\)$", "\\1", z)), use.names = FALSE)
  local <- targets[!grepl("^[a-z]+://", targets)]
  if (any(!file.exists(file.path(out, local)))) {
    stop("Pooled report references missing artifacts: ", paste(local[!file.exists(file.path(out, local))], collapse = ", "))
  }
  writeLines(report, file.path(out, "README.md"), useBytes = TRUE)
  invisible(file.path(out, "README.md"))
}
