# WP3 scientific figures and output-derived report. Source support functions first.
# Figure selection is fixed by the original FCG network, never CF453/CF459 results.
wp3_label <- function(x) {
  labels <- c(male_periportal = "Male periportal", male_pericentral = "Male pericentral")
  result <- unname(labels[x]); result[is.na(result)] <- x[is.na(result)]; result
}
wp3_number <- function(x, digits = 3L) vapply(x, function(v) {
  if (!is.finite(v)) "NA" else formatC(v, format = "f", digits = digits)
}, character(1L))
wp3_p <- function(x) vapply(x, function(v) {
  if (!is.finite(v)) return("NA")
  if (v == 0) return("0 (numerical underflow)")
  if (v < .001) formatC(v, format = "e", digits = 2L) else formatC(v, format = "f", digits = 4L)
}, character(1L))
wp3_inference_row <- function(analysis, context, universe = "full", gene = NULL) {
  x <- analysis$inference
  keep <- x$context == context & x$universe == universe &
    if (is.null(gene)) x$scope == "network" else x$scope == "gene" & !is.na(x$gene) & x$gene == gene
  row <- x[which(keep), , drop = FALSE]
  if (nrow(row) > 1L) stop("Duplicate WP3 inference row.")
  row
}
wp3_scatter <- function(analysis, context, universe = "full", gene = NULL) {
  e <- analysis$edges
  e <- e[e$context == context & e$eligible & (universe == "full" | e$primary_edge), , drop = FALSE]
  if (!is.null(gene)) e <- e[e$gene == gene, , drop = FALSE]
  title <- paste(wp3_label(context), if (is.null(gene)) "" else paste0(" | ", gene))
  if (!nrow(e)) {graphics::plot.new(); graphics::title(main = title); return(invisible(NULL))}
  colour <- if (universe == "full") "#306A96" else "#AD682B"
  graphics::plot(e$r_CF453_WD, e$r_CF459_WD, xlim = c(-1, 1), ylim = c(-1, 1),
    xlab = "CF453 male WD gene-lipid Pearson r (7 mice)", ylab = "CF459 male WD gene-lipid Pearson r (4 mice)",
    main = title, pch = 16, cex = if (nrow(e) > 1000L) .25 else .65,
    col = grDevices::adjustcolor(colour, alpha.f = if (nrow(e) > 1000L) .15 else .65),
    cex.main = .97, cex.lab = .84, axes = FALSE, xaxs = "i", yaxs = "i")
  graphics::axis(1, seq(-1, 1, .5)); graphics::axis(2, seq(-1, 1, .5), las = 1)
  graphics::box(col = "#697780"); graphics::abline(h = 0, v = 0, col = "#D7DFE4", lty = 3)
  graphics::abline(0, 1, col = "#81909A", lty = 2)
  if (stats::sd(e$r_CF453_WD) > 0) graphics::abline(stats::lm(e$r_CF459_WD ~ e$r_CF453_WD), col = colour, lwd = 1.7)
  row <- wp3_inference_row(analysis, context, universe, gene)
  if (nrow(row)) {
    annotation <- paste0("Coefficient-vector r = ", wp3_number(row$coefficient_vector_r),
      "\nULM = ", wp3_number(row$ULM_score), "\n", nrow(e), " fixed pairs")
    graphics::rect(-.98, .96 - graphics::strheight(annotation, cex = .79) - .04,
      -.95 + graphics::strwidth(annotation, cex = .79) + .03, .99,
      col = grDevices::adjustcolor("white", alpha.f = .94), border = NA)
    graphics::text(-.95, .96, annotation, adj = c(0, 1), cex = .79, col = "#172F41")
  }
  graphics::mtext(if (universe == "full") "Full fixed universe (primary)" else "Original displayed subset (secondary)",
                  side = 3, line = .25, cex = .73)
}
wp3_hist <- function(draws, observed, title, subtitle, colour = "#89AFCC") {
  valid <- is.finite(draws)
  if (!any(valid)) {graphics::plot.new(); graphics::title(main = title); graphics::text(.5, .5, "No defined draws"); return(invisible(NULL))}
  limits <- range(c(draws[valid], observed, 0)); padding <- max(diff(limits) * .10, .015)
  limits <- limits + c(-padding, padding)
  graphics::hist(draws[valid], breaks = pretty(limits, n = 18L),
    xlim = limits, col = colour, border = "white", main = title,
    xlab = "Coefficient-vector r", ylab = "Draws", cex.main = .9, cex.lab = .87)
  graphics::abline(v = observed, col = "#A53738", lwd = 2)
  graphics::mtext(subtitle, side = 3, line = .22, cex = .67)
  graphics::legend("topright", c(paste0("Observed ", wp3_number(observed)),
    paste0(sum(valid), " defined / ", length(draws)), paste0(sum(!valid), " undefined")),
    bty = "n", cex = .64, text.col = c("#A53738", "#172F41", "#172F41"))
}
wp3_pdf <- function(path, width, height, code) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE,
    title = paste("WP3", basename(path)), onefile = TRUE)
  on.exit(grDevices::dev.off())
  force(code)
}
# The original-cohort extension changes the comparison anchor, not edge selection.
wp3_fcg_network_rows <- function(fcg_comparison) {
  need <- c("edges", "inference", "paired_draws", "omission_draws", "omission_summary", "metadata")
  if (!is.list(fcg_comparison) || !all(need %in% names(fcg_comparison))) stop("Incomplete WP3 original-cohort comparison for presentation.")
  x <- fcg_comparison$inference
  x <- x[x$scope == "network", , drop = FALSE]
  keys <- paste(x$context, x$universe, sep = "__")
  expected <- c("male_periportal__full", "male_periportal__displayed", "male_pericentral__full", "male_pericentral__displayed")
  if (nrow(x) != 4L || anyDuplicated(keys) || !setequal(keys, expected)) stop("Expected two zones and identical full/displayed original-cohort comparisons.")
  x[match(expected, keys), , drop = FALSE]
}
wp3_age_report_note <- function(fcg_comparison) {
  supplied <- fcg_comparison$metadata$age_report_note
  if (is.character(supplied) && length(supplied) == 1L && nzchar(supplied)) return(supplied)
  "The user clarified the ages on 14 October 2026: **CF459 is 3-month WD, the original FCG cohort is 3-month SD, and CF453 is 10-month WD**. The requested comparison therefore contrasts the age-matched original SD cohort with an older WD cohort. These are ages at sampling; they do not establish diet duration, strain, litter or processing batch."
}
wp3_make_fcg_plots <- function(out, fcg_comparison) {
  rows <- wp3_fcg_network_rows(fcg_comparison)
  plot_dir <- file.path(out, "plots"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c("FCG_age_diet_comparison.pdf", "FCG_closeness_stability.pdf")
  age_note <- fcg_comparison$metadata$age_figure_note
  if (is.null(age_note)) age_note <- "Ages clarified by the user, 14 October 2026: CF459 3-month WD; original FCG 3-month SD; CF453 10-month WD."
  wp3_pdf(file.path(plot_dir, files[1L]), 11.0, 9.6, {
    graphics::par(mfrow = c(2, 2), mar = c(5.3, 5.2, 4.0, 1.2), oma = c(3.0, 0, 1.7, 0))
    for (context in c("male_periportal", "male_pericentral")) for (anchor in c("FCG", "CF453")) {
      e <- fcg_comparison$edges[fcg_comparison$edges$context == context & fcg_comparison$edges$eligible, , drop = FALSE]
      row <- rows[rows$context == context & rows$universe == "full", , drop = FALSE]
      xx <- if (anchor == "FCG") e$r_FCG_SD else e$r_CF453_WD
      yy <- e$r_CF459_WD
      r <- if (anchor == "FCG") row$r_FCG_CF459 else row$r_CF453_CF459
      colour <- if (anchor == "FCG") "#557A5C" else "#306A96"
      axis_label <- if (anchor == "FCG") "Original FCG male SD (3 months)\nGene-lipid Pearson r; 22 mice (11 XXT + 11 XYT)" else "CF453 male WD (10 months)\nGene-lipid Pearson r; 7 mice"
      graphics::plot(xx, yy, xlim = c(-1, 1), ylim = c(-1, 1),
        xlab = axis_label, ylab = "CF459 male WD (3 months; 4 mice)\nGene-lipid Pearson r",
        main = paste(wp3_label(context), "|", anchor, "anchor"), pch = 16, cex = .25,
        col = grDevices::adjustcolor(colour, alpha.f = .15), cex.main = .96, cex.lab = .80,
        axes = FALSE, xaxs = "i", yaxs = "i")
      graphics::axis(1, seq(-1, 1, .5)); graphics::axis(2, seq(-1, 1, .5), las = 1)
      graphics::box(col = "#697780"); graphics::abline(h = 0, v = 0, col = "#D7DFE4", lty = 3)
      graphics::abline(0, 1, col = "#81909A", lty = 2)
      graphics::abline(stats::lm(yy ~ xx), col = colour, lwd = 1.7)
      annotation <- paste0("Coefficient-vector r = ", wp3_number(r), "\n", nrow(e), " identical fixed pairs")
      graphics::rect(-.98, .96 - graphics::strheight(annotation, cex = .78) - .04,
        -.95 + graphics::strwidth(annotation, cex = .78) + .03, .99,
        col = grDevices::adjustcolor("white", alpha.f = .95), border = NA)
      graphics::text(-.95, .96, annotation, adj = c(0, 1), cex = .78, col = "#172F41")
      graphics::mtext("Full frozen universe (primary); no new correlation threshold", side = 3, line = .23, cex = .65)
    }
    graphics::mtext("WP3 | Is CF459 closer to original FCG SD or CF453 WD?", side = 3, outer = TRUE, line = .3, font = 2, cex = .97)
    graphics::mtext("Matched genes and 610 strict lipid identities in all panels; edge counts are not biological sample sizes.",
      side = 1, outer = TRUE, line = .45, cex = .72)
    graphics::mtext(age_note, side = 1, outer = TRUE, line = 1.65, cex = .67)
  })
  wp3_pdf(file.path(plot_dir, files[2L]), 11.0, 9.2, {
    graphics::par(mfrow = c(2, 2), mar = c(4.8, 4.5, 4.1, 1.0), oma = c(3.6, 0, 1.7, 0))
    for (context in c("male_periportal", "male_pericentral")) {
      row <- rows[rows$context == context & rows$universe == "full", , drop = FALSE]
      x <- fcg_comparison$paired_draws[fcg_comparison$paired_draws$specification == row$specification, , drop = FALSE]
      draws <- x$delta_r; valid <- is.finite(draws)
      limits <- range(c(draws[valid], row$delta_r, 0)); padding <- max(diff(limits) * .10, .015)
      limits <- limits + c(-padding, padding)
      graphics::hist(draws[valid], breaks = pretty(limits, n = 18L), xlim = limits, col = "#A8C6B2", border = "white",
        main = paste(wp3_label(context), "| three-cohort bootstrap"),
        xlab = "Delta r (original FCG anchor minus CF453 anchor)", ylab = "Defined draws", cex.main = .86, cex.lab = .78)
      graphics::abline(v = 0, col = "#697780", lty = 3)
      graphics::abline(v = row$delta_r, col = "#A53738", lwd = 2)
      graphics::mtext("Same paired CF459 draw for both anchor comparisons", side = 3, line = .22, cex = .65)
      graphics::legend("topright", c(paste0("Observed delta ", wp3_number(row$delta_r)),
        paste0(sum(valid), " defined / ", length(draws)), paste0(sum(!valid), " undefined"),
        paste0(wp3_number(100 * mean(draws[valid] > 0), 1), "% positive among defined")),
        bty = "n", cex = .63, text.col = c("#A53738", rep("#172F41", 3L)))
      z <- fcg_comparison$omission_draws[fcg_comparison$omission_draws$specification == row$specification, , drop = FALSE]
      colour <- c(joint_omission = "#306A96", CF453_only = "#AD682B", CF459_only = "#557A5C", FCG_only = "#8C6396")
      limits <- range(c(z$delta_r, row$delta_r, 0), na.rm = TRUE) + c(-.015, .015)
      graphics::plot(seq_len(nrow(z)), z$delta_r, ylim = limits, pch = 16, col = colour[z$method],
        main = paste(wp3_label(context), "| omissions"), xlab = "Prespecified omission configuration",
        ylab = "Delta r", cex.main = .89, cex.lab = .82)
      graphics::abline(h = row$delta_r, col = "#A53738", lwd = 1.7, lty = 2)
      graphics::abline(h = 0, col = "#CDD8DE", lty = 3)
      graphics::legend("bottomleft", c("28 joint WD omissions", "7 CF453 only", "4 CF459 only", "22 original FCG only"),
        col = colour, pch = 16, bty = "n", cex = .63)
      graphics::mtext(paste0(sum(is.finite(z$delta_r)), "/", nrow(z), " omission deltas defined; observed delta dashed"), side = 3, line = .22, cex = .64)
    }
    graphics::mtext("Closeness stability | Positive delta means closer to original FCG SD", side = 3, outer = TRUE, line = .3, font = 2, cex = .97)
    graphics::mtext("Paired whole mice: FCG stratified 11 XXT / 11 XYT; independent CF453 and FCG draws; shared CF459 draw and zone indices.",
      side = 1, outer = TRUE, line = .45, cex = .69)
    graphics::mtext("Positive-delta fraction is descriptive conditional on defined draws, not a null P or a probability of an age effect.",
      side = 1, outer = TRUE, line = 1.65, cex = .69)
    graphics::mtext("Any undefined frozen edge invalidates the complete comparison; no draw-specific feature deletion.",
      side = 1, outer = TRUE, line = 2.80, cex = .69)
  })
  invisible(file.path("plots", files))
}

# Refresh the intake figure from saved intake ledgers without statistical refitting.
wp3_refresh_cohort_QC <- function(out, cf453 = NULL, cf459 = NULL, setup = NULL, analysis = NULL) {
  if (is.null(cf453)) cf453 <- list(rna_profile_qc = utils::read.csv(file.path(out, "CF453_rna_profile_qc.csv"), check.names = FALSE))
  if (is.null(cf459)) cf459 <- list(rna_profile_qc = utils::read.csv(file.path(out, "CF459_rna_profile_qc.csv"), check.names = FALSE),
    zonation_marker_diagnostic = utils::read.csv(file.path(out, "CF459_zonation_marker_diagnostic.csv"), check.names = FALSE))
  if (is.null(setup)) setup <- list(coverage = utils::read.csv(file.path(out, "coverage.csv"), check.names = FALSE))
  plot_dir <- file.path(out, "plots"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  wp3_pdf(file.path(plot_dir, "QC_and_coverage.pdf"), 10.8, 8.9, {
    graphics::par(mfrow = c(2, 2), mar = c(5.5, 4.8, 3.8, 1.2), oma = c(3.0, 0, 1.5, 0))
    c <- setup$coverage
    counts <- rbind(Original = c$original_full_edges, Shared = c$eligible_full_edges)
    graphics::barplot(counts, beside = TRUE, names.arg = c("Periportal", "Pericentral"), col = c("#BECBD3", "#306A96"),
      ylab = "Fixed full gene-lipid pairs", main = "Strict common feature coverage", cex.main = .95, cex.names = .83)
    graphics::legend("topleft", c("Original FCG full", "Shared eligible"), fill = c("#BECBD3", "#306A96"), bty = "n", cex = .70)
    graphics::mtext("Full universe primary; fixed displayed subsets secondary", side = 1, line = 3.3, cex = .7)
    profiles <- rbind(transform(cf453$rna_profile_qc, cohort = "CF453"),
      transform(cf459$rna_profile_qc[, names(cf453$rna_profile_qc), drop = FALSE], cohort = "CF459"))
    profiles <- profiles[profiles$object_kind == "zonation" & profiles$sex == "male" & profiles$diet == "WD", ]
    groups <- paste(profiles$cohort, profiles$zonation_type, sep = "\n")
    graphics::boxplot(log10(profiles$supplied_count_library_size) ~ groups,
      col = c("#8DB2CC", "#8DB2CC", "#D5AD84", "#D5AD84"),
      ylab = "log10 supplied aggregate counts", xlab = "", main = "Primary male WD RNA profiles", cex.main = .93, cex.axis = .72)
    graphics::mtext("Count totals are not nucleus contributions", side = 1, line = 4.2, cex = .69)
    markers <- cf459$zonation_marker_diagnostic
    markers <- markers[markers$sex == "male", ]
    colour <- ifelse(markers$zonation_type == "Periportal", "#557A5C", "#8C6396")
    graphics::plot(markers$Cyp2f2_data, markers$Cyp2e1_data, pch = 16, col = colour,
      xlab = "CF459 Cyp2f2 stored data", ylab = "CF459 Cyp2e1 stored data", main = "Aggregate zonation marker diagnostic",
      cex.main = .9, cex.lab = .80)
    graphics::abline(0, 1, col = "#9AAAB4", lty = 2)
    graphics::legend("topright", c("Periportal", "Pericentral"), pch = 16, col = c("#557A5C", "#8C6396"), bty = "n", cex = .70)
    graphics::mtext("Aggregate means cannot validate nucleus-level gates", side = 1, line = 3.5, cex = .66)
    graphics::plot.new(); graphics::title(main = "Design and unavailable information", cex.main = .95)
    notes <- c("CF453: 7 male WD mice", "CF459: 4 male WD mice", "Each mouse: 2 RNA zones, 1 bulk lipid profile", "Primary RNA: stored scale.data, unchanged", "CF459: no analytical QC or Fill export", "Nucleus counts and processing provenance unknown", "Export lacks age/background/diet-duration fields", "CF459 ctrl labels mean WD, not SD")
    graphics::text(.02, .94, paste(notes, collapse = "\n\n"), adj = c(0, 1), cex = .80, col = "#172F41")
    graphics::mtext("WP3 | Intake and fixed-universe coverage", side = 3, outer = TRUE, line = .2, font = 2)
    graphics::mtext("No cohort or lipid was selected on observed concordance. Strict matching preserves chain order, ether prefixes and oxygenation.",
      side = 1, outer = TRUE, line = .5, cex = .72)
  })
  invisible(file.path("plots", "QC_and_coverage.pdf"))
}

wp3_make_plots <- function(out, cf453, cf459, setup, analysis, sensitivities = NULL, fcg_comparison = NULL) {
  if (!all(c("edges", "inference", "metadata", "loco_draws", "loco_summary") %in% names(analysis))) stop("Incomplete WP3 analysis for presentation.")
  contexts <- setup$contexts
  names_contexts <- c("male_periportal", "male_pericentral")
  if (!identical(names(contexts), names_contexts)) stop("Expected ordered male zone contexts.")
  plot_dir <- file.path(out, "plots"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  wp3_pdf(file.path(plot_dir, "network_conservation.pdf"), 10.8, 9.3, {
    graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.5, 4.0, 1.2), oma = c(2.5, 0, 1.5, 0))
    for (context in names_contexts) for (universe in c("full", "displayed")) wp3_scatter(analysis, context, universe)
    graphics::mtext("WP3 | Independent male Western-diet cohorts", side = 3, outer = TRUE, line = .2, font = 2)
    graphics::mtext("Dotted grey lines: zero; dashed: identity. Lines fit coefficients, not independent biological replicates.",
      side = 1, outer = TRUE, line = .5, cex = .75)
  })
  wp3_pdf(file.path(plot_dir, "null_checks.pdf"), 11.6, 8.5, {
    graphics::par(mfrow = c(2, 3), mar = c(4.2, 3.9, 4.1, .9), oma = c(3.0, 0, 1.5, 0))
    for (context in names_contexts) {
      row <- wp3_inference_row(analysis, context); key <- row$specification
      wp3_hist(analysis$permutation_CF459_draws[[key]], row$coefficient_vector_r,
        paste(wp3_label(context), "| 4! exact"), "CF459 pairing; fixed observed CF453")
      wp3_hist(analysis$permutation_CF453_draws[[key]], row$coefficient_vector_r,
        paste(wp3_label(context), "| 7! exact"), "CF453 pairing; fixed observed CF459")
      wp3_hist(analysis$bootstrap_null_draws[[key]], row$coefficient_vector_r,
        paste(wp3_label(context), "| null bootstrap"), "Independent CF459 margins; fixed CF453", "#D7B28C")
    }
    graphics::mtext("Whole-mouse pairing nulls | Full frozen universe", side = 3, outer = TRUE, line = .2, font = 2)
    graphics::mtext("Two-sided tails use |r|. The 24-configuration test is primary; reciprocal 5,040 configurations are complementary.",
      side = 1, outer = TRUE, line = .45, cex = .72)
    graphics::mtext("Approximate null-bootstrap P is conditional on defined frozen vectors; undefined draws are retained in audit counts.",
      side = 1, outer = TRUE, line = 1.65, cex = .71)
  })
  wp3_pdf(file.path(plot_dir, "stability_and_omissions.pdf"), 10.8, 8.9, {
    graphics::par(mfrow = c(2, 2), mar = c(4.6, 4.4, 4.0, 1), oma = c(2.8, 0, 1.5, 0))
    for (context in names_contexts) {
      row <- wp3_inference_row(analysis, context); key <- row$specification
      wp3_hist(analysis$paired_stability_draws[[key]], row$coefficient_vector_r,
        paste(wp3_label(context), "| paired bootstrap"), "Paired mice resampled independently per cohort", "#A8C6B2")
      x <- analysis$loco_draws[analysis$loco_draws$specification == key, , drop = FALSE]
      groups <- c("joint_omission", "CF453_only", "CF459_only")
      colour <- c(joint_omission = "#306A96", CF453_only = "#AD682B", CF459_only = "#557A5C")
      limits <- range(c(x$coefficient_vector_r, row$coefficient_vector_r, 0), na.rm = TRUE) + c(-.015, .015)
      graphics::plot(seq_len(nrow(x)), x$coefficient_vector_r, ylim = limits,
        pch = 16, col = colour[x$method], main = paste(wp3_label(context), "| omissions"),
        xlab = "Prespecified omission configuration", ylab = "Coefficient-vector r", cex.main = .92, cex.lab = .84)
      graphics::abline(h = row$coefficient_vector_r, col = "#A53738", lwd = 1.7, lty = 2)
      graphics::abline(h = 0, col = "#CDD8DE", lty = 3)
      graphics::legend("bottomleft", c("28 joint (6 + 3 mice)", "7 CF453-only (6 + 4)", "4 CF459-only (7 + 3)"),
        col = colour[groups], pch = 16, bty = "n", cex = .65)
      graphics::mtext(paste0("Observed r dashed; ", sum(is.finite(x$coefficient_vector_r)), "/", nrow(x), " omission vectors defined"), side = 3, line = .2, cex = .65)
    }
    graphics::mtext("Stability diagnostics | Whole mice, shared indices across zones", side = 3, outer = TRUE, line = .2, font = 2)
    graphics::mtext("Stability is descriptive. Resampling four mice often gives low diversity; an undefined edge invalidates the entire fixed vector.",
      side = 1, outer = TRUE, line = .55, cex = .73)
  })
  wp3_pdf(file.path(plot_dir, "Mtmr11_coefficients.pdf"), 10.8, 5.4, {
    graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.6, 4.1, 1.0), oma = c(2.1, 0, 1.3, 0))
    for (universe in c("full", "displayed")) wp3_scatter(analysis, "male_periportal", universe, "Mtmr11")
    graphics::mtext("Mtmr11 | Prespecified periportal example", side = 3, outer = TRUE, font = 2, line = .1)
    graphics::mtext("Species membership follows the original FCG network. This gene was not selected from CF459 results.",
      side = 1, outer = TRUE, cex = .75, line = .5)
  })
  e <- analysis$edges
  selected <- e[e$context == "male_periportal" & e$gene == "Mtmr11" & e$primary_edge & e$eligible, , drop = FALSE]
  selected <- selected[order(-abs(selected$r_reference), selected$reference_lipid), , drop = FALSE]
  selected <- utils::head(selected, 3L)
  selected$selection_rule <- "top three absolute frozen FCG displayed Mtmr11 periportal coefficients; no new-cohort selection"
  utils::write.csv(selected, file.path(out, "predefined_pairs.csv"), row.names = FALSE)
  wp3_pdf(file.path(plot_dir, "Mtmr11_selected_pairs.pdf"), 8.3, 11.5, {
    graphics::par(mfrow = c(max(1L, nrow(selected)), 2), mar = c(4.4, 4.8, 5.0, 1.1), oma = c(3.0, 0, 1.8, 0))
    for (i in seq_len(nrow(selected))) for (cohort in c("CF453", "CF459")) {
      edge <- selected[i, ]; x <- contexts$male_periportal
      xx <- x[[paste0("rna_", cohort)]][edge$gene, ]; yy <- x[[paste0("lipid_", cohort)]][edge$feature_id, ]
      expand <- function(z) range(z) + c(-.15, .15) * max(diff(range(z)), .2)
      colour <- if (cohort == "CF453") "#306A96" else "#AD682B"
      title <- paste0(cohort, " WD | ", edge$gene, "\n", paste(strwrap(edge$reference_lipid, width = 36L), collapse = "\n"))
      graphics::plot(xx, yy, xlim = expand(xx), ylim = expand(yy), pch = if (cohort == "CF453") 16 else 17,
        col = colour, cex = .95, xlab = "Stored scaled RNA", ylab = "Supplied lipid intensity (log2)",
        main = title, cex.main = .83, cex.lab = .78)
      graphics::abline(stats::lm(yy ~ xx), col = colour, lwd = 1.5)
      graphics::text(xx, yy, labels = sub("^cf[0-9]+-", "", colnames(x[[paste0("rna_", cohort)]])),
        pos = 3, offset = .3, cex = .62, col = colour)
      graphics::mtext(paste0("n=", length(xx), " mice; within-cohort r=", wp3_number(stats::cor(xx, yy))),
        side = 3, line = .3, cex = .74)
    }
    graphics::mtext("Mtmr11 | Original-reference-selected mouse-level examples", side = 3, outer = TRUE, line = .3, font = 2, cex = .94)
    graphics::mtext("Point labels are mouse suffixes. Each panel uses its own axes; absolute RNA/lipid levels are not compared across exports.",
      side = 1, outer = TRUE, line = .5, cex = .64)
    graphics::mtext("All examples use the male periportal network and the same fixed molecular annotation in both cohorts.",
      side = 1, outer = TRUE, line = 1.7, cex = .68)
  })
  wp3_refresh_cohort_QC(out, cf453, cf459, setup, analysis)
  # Summarize a concrete source of bootstrap degeneracy without changing selection.
  low_diversity <- do.call(rbind, lapply(names_contexts, function(context) {
    rn <- contexts[[context]]$rna_CF459
    rows <- rownames(rn)[apply(rn, 1L, function(v) length(unique(v))) <= 3L]
    if (!length(rows)) return(NULL)
    do.call(rbind, lapply(rows, function(gene) {
      v <- rn[gene, ]; indices <- analysis$resampling_index_matrices$CF459_null_rna
      data.frame(context = context, gene = gene, unique_expression_values = length(unique(v)),
        largest_identical_group = max(table(v)), n_mice = length(v),
        constant_null_rna_draws = sum(apply(indices, 1L, function(index) length(unique(v[index])) == 1L)),
        requested_null_draws = nrow(indices), stringsAsFactors = FALSE)
    }))
  }))
  bundle <- list(low_diversity = low_diversity, setup = setup[c("coverage", "lipid_crosswalk", "WP2_replay")],
    analysis = analysis[c("inference", "metadata", "loco_summary", "resampling_draw_summary")],
    sensitivities = sensitivities, fcg_comparison = fcg_comparison, predefined_pairs = selected,
    cf453_qc = cf453$qc_summary, cf459_qc = cf459$qc_summary,
    cf453_layers = cf453$layer_inventory, cf459_layers = cf459$layer_inventory,
    cf453_metadata = cf453$metadata, cf459_metadata = cf459$metadata,
    cf453_limitations = cf453$limitations, cf459_limitations = cf459$limitations,
    plot_files = file.path("plots", c("network_conservation.pdf", "null_checks.pdf", "stability_and_omissions.pdf",
      "Mtmr11_coefficients.pdf", "Mtmr11_selected_pairs.pdf", "QC_and_coverage.pdf")))
  if (!is.null(fcg_comparison)) bundle$plot_files <- c(wp3_make_fcg_plots(out, fcg_comparison), bundle$plot_files)
  saveRDS(bundle, file.path(out, "report_bundle.rds"))
  invisible(bundle$plot_files)
}

wp3_write_fcg_section <- function(add, table_md, fcg_comparison) {
  rows <- wp3_fcg_network_rows(fcg_comparison)
  primary <- rows[rows$universe == "full", , drop = FALSE]
  label <- function(x) paste(wp3_label(x$context), x$universe)
  add("## 1. Requested comparison: CF459 versus original FCG and CF453", "",
    "The new question is whether the CF459 WD network is closer to the original paper's male SD network than to the CF453 WD network. **All three correlations use exactly the same frozen, strictly shared 610 lipid identities and original context-specific genes**: 15,250 periportal and 18,300 pericentral full pairs. Original displayed pairs remain secondary (69 and 202). No cohort gets a larger feature universe in this comparison.", "",
    wp3_age_report_note(fcg_comparison), "",
    "The original male cohort contains **22 testes-bearing FCG mice: 11 XXT and 11 XYT**. This is the cohort used to construct the original male Figure 5 networks; it is not the five additional wild-type males. CF453 has seven WD males and CF459 has four. The original FCG SD coefficients are recalculated from paired original RNA and log2-lipid matrices and replay the archived network before matching.", "",
    "Define **delta r = correlation(CF459, original FCG SD) - correlation(CF459, CF453 WD)**. Positive delta means CF459 is closer to the original FCG SD coefficient vector; negative delta means it is closer to CF453 WD. This compares coefficient-vector alignment, not equality of gene-lipid correlations or absolute molecular abundance.", "")
  table_md(data.frame(Comparison = label(rows), Pairs = rows$n_edges,
    `Original FCG / CF459 r` = wp3_number(rows$r_FCG_CF459),
    `CF453 / CF459 r` = wp3_number(rows$r_CF453_CF459),
    `Original FCG / CF453 r` = wp3_number(rows$r_FCG_CF453),
    `Delta r` = wp3_number(rows$delta_r), check.names = FALSE))
  direction <- vapply(seq_len(nrow(primary)), function(i) paste0(wp3_label(primary$context[i]), ": CF459 is ",
    if (primary$delta_r[i] > 0) "closer to original FCG SD" else if (primary$delta_r[i] < 0) "closer to CF453 WD" else "equally aligned to the two anchors",
    " in the observed full network (delta r = ", wp3_number(primary$delta_r[i]), ")."), character(1L))
  add(paste0("**", paste(direction, collapse = " "), "**"), "",
    "Absolute full-network alignment with original FCG remains weak despite being higher than with CF453. This is an observed, cohort-specific comparison. Original feature selection, different sample sizes and unmeasured cohort/processing differences can affect alignment. Even a positive delta does not identify an age effect or show that diet is unimportant. Age, diet, genetic background, diet duration and batch are not independently varied here; there is no within-CF459 SD arm.", "",
    "[Matched-anchor coefficient figure](plots/FCG_age_diet_comparison.pdf) shows both full-network anchors beside one another, with CF459 coefficients on the same y-axis and exactly the same pairs in each zone.", "",
    "### 1.1 ULM and whole-mouse evidence for the original-FCG anchor", "",
    "The native ULM below regresses CF459 coefficients on the original FCG SD coefficients with an intercept. Its nominal edge-level p-value is descriptive because edges share mice, genes and lipids. The exact null enumerates the same 24 CF459 RNA-lipid pairings as the original WP3 comparison, now conditional on the fixed original-FCG anchor. The approximate null bootstrap reuses the same complete CF459 RNA and lipid index banks; no smaller p-value is selected across anchors.", "")
  table_md(data.frame(Comparison = label(rows), `ULM score` = wp3_number(rows$ULM_FCG_CF459),
    `Nominal ULM p` = wp3_p(rows$ULM_p_FCG_CF459),
    `Exact CF459 p (24)` = wp3_p(rows$pairing_CF459_given_FCG_p),
    `Null bootstrap p (defined only)` = wp3_p(rows$null_bootstrap_p), check.names = FALSE))
  if (is.data.frame(fcg_comparison$exact_zero_null_draws)) {
    certified <- sum(fcg_comparison$exact_zero_null_draws$exact_zero_source_contingency)
    add(paste0("A numerical-degeneracy audit identifies **", certified, " of ", fcg_comparison$metadata$B,
      " CF459 null index configurations** whose exact integer source-identity contingency makes every feature covariance mathematically zero. Tiny floating-point residuals cannot supply a defined coefficient-vector correlation. The new original-FCG-anchor calculation therefore records these vectors as undefined, without using an arbitrary coefficient cutoff. Its conditional denominator can differ from the original WP3 tables retained below. The [degeneracy ledger](fcg_comparison/null_degeneracy_ledger.csv) records every affected specification and historical/raw value; observed correlations and delta are unaffected by this null-only correction."), "")
  }
  audit <- fcg_comparison$historical_null_correction_summary
  if (is.data.frame(audit) && nrow(audit)) {
    audit <- audit[audit$scope == "network", , drop = FALSE]
    add("For transparency, the same mathematical-zero correction also gives the following audit of the original CF453-anchor null bootstrap. Source files and the historical figures are preserved; the corrected values below supersede the historical bootstrap p-values for interpretation. This correction does not alter the observed WD-versus-WD correlations or exact pairing tests.", "")
    table_md(data.frame(Comparison = label(audit),
      `Historical bootstrap p` = wp3_p(audit$original_null_bootstrap_p),
      `Corrected bootstrap p` = wp3_p(audit$corrected_null_bootstrap_p),
      `Original defined draws` = audit$original_valid_draws,
      `Corrected defined draws` = audit$corrected_valid_draws, check.names = FALSE))
    add("Complete corrected draws, completion bounds and predefined-gene checks are saved in the [historical-null correction audit](fcg_comparison/historical_null_correction_summary.csv). The next table gives the corrected original-FCG-anchor denominator and the separate paired-delta stability counts.", "")
  }
  table_md(data.frame(Comparison = label(rows),
    `Null defined / requested` = paste(rows$bootstrap_valid_draws, "/", rows$null_bootstrap_draws),
    `Null completion range` = paste(wp3_p(rows$bootstrap_completion_lower), wp3_p(rows$bootstrap_completion_upper), sep = " to "),
    `Delta defined / requested` = paste(rows$paired_delta_valid_draws, "/", rows$paired_draws),
    `Positive delta among defined` = paste0(wp3_number(100 * rows$paired_delta_positive_fraction, 1), "%"), check.names = FALSE))
  add(paste0("The closeness comparison uses **", fcg_comparison$metadata$paired_B,
    " three-cohort paired whole-mouse bootstrap draws**. Each draw resamples original FCG RNA/lipids together within the 11 XXT and 11 XYT strata; CF453 and FCG are independently sampled, and the **same CF459 paired draw** enters both anchor correlations. Both zones also share indices. The positive-delta fraction describes stability conditional on defined draws; **it is not a null p-value, posterior probability or probability of an age effect**."), "",
    if (!is.null(fcg_comparison$metadata$FCG_paired_seed)) paste0("The original WP3 whole-mouse index banks are reused; the additional independent FCG stratified paired bank uses recorded seed ", fcg_comparison$metadata$FCG_paired_seed, ".") else "",
    "The full frozen vector is retained in every draw. If a required gene or lipid is constant, or the complete coefficient vector is constant, the affected comparison and delta are undefined. No draw-specific edge deletion is used. Null completion ranges assign every undefined null draw to either tail extreme; these are mathematical missing-draw bounds, not confidence intervals. The n=4 CF459 limitation therefore remains visible in the new comparison.", "")
  om <- fcg_comparison$omission_summary
  om <- om[om$scope == "network" & om$universe == "full", , drop = FALSE]
  table_md(data.frame(Zone = wp3_label(om$context), Omissions = om$method,
    `Defined / planned` = paste(om$valid_draws, "/", om$draws),
    `Minimum delta r` = wp3_number(om$minimum_delta_r),
    `Maximum delta r` = wp3_number(om$maximum_delta_r),
    `Positive delta among defined` = paste0(wp3_number(100 * om$delta_positive_fraction, 1), "%"), check.names = FALSE))
  add("The omission audit preserves the original 28 joint WD omissions, seven CF453-only and four CF459-only omissions, and adds all 22 original-FCG single-mouse omissions. These 61 overlapping configurations are descriptive influence checks, not independent experiments. [Closeness stability figure](plots/FCG_closeness_stability.pdf) displays full-network deltas and undefined-draw counts. The [extension inference](fcg_comparison/inference.csv), [paired-draw ledger](fcg_comparison/paired_draws.csv) and [omission summary](fcg_comparison/omission_summary.csv) retain the complete full/displayed and predefined-gene results.", "")
  focus <- fcg_comparison$inference
  focus <- focus[focus$scope == "gene" & !is.na(focus$gene) & focus$gene == "Mtmr11", , drop = FALSE]
  if (nrow(focus)) {
    add("### 1.2 Prespecified Mtmr11 does not mirror the full-network pattern", "",
      "The full-network comparison does not imply that every gene is closer to the original FCG network. The predefined periportal Mtmr11 example has a negative observed delta in both full and displayed universes; its paired delta direction is unstable. This is retained alongside the global network result rather than selecting only concordant genes.", "")
    table_md(data.frame(Example = paste("Mtmr11", focus$universe), Pairs = focus$n_edges,
      `Original FCG / CF459 r` = wp3_number(focus$r_FCG_CF459), `CF453 / CF459 r` = wp3_number(focus$r_CF453_CF459),
      `Delta r` = wp3_number(focus$delta_r), `Positive delta among defined` = paste0(wp3_number(100 * focus$paired_delta_positive_fraction, 1), "%"), check.names = FALSE))
    table_md(data.frame(Comparison = paste("Mtmr11", focus$universe), `ULM score` = wp3_number(focus$ULM_FCG_CF459),
      `Nominal ULM p` = wp3_p(focus$ULM_p_FCG_CF459), `Exact CF459 p (24)` = wp3_p(focus$pairing_CF459_given_FCG_p),
      `Null bootstrap p (defined only)` = wp3_p(focus$null_bootstrap_p), check.names = FALSE))
  }
  invisible(NULL)
}

wp3_write_report <- function(out) {
  bundle <- readRDS(file.path(out, "report_bundle.rds"))
  fcg_path <- file.path(out, "fcg_comparison", "analysis.rds")
  fcg_comparison <- if (file.exists(fcg_path)) readRDS(fcg_path) else bundle$fcg_comparison
  if (!is.null(fcg_comparison)) {
    wp3_fcg_network_rows(fcg_comparison)
    # Figures depend only on saved numbers; report-only never recalculates fits.
    new_plots <- wp3_make_fcg_plots(out, fcg_comparison)
    bundle$plot_files <- unique(c(new_plots, bundle$plot_files))
  }
  analysis <- bundle$analysis; inference <- analysis$inference
  network <- inference[inference$scope == "network", , drop = FALSE]
  primary <- network[network$universe == "full", , drop = FALSE]
  if (nrow(primary) != 2L) stop("Report requires two primary male zone inference rows.")
  text <- character()
  add <- function(...) {text <<- c(text, ...); invisible(NULL)}
  table_md <- function(x) {
    if (!nrow(x)) return(invisible(NULL))
    values <- as.matrix(data.frame(lapply(x, function(z) {
      v <- as.character(z); v[is.na(v)] <- "NA"; gsub("|", "/", v, fixed = TRUE)
    }), check.names = FALSE))
    add(paste0("| ", paste(names(x), collapse = " | "), " |"), paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"))
    for (j in seq_len(nrow(values))) add(paste0("| ", paste(values[j, ], collapse = " | "), " |"))
    add("")
  }
  label <- function(x) paste(wp3_label(x$context), ifelse(x$universe == "full", "full", "displayed"))
  add(if (is.null(fcg_comparison)) "# WP3 - Male Western-diet cross-cohort network comparison" else "# WP3 - CF459 WD comparison with original FCG SD and CF453 WD", "",
      if (is.null(fcg_comparison)) "Report date: 14 October 2026. CF453 WDSD Western-diet males versus the independent CF459 WD-only cohort." else
        "Report date: 14 October 2026. Amended comparison of the original male FCG SD cohort, CF453 WD and independent CF459 WD.", "")
  if (!is.null(fcg_comparison)) wp3_write_fcg_section(add, table_md, fcg_comparison)
  section <- function(number, title) paste0("## ", number + as.integer(!is.null(fcg_comparison)), ". ", title)
  add(section(1L, "Original WP3 scope and observed WD network agreement"), "",
      "This secondary robustness analysis compares gene-lipid correlation coefficients in **seven CF453 male WD mice and four CF459 male WD mice**, separately in periportal and pericentral hepatocytes. The primary universe is every measurable common pair from the fixed original FCG full networks; the original displayed Figure 5 pairs are a secondary subset. Selection is frozen before examining either new cohort. No female network or within-CF459 diet contrast is included.", "",
      paste0("The observed full-network coefficient-vector correlations are **r = ", wp3_number(primary$coefficient_vector_r[primary$context == "male_periportal"]),
        "** in periportal and **r = ", wp3_number(primary$coefficient_vector_r[primary$context == "male_pericentral"]),
        "** in pericentral hepatocytes. These quantify alignment of the two edge-coefficient vectors, not correlation of independent biological replicates. Positive values denote similarly ordered coefficients; negative values denote opposing ordering."), "")
  if (all(abs(primary$coefficient_vector_r) < .10)) add(
    "**Neither primary network shows substantial positive observed alignment. This comparison does not demonstrate robust positive cross-cohort conservation in the measured male WD features.** The near-zero coefficients, four-mouse resolution and unknown cohort differences limit interpretation; they do not establish absence of shared regulation.", "")
  table_md(data.frame(Comparison = label(network), Pairs = network$n_edges,
    `Vector r` = wp3_number(network$coefficient_vector_r), `ULM score` = wp3_number(network$ULM_score),
    `Same coefficient sign` = paste0(wp3_number(100 * network$sign_agreement_fraction, 1), "%"), check.names = FALSE))
  add("The native ULM score is the intercept-OLS slope t-statistic for CF459 coefficients regressed on CF453 coefficients. Its sign follows coefficient-vector r. Thousands of overlapping gene-lipid edges share mice, genes and lipids; therefore its nominal p-value is descriptive and does not turn edge count into biological replication. No FDR, signed concordance score or confidence interval is reported.", "",
      "[Network coefficient figure](plots/network_conservation.pdf) shows full networks first and original displayed subsets second. Axes use CF453 WD coefficients on x and CF459 WD coefficients on y, each within [-1, 1]. The dashed identity line is a visual reference; coefficient equality is not tested.", "",
      section(2L, "Common-feature coverage and exact pairing"), "")
  coverage <- bundle$setup$coverage
  table_md(data.frame(Zone = wp3_label(coverage$context), `Original full pairs` = coverage$original_full_edges,
    `Shared full pairs` = coverage$eligible_full_edges, Genes = coverage$eligible_genes, Lipids = coverage$eligible_lipids,
    `CF453 / CF459 mice` = paste(coverage$n_CF453, coverage$n_CF459, sep = " / "), check.names = FALSE))
  table_md(data.frame(Zone = wp3_label(coverage$context), `Original displayed pairs` = coverage$original_displayed_edges,
    `Shared displayed pairs` = coverage$eligible_displayed_edges, Genes = coverage$eligible_displayed_genes,
    Lipids = coverage$eligible_displayed_lipids, check.names = FALSE))
  cw <- bundle$setup$lipid_crosswalk
  add(paste0("Strict annotation matching retains **", sum(cw$shared_primary), " shared eligible FCG lipids** across the two cohorts. Unique source feature IDs remain separate; matching uses the final pipe annotation and class formatting only. Chain order, slash/underscore, ether prefixes and oxygenation are preserved. Alignment IDs are never cross-cohort identifiers and no alias or coarser resolution is forced into the primary analysis."), "",
    "The full [lipid crosswalk](lipid_crosswalk.csv) and [edge eligibility ledger](edge_eligibility_ledger.csv) document detected/excluded, missing, ambiguous, nonfinite and invariant features. No SD variance requirement is imposed on this WD-only comparison. Both zones use the same ordered mouse IDs within each cohort; each animal supplies two RNA profiles but one bulk liver lipid profile.", "")
  replay <- bundle$setup$WP2_replay
  if (nrow(replay) && all(replay$WP2_available)) add(paste0("As a reproducibility check, CF453 WD coefficients reproduce the overlapping saved WP2 results with maximum absolute difference **", format(max(replay$max_abs_difference), scientific = TRUE, digits = 3L), "**. WP3 reads original inputs and frozen reference selections; it does not need prior WP2 inference to calculate new CF459 coefficients."), "")
  add(section(3L, "Whole-mouse null tests and their limits"), "",
    if (!is.null(fcg_comparison) && is.data.frame(fcg_comparison$exact_zero_null_draws))
      "The following two-WD-cohort null tables preserve the original saved WP3 calculations. They predate the exact-zero null-degeneracy correction documented in Section 1.1 and are shown as the historical analysis, not newly corrected bootstrap estimates." else "",
    "The primary exact test keeps observed CF453 coefficients fixed and enumerates all **4! = 24** complete CF459 lipid-to-RNA mouse pairings. The complementary test keeps observed CF459 coefficients fixed and enumerates all **7! = 5,040** CF453 pairings. The identity pairing and ties are included. Both tests use the two-sided absolute coefficient-vector r, equivalent to absolute ULM at a fixed edge count. These are separate conditional pairing-independence nulls, not a test of equal cohort correlations and not two independent replications.", "")
  table_md(data.frame(Comparison = label(network), `Nominal ULM p` = wp3_p(network$ULM_p_value),
    `Exact CF459 p (24)` = wp3_p(network$pairing_CF459_given_CF453_p),
    `Exact CF453 p (5040)` = wp3_p(network$pairing_CF453_given_CF459_p),
    `Null bootstrap p (defined only)` = wp3_p(network$null_bootstrap_p), check.names = FALSE))
  add(paste0("The approximate null bootstrap independently samples complete CF459 RNA and lipid mouse vectors from their empirical margins, holding CF453 coefficients fixed (**", analysis$metadata$B, "** requested draws; seed **", analysis$metadata$seed, "**). It tests the specified RNA-lipid independence null conditional on the empirical four-mouse cohort. It is not an ordinary paired bootstrap and does not estimate a probability that the observed biological network is true."), "",
    "With four mice, null resampling often has very little animal diversity. If any edge in the fixed vector becomes undefined, the entire vector draw is undefined; the code never deletes a different set of edges in each draw. The reported bootstrap p-value is conditional on defined draws. The completion range below assigns all undefined draws to either nonexceedance or exceedance; it is a mathematical missing-draw bound, **not a confidence interval or calibrated correction**.", "")
  table_md(data.frame(Comparison = label(network), `Exact 24 defined` = paste(network$permutation_CF459_valid_draws, "/", network$permutation_CF459_draws),
    `Exact 5040 defined` = paste(network$permutation_CF453_valid_draws, "/", network$permutation_CF453_draws),
    `Bootstrap defined / requested` = paste(network$bootstrap_valid_draws, "/", network$null_bootstrap_draws),
    `Completion range` = paste(wp3_p(network$bootstrap_completion_lower), wp3_p(network$bootstrap_completion_upper), sep = " to "), check.names = FALSE))
  low <- bundle$low_diversity
  if (!is.null(low) && nrow(low)) {
    low <- low[low$context == "male_pericentral", , drop = FALSE]
    if (nrow(low)) {
      details <- paste0(low$gene, " (", low$unique_expression_values, " observed expression levels; ",
        low$constant_null_rna_draws, "/", low$requested_null_draws, " constant RNA resamples)")
      add(paste0("A concrete source of pericentral undefined vectors is limited CF459 expression diversity: ",
        paste(details, collapse = "; "), ". A resample can omit an animal carrying the distinguishing expression level and make the gene constant. Complete exact pairing permutations preserve the original expression margins and remain defined. These genes and their fixed edges are retained; no post hoc rare-gene exclusion is used to improve the result."), "")
    }
  }
  ds <- analysis$resampling_draw_summary
  diversity <- do.call(rbind, lapply(unique(ds$method), function(method) {
    x <- ds[ds$method == method, ]; data.frame(Method = method, `Draws` = nrow(x),
      `One unique mouse` = sum(x$n_unique_mice == 1L), `Two or fewer unique mice` = sum(x$n_unique_mice <= 2L),
      `Median unique mice` = wp3_number(stats::median(x$n_unique_mice), 1L), check.names = FALSE)
  }))
  table_md(diversity)
  add("The four-mouse exact test is intrinsically coarse: with 24 defined configurations its smallest possible p-value is 1/24 (about 0.0417), before accounting for tied statistics. A weak or inconclusive animal-level result does not establish absence of shared regulation. Conversely, a strong nominal ULM result alone cannot overcome four biological replicates. [Null-check figures](plots/null_checks.pdf) display the full-vector nulls and defined-draw counts.", "",
    section(4L, "Stability and prespecified Mtmr11 example"), "",
    paste0("Ordinary paired bootstrap resamples whole paired RNA/lipid mice independently within each cohort (", analysis$metadata$paired_B, " requested draws). The same animal indices are used in both zones. This describes sensitivity to animal composition, without supplying an additional null-test p-value. Frozen-vector undefined draws remain visible."), "")
  table_md(data.frame(Comparison = label(network), `Paired defined / requested` = paste(network$paired_valid_draws, "/", network$paired_stability_draws),
    `Positive vector r among defined` = paste0(wp3_number(100 * network$paired_positive_fraction, 1), "%"),
    `Observed direction among defined` = paste0(wp3_number(100 * network$paired_same_direction_fraction, 1), "%"), check.names = FALSE))
  loco <- analysis$loco_summary
  loco <- loco[loco$scope == "network" & loco$universe == "full", , drop = FALSE]
  table_md(data.frame(Zone = wp3_label(loco$context), Omissions = loco$method,
    `Defined / planned` = paste(loco$valid_draws, "/", loco$draws),
    `Minimum r` = wp3_number(loco$minimum_r), `Maximum r` = wp3_number(loco$maximum_r),
    `Observed direction` = paste0(wp3_number(100 * loco$same_direction_fraction, 1), "%"), check.names = FALSE))
  add("There are **28 joint** one-CF453/one-CF459 omissions (six plus three retained mice), **seven CF453-only** omissions and **four CF459-only** omissions. These prespecified configurations are descriptive and share animals; they are not 39 independent validation experiments. [Stability figures](plots/stability_and_omissions.pdf) show the complete bank.", "")
  focus <- inference[inference$scope == "gene" & inference$gene == "Mtmr11", , drop = FALSE]
  corrected_focus <- NULL
  audit <- fcg_comparison$historical_null_correction_summary
  if (!is.null(fcg_comparison) && is.data.frame(audit) && nrow(focus)) {
    ai <- match(focus$specification, audit$specification)
    if (anyNA(ai)) stop("Missing historical-null correction for the predefined Mtmr11 report.")
    # Only the rendered table is updated; historical inference files are intact.
    focus$null_bootstrap_p <- audit$corrected_null_bootstrap_p[ai]
    corrected_focus <- TRUE
  }
  if (nrow(focus)) table_md(data.frame(Example = label(focus), Pairs = focus$n_edges,
    `Vector r` = wp3_number(focus$coefficient_vector_r), `ULM score` = wp3_number(focus$ULM_score),
    `Nominal ULM p` = wp3_p(focus$ULM_p_value), `Null bootstrap p (defined only)` = wp3_p(focus$null_bootstrap_p), check.names = FALSE))
  if (isTRUE(corrected_focus)) add("The Mtmr11 null-bootstrap p-values in this table use the corrected CF453-anchor null bank documented in Section 1.1. Observed coefficients, native ULM and nominal ULM p are unchanged; the original saved inference and scientific figures remain historical records.", "")
  add("Mtmr11 periportal was chosen from the original Figure 5 context before this comparison. The [coefficient figure](plots/Mtmr11_coefficients.pdf) uses its full and displayed shared species. The three [individual-mouse examples](plots/Mtmr11_selected_pairs.pdf) were selected by the largest absolute **original FCG displayed** coefficients, with a lipid-name tie break; selection never uses CF453/CF459 agreement. They remain illustrative associations. Cohort panels have separate axes because absolute export scales are not assumed comparable.", "")
  selected <- bundle$predefined_pairs
  if (nrow(selected)) table_md(data.frame(Lipid = selected$reference_lipid,
    `Frozen FCG r` = wp3_number(selected$r_reference), `CF453 WD r` = wp3_number(selected$r_CF453_WD),
    `CF459 WD r` = wp3_number(selected$r_CF459_WD), check.names = FALSE))
  add(section(5L, "Prespecified sensitivities and intake QC"), "")
  sens <- bundle$sensitivities
  if (is.data.frame(sens) && nrow(sens) && all(c("context", "universe", "coefficient_vector_r") %in% names(sens))) {
    if ("scope" %in% names(sens)) sens <- sens[sens$scope == "network", , drop = FALSE]
    name_col <- intersect(c("sensitivity", "type", "sensitivity_name"), names(sens))
    stype <- if (length(name_col)) as.character(sens[[name_col[1L]]]) else rep("Prespecified sensitivity", nrow(sens))
    table_md(data.frame(Sensitivity = stype, Comparison = label(sens),
      Pairs = sens$n_edges, `Vector r` = wp3_number(sens$coefficient_vector_r),
      `ULM score` = wp3_number(sens$ULM_score), check.names = FALSE))
  }
  add("CF453 QC-CV <=30% and supplied Fill >=0.5 sensitivities restrict only the CF453 annotation/QC ledger; no equivalent CF459 fields were supplied. The RNA data-layer sensitivity uses stored data expression while preserving fixed molecular membership, animal pairing and original selection. Pearson correlations are invariant to positive row scaling within a cohort, so agreement of stored data and scale.data results is expected and does not independently validate upstream normalization. Full sensitivity tables and resampling remain saved separately.", "")
  q <- function(qc, key) qc$value[match(key, qc$metric)]
  table_md(data.frame(Cohort = c("CF453", "CF459"), `Male WD mice` = c("7", "4"),
    `Supplied lipid rows` = c(q(bundle$cf453_qc, "supplied_lipid_feature_rows"), q(bundle$cf459_qc, "supplied_lipid_feature_rows")),
    `Eligible lipids` = c(q(bundle$cf453_qc, "primary_lipid_features"), q(bundle$cf459_qc, "primary_lipid_features")),
    `Analytical QC injections` = c(q(bundle$cf453_qc, "analytical_QC_injections"), "Unavailable in export"), check.names = FALSE))
  add("CF459 import correctly handles a UTF-8 BOM, semicolon delimiter, decimal commas, scientific notation and quoted annotation semicolons. Its 1,463 rows retain 1,395 eligible lipids after 68 repository-class exclusions, with no supplied internal-standard rows. All source annotations and 244 duplicate Alignment-ID excess rows remain auditable. No CF459 QC CV, Fill, adduct or RT/m/z evidence is invented.", "",
    "Both cohorts supply named Seurat 5.3.0 RNA Assay5 counts, data and scale.data layers with 20,123 genes. Stored scale.data is primary and unchanged. CF459 scaling diagnostics use all 40 original cell-type profiles or all 16 zone profiles before hepatocyte/male subsetting. Existing data replays stored scaling to approximately 1.25 x 10^-14. Counts are finite, nonnegative integers; aggregate library sizes and marker ordering are diagnostics, not nucleus counts or proof that nucleus-level QC/gates passed.", "",
    paste0("The CF459 RNA labels f-ctrl/m-ctrl identify recorded sex; the meeting establishes WD-only exposure. Four female mice appear only in intake and scaling-roster QC, not new female networks. The supplied assay exports lack age, genetic-background, diet-duration, litter and batch fields. ",
      if (is.null(fcg_comparison)) "The user clarified CF459 age as 3 months and CF453 age as 10 months on 14 October 2026; these are ages at sampling, not diet durations. " else "The age/diet descriptions and their provenance are recorded in the requested-comparison section; none establishes controlled matching or equal diet exposure. ",
      "[QC and coverage figures](plots/QC_and_coverage.pdf) show the design and aggregate diagnostics."), "",
    section(6L, "Interpretation and reproducibility"), "",
    "Cross-cohort coefficient agreement supports robustness only to the extent observed in these measured, strictly matched male WD features. It does not establish equal cohort effects, eliminate age/background/batch confounding, prove causal gene regulation, or identify WD-induced change. CF459 has no SD arm. Conversely, low or unstable agreement with four mice cannot establish absence of shared biology. Report the full-network r and the displayed subset together, alongside distinct conditional null tests and animal-composition stability.", "",
    "Primary results, full exclusions, bootstrap/permutation indices and omissions are retained under this cohort-specific output directory. Original FCG, siRNA and completed WP1/WP2 numerical results are preserved. Input byte hashes and remaining metadata limits are documented in data/WP3_INPUTS.md and data/wp3_input_migration_manifest.csv. The report is regenerated from saved numeric results and report_bundle.rds; no statistical rerun is needed for --report-only.", "")
  artifacts <- c("crosscohort_inference.csv", "network_edges.csv", "coverage.csv", "lipid_crosswalk.csv", "edge_eligibility_ledger.csv",
    "sensitivities.csv", "loco_summary.csv", "resampling_draw_summary.csv", "predefined_pairs.csv", "analysis.rds", "sensitivity_analyses.rds")
  if (!is.null(fcg_comparison)) artifacts <- c("fcg_comparison/inference.csv", "fcg_comparison/replay.csv",
    "fcg_comparison/paired_summary.csv", "fcg_comparison/omission_summary.csv", "fcg_comparison/analysis.rds", artifacts)
  for (file in artifacts[file.exists(file.path(out, artifacts))]) add(paste0("- [", file, "](", file, ")"))
  add("", "Reproduce from the repository root, using a fresh output directory for a full run:", "", "```bash",
    paste0("Rscript scripts/workflows/13_wd_crosscohort_robustness.R --output=results/wd_crosscohort_reproduction --bootstraps=", analysis$metadata$B,
      " --paired-bootstraps=", analysis$metadata$paired_B, " --seed=", analysis$metadata$seed),
    "Rscript scripts/workflows/13_wd_crosscohort_robustness.R --report-only --output=results/wd_crosscohort_reproduction",
    "python3 scripts/reports/build_wp3_report_pdf.py --input results/wd_crosscohort_reproduction/README.md", "```", "",
    if (is.null(fcg_comparison)) "The PDF appendix preserves six scientific vector pages in report order: network coefficients; exact and null checks; stability/omissions; prespecified Mtmr11 coefficients; original-reference-selected individual pairs; QC/coverage." else
      "The PDF appendix begins with two new scientific vector pages for the matched original-FCG/CF453 anchors and direct-closeness stability, followed by the six original WP3 figure pages: WD network coefficients; historical exact and null checks; stability/omissions; prespecified Mtmr11 coefficients; original-reference-selected individual pairs; QC/coverage.")
  for (path in bundle$plot_files) if (!file.exists(file.path(out, path))) stop("Missing report figure: ", path)
  utils::write.table(text, file.path(out, "README.md"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(file.path(out, "README.md"))
}
