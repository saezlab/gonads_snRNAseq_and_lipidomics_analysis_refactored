# Small read-only numerical/export QA for the two Figure 6A layout candidates.
# Run: Rscript tests/qa_cf453_figure6a_layouts.R OUTPUT_DIR [ANALYSIS_RDS]
# No statistical model, selection workflow, correlation or resampling is rerun.

qa_cf453_figure6a_layouts <- function(output,
    analysis_file = file.path(dirname(output), "rewiring_analysis.rds")) {
  if (dir.exists(".local/R-library")) .libPaths(c(".local/R-library", .libPaths()))
  a <- readRDS(analysis_file); anchors <- a$anchors
  checks <- list()
  record <- function(label, difference = 0, detail = "") {
    if (!is.finite(difference) || difference > 1e-10) stop("Independent Figure6A layout QA failed: ", label, "; difference=", difference)
    checks[[length(checks) + 1L]] <<- data.frame(check = label,
      max_abs_difference = difference, detail = detail, stringsAsFactors = FALSE)
  }
  same <- function(actual, expected, label) {
    stopifnot(length(actual) == length(expected), identical(is.na(actual), is.na(expected)))
    if (is.numeric(expected) || is.logical(expected)) {
      actual <- as.numeric(actual); expected <- as.numeric(expected)
      good <- is.finite(actual) & is.finite(expected)
      stopifnot(all(good | (is.na(actual) & is.na(expected))))
      record(label, if (any(good)) max(abs(actual[good] - expected[good])) else 0)
    } else {
      stopifnot(identical(as.character(actual), as.character(expected)))
      record(label)
    }
  }
  key <- function(d, fields) do.call(paste, c(d[, fields, drop = FALSE], sep = "|"))
  contexts <- c("female_PP", "female_PC", "male_PP", "male_PC")
  effects <- c("effect_female_SD", "effect_female_WD", "effect_male_SD", "effect_male_WD")
  invariant <- c("gene", "context", "sex", "zone", "threeway_effect", "threeway_p", "threeway_q", effects)
  stopifnot(nrow(anchors) > 0L)
  expected <- lapply(unique(anchors$gene), function(gene) {
    z <- anchors[anchors$gene == gene, , drop = FALSE]
    stopifnot(nrow(unique(z[, invariant, drop = FALSE])) == 1L,
      !anyDuplicated(z$anchor_diet), all(z$anchor_diet %in% c("SD", "WD")))
    row <- z[1L, invariant, drop = FALSE]
    row$selection <- if (setequal(z$anchor_diet, c("SD", "WD"))) "Both" else z$anchor_diet
    row$n_anchors <- nrow(z); row$direct_supported <- row$threeway_q <= .05
    row
  })
  expected <- do.call(rbind, expected)
  expected <- expected[order(match(expected$context, contexts), -abs(expected$threeway_effect), expected$gene), , drop = FALSE]
  expected$display_rank <- as.integer(ave(seq_len(nrow(expected)), expected$context, FUN = seq_along))
  expected$global_rank <- seq_len(nrow(expected)); rownames(expected) <- NULL
  ledger <- utils::read.csv(file.path(output, "common_gene_ledger.csv"), stringsAsFactors = FALSE)
  stopifnot(identical(ledger$gene, expected$gene), !anyDuplicated(ledger$gene))
  for (field in names(expected)) same(ledger[[field]], expected[[field]], paste("Common unique-gene ledger", field))
  stopifnot(sum(ledger$selection %in% c("SD", "Both")) == sum(anchors$anchor_diet=="SD"),
    sum(ledger$selection %in% c("WD", "Both")) == sum(anchors$anchor_diet=="WD"),
    sum(ledger$selection == "Both") == sum(table(anchors$gene)==2L),
    sum(ledger$direct_supported) == sum(expected$direct_supported))
  record("Saved selections deduplicate without losing anchor membership or direct-test markers",
    detail=paste(nrow(anchors),"selections;",nrow(expected),"genes;",sum(expected$direct_supported),"markers"))
  expected_long <- expand.grid(gene = expected$gene, sex = c("female", "male"), diet = c("SD", "WD"), stringsAsFactors = FALSE)
  expected_long$effect <- vapply(seq_len(nrow(expected_long)), function(j) {
    expected[match(expected_long$gene[j], expected$gene),
             paste("effect", expected_long$sex[j], expected_long$diet[j], sep = "_")]
  }, numeric(1L))
  point_key <- function(d) key(d, c("gene", "sex", "diet"))
  effect_range <- range(expected_long$effect)
  manifest <- utils::read.csv(file.path(output, "panel_manifest.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(manifest) == 6L, length(unique(manifest$panel)) == 2L,
    length(unique(manifest$source_table)) == 2L, all(table(manifest$panel, manifest$format) == 1L),
    all(file.exists(manifest$file)), all(file.exists(manifest$source_table)))
  saved_sources <- list()
  for (file in unique(manifest$source_table)) {
    s <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    label <- basename(file)
    stopifnot(nrow(s) == 4L*nrow(expected), !anyDuplicated(point_key(s)),
              setequal(point_key(s), point_key(expected_long)),
              all(table(s$gene) == 4L), all(table(s$sex, s$diet) == nrow(expected)))
    idx <- match(point_key(s), point_key(expected_long)); gi <- match(s$gene, expected$gene)
    same(s$effect, expected_long$effect[idx], paste(label, "all saved raw PP-minus-PC effects"))
    same(s$anchor_sex, expected$sex[gi], paste(label, "anchor sex remains distinct from plotted sex"))
    same(s$state, paste(s$sex, s$diet, sep = "_"), paste(label, "four state identities"))
    for (field in setdiff(names(expected), "sex")) {
      same(s[[field]], expected[[field]][gi], paste(label, "shared saved metadata", field))
    }
    if ("colour_limit" %in% names(s)) same(s$colour_limit, rep(max(abs(effect_range)), nrow(s)),
      "Heatmap source records one global symmetric raw-effect colour range")
    if ("common_effect_min" %in% names(s)) {
      same(s$common_effect_min, rep(effect_range[1L], nrow(s)), "Dot source common minimum is the raw-effect minimum")
      same(s$common_effect_max, rep(effect_range[2L], nrow(s)), "Dot source common maximum is the raw-effect maximum")
    }
    saved_sources[[label]] <- s
  }
  record("Both alternatives retain every untransformed estimate and common per-programme display order")

  # Inspect actual plot objects without rendering/exporting files. Source data
  # above are independently reconstructed from the historical saved anchors.
  env <- new.env(parent = globalenv())
  for (file in c("scripts/core/cf453_figure6a_layouts.R", "scripts/core/cf453_figure6a_heatmap.R", "scripts/core/cf453_figure6a_dots.R")) {
    sys.source(file, envir = env)
  }
  heatmap <- env$cf6a_heatmap(expected); dots <- env$cf6a_dots(expected)
  hb <- ggplot2::ggplot_build(heatmap$plot); db <- ggplot2::ggplot_build(dots$plot)
  fill <- hb$plot$scales$get_scales("fill")
  colour_limits <- as.numeric(fill$get_limits())
  stopifnot(length(colour_limits) == 2L, all(is.finite(colour_limits)),
            abs(sum(colour_limits)) < 1e-12,
            colour_limits[1L] <= effect_range[1L], colour_limits[2L] >= effect_range[2L])
  trans_name <- function(scale) {
    transformation <- if (is.function(scale$get_transformation)) scale$get_transformation() else scale$trans
    transformation$name
  }
  stopifnot(identical(trans_name(fill), "identity"))
  tile_layer <- which(vapply(heatmap$plot$layers, function(z) inherits(z$geom, "GeomTile"), logical(1L)))
  stopifnot(length(tile_layer) == 1L, nrow(hb$data[[tile_layer]]) == 4L*nrow(expected))
  same(hb$data[[tile_layer]]$fill, fill$map(heatmap$source$effect),
       "All heatmap tile colours directly map raw effects without row scaling")
  support_layer <- which(vapply(heatmap$plot$layers, function(z) inherits(z$geom, "GeomPoint"), logical(1L)))
  stopifnot(length(support_layer) == 1L, nrow(hb$data[[support_layer]]) == sum(expected$direct_supported))
  record("Heatmap uses one identity colour scale, complete coverage and the saved gene-level test markers",
         detail = paste("Colour limits:", paste(colour_limits, collapse = ", ")))
  stopifnot(length(db$layout$panel_params) == 4L)
  x_domains <- lapply(db$layout$panel_params, function(p) p$x.range)
  stopifnot(all(vapply(x_domains, function(v) identical(v, x_domains[[1L]]), logical(1L))),
            x_domains[[1L]][1L] <= effect_range[1L], x_domains[[1L]][2L] >= effect_range[2L],
            all(vapply(db$layout$panel_scales_x, function(z) identical(trans_name(z), "identity"), logical(1L))))
  point_layer <- which(vapply(dots$plot$layers, function(z) inherits(z$geom, "GeomPoint"), logical(1L)))
  stopifnot(length(point_layer) == 1L, nrow(db$data[[point_layer]]) == 4L*nrow(expected),
            all(is.finite(db$data[[point_layer]]$x)))
  points <- db$data[[point_layer]]
  same(sort(points$x[points$shape == 1]), sort(expected_long$effect[expected_long$diet == "SD"]),
       "All open dot positions are exact SD raw effects")
  same(sort(points$x[points$shape == 16]), sort(expected_long$effect[expected_long$diet == "WD"]),
       "All filled dot positions are exact WD raw effects")
  segment_layer <- which(vapply(dots$plot$layers, function(z) inherits(z$geom, "GeomSegment"), logical(1L)))
  stopifnot(length(segment_layer) == 1L)
  segments <- dots$plot$layers[[segment_layer]]$data
  stopifnot(nrow(segments) == 2L*nrow(expected), !anyDuplicated(key(segments, c("gene", "sex"))), all(segments$diet == "SD"))
  si <- match(segments$gene, expected$gene)
  expected_start <- expected_end <- numeric(nrow(segments))
  for (j in seq_len(nrow(segments))) {
    expected_start[j] <- expected[si[j], paste("effect", segments$sex[j], "SD", sep = "_")]
    expected_end[j] <- expected[si[j], paste("effect", segments$sex[j], "WD", sep = "_")]
  }
  same(segments$effect, expected_start, "Every dot segment starts at the same gene/sex SD effect")
  same(segments$effect_end, expected_end, "Every dot segment ends at the same gene/sex WD effect")
  for (gene in unique(dots$source$gene)) for (sex in c("female", "male")) {
    rows <- dots$source$gene == gene & dots$source$sex == sex
    stopifnot(sum(rows) == 2L, length(unique(dots$source$y[rows])) == 1L)
  }
  record("Paired dots use one identity x-domain across four blocks, separate sex tracks and correctly matched SD-to-WD segments",
         detail = paste("Axis limits including label gutter:", paste(x_domains[[1L]], collapse = ", ")))

  for (file in manifest$file[manifest$format == "SVG"]) {
    text <- paste(readLines(file, warn = FALSE), collapse = "\n")
    stopifnot(!grepl("<image\\b", text, perl = TRUE))
    record(paste("Vector-only SVG", basename(file)))
  }
  for (j in which(manifest$format == "PNG")) {
    con <- file(manifest$file[j], "rb"); raw <- readBin(con, "raw", 24L); close(con)
    stopifnot(identical(as.integer(raw[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
    decode <- function(v) sum(as.integer(v) * 256^(3:0))
    stopifnot(abs(decode(raw[17:20]) - 300 * manifest$width_in[j]) <= 1,
              abs(decode(raw[21:24]) - 300 * manifest$height_in[j]) <= 1)
    record(paste("300-dpi PNG dimensions", basename(manifest$file[j])))
  }
  retained <- utils::read.csv(file.path(output, "retained_original_sha256.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(retained) > 0L, !anyDuplicated(retained$path), all(file.exists(retained$path)))
  hashes <- vapply(retained$path, digest::digest, character(1L), file = TRUE, algo = "sha256")
  stopifnot(identical(unname(hashes), unname(retained$sha256)))
  replay <- file.path(dirname(analysis_file), "QA", "cli_replay", "rewiring_analysis.rds")
  if(file.exists(replay)) stopifnot(digest::digest(analysis_file, file = TRUE, algo = "sha256") == digest::digest(replay, file = TRUE, algo = "sha256"))
  record("Protected parent files retain before-export hashes; optional saved CLI replay checked when present",
    detail=paste(nrow(retained),"protected files; CLI replay present:",file.exists(replay)))
  provenance <- utils::read.csv(file.path(output, "provenance_sha256.csv"), stringsAsFactors = FALSE)
  hashes <- vapply(provenance$path, digest::digest, character(1L), file = TRUE, algo = "sha256")
  stopifnot(identical(unname(hashes), unname(provenance$sha256)))
  record("All recorded source, saved-analysis, test and plan SHA256 values match")
  qa_dir <- file.path(output, "QA")
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- do.call(rbind, checks)
  utils::write.csv(summary, file.path(qa_dir, "independent_checks.csv"), row.names = FALSE)
  saveRDS(list(checks = summary, colour_limits = colour_limits, dot_axis_limits = x_domains[[1L]],
    raw_effect_range = effect_range, session = utils::sessionInfo()), file.path(qa_dir, "independent_qa.rds"))
  cat("Independent Figure6A layout QA passed:", nrow(summary), "check groups.\n")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || identical(args, "--help")) {
    cat("Usage: Rscript tests/qa_cf453_figure6a_layouts.R OUTPUT_DIR [ANALYSIS_RDS]\n")
  } else {
    qa_cf453_figure6a_layouts(args[1L], if (length(args) >= 2L) args[2L] else
      file.path(dirname(args[1L]), "rewiring_analysis.rds"))
  }
}
