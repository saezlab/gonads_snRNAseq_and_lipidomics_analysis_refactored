# Independent actual-data replay of the selected CF453 programmes and displays.
# Run: Rscript tests/qa_cf453_diet_rewiring.R REWIRING_RDS [QA_DIR]
# Final exports: Rscript tests/qa_cf453_diet_rewiring.R --exports-only REWIRING_RDS [QA_DIR]
# No production analysis or plotting helper is sourced.

qa_cf453_diet_rewiring <- function(analysis_file,
    qa_dir = "tests/outputs/cf453_diet_rewiring") {
  if (dir.exists(".local/R-library")) .libPaths(c(".local/R-library", .libPaths()))
  a <- readRDS(analysis_file)
  parent <- readRDS(a$source_paths[["parent"]])
  direct <- readRDS(a$source_paths[["direct"]])
  checks <- list()
  record <- function(label, difference = 0, detail = "") {
    if (!is.finite(difference) || difference > 1e-8) stop("Independent CF453 programme QA failed: ", label, "; difference=", difference)
    checks[[length(checks) + 1L]] <<- data.frame(check = label,
      max_abs_difference = difference, detail = detail, stringsAsFactors = FALSE)
  }
  same <- function(actual, expected, label) {
    actual <- as.numeric(actual); expected <- as.numeric(expected)
    if (length(actual) != length(expected) || !identical(is.na(actual), is.na(expected))) {
      stop("Independent programme QA length/missingness mismatch: ", label)
    }
    valid <- is.finite(actual) & is.finite(expected)
    stopifnot(all(valid | (is.na(actual) & is.na(expected)) | (is.infinite(actual) & actual == expected)))
    record(label, if (any(valid)) max(abs(actual[valid] - expected[valid])) else 0)
  }
  key <- function(d, fields) do.call(paste, c(d[, fields, drop = FALSE], sep = "|"))
  anchored_key <- function(d) key(d, c("gene", "anchor_diet", "sex", "zone"))
  programme_key <- function(d) key(d, c("context", "gene"))
  edge_key <- function(d) key(d, c("context", "gene", "lipid"))
  md <- parent$metadata; pm <- parent$profile_metadata
  stopifnot(identical(a$metadata, md), nrow(md) == 28L, nrow(pm) == 56L,
    all(table(md$sex, md$diet) == 7L), all(table(pm$mouse_id, pm$zone) == 1L),
    identical(parent$modeled_rna, direct$modeled_rna), identical(md, direct$metadata))
  record("Exactly 28 independent matched mice, seven per sex/diet, with paired RNA zones")

  # Rebuild every anchor eligibility row and its within-anchor-diet ranking.
  rr <- parent$rna_results
  eligibility <- selections <- list()
  for (diet in c("SD", "WD")) for (sex in c("female", "male")) {
    inter <- rr[rr$contrast == paste0("sex_zone_", diet), , drop = FALSE]
    own <- rr[rr$contrast == paste("zone", sex, diet, sep = "_"), , drop = FALSE]
    other <- rr[rr$contrast == paste("zone", if (sex == "female") "male" else "female", diet, sep = "_"), , drop = FALSE]
    own <- own[match(inter$gene, own$gene), , drop = FALSE]
    other <- other[match(inter$gene, other$gene), , drop = FALSE]
    d <- data.frame(gene = inter$gene, anchor_diet = diet, sex = sex,
      zone = ifelse(own$effect > 0, "PP", "PC"), anchor_effect = own$effect,
      anchor_q = inter$q, own_q = own$q, other_t = other$t,
      eligible = inter$q <= .05 & own$q <= .05 & abs(other$t) < 1,
      rank = NA_integer_, stringsAsFactors = FALSE)
    d$context <- paste(d$sex, d$zone, sep = "_")
    for (zone in c("PP", "PC")) {
      rows <- which(d$eligible & d$zone == zone)
      rows <- rows[order(-abs(d$anchor_effect[rows]), d$anchor_q[rows], d$own_q[rows], d$gene[rows])]
      d$rank[rows] <- seq_along(rows)
      selections[[length(selections) + 1L]] <- d[head(rows, 5L), , drop = FALSE]
      record(paste(diet, sex, zone, "strict RNA anchor eligibility and magnitude ranking"),
             detail = paste(length(rows), "eligible;", min(5L, length(rows)), "selected"))
    }
    eligibility[[length(eligibility) + 1L]] <- d
  }
  eligibility <- do.call(rbind, eligibility); anchors <- do.call(rbind, selections)
  stopifnot(!anyDuplicated(anchored_key(a$eligibility)),
            setequal(anchored_key(eligibility), anchored_key(a$eligibility)),
            identical(anchored_key(anchors), anchored_key(a$anchors)),
            all(table(a$anchors$anchor_diet, a$anchors$sex, a$anchors$zone) == 5L))
  actual <- a$eligibility[match(anchored_key(eligibility), anchored_key(a$eligibility)), , drop = FALSE]
  for (field in c("anchor_effect", "anchor_q", "own_q", "other_t", "eligible", "rank")) {
    same(actual[[field]], eligibility[[field]], paste("All eligibility rows", field))
    same(a$anchors[[field]], anchors[[field]], paste("All selected anchor rows", field))
  }
  for (diet in c("SD", "WD")) for (sex in c("female", "male")) {
    r <- rr[rr$contrast == paste("zone", sex, diet, sep = "_"), , drop = FALSE]
    ii <- match(a$anchors$gene, r$gene)
    coordinate <- paste0(if (sex == "female") "x_" else "y_", diet)
    same(a$anchors[[coordinate]], r$t[ii], paste("Tracked anchor coordinates", coordinate))
    same(a$anchors[[paste("effect", sex, diet, sep = "_")]], r$effect[ii],
         paste("Tracked anchor effects", sex, diet))
  }
  three <- direct$overview[direct$overview$contrast == "diet_sex_zone", , drop = FALSE]
  for (field in c("effect", "p", "q")) {
    same(a$anchors[[paste0("threeway_", field)]], three[[field]][match(a$anchors$gene, three$gene)],
         paste("Arrow direct Workflow15 three-way annotation", field))
  }
  conditional_q <- vapply(seq_len(nrow(a$anchors)), function(j) {
    rows <- direct$conditional$gene == a$anchors$gene[j] &
      direct$conditional$contrast == paste0("diet_zone_", a$anchors$sex[j])
    stopifnot(sum(rows) == 1L)
    direct$conditional$q[rows]
  }, numeric(1L))
  same(a$anchors$conditional_q, conditional_q, "All anchor own-sex direct zonation-change q annotations")
  programmes <- unique(anchors[, c("context", "sex", "zone", "gene"), drop = FALSE])
  stopifnot(setequal(programme_key(programmes), programme_key(a$programmes)), nrow(programmes) == nrow(a$programmes))
  for (j in seq_len(nrow(a$programmes))) {
    p <- a$programmes[j, ]; rows <- anchors$context == p$context & anchors$gene == p$gene
    stopifnot(p$anchor_diets == paste(sort(unique(anchors$anchor_diet[rows])), collapse = ";"))
    same(p$best_anchor_q, min(anchors$anchor_q[rows]), paste(p$context, p$gene, "best anchor evidence"))
  }
  record("Primary programmes are exactly the selected-anchor union, without lipid or destination-diet selection")

  lipid <- parent$lipid_matrix
  ledger <- parent$lipid_ledger[match(rownames(lipid), parent$lipid_ledger$feature_id), , drop = FALSE]
  classes <- ledger$lipid_class; class_names <- sort(unique(classes))
  stopifnot(nrow(lipid) == 1129L, identical(colnames(lipid), md$mouse_id), all(ledger$keep_primary))
  inputs <- function(genes, sex, diet, zone) {
    mice <- md$mouse_id[md$sex == sex & md$diet == diet]
    profiles <- pm$profile_id[match(paste(mice, zone), paste(pm$mouse_id, pm$zone))]
    stopifnot(length(mice) == 7L, !anyNA(profiles))
    list(x = parent$modeled_rna[genes, profiles, drop = FALSE],
         y = lipid[, mice, drop = FALSE], mice = mice, profiles = profiles)
  }
  metrics <- function(r, threshold = .65) {
    result <- lapply(class_names, function(class) {
      inside <- classes == class; outside <- !inside
      nr <- nrow(r); nc <- sum(inside); no <- sum(outside)
      hits <- sum(abs(r[, inside, drop = FALSE]) >= threshold)
      other <- sum(abs(r[, outside, drop = FALSE]) >= threshold)
      data.frame(lipid_class = class, n_genes = nr, n_lipids = nc,
        n_pairs_class = nr * nc, n_pairs_rest = nr * no,
        hits = hits, hits_rest = other, rate = hits / (nr * nc), rest = other / (nr * no),
        E = hits / (nr * nc) - other / (nr * no),
        mean_abs_excess = mean(abs(r[, inside, drop = FALSE])) - mean(abs(r[, outside, drop = FALSE])))
    })
    do.call(rbind, result)
  }
  bank <- a$bootstrap_indices; B <- a$settings$B
  groups <- c("female_SD", "female_WD", "male_SD", "male_WD")
  stopifnot(identical(dim(bank), c(7L, as.integer(B), 4L)),
            identical(dimnames(bank)[[3L]], groups), !anyNA(bank),
            all(bank == trunc(bank)), all(bank >= 1L & bank <= 7L))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(a$settings$seed)
  expected_bank <- array(sample.int(7L, 7L * B * 4L, TRUE), c(7L, B, 4L), dimnames = dimnames(bank))
  if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv)
  stopifnot(identical(bank, expected_bank))
  record("All complete-mouse bootstrap banks exactly replay seed and shared zone pairing")
  all_expected_edges <- list(); boot_replay_count <- 0L
  context_order <- c("female_PP", "female_PC", "male_PP", "male_PC")
  for (context in context_order) {
    sex <- sub("_.*$", "", context); zone <- sub("^.*_", "", context)
    genes <- sort(programmes$gene[programmes$context == context])
    e <- a$edges[a$edges$context == context, , drop = FALSE]
    expected_grid <- expand.grid(gene = genes, lipid = rownames(lipid), stringsAsFactors = FALSE)
    stopifnot(identical(e$gene, expected_grid$gene), identical(e$lipid, expected_grid$lipid),
              all(e$sex == sex), all(e$zone == zone),
              all(e$lipid_class == classes[match(e$lipid, rownames(lipid))]))
    inp <- obs <- loo <- met <- list()
    for (diet in c("SD", "WD")) {
      v <- inputs(genes, sex, diet, zone); inp[[diet]] <- v
      order <- a$sample_order[a$sample_order$context == context & a$sample_order$diet == diet, , drop = FALSE]
      stopifnot(identical(order$mouse_id, v$mice), identical(order$profile_id, v$profiles),
                identical(as.integer(order$index), 1:7))
      obs[[diet]] <- stats::cor(t(v$x), t(v$y)); met[[diet]] <- metrics(obs[[diet]])
      same(e[[paste0("r_", diet)]], obs[[diet]], paste(context, diet, "every independently mouse-paired correlation"))
      loo[[diet]] <- vapply(1:7, function(j)
        as.vector(suppressWarnings(stats::cor(t(v$x[, -j, drop = FALSE]), t(v$y[, -j, drop = FALSE])))), numeric(nrow(e)))
      same(a$loo[[context]][[diet]], loo[[diet]], paste(context, diet, "all seven complete-mouse omissions"))
      for (threshold in c(.60, .65, .70, .80)) {
        tm <- metrics(obs[[diet]], threshold)
        saved <- a$threshold_sensitivity[a$threshold_sensitivity$context == context &
          a$threshold_sensitivity$diet == diet & a$threshold_sensitivity$threshold == threshold, , drop = FALSE]
        saved <- saved[match(tm$lipid_class, saved$lipid_class), , drop = FALSE]
        for (field in setdiff(names(tm), "lipid_class")) same(saved[[field]], tm[[field]],
          paste(context, diet, threshold, "threshold sensitivity", field))
      }
    }
    delta_r <- as.vector(obs$WD - obs$SD)
    delta_abs <- as.vector(abs(obs$WD) - abs(obs$SD))
    same(e$delta_r, delta_r, paste(context, "all signed correlation changes"))
    same(e$delta_abs, delta_abs, paste(context, "all absolute-strength changes"))
    gain <- abs(as.vector(obs$WD)) >= .8 & abs(as.vector(obs$SD)) <= .4 & delta_abs >= .5
    loss <- abs(as.vector(obs$SD)) >= .8 & abs(as.vector(obs$WD)) <= .4 & delta_abs <= -.5
    direction <- ifelse(gain, "gained", ifelse(loss, "lost", "other"))
    stopifnot(identical(e$direction, direction))
    robust <- rep(FALSE, nrow(e)); separation <- rep(NA_real_, nrow(e))
    for (dir in c("gained", "lost")) {
      rows <- which(direction == dir)
      high <- if (dir == "gained") "WD" else "SD"; low <- if (dir == "gained") "SD" else "WD"
      for (j in rows) {
        hi <- loo[[high]][j, ]; lo <- loo[[low]][j, ]
        if (!all(is.finite(c(hi, lo)))) next
        separation[j] <- min(as.vector(outer(abs(hi), abs(lo), "-")))
        robust[j] <- all(abs(hi) >= .65 & sign(hi) == sign(e[[paste0("r_", high)]][j])) && all(abs(lo) <= .65)
      }
    }
    same(e$loo_robust, robust, paste(context, "every independent-arm high/low omission gate"))
    same(e$loo_separation_min, separation, paste(context, "minimum strength separation across all49 omission pairs"))
    candidates <- which(robust)
    stopifnot(identical(as.integer(a$bootstrap[[context]]$candidate_rows), as.integer(candidates)))
    expected_E <- expected_candidate_r <- valid <- list()
    for (diet in c("SD", "WD")) {
      expected_E[[diet]] <- matrix(NA_real_, length(class_names), B, dimnames = list(class_names, NULL))
      expected_candidate_r[[diet]] <- matrix(NA_real_, length(candidates), B)
      valid[[diet]] <- logical(B); v <- inp[[diet]]
      index <- bank[, , paste(sex, diet, sep = "_")]
      for (b in seq_len(B)) {
        ii <- index[, b]
        r <- suppressWarnings(stats::cor(t(v$x[, ii, drop = FALSE]), t(v$y[, ii, drop = FALSE])))
        if (length(candidates)) expected_candidate_r[[diet]][, b] <- as.vector(r)[candidates]
        valid[[diet]][b] <- all(is.finite(r))
        if (valid[[diet]][b]) expected_E[[diet]][, b] <- metrics(r)$E
      }
      same(a$bootstrap[[context]]$valid[[diet]], valid[[diet]], paste(context, diet, "context-wide all-features-defined bootstrap mask"))
      same(a$bootstrap[[context]]$E[[diet]], expected_E[[diet]], paste(context, diet, "all fixed-denominator class bootstrap excesses"))
      same(a$bootstrap[[context]]$candidate_r[[diet]], expected_candidate_r[[diet]], paste(context, diet, "every edge-specific candidate bootstrap correlation"))
      boot_replay_count <- boot_replay_count + length(expected_candidate_r[[diet]])
    }
    sm <- a$class_summary[a$class_summary$context == context, , drop = FALSE]
    sm <- sm[match(class_names, sm$lipid_class), , drop = FALSE]
    for (field in c("n_genes", "n_lipids", "n_pairs_class", "n_pairs_rest")) {
      same(sm[[field]], met$SD[[field]], paste(context, "fixed class/rest denominator", field))
      same(met$SD[[field]], met$WD[[field]], paste(context, "same SD/WD denominator", field))
    }
    for (diet in c("SD", "WD")) for (field in c("hits", "hits_rest", "rate", "rest", "E", "mean_abs_excess")) {
      same(sm[[paste0(field, "_", diet)]], met[[diet]][[field]], paste(context, diet, "observed class/rest metric", field))
    }
    expected_delta_E <- met$WD$E - met$SD$E
    same(sm$delta_E, expected_delta_E, paste(context, "difference in class-versus-rest excess rates"))
    same(sm$delta_rate, met$WD$rate - met$SD$rate, paste(context, "unadjusted class pair-rate changes"))
    same(sm$delta_mean_abs_excess, met$WD$mean_abs_excess - met$SD$mean_abs_excess,
         paste(context, "continuous absolute-correlation excess changes"))
    bd <- expected_E$WD - expected_E$SD
    qfun <- function(v, p) if (any(is.finite(v))) as.numeric(stats::quantile(v, p, na.rm = TRUE)) else NA_real_
    same(sm$bootstrap_low, apply(bd, 1L, qfun, p = .025), paste(context, "class change percentile lower limits"))
    same(sm$bootstrap_high, apply(bd, 1L, qfun, p = .975), paste(context, "class change percentile upper limits"))
    same(sm$bootstrap_valid, rowSums(is.finite(bd)), paste(context, "class change valid complete-context draws"))
    stopifnot(all(sm$bootstrap_requested == B))
    bd_direction <- vapply(seq_len(nrow(bd)), function(j) {
      v <- bd[j, is.finite(bd[j, ])]
      if (!length(v) || expected_delta_E[j] == 0) NA_real_ else mean(sign(v) == sign(expected_delta_E[j]))
    }, numeric(1L))
    same(sm$bootstrap_direction, bd_direction, paste(context, "class directional stability fractions"))
    expected_valid <- expected_direction <- expected_low <- expected_high <- rep(NA_real_, nrow(e))
    if (length(candidates)) {
      dd <- abs(expected_candidate_r$WD) - abs(expected_candidate_r$SD)
      expected_valid[candidates] <- rowSums(is.finite(dd))
      expected_low[candidates] <- apply(dd, 1L, qfun, p = .025)
      expected_high[candidates] <- apply(dd, 1L, qfun, p = .975)
      expected_direction[candidates] <- vapply(seq_len(nrow(dd)), function(j) {
        v <- dd[j, is.finite(dd[j, ])]
        if (!length(v)) NA_real_ else mean(sign(v) == sign(delta_abs[candidates[j]]))
      }, numeric(1L))
    }
    same(e$bootstrap_valid, expected_valid, paste(context, "candidate edge-specific paired-arm valid draw counts"))
    same(e$bootstrap_direction, expected_direction, paste(context, "candidate absolute-strength change directional fractions"))
    same(e$bootstrap_low, expected_low, paste(context, "candidate absolute-strength lower percentiles"))
    same(e$bootstrap_high, expected_high, paste(context, "candidate absolute-strength upper percentiles"))
    eligible <- robust & !is.na(expected_direction) & expected_direction >= .90 & expected_valid >= ceiling(.90 * B)
    same(e$network_eligible, eligible, paste(context, "every bootstrap-supported exploratory network gate"))
    for (diet in c("SD", "WD")) {
      rv <- as.vector(obs[[diet]])
      p <- 2 * stats::pt(-abs(rv) * sqrt(5 / pmax(1 - rv^2, .Machine$double.eps)), df = 5)
      same(e[[paste0("p_", diet)]], p, paste(context, diet, "all nominal Pearson probabilities"))
    }
    broad_genes <- sort(unique(eligibility$gene[eligibility$eligible & eligibility$context == context]))
    for (diet in c("SD", "WD")) {
      v <- inputs(broad_genes, sex, diet, zone)
      expected_broad <- metrics(stats::cor(t(v$x), t(v$y)))
      saved <- a$broad_sensitivity[a$broad_sensitivity$context == context & a$broad_sensitivity$diet == diet, , drop = FALSE]
      saved <- saved[match(expected_broad$lipid_class, saved$lipid_class), , drop = FALSE]
      for (field in setdiff(names(expected_broad), "lipid_class")) same(saved[[field]], expected_broad[[field]],
        paste(context, diet, "all-strict-programme sensitivity", field))
    }
    all_expected_edges[[context]] <- e
    record(paste(context, "independent observed, omission, bootstrap and class replay completed"),
           detail = paste(nrow(e), "fixed pairs;", length(candidates), "LOO candidates;", sum(eligible), "display-eligible"))
  }
  edges <- do.call(rbind, all_expected_edges)
  stopifnot(identical(edge_key(edges), edge_key(a$edges)))
  q <- stats::p.adjust(c(edges$p_SD, edges$p_WD), "BH"); n <- nrow(edges)
  same(edges$q_SD, q[seq_len(n)], "Global selected-universe correlation BH family SD component")
  same(edges$q_WD, q[n + seq_len(n)], "Global selected-universe correlation BH family WD component")

  for (context in context_order) {
    e <- edges[edges$context == context & edges$network_eligible, , drop = FALSE]
    g <- a$network_genes[a$network_genes$context == context, , drop = FALSE]
    actual <- a$network_edges[a$network_edges$context == context, , drop = FALSE]
    stopifnot(nrow(g) == 1L)
    if (!nrow(e)) {
      stopifnot(is.na(g$gene), g$n_edges == 0L, nrow(actual) == 0L)
      next
    }
    genes <- unique(e$gene)
    ranking <- data.frame(gene = genes,
      n = vapply(genes, function(gene) sum(e$gene == gene), integer(1L)),
      max_change = vapply(genes, function(gene) max(abs(e$delta_abs[e$gene == gene])), numeric(1L)))
    pr <- a$programmes[a$programmes$context == context, , drop = FALSE]
    ranking$q <- pr$best_anchor_q[match(ranking$gene, pr$gene)]
    ranking <- ranking[order(-ranking$n, -ranking$max_change, ranking$q, ranking$gene), , drop = FALSE]
    stopifnot(g$gene == ranking$gene[1L])
    selected <- e[e$gene == ranking$gene[1L], , drop = FALSE]
    selected <- selected[order(-abs(selected$delta_abs), selected$lipid), , drop = FALSE]
    first <- c(head(which(selected$direction == "gained"), 3L), head(which(selected$direction == "lost"), 3L))
    take <- head(c(first, setdiff(seq_len(nrow(selected)), first)), 6L)
    selected <- selected[take, , drop = FALSE]
    stopifnot(identical(edge_key(actual), edge_key(selected)), g$n_edges == nrow(selected))
    same(actual$display_rank, seq_len(nrow(selected)), paste(context, "deterministic gained/lost partner display order"))
    record(paste(context, "gene rank and up-to-three-per-direction plus fill partner selection"))
  }
  prior <- utils::read.csv(file.path(dirname(a$source_paths[["parent"]]), "provenance_sha256.csv"), stringsAsFactors = FALSE)
  paths <- unname(parent$input_paths)
  expected <- prior$sha256[match(paths, prior$path)]
  actual <- vapply(paths, digest::digest, character(1L), file = TRUE, algo = "sha256")
  stopifnot(!anyNA(expected), identical(unname(expected), unname(actual)),
            all(paths %in% a$source_paths), !any(grepl("CF459|scaled_RNA|results/zonation/|results/mofa/", a$source_paths)))
  record("All original CF453 input hashes unchanged; no external cohort/network reference dependencies")
  record("Every stored edge-specific bootstrap correlation independently reconstructed",
         detail = paste(boot_replay_count, "candidate edge/draw/diet values; complete class-draw arrays also replayed"))
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- do.call(rbind, checks)
  utils::write.csv(summary, file.path(qa_dir, "independent_checks.csv"), row.names = FALSE)
  saveRDS(list(checks = summary, analysis = normalizePath(analysis_file),
               bootstrap_values_checked = boot_replay_count, session = utils::sessionInfo()),
          file.path(qa_dir, "independent_qa.rds"))
  cat("Independent CF453 selected-programme QA passed:", nrow(summary), "check groups.\n")
  invisible(summary)
}

qa_cf453_rewiring_exports <- function(analysis_file,
    qa_dir = "tests/outputs/cf453_diet_rewiring") {
  if (dir.exists(".local/R-library")) .libPaths(c(".local/R-library", .libPaths()))
  a <- readRDS(analysis_file); output <- dirname(analysis_file)
  checks <- list()
  record <- function(label, difference = 0) {
    if (!is.finite(difference) || difference > 1e-8) stop("Independent export QA failed: ", label, "; difference=", difference)
    checks[[length(checks) + 1L]] <<- data.frame(check = label, max_abs_difference = difference)
  }
  same <- function(actual, expected, label) {
    stopifnot(length(actual) == length(expected), identical(is.na(actual), is.na(expected)))
    if (is.numeric(expected) || is.logical(expected)) {
      good <- is.finite(as.numeric(actual)) & is.finite(as.numeric(expected))
      difference <- if (any(good)) max(abs(as.numeric(actual)[good] - as.numeric(expected)[good])) else 0
      record(label, difference)
    } else {
      stopifnot(identical(as.character(actual), as.character(expected)))
      record(label)
    }
  }
  key <- function(x, fields) do.call(paste, c(x[, fields, drop = FALSE], sep = "|"))
  replay_file <- file.path(output, "QA", "cli_replay", "rewiring_analysis.rds")
  stopifnot(file.exists(replay_file), identical(a, readRDS(replay_file)))
  record("Complete fresh models-only CLI object is exactly identical to production")
  mf <- utils::read.csv(file.path(output, "panel_manifest.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(mf) == 9L, identical(sort(unique(mf$family)), c("A", "B", "C")),
            all(table(mf$family, mf$format) == 1L),
            all(file.exists(mf$file)), all(file.exists(mf$source_table)),
            length(unique(mf$source_table)) == 3L)
  graphics <- list.files(file.path(output, "panels"), pattern = "\\.(pdf|svg|png)$", full.names = TRUE)
  stopifnot(length(graphics) == 9L)
  record("Exactly three independent panel groups and nine alternative graphics exports")
  read_source <- function(family) {
    files <- unique(mf$source_table[mf$family == family]); stopifnot(length(files) == 1L)
    utils::read.csv(files, stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = if (family == "B") c(label = "character") else NA)
  }
  sa <- read_source("A"); akey <- c("gene", "anchor_diet", "sex", "zone")
  stopifnot(nrow(sa) == nrow(a$anchors), !anyDuplicated(key(sa, akey)),
            setequal(key(sa, akey), key(a$anchors, akey)))
  aa <- a$anchors[match(key(sa, akey), key(a$anchors, akey)), , drop = FALSE]
  for (field in names(aa)) same(sa[[field]], aa[[field]], paste("A exact anchor source join", field))
  same(sa$x_start, ifelse(sa$anchor_diet == "SD", sa$x_SD, sa$x_WD), "A arrows start in their selection diet, female coordinate")
  same(sa$y_start, ifelse(sa$anchor_diet == "SD", sa$y_SD, sa$y_WD), "A arrows start in their selection diet, male coordinate")
  same(sa$x_end, ifelse(sa$anchor_diet == "SD", sa$x_WD, sa$x_SD), "A arrows end in the other diet, female coordinate")
  same(sa$y_end, ifelse(sa$anchor_diet == "SD", sa$y_WD, sa$y_SD), "A arrows end in the other diet, male coordinate")
  same(sa$direct_test, ifelse(sa$threeway_q <= .05, "supported", "other"),
       "A line styles use the direct three-way adjusted test")
  sb <- read_source("B"); bkey <- c("context", "lipid_class")
  stopifnot(nrow(sb) == nrow(a$class_summary), !anyDuplicated(key(sb, bkey)),
            setequal(key(sb, bkey), key(a$class_summary, bkey)))
  bb <- a$class_summary[match(key(sb, bkey), key(a$class_summary, bkey)), , drop = FALSE]
  for (field in names(bb)) same(sb[[field]], bb[[field]], paste("B exact fixed-denominator source join", field))
  same(sb$delta_percentage_points, 100 * sb$delta_E, "B displayed units are excess-rate percentage points")
  same(sb$stability_marker, is.finite(sb$bootstrap_direction) & sb$bootstrap_direction >= .90 & sb$bootstrap_valid > 0,
       "B markers mean conditional directional stability, not p values")
  same(sb$label, ifelse(is.finite(sb$delta_percentage_points), sprintf("%+.1f", sb$delta_percentage_points), "NA"),
       "B printed cell values match the displayed percentage-point scale")
  sc <- read_source("C"); ec <- sc[sc$row_type == "selected_edge", , drop = FALSE]
  gc <- sc[sc$row_type == "context_representative", , drop = FALSE]
  ekey <- c("context", "gene", "lipid")
  stopifnot(nrow(ec) == nrow(a$network_edges), !anyDuplicated(key(ec, ekey)),
            setequal(key(ec, ekey), key(a$network_edges, ekey)))
  ee <- a$network_edges[match(key(ec, ekey), key(a$network_edges, ekey)), , drop = FALSE]
  for (field in names(ee)) same(ec[[field]], ee[[field]], paste("C exact selected-edge source join", field))
  same(ec$plot_direction, ec$direction, "C gain/loss edge colours follow the saved strength-change category")
  same(ec$edge_width, abs(ec$delta_abs), "C edge width shows absolute change in correlation strength")
  sd_label <- formatC(ec$r_SD, digits = 2, format = "f")
  wd_label <- formatC(ec$r_WD, digits = 2, format = "f")
  sd_label[sd_label == "-0.00"] <- "0.00"
  wd_label[wd_label == "-0.00"] <- "0.00"
  same(ec$edge_label, paste0(sd_label, " -> ", wd_label),
       "C printed SD-to-WD coefficients exactly match the selected pairs")
  gg <- a$network_genes[match(gc$context, a$network_genes$context), , drop = FALSE]
  stopifnot(nrow(gc) == nrow(a$network_genes), !anyDuplicated(gc$context))
  for (field in names(gg)) same(gc[[field]], gg[[field]], paste("C representative source join", field))
  for (context in unique(ec$context)) {
    rows <- ec$context == context
    stopifnot(sum(rows) <= 6L, !anyDuplicated(ec$slot[rows]), all(ec$slot[rows] %in% 1:6),
              length(unique(ec$gene[rows])) == 1L)
  }
  for (file in mf$file[mf$format == "SVG"]) {
    text <- paste(readLines(file, warn = FALSE), collapse = "\n")
    stopifnot(!grepl("<image\\b", text, perl = TRUE))
    record(paste("Pure-vector SVG", basename(file)))
  }
  for (j in which(mf$format == "PNG")) {
    con <- file(mf$file[j], "rb"); bytes <- readBin(con, "raw", n = 24L); close(con)
    stopifnot(identical(as.integer(bytes[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
    decode <- function(raw) sum(as.integer(raw) * 256^(3:0))
    stopifnot(abs(decode(bytes[17:20]) - mf$width_in[j] * 300) <= 1,
              abs(decode(bytes[21:24]) - mf$height_in[j] * 300) <= 1)
    record(paste("300-dpi PNG dimensions", basename(mf$file[j])))
  }
  provenance <- utils::read.csv(file.path(output, "provenance_sha256.csv"), stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(provenance$path), all(file.exists(provenance$path)))
  hashes <- vapply(provenance$path, digest::digest, character(1L), file = TRUE, algo = "sha256")
  stopifnot(identical(unname(hashes), unname(provenance$sha256)))
  record("Every recorded input, parent, code and plan SHA256 matches the final saved file")
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- do.call(rbind, checks)
  utils::write.csv(summary, file.path(qa_dir, "independent_export_checks.csv"), row.names = FALSE)
  cat("Independent CF453 programme export QA passed:", nrow(summary), "check groups.\n")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || identical(args, "--help")) {
    cat("Usage: Rscript tests/qa_cf453_diet_rewiring.R REWIRING_RDS [QA_DIR]\n")
    cat("       Rscript tests/qa_cf453_diet_rewiring.R --exports-only REWIRING_RDS [QA_DIR]\n")
  } else if (args[1L] == "--exports-only") {
    stopifnot(length(args) >= 2L)
    qa_cf453_rewiring_exports(args[2L], if (length(args) >= 3L) args[3L] else "tests/outputs/cf453_diet_rewiring")
  } else {
    qa_cf453_diet_rewiring(args[1L], if (length(args) >= 2L) args[2L] else "tests/outputs/cf453_diet_rewiring")
  }
}
