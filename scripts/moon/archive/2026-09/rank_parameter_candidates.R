# Summary-only validation/ranking; does not score or change a network.
# A run_id identifies ONE fixed, fully scored/coherence-filtered PKN.
# Run this file directly to execute synthetic tests.

moon_validate_parameter_summary <- function(x) {
  required <- c("run_id", "depth", "compression", "primary", "secondary", "edges",
                "upstream_tfs", "downstream_lipids", "converged", "audit_passed")
  stopifnot(is.data.frame(x), all(required %in% names(x)), !anyNA(x[, required]),
            is.logical(x$converged), is.logical(x$audit_passed),
            all(nzchar(as.character(x$run_id))))
  numeric_names <- c("depth", "primary", "secondary", "edges", "upstream_tfs", "downstream_lipids")
  stopifnot(all(vapply(x[, numeric_names, drop = FALSE], is.numeric, logical(1))),
            all(is.finite(as.matrix(x[, numeric_names, drop = FALSE]))),
            all(x$primary >= x$secondary), all(x$secondary >= 0), all(x$depth > 0),
            all(x$edges >= 0), all(x$upstream_tfs >= 0), all(x$downstream_lipids >= 0))
  count_names <- c("depth", "edges", "upstream_tfs", "downstream_lipids")
  stopifnot(all(as.matrix(x[, count_names]) == floor(as.matrix(x[, count_names]))))
  # A stricter pair of cutoffs can only remove edges/nodes on the same scored PKN.
  # Compare no counts across different depth, compression or coherence runs.
  violations <- list()
  for (idx in split(seq_len(nrow(x)), x$run_id)) {
    stopifnot(length(unique(x$depth[idx])) == 1L,
              length(unique(x$compression[idx])) == 1L)
    for (i in idx) {
      stricter <- idx[x$primary[idx] >= x$primary[i] & x$secondary[idx] >= x$secondary[i] &
                       (x$primary[idx] > x$primary[i] | x$secondary[idx] > x$secondary[i])]
      for (j in stricter) {
        failed <- c("edges", "upstream_tfs", "downstream_lipids")
        failed <- failed[vapply(failed, function(n) x[[n]][j] > x[[n]][i], logical(1))]
        if (length(failed)) violations[[length(violations) + 1L]] <- data.frame(
          run_id = x$run_id[i], lower_row = i, stricter_row = j,
          increased_counts = paste(failed, collapse = ","))
      }
    }
  }
  bad <- if (length(violations)) do.call(rbind, violations) else
    data.frame(run_id = character(), lower_row = integer(), stricter_row = integer(),
               increased_counts = character())
  list(valid = nrow(bad) == 0L, monotonicity_violations = bad)
}

moon_rank_parameter_candidates <- function(x, max_edges = 500L,
                                           min_tfs = 1L, min_lipids = 1L) {
  check <- moon_validate_parameter_summary(x)
  stopifnot(check$valid, length(max_edges) == 1L, is.finite(max_edges), max_edges > 0,
            min_tfs >= 1, min_lipids >= 1)
  x$eligible <- x$converged & x$audit_passed & x$edges > 0 & x$edges < max_edges &
    x$upstream_tfs >= min_tfs & x$downstream_lipids >= min_lipids
  # All eligible results remain available; a one-TF/one-lipid graph is not hidden.
  x$multiple_tfs_and_lipids <- x$upstream_tfs >= 2 & x$downstream_lipids >= 2
  x$threshold_frontier_within_run <- FALSE
  for (idx in split(which(x$eligible), x$run_id[x$eligible])) {
    for (i in idx) {
      dominated <- x$primary[idx] >= x$primary[i] & x$secondary[idx] >= x$secondary[i] &
        (x$primary[idx] > x$primary[i] | x$secondary[idx] > x$secondary[i])
      x$threshold_frontier_within_run[i] <- !any(dominated)
    }
  }
  # Declared preference: maximize the weakest inferred-node cutoff first, then
  # the primary cutoff; retain coverage and shallower explanations as tie-breaks.
  # These are display/search preferences, not statistical significance ranks.
  idx <- order(!x$eligible, -x$secondary, -x$primary, -x$downstream_lipids,
               -x$upstream_tfs, x$depth, x$edges, x$run_id)
  x <- x[idx, , drop = FALSE]
  x$eligible_rank <- NA_integer_
  x$eligible_rank[x$eligible] <- seq_len(sum(x$eligible))
  rownames(x) <- NULL
  x
}

moon_parameter_rank_self_test <- function() {
  x <- data.frame(run_id = "toy", depth = 8L, compression = "safe",
                  primary = c(1, 2, 3, 4, 2, 5), secondary = c(0, 1, 1, 1, 2, 2),
                  edges = c(500L, 200L, 100L, 80L, 150L, 0L),
                  upstream_tfs = c(4L, 3L, 2L, 1L, 2L, 0L),
                  downstream_lipids = c(4L, 3L, 2L, 1L, 2L, 0L),
                  converged = TRUE, audit_passed = TRUE)
  stopifnot(moon_validate_parameter_summary(x)$valid)
  r <- moon_rank_parameter_candidates(x)
  stopifnot(sum(r$eligible) == 4L, r$primary[1] == 2, r$secondary[1] == 2,
            !any(r$eligible[r$edges %in% c(0, 500)]),
            r$threshold_frontier_within_run[r$primary == 4 & r$secondary == 1],
            !r$threshold_frontier_within_run[r$primary == 3 & r$secondary == 1])
  bad <- x; bad$edges[3] <- 250L
  stopifnot(!moon_validate_parameter_summary(bad)$valid)
  failed <- x; failed$converged[5] <- FALSE
  rf <- moon_rank_parameter_candidates(failed)
  stopifnot(!rf$eligible[rf$primary == 2 & rf$secondary == 2])
  message("PASS: strict <500 constraint, empty/failure exclusion, declared ranking, ",
          "within-run frontier, and threshold monotonicity checks.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) moon_parameter_rank_self_test()
