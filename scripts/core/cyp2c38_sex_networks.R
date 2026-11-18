# Full Cyp2c38-to-lipid FCG reference networks for Workflow 07c.
# Source scripts/support_functions.R and cyp2c38_reference.R before calling.
# Sourcing this helper defines functions without generating result files.

#' Regenerate pooled and gonad-specific Cyp2c38 lipid associations.
#'
#' All 1,254 measured original lipid species enter the regenerated networks.
#' Scaled hepatocyte RNA and log2 lipid intensities retain the original 43-mouse
#' pairing. No MOFA selection, KO outcome, correlation cutoff or p-value threshold
#' changes this complete feature universe. Ovaries and Testes denote FCG gonadal
#' groups rather than sex-chromosome groups.
#'
#' @param root Existing analysis repository root.
#' @return A named list with gene, lipids, metadata, coefficients, a named logical
#'   fixed_reference mask, and historical-subset numerical replay diagnostics.
wp1_full_reference <- function(root = ".") {
  if (!exists("wp1_load_reference", mode = "function")) {
    stop("Source scripts/core/cyp2c38_reference.R before regenerating full networks.")
  }
  if (length(root) != 1L || is.na(root) || !dir.exists(root)) {
    stop("root must identify one existing analysis repository directory.")
  }
  root <- normalizePath(root, mustWork = TRUE)
  fixed <- wp1_load_reference(root)
  lipid_file <- file.path(root, "data/metabolomics/Lipid_data_file_JPC_cleaned.csv")
  original <- utils::read.csv(lipid_file, check.names = FALSE,
                             stringsAsFactors = FALSE)
  if (!identical(names(original)[1L], "Unnamed: 0") ||
      anyNA(names(original)) || any(!nzchar(names(original))) ||
      anyDuplicated(names(original)) || ncol(original) - 1L != 1254L) {
    stop("Expected 1,254 unique original lipid columns and a first Unnamed: 0 mouse-ID column.")
  }
  mouse_ids <- gsub("-", "_", gsub("-0", "-", original[["Unnamed: 0"]],
                                  fixed = TRUE))
  if (anyNA(mouse_ids) || any(!nzchar(mouse_ids)) || anyDuplicated(mouse_ids)) {
    stop("Original lipid mouse identifiers must remain unique after normalization.")
  }
  paired_index <- match(names(fixed$gene), mouse_ids)
  if (anyNA(paired_index)) stop("Every FCG RNA mouse requires its original paired lipid sample.")
  abundance <- as.matrix(original[paired_index, -1L, drop = FALSE])
  if (!is.numeric(abundance) || any(!is.finite(abundance)) || any(abundance <= 0)) {
    stop("All original paired lipid intensities must be finite positive numeric values.")
  }
  lipids <- t(log2(abundance))
  dimnames(lipids) <- list(names(original)[-1L], names(fixed$gene))
  metadata <- fixed$metadata
  gene <- fixed$gene
  if (!identical(metadata$sample_id, colnames(lipids)) ||
      !identical(names(gene), colnames(lipids))) {
    stop("Full-network gene, lipid and metadata mouse order must agree exactly.")
  }

  strata <- list(pooled = rep(TRUE, length(gene)),
                 female = metadata$gonadal_type == "Ovaries",
                 male = metadata$gonadal_type == "Testes")
  coefficients <- data.frame(reference_lipid = rownames(lipids),
                             stringsAsFactors = FALSE)
  for (stratum in names(strata)) {
    keep <- strata[[stratum]]
    # Complete paired data: no feature-specific exclusions or imputation.
    if (sum(keep) < 3L || !is.finite(stats::var(gene[keep])) ||
        stats::var(gene[keep]) <= 0) {
      stop("Cyp2c38 expression must vary within the ", stratum, " group.")
    }
    lipid_variance <- apply(lipids[, keep, drop = FALSE], 1L, stats::var)
    if (any(!is.finite(lipid_variance)) || any(lipid_variance <= 0)) {
      stop("Undefined Pearson association for constant lipid feature(s) in ", stratum, ".")
    }
    correlation <- as.numeric(stats::cor(t(lipids[, keep, drop = FALSE]),
                                         gene[keep], method = "pearson"))
    if (length(correlation) != nrow(lipids) || any(!is.finite(correlation))) {
      stop("Full-network correlations must be finite for every measured species in ", stratum, ".")
    }
    coefficients[[paste0("r_", stratum)]] <- correlation
    coefficients[[paste0("n_", stratum)]] <- sum(keep)
  }
  fixed_mask <- coefficients$reference_lipid %in% fixed$coefficients$reference_lipid
  if (sum(fixed_mask) != 828L) stop("The full network must contain all 828 archived reference species.")
  names(fixed_mask) <- coefficients$reference_lipid
  coefficients$fixed_reference <- unname(fixed_mask)
  fixed_index <- match(fixed$coefficients$reference_lipid, coefficients$reference_lipid)
  replay_differences <- vapply(names(strata), function(stratum) {
    column <- paste0("r_", stratum)
    max(abs(coefficients[[column]][fixed_index] - fixed$coefficients[[column]]))
  }, numeric(1L))
  if (any(!is.finite(replay_differences)) || any(replay_differences > 1e-12)) {
    stop("The regenerated full-network subset must reproduce all prior pooled and sex-specific coefficients.")
  }
  archive_difference <- max(abs(coefficients$r_pooled[fixed_index] -
                                  fixed$coefficients$r_archived))
  if (!is.finite(archive_difference) || archive_difference > 1e-12) {
    stop("Full-network pooled coefficients fail the archived Cyp2c38 replay check.")
  }
  cat("Full Cyp2c38 networks: ", nrow(lipids), " measured species x ",
      ncol(lipids), " paired FCG mice; historical subset maximum difference = ",
      format(max(replay_differences), digits = 4L), "\n", sep = "")
  list(gene = gene, lipids = lipids, metadata = metadata,
       coefficients = coefficients, fixed_reference = fixed_mask,
       fixed_subset_max_difference = replay_differences,
       max_archive_difference = archive_difference)
}

#' Export complete and fixed-historical Cyp2c38 reference networks.
#'
#' Every SIF row retains one specific molecular annotation. There is no collapse
#' of chain order, slash/underscore connectivity or sum/species resolution. mor is
#' the original Pearson correlation coefficient, suitable for a signed ULM input.
#' The first three columns are node1, node2 and mor; remaining columns record the
#' contributing mouse count, FCG gonadal group and historical-feature membership.
#'
#' @param out Workflow output root; networks are written to its networks folder.
#' @param reference Result of wp1_full_reference().
#' @return A manifest of six CSV files with network stratum, universe and edge count.
wp1_export_sex_networks <- function(out, reference) {
  if (length(out) != 1L || is.na(out) || !nzchar(out)) {
    stop("out must identify one output directory.")
  }
  required <- c("reference_lipid", "r_pooled", "r_female", "r_male",
                "n_pooled", "n_female", "n_male", "fixed_reference")
  if (!is.list(reference) || !is.data.frame(reference$coefficients) ||
      !all(required %in% names(reference$coefficients))) {
    stop("reference must come from wp1_full_reference().")
  }
  coefficients <- reference$coefficients
  if (nrow(coefficients) != 1254L ||
      anyDuplicated(coefficients$reference_lipid) ||
      anyNA(coefficients$fixed_reference) ||
      !is.logical(coefficients$fixed_reference) ||
      sum(coefficients$fixed_reference) != 828L) {
    stop("Expected unique full-network species and an intact 828-feature historical mask.")
  }
  network_dir <- file.path(out, "networks")
  if (!dir.exists(network_dir) && !dir.create(network_dir, recursive = TRUE)) {
    stop("Unable to create network output directory: ", network_dir)
  }
  strata <- data.frame(stratum = c("pooled", "female", "male"),
                       label = c("pooled", "ovaries", "testes"),
                       gonadal_type = c("Ovaries_and_Testes", "Ovaries", "Testes"),
                       stringsAsFactors = FALSE)
  manifest <- list()
  for (i in seq_len(nrow(strata))) {
    stratum <- strata$stratum[i]
    network <- data.frame(node1 = "Cyp2c38",
                          node2 = coefficients$reference_lipid,
                          mor = coefficients[[paste0("r_", stratum)]],
                          n_mice = coefficients[[paste0("n_", stratum)]],
                          gonadal_type = strata$gonadal_type[i],
                          fixed_reference = coefficients$fixed_reference,
                          stringsAsFactors = FALSE)
    if (any(!is.finite(network$mor)) || any(network$n_mice < 3L)) {
      stop("Only finite Pearson correlations with known mouse counts can be exported.")
    }
    for (universe in c("full", "fixed828")) {
      selected <- if (universe == "full") rep(TRUE, nrow(network)) else network$fixed_reference
      path <- file.path(network_dir,
                        paste0("Cyp2c38_", strata$label[i], "_", universe, "_SIF.csv"))
      utils::write.csv(network[selected, , drop = FALSE], path, row.names = FALSE)
      manifest[[length(manifest) + 1L]] <- data.frame(
        stratum = stratum, gonadal_type = strata$gonadal_type[i],
        universe = universe, n_edges = sum(selected),
        n_mice = unique(network$n_mice[selected]),
        file = normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, manifest)
}
