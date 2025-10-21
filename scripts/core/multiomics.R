#' Core Multi-Omics Integration Functions (MOFA)
#'
#' This module provides functions for MOFA (Multi-Omics Factor Analysis) integration.
#'
#' ⚠️ **CRITICAL WARNING** ⚠️
#' The MOFA input file (data/mofa/MOFA_RNA_data.csv) is IMMUTABLE.
#' Pre-trained models depend on this exact input. Any change requires hours of retraining.
#' All functions here are tested for EXACT numerical identity with original code.
#'
#' Workflow:
#' 1. Prepare MOFA input (RNA + lipidomics → long format)
#' 2. Load pre-trained models
#' 3. Analyze Z matrix (sample scores)
#' 4. Analyze W matrix (feature weights)
#' 5. Build correlation networks (TF → Gene → Lipid)

library(MOFA2)
library(reshape2)
library(decoupleR)
library(dplyr)
library(readr)

# Source support functions
source("scripts/support_functions.R")


#' Prepare MOFA Input Data
#'
#' Combines RNA and lipidomic data into MOFA long format.
#'
#' ⚠️ **CRITICAL**: This function MUST produce IDENTICAL output to original code.
#' Use test_mofa_input() to verify before using in production.
#'
#' @param rna_data Matrix or data frame. RNA expression data (genes x samples).
#'                Sample names should be: "MouseID_CellType_GonadType"
#' @param lipidomic_data Matrix or data frame. Lipidomic data (lipids x samples).
#'                       Sample names should be: "MouseID"
#' @param remove_gonad_suffix Logical. Remove "_Testes" and "_Ovaries" from RNA sample names.
#'                            Default: TRUE (REQUIRED for MOFA)
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Data frame in MOFA long format with columns: sample, feature, view, value
#'
#' @details
#' Processing steps (MUST match original exactly):
#' 1. RNA data:
#'    - Extract cell type (view) from sample name (text after last "_")
#'    - Remove gonad suffix from sample names
#'    - Filter to samples present in lipidomic data
#'    - Melt to long format
#' 2. Lipidomic data:
#'    - Add view = "lipidomic"
#'    - Melt to long format
#' 3. Combine:
#'    - Rbind RNA and lipidomic
#'    - Append view to feature names: "FeatureName_ViewName"
#'
#' @examples
#' mofa_input <- prepare_mofa_input(rna_data, lipidomic_data)
#' write_csv(mofa_input, "data/mofa/MOFA_RNA_data.csv")
prepare_mofa_input <- function(rna_data,
                               lipidomic_data,
                               remove_gonad_suffix = TRUE,
                               verbose = TRUE) {

  if (verbose) {
    cat("Preparing MOFA input data...\n")
    cat("  RNA data:", nrow(rna_data), "genes x", ncol(rna_data), "samples\n")
    cat("  Lipidomic data:", nrow(lipidomic_data), "lipids x",
        ncol(lipidomic_data), "samples\n")
  }

  # ============================================================================
  # Step 1: Process RNA data
  # ============================================================================

  # Make a copy to avoid modifying input
  rna_mofa <- as.data.frame(rna_data)

  # Remove gonad suffixes from column names (CRITICAL: must match original)
  if (remove_gonad_suffix) {
    names(rna_mofa) <- gsub("_Testes", "", names(rna_mofa))
    names(rna_mofa) <- gsub("_Ovaries", "", names(rna_mofa))
  }

  # Add feature column
  rna_mofa$feature <- row.names(rna_mofa)

  # Melt to long format
  rna_mofa_long <- melt(rna_mofa)

  # Extract view (cell type) from sample name - text after last "_"
  rna_mofa_long$variable <- as.character(rna_mofa_long$variable)
  rna_mofa_long$view <- gsub(".*_", "", rna_mofa_long$variable)

  # Extract sample base name - text before cell type
  rna_mofa_long$variable <- gsub("_[A-Z].*", "", rna_mofa_long$variable)

  # Rename columns to MOFA format
  names(rna_mofa_long) <- c("feature", "sample", "value", "view")

  # Reorder columns: sample, feature, view, value (CRITICAL: exact order)
  rna_mofa_long <- rna_mofa_long[, c("sample", "feature", "view", "value")]

  # Filter to samples present in lipidomic data (CRITICAL: exact filtering)
  samples_in_lipid <- names(lipidomic_data)
  rna_mofa_long <- rna_mofa_long[rna_mofa_long$sample %in% samples_in_lipid, ]

  if (verbose) {
    cat("  RNA long format:", nrow(rna_mofa_long), "rows\n")
    cat("  RNA views:", paste(unique(rna_mofa_long$view), collapse = ", "), "\n")
  }

  # ============================================================================
  # Step 2: Process lipidomic data
  # ============================================================================

  # Make a copy
  lipid_mofa <- as.data.frame(lipidomic_data)

  # Add feature column
  lipid_mofa$feature <- row.names(lipid_mofa)

  # Melt to long format
  lipid_mofa_long <- melt(lipid_mofa)
  lipid_mofa_long$variable <- as.character(lipid_mofa_long$variable)

  # Add view label
  lipid_mofa_long$view <- "lipidomic"

  # Rename columns
  names(lipid_mofa_long) <- c("feature", "sample", "value", "view")

  if (verbose) {
    cat("  Lipidomic long format:", nrow(lipid_mofa_long), "rows\n")
  }

  # ============================================================================
  # Step 3: Combine RNA and lipidomic
  # ============================================================================

  mofa_data <- as.data.frame(rbind(rna_mofa_long, lipid_mofa_long))

  # Append view to feature names: "FeatureName_ViewName" (CRITICAL)
  mofa_data$feature <- paste(mofa_data$feature, mofa_data$view, sep = "_")

  if (verbose) {
    cat("  Combined MOFA data:", nrow(mofa_data), "rows\n")
    cat("  Samples:", length(unique(mofa_data$sample)), "\n")
    cat("  Features:", length(unique(mofa_data$feature)), "\n")
    cat("  Views:", paste(unique(mofa_data$view), collapse = ", "), "\n")
  }

  return(mofa_data)
}


#' Load Pre-Trained MOFA Models
#'
#' Loads multiple pre-trained MOFA models for comparison.
#'
#' @param model_dir Character. Directory containing model files.
#'                  Default: "results/mofa/RNA"
#' @param factor_range Numeric vector. Factor numbers to load. Default: 5:15
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Named list of MOFA models (model5, model6, ..., model15)
#'
#' @examples
#' models <- load_mofa_models()
#' selected_model <- models$model12
load_mofa_models <- function(model_dir = "results/mofa/RNA",
                             factor_range = 5:15,
                             verbose = TRUE) {

  if (verbose) {
    cat("Loading pre-trained MOFA models...\n")
    cat("  Directory:", model_dir, "\n")
    cat("  Factors:", paste(factor_range, collapse = ", "), "\n")
  }

  models <- list()

  for (i in factor_range) {
    filepath <- file.path(model_dir, paste0("mofa_res_", i, "factor.hdf5"))

    if (!file.exists(filepath)) {
      warning("Model file not found: ", filepath)
      next
    }

    if (verbose) cat("  Loading model", i, "...\n")
    model <- load_model(filepath)
    models[[paste0("model", i)]] <- model
  }

  if (verbose) cat("  ✓ Loaded", length(models), "models\n")

  return(models)
}


#' Enrich MOFA Factors for Experimental Conditions
#'
#' Tests which MOFA factors are enriched for specific experimental conditions
#' (e.g., gonadal type, allosome type).
#'
#' @param z_matrix Matrix. Z matrix from MOFA model (samples x factors).
#' @param sample_info Data frame. Sample information with columns "sample" and condition columns.
#' @param condition_column Character. Column name for condition to test.
#'                         Default: "gonadal_type"
#' @param method Character. decoupleR method. Default: "ulm"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Data frame with enrichment scores (long format: factor, condition, score)
#'
#' @details
#' Uses decoupleR to test if factor scores differ between conditions.
#' Positive score = factor higher in that condition.
#'
#' @examples
#' z_matrix <- get_factors(model)
#' enrichment <- enrich_factors_for_condition(z_matrix, metadata, "gonadal_type")
enrich_factors_for_condition <- function(z_matrix,
                                        sample_info,
                                        condition_column = "gonadal_type",
                                        method = "ulm",
                                        verbose = TRUE) {

  if (verbose) {
    cat("Enriching MOFA factors for:", condition_column, "\n")
  }

  # Prepare network (sample → condition mapping)
  network <- sample_info[, c("sample", condition_column)]
  names(network) <- c("target", "source")
  network$mor <- 1
  network <- unique(network)

  # Run enrichment
  enrichment <- decoupleR::run_ulm(z_matrix, network)

  # Convert to wide format
  enrichment_wide <- dcast(enrichment, source ~ condition, value.var = "score")

  if (verbose) {
    cat("  ✓ Enrichment complete\n")
    cat("  Conditions tested:", paste(unique(enrichment$source), collapse = ", "), "\n")
  }

  return(list(
    long = enrichment,
    wide = enrichment_wide
  ))
}


#' Extract and Filter MOFA Weights
#'
#' Extracts feature weights from MOFA model and filters to top features.
#'
#' @param model MOFA model object
#' @param top_percent Numeric. Percentage of top features to keep (each tail).
#'                    Default: 0.025 (top 2.5% positive + top 2.5% negative = 5% total)
#' @param views Character vector. Views to extract. Default: "all"
#' @param factors Character vector or numeric. Factors to extract. Default: "all"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return List of weight matrices (one per view), filtered to top features
#'
#' @details
#' Filtering keeps top X% of features by absolute weight (both positive and negative tails).
#' This focuses downstream analysis on features that most strongly drive each factor.
#'
#' @examples
#' weights <- extract_mofa_weights(model, top_percent = 0.025)
#' hepatocyte_weights <- weights$Hepatocytes
extract_mofa_weights <- function(model,
                                 top_percent = 0.025,
                                 views = "all",
                                 factors = "all",
                                 verbose = TRUE) {

  if (verbose) {
    cat("Extracting MOFA weights...\n")
    cat("  Top percent (each tail):", top_percent * 100, "%\n")
  }

  # Get weights
  weights <- get_weights(model, views = views, factors = factors)

  # Filter each view's weights to top features
  weights_filtered <- lapply(weights, function(factor_weights) {
    apply(factor_weights, 2, function(x) {
      # Sort by value (not absolute)
      x_sorted <- x[order(x, decreasing = TRUE)]

      # Calculate number to keep from each tail
      top_n <- round(length(x) * top_percent)

      # Set middle values to NA
      x_sorted[top_n:(length(x) - top_n)] <- NA

      # Reorder back to original
      x[names(x_sorted)] <- x_sorted

      return(x)
    }, simplify = FALSE)
  })

  if (verbose) {
    for (view_name in names(weights_filtered)) {
      n_features <- sum(!is.na(weights_filtered[[view_name]][[1]]))
      cat("  ", view_name, ":", n_features, "features kept\n")
    }
  }

  return(weights_filtered)
}


#' Run TF Enrichment on MOFA Weights
#'
#' Infers which transcription factors are enriched in MOFA factor weights.
#'
#' @param weights Matrix. MOFA weights for one view (e.g., Hepatocytes).
#' @param regulon_file Character. Path to TF regulons RData. Default: "support/TF_regulons.RData"
#' @param hormone_gene_file Character. Path to hormone target genes CSV.
#'                          Default: "support/Gonadal_Hormone_Genes.csv"
#' @param remove_view_suffix Logical. Remove view suffix from feature names. Default: TRUE
#' @param minsize Integer. Minimum regulon size. Default: 20
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Matrix of TF scores (TFs x factors)
#'
#' @details
#' Combines CollecTRI regulons with gonadal hormone target genes (ESR1, AR)
#' to comprehensively test TF activity.
#'
#' @examples
#' hepatocyte_weights <- weights$Hepatocytes
#' tf_scores <- run_tf_enrichment_on_weights(hepatocyte_weights)
run_tf_enrichment_on_weights <- function(weights,
                                         regulon_file = "support/TF_regulons.RData",
                                         hormone_gene_file = "support/Gonadal_Hormone_Genes.csv",
                                         remove_view_suffix = TRUE,
                                         minsize = 20,
                                         verbose = TRUE) {

  if (verbose) {
    cat("Running TF enrichment on MOFA weights...\n")
  }

  # Load TF regulons
  load(regulon_file)  # Loads TF_regulons
  TF_regulons <- TF_regulons[complete.cases(TF_regulons), ]

  # Load hormone target genes
  hormone_genes <- as.data.frame(read_csv(hormone_gene_file, show_col_types = FALSE))
  hormone_genes$source <- gsub("_spec", "", hormone_genes$source)
  hormone_genes$source <- gsub("ESR1", "Esr1", hormone_genes$source)
  hormone_genes$source <- gsub("AR", "Ar", hormone_genes$source)

  # Remove Ar and Esr1 from CollecTRI (will use hormone-specific targets instead)
  TF_regulons <- TF_regulons[!TF_regulons$source %in% c("Ar", "Esr1"), ]

  # Combine regulons
  TF_regulons <- rbind(TF_regulons, hormone_genes)

  if (verbose) {
    cat("  TF regulons:", length(unique(TF_regulons$source)), "TFs\n")
  }

  # Prepare weight matrix
  weight_matrix <- as.matrix(weights)

  # Remove view suffix from row names if requested
  if (remove_view_suffix) {
    rownames(weight_matrix) <- gsub("_.*", "", rownames(weight_matrix))
  }

  # Run TF enrichment
  tf_scores <- decoupleR::run_ulm(
    mat = weight_matrix,
    network = TF_regulons,
    minsize = minsize
  )

  # Convert to wide format
  tf_scores_wide <- dcast(tf_scores, source ~ condition, value.var = "score")
  rownames(tf_scores_wide) <- tf_scores_wide$source
  tf_scores_wide <- tf_scores_wide[, -1]

  if (verbose) {
    cat("  ✓ TF enrichment complete:", nrow(tf_scores_wide), "TFs\n")
  }

  return(as.matrix(tf_scores_wide))
}


#' Save MOFA Weights with TF Scores
#'
#' Combines MOFA weights with TF enrichment scores and saves.
#'
#' @param weights List. MOFA weights (from extract_mofa_weights())
#' @param tf_scores Matrix. TF scores (from run_tf_enrichment_on_weights())
#' @param tf_view_name Character. Name for TF view. Default: "Hepatocytes_TFs"
#' @param output_file Character. Output file path. Default: "results/mofa/main_model_RNA_lipid_weight.RData"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @examples
#' save_mofa_weights_with_tfs(weights, tf_scores)
save_mofa_weights_with_tfs <- function(weights,
                                       tf_scores,
                                       tf_view_name = "Hepatocytes_TFs",
                                       output_file = "results/mofa/main_model_RNA_lipid_weight.RData",
                                       verbose = TRUE) {

  # Add TF view to weights
  # Add TF suffix to row names
  rownames(tf_scores) <- paste0(rownames(tf_scores), "_", tf_view_name)

  # Add to weights list
  weights[[tf_view_name]] <- tf_scores

  # Save
  save(weights, file = output_file)

  if (verbose) {
    cat("✓ Saved weights with TF scores to:", output_file, "\n")
    cat("  Views:", paste(names(weights), collapse = ", "), "\n")
  }
}


#' Compute Gene-Lipid Correlation Matrix
#'
#' Computes correlations between gene expression and lipid abundances
#' for features selected by MOFA.
#'
#' @param rna_data Matrix. RNA expression data (genes x samples)
#' @param lipid_data Matrix. Lipidomic data (lipids x samples)
#' @param mofa_features Character vector. Feature names from MOFA (with view suffix)
#' @param cell_type Character. Cell type to analyze. Default: "Hepatocytes"
#' @param filter_to_lipid_genes Logical. Filter genes to lipid metabolism genes only.
#'                              Default: TRUE
#' @param pathway_file Character. GMT file for lipid metabolism genes.
#'                     Default: "support/m5.all.v2023.2.Mm.symbols.gmt"
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Matrix of correlations (lipids x genes)
#'
#' @examples
#' cor_matrix <- compute_gene_lipid_correlations(rna_data, lipid_data, mofa_features)
compute_gene_lipid_correlations <- function(rna_data,
                                           lipid_data,
                                           mofa_features,
                                           cell_type = "Hepatocytes",
                                           filter_to_lipid_genes = TRUE,
                                           pathway_file = "support/m5.all.v2023.2.Mm.symbols.gmt",
                                           verbose = TRUE) {

  if (verbose) {
    cat("Computing gene-lipid correlations...\n")
    cat("  Cell type:", cell_type, "\n")
  }

  # Extract gene and lipid features from MOFA
  cell_suffix <- paste0("_", cell_type)
  genes_in_mofa <- mofa_features[grepl(cell_suffix, mofa_features)]
  genes_in_mofa <- gsub(cell_suffix, "", genes_in_mofa)

  lipids_in_mofa <- mofa_features[grepl("_lipidomic", mofa_features)]
  lipids_in_mofa <- gsub("_lipidomic", "", lipids_in_mofa)

  # Subset RNA to cell type and MOFA genes
  rna_celltype <- rna_data[, grepl(cell_suffix, colnames(rna_data))]
  colnames(rna_celltype) <- gsub(cell_suffix, "", colnames(rna_celltype))
  rna_celltype <- rna_celltype[rownames(rna_celltype) %in% genes_in_mofa, ]

  # Subset lipidomic to MOFA lipids
  lipid_subset <- lipid_data[rownames(lipid_data) %in% lipids_in_mofa, ]

  # Match samples
  common_samples <- intersect(colnames(rna_celltype), colnames(lipid_subset))
  rna_celltype <- rna_celltype[, common_samples]
  lipid_subset <- lipid_subset[, common_samples]

  if (verbose) {
    cat("  Genes:", nrow(rna_celltype), "\n")
    cat("  Lipids:", nrow(lipid_subset), "\n")
    cat("  Samples:", ncol(rna_celltype), "\n")
  }

  # Compute correlations
  cor_matrix <- as.data.frame(cor(t(rbind(lipid_subset, rna_celltype))))
  cor_matrix <- cor_matrix[rownames(lipid_subset), rownames(rna_celltype)]

  # Filter to lipid metabolism genes if requested
  if (filter_to_lipid_genes) {
    if (verbose) cat("  Filtering to lipid metabolism genes...\n")

    pathways <- import_gmt(pathway_file)
    names(pathways) <- c("target", "source")
    lipid_genes <- pathways[pathways$source == "GOBP_CELLULAR_LIPID_METABOLIC_PROCESS", "target"]

    cor_matrix <- cor_matrix[, colnames(cor_matrix) %in% lipid_genes]

    if (verbose) cat("  Filtered to", ncol(cor_matrix), "lipid metabolism genes\n")
  }

  if (verbose) cat("  ✓ Correlation matrix complete\n")

  return(cor_matrix)
}
