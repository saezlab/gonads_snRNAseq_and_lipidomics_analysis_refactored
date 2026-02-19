#' Lipid Structural Analysis Module
#'
#' This module provides functions for analyzing how lipid structural features
#' (chain length, double bond position) predict gene-lipid correlation patterns.
#'
#' This analysis reveals substrate specificity - which genes/enzymes preferentially
#' associate with lipids having specific structural characteristics.
#'
#' Key workflow:
#' 1. Compute gene-lipid correlations from MOFA-selected features
#' 2. Load lipid structural data (chain length C, double bond position DB)
#' 3. Run regression models testing structural predictors
#' 4. Identify significant associations (p < 0.05)
#' 5. Summarize genes with structural preferences
#'
#' Regression model:
#'   gene_correlation ~ lipid_class + lipid_abundance + chain_length +
#'                      db_position + lipid_abundance:chain_length +
#'                      lipid_abundance:db_position + chain_length:db_position
#'
#' Consolidates code from:
#' - lipid_rna_correaltion_for_structural_analysis.R
#' - lipid_structure_analysis.R
#' - lipid_structure_visualisation.R (data prep)

library(dplyr)
library(tidyr)
library(readr)


#' ==============================================================================
#' LIPID NAME PARSING FUNCTIONS
#' ==============================================================================
#'
#' Functions to extract structural features (chain length, double bonds) from
#' lipid names following standard nomenclature.
#'
#' Examples of lipid name formats:
#'   - TG(16:0_18:1_22:6) -> C = 56, DB = 7
#'   - PC(16:0/18:1) -> C = 34, DB = 1
#'   - Cer(18:1;O2/18:0) -> C = 36, DB = 1
#'   - LPC(18:0) -> C = 18, DB = 0


#' Parse Lipid Name to Extract Chain Length and Double Bonds
#'
#' Extracts total carbon chain length (C) and total double bonds (DB) from
#' lipid names following standard nomenclature (e.g., "PC(16:0/18:1)").
#'
#' @param lipid_name Character. A lipid name containing fatty acid specifications
#'   in the format "X:Y" where X is carbon count and Y is double bond count.
#'
#' @return Named list with elements:
#'   - chain_length: Total carbon atoms across all fatty acid chains
#'   - db_position: Total double bonds across all fatty acid chains
#'   - lipid_class: Extracted lipid class (e.g., "PC", "TG", "Cer")
#'   Returns NA values if parsing fails.
#'
#' @examples
#' parse_lipid_name("TG(16:0_18:1_22:6)")
#' # Returns: list(chain_length = 56, db_position = 7, lipid_class = "TG")
#'
#' parse_lipid_name("PC(16:0/18:1)")
#' # Returns: list(chain_length = 34, db_position = 1, lipid_class = "PC")
#'
#' @export
parse_lipid_name <- function(lipid_name) {

  # Default return values
  result <- list(
    chain_length = NA_real_,
    db_position = NA_real_,
    lipid_class = NA_character_
  )

  if (is.na(lipid_name) || lipid_name == "") {
    return(result)
  }

  # Extract lipid class (everything before the first parenthesis)

  # This matches the lipid class extraction in run_lipid_class_enrichment()
  # from scripts/core/enrichment_analysis.R for consistency across workflows
  lipid_class_extracted <- gsub("\\(.*", "", lipid_name)
  if (nchar(lipid_class_extracted) > 0) {
    result$lipid_class <- lipid_class_extracted
  }

  # Find all fatty acid specifications in format "X:Y"
  # Pattern: one or more digits, colon, one or more digits
  fa_pattern <- "([0-9]+):([0-9]+)"
  fa_matches <- gregexpr(fa_pattern, lipid_name, perl = TRUE)
  matched_strings <- regmatches(lipid_name, fa_matches)[[1]]

  if (length(matched_strings) == 0) {
    return(result)
  }

  # Parse each fatty acid specification
  total_c <- 0
  total_db <- 0

  for (fa in matched_strings) {
    parts <- strsplit(fa, ":")[[1]]
    if (length(parts) == 2) {
      c_val <- as.numeric(parts[1])
      db_val <- as.numeric(parts[2])
      if (!is.na(c_val)) total_c <- total_c + c_val
      if (!is.na(db_val)) total_db <- total_db + db_val
    }
  }

  result$chain_length <- total_c
  result$db_position <- total_db

  return(result)
}


#' Parse Multiple Lipid Names to Data Frame
#'
#' Applies parse_lipid_name() to a vector of lipid names and returns a data frame.
#'
#' @param lipid_names Character vector. Lipid names to parse.
#' @param verbose Logical. Print progress messages?
#'
#' @return Data frame with columns:
#'   - lipid_name: Original lipid name
#'   - lipid_class: Extracted lipid class
#'   - chain_length: Total carbon atoms
#'   - db_position: Total double bonds
#'
#' @examples
#' lipids <- c("TG(16:0_18:1_22:6)", "PC(16:0/18:1)", "Cer(18:1;O2/18:0)")
#' parse_lipid_names_to_df(lipids)
#'
#' @export
parse_lipid_names_to_df <- function(lipid_names, verbose = TRUE) {

  if (verbose) {
    cat("Parsing", length(lipid_names), "lipid names...\n")
  }

  # Parse each lipid name
  parsed_list <- lapply(lipid_names, parse_lipid_name)

  # Convert to data frame
  result_df <- data.frame(
    lipid_name = lipid_names,
    lipid_class = sapply(parsed_list, function(x) x$lipid_class),
    chain_length = sapply(parsed_list, function(x) x$chain_length),
    db_position = sapply(parsed_list, function(x) x$db_position),
    stringsAsFactors = FALSE
  )

  # Report parsing success
  n_parsed <- sum(!is.na(result_df$chain_length))
  n_failed <- sum(is.na(result_df$chain_length))

  if (verbose) {
    cat("  Successfully parsed:", n_parsed, "lipids\n")
    if (n_failed > 0) {
      cat("  Failed to parse:", n_failed, "lipids\n")
      failed_names <- head(lipid_names[is.na(result_df$chain_length)], 5)
      cat("  Examples of failed names:", paste(failed_names, collapse = ", "), "\n")
    }
  }

  return(result_df)
}


#' Generate Lipid Structural Data from Limma Results
#'
#' Extracts structural features (C, DB) from lipid names in limma results
#' and combines them with fold change values. This generates the full dataset
#' needed for structural analysis (not filtered by |logFC| threshold).
#'
#' @param ttop Data frame. Limma topTable results with columns:
#'   - ID or rownames: Lipid names
#'   - logFC: Log fold change values
#' @param verbose Logical. Print progress messages?
#'
#' @return Data frame with columns:
#'   - lipid_name: Lipid identifier
#'   - lipid_class: Extracted lipid class
#'   - chain_length: Total carbon atoms
#'   - db_position: Total double bonds
#'   - logfc: Log fold change from limma
#'
#' @examples
#' load("results/lipidomic/ttop_full_list.RData")
#' ttop <- ttop_full_list[["Testes vs Ovaries"]][["Testes vs Ovaries"]]
#' lipid_structure <- generate_lipid_structural_data(ttop)
#'
#' @export
generate_lipid_structural_data <- function(ttop, verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Generating Lipid Structural Data from Limma Results\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  # Extract lipid names
  if ("ID" %in% colnames(ttop)) {
    lipid_names <- ttop$ID
  } else {
    lipid_names <- rownames(ttop)
  }

  if (verbose) {
    cat("Total lipids in limma results:", length(lipid_names), "\n")
  }

  # Parse lipid names
  structural_df <- parse_lipid_names_to_df(lipid_names, verbose = verbose)

  # Add logFC values
  structural_df$logfc <- ttop$logFC

  # Remove lipids that failed parsing
  n_before <- nrow(structural_df)
  structural_df <- structural_df[!is.na(structural_df$chain_length), ]
  n_after <- nrow(structural_df)

  if (verbose) {
    cat("\nLipids with valid structural data:", n_after, "/", n_before, "\n")
    cat("LogFC range:", round(min(structural_df$logfc), 3), "to",
        round(max(structural_df$logfc), 3), "\n")
    cat("Lipid classes:", length(unique(structural_df$lipid_class)), "\n")
    cat("  Classes:", paste(head(sort(unique(structural_df$lipid_class)), 10), collapse = ", "))
    if (length(unique(structural_df$lipid_class)) > 10) cat(", ...")
    cat("\n\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  return(structural_df)
}


#' ==============================================================================
#' LIPID-ONLY STRUCTURAL ANALYSIS FUNCTIONS
#' ==============================================================================
#'
#' These functions analyze how lipid structural features (chain length, double
#' bond position) predict lipid differential abundance (logFC from gonad comparison)
#' WITHOUT involving gene expression data.
#'
#' This is a simpler analysis that characterizes lipid properties before the
#' more complex gene-lipid correlation analysis.
#'
#' PRIMARY APPROACH: Class-specific models
#' For each lipid class: logfc ~ chain_length + db_position + chain_length:db_position
#'
#' This approach allows each lipid class to have its own structural effects,
#' reflecting that different lipid classes may have distinct structure-abundance
#' relationships.


#' Run Lipid-Only Class-Specific Structural Regression
#'
#' Tests how structural features predict differential abundance within each lipid class.
#'
#' Class-specific regression model:
#'   logfc ~ chain_length + db_position + chain_length:db_position
#'
#' @param lipid_structure Data frame. Lipid structural features.
#' @param min_lipids_per_class Numeric. Minimum lipids per class (default: 10).
#' @param verbose Logical. Print progress messages?
#'
#' @return Data frame with columns:
#'   - lipid_class: Lipid class analyzed
#'   - term: Regression term
#'   - estimate: Coefficient
#'   - std_error: Standard error
#'   - t_value: t-statistic
#'   - p_value: P-value
#'   - significant: Logical, TRUE if p < 0.05
#'   - r_squared: Model R-squared for that class
#'   - n_lipids: Number of lipids in the class
#'
#' @examples
#' class_results <- run_lipid_only_class_regression(lipid_structure)
#'
#' @export
run_lipid_only_class_regression <- function(lipid_structure,
                                            min_lipids_per_class = 10,
                                            verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Running Lipid-Only Class-Specific Structural Regression\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
    cat("Model: logfc ~ chain_length + db_position + chain_length:db_position\n\n")
  }

  # Get valid lipid classes
  class_counts <- table(lipid_structure$lipid_class)
  valid_classes <- names(class_counts)[class_counts >= min_lipids_per_class]

  if (verbose) {
    cat("Valid lipid classes (n >=", min_lipids_per_class, "):\n")
    for (cls in valid_classes) {
      cat("  ", cls, ":", class_counts[cls], "lipids\n")
    }
    cat("\n")
  }

  # Run analysis for each class
  all_class_results <- list()
  class_summaries <- list()

  for (cls in valid_classes) {
    if (verbose) cat("Processing", cls, "...\n")

    # Filter to this class
    class_data <- lipid_structure %>%
      dplyr::filter(lipid_class == cls) %>%
      dplyr::filter(!is.na(chain_length) & !is.na(db_position) & !is.na(logfc))

    n_lipids <- nrow(class_data)

    if (n_lipids < 5) {
      if (verbose) cat("  Skipping - insufficient data after NA removal\n")
      next
    }

    # Check for variance in predictors
    if (sd(class_data$chain_length) == 0 || sd(class_data$db_position) == 0) {
      if (verbose) cat("  Skipping - no variance in structural features\n")
      next
    }

    # Run regression
    tryCatch({
      model <- lm(logfc ~ chain_length + db_position + chain_length:db_position,
                  data = class_data)

      model_summary <- summary(model)

      # Extract coefficients
      coef_df <- as.data.frame(model_summary$coefficients)
      coef_df$term <- rownames(coef_df)

      terms_of_interest <- c("chain_length", "db_position", "chain_length:db_position")

      coef_df <- coef_df %>%
        dplyr::filter(term %in% terms_of_interest) %>%
        dplyr::mutate(
          lipid_class = cls,
          estimate = Estimate,
          std_error = `Std. Error`,
          t_value = `t value`,
          p_value = `Pr(>|t|)`,
          significant = p_value < 0.05,
          r_squared = model_summary$r.squared,
          n_lipids = n_lipids
        ) %>%
        dplyr::select(lipid_class, term, estimate, std_error, t_value, p_value,
                      significant, r_squared, n_lipids)

      all_class_results[[cls]] <- coef_df

      # Store class summary
      class_summaries[[cls]] <- data.frame(
        lipid_class = cls,
        n_lipids = n_lipids,
        r_squared = model_summary$r.squared,
        adj_r_squared = model_summary$adj.r.squared,
        model_p_value = pf(model_summary$fstatistic[1],
                           model_summary$fstatistic[2],
                           model_summary$fstatistic[3],
                           lower.tail = FALSE)
      )

    }, error = function(e) {
      if (verbose) cat("  Error fitting model:", conditionMessage(e), "\n")
    })
  }

  # Combine results
  combined_results <- dplyr::bind_rows(all_class_results)
  class_summary_df <- dplyr::bind_rows(class_summaries)

  if (verbose) {
    cat("\n")
    cat("Class-specific results summary:\n")
    cat("  Total class-term combinations:", nrow(combined_results), "\n")
    cat("  Significant associations (p < 0.05):",
        sum(combined_results$significant, na.rm = TRUE), "\n\n")

    cat("Classes with significant structural effects:\n")
    sig_summary <- combined_results %>%
      dplyr::filter(significant) %>%
      dplyr::arrange(lipid_class, term)

    if (nrow(sig_summary) > 0) {
      print(sig_summary %>% dplyr::select(lipid_class, term, estimate, p_value))
    } else {
      cat("  None\n")
    }

    cat("\nModel fit by class:\n")
    print(class_summary_df %>%
            dplyr::arrange(dplyr::desc(r_squared)) %>%
            dplyr::mutate(r_squared = round(r_squared, 3),
                          adj_r_squared = round(adj_r_squared, 3),
                          model_p_value = format.pval(model_p_value)))

    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  return(list(
    coefficients = combined_results,
    class_summaries = class_summary_df
  ))
}


#' Extract Significant Lipid Structural Associations
#'
#' Summarizes which lipid classes show significant structural effects.
#'
#' @param class_results Output from run_lipid_only_class_regression().
#' @param p_cutoff Numeric. P-value cutoff (default: 0.05).
#'
#' @return List with elements:
#'   - classes_with_chain_effect: Classes where chain length predicts logFC
#'   - classes_with_db_effect: Classes where DB position predicts logFC
#'   - classes_with_interaction: Classes with chain_length:db_position interaction
#'   - summary_table: Summary of effects by class
#'
#' @export
extract_lipid_structural_associations <- function(class_results,
                                                  p_cutoff = 0.05) {

  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Extracting Significant Lipid Structural Associations\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  coef_df <- class_results$coefficients

  if (nrow(coef_df) == 0) {
    cat("WARNING: No results to extract.\n\n")
    return(list(
      classes_with_chain_effect = character(),
      classes_with_db_effect = character(),
      classes_with_interaction = character(),
      summary_table = data.frame()
    ))
  }

  # Filter to significant
  sig_results <- coef_df %>%
    dplyr::filter(p_value < p_cutoff)

  # Classes with chain length effect
  chain_classes <- sig_results %>%
    dplyr::filter(term == "chain_length") %>%
    dplyr::select(lipid_class, estimate, p_value) %>%
    dplyr::arrange(p_value)

  # Classes with DB position effect
  db_classes <- sig_results %>%
    dplyr::filter(term == "db_position") %>%
    dplyr::select(lipid_class, estimate, p_value) %>%
    dplyr::arrange(p_value)

  # Classes with interaction effect
  interaction_classes <- sig_results %>%
    dplyr::filter(term == "chain_length:db_position") %>%
    dplyr::select(lipid_class, estimate, p_value) %>%
    dplyr::arrange(p_value)

  # Create summary table
  summary_table <- coef_df %>%
    dplyr::select(lipid_class, term, estimate, p_value, significant, r_squared) %>%
    tidyr::pivot_wider(
      id_cols = c(lipid_class, r_squared),
      names_from = term,
      values_from = c(estimate, p_value, significant),
      names_glue = "{term}_{.value}"
    )

  cat("Classes with chain length effect (p <", p_cutoff, "):", nrow(chain_classes), "\n")
  if (nrow(chain_classes) > 0) {
    cat("  Positive effect (longer chains → higher logFC):\n")
    pos_chain <- chain_classes %>% dplyr::filter(estimate > 0)
    if (nrow(pos_chain) > 0) {
      for (i in seq_len(nrow(pos_chain))) {
        cat("    ", pos_chain$lipid_class[i], " (est = ", round(pos_chain$estimate[i], 4),
            ", p = ", format.pval(pos_chain$p_value[i]), ")\n", sep = "")
      }
    } else {
      cat("    None\n")
    }
    cat("  Negative effect (longer chains → lower logFC):\n")
    neg_chain <- chain_classes %>% dplyr::filter(estimate < 0)
    if (nrow(neg_chain) > 0) {
      for (i in seq_len(nrow(neg_chain))) {
        cat("    ", neg_chain$lipid_class[i], " (est = ", round(neg_chain$estimate[i], 4),
            ", p = ", format.pval(neg_chain$p_value[i]), ")\n", sep = "")
      }
    } else {
      cat("    None\n")
    }
  }

  cat("\nClasses with double bond effect (p <", p_cutoff, "):", nrow(db_classes), "\n")
  if (nrow(db_classes) > 0) {
    cat("  Positive effect (more unsaturation → higher logFC):\n")
    pos_db <- db_classes %>% dplyr::filter(estimate > 0)
    if (nrow(pos_db) > 0) {
      for (i in seq_len(nrow(pos_db))) {
        cat("    ", pos_db$lipid_class[i], " (est = ", round(pos_db$estimate[i], 4),
            ", p = ", format.pval(pos_db$p_value[i]), ")\n", sep = "")
      }
    } else {
      cat("    None\n")
    }
    cat("  Negative effect (more unsaturation → lower logFC):\n")
    neg_db <- db_classes %>% dplyr::filter(estimate < 0)
    if (nrow(neg_db) > 0) {
      for (i in seq_len(nrow(neg_db))) {
        cat("    ", neg_db$lipid_class[i], " (est = ", round(neg_db$estimate[i], 4),
            ", p = ", format.pval(neg_db$p_value[i]), ")\n", sep = "")
      }
    } else {
      cat("    None\n")
    }
  }

  cat("\nClasses with interaction effect (p <", p_cutoff, "):", nrow(interaction_classes), "\n")
  if (nrow(interaction_classes) > 0) {
    for (i in seq_len(nrow(interaction_classes))) {
      cat("  ", interaction_classes$lipid_class[i], " (est = ", round(interaction_classes$estimate[i], 4),
          ", p = ", format.pval(interaction_classes$p_value[i]), ")\n", sep = "")
    }
  }

  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  return(list(
    classes_with_chain_effect = chain_classes,
    classes_with_db_effect = db_classes,
    classes_with_interaction = interaction_classes,
    summary_table = summary_table
  ))
}


#' ==============================================================================
#' GENE-LIPID CORRELATION ANALYSIS FUNCTIONS
#' ==============================================================================


#' Compute Gene-Lipid Correlations for Structural Analysis
#'
#' Computes Pearson correlations between gene expression and lipid abundances,
#' filtered to MOFA-selected features and lipid metabolism genes.
#'
#' This function replicates the exact workflow from
#' lipid_rna_correaltion_for_structural_analysis.R
#'
#' @param mofa_weights Data frame with MOFA weights. Must have columns:
#'   - feature: Feature names (genes or lipids with view suffix)
#'   - view: View name (cell type or "lipidomic")
#'   - factor: Factor number
#'   - weight: MOFA weight value
#' @param rna_data Matrix or data frame. Gene expression data (genes x samples).
#'   Row names = gene symbols, column names = sample IDs.
#' @param lipidomic_data Matrix or data frame. Lipid abundance data (lipids x samples).
#' @param cell_type Character. Which RNA view to use (e.g., "Hepatocytes").
#' @param factor_number Numeric. Which MOFA factor to use (default: 1).
#' @param top_percent Numeric. Top percent of features by absolute weight (default: 0.5 = top 50%).
#' @param lipid_metabolism_genes Character vector. Gene symbols for lipid metabolism genes.
#'   If NULL, loads from GO:BP "GOBP_CELLULAR_LIPID_METABOLIC_PROCESS".
#' @param verbose Logical. Print progress messages?
#'
#' @return Matrix of correlations (genes x lipids)
#'
#' @examples
#' # Load MOFA weights
#' load("results/mofa/main_model_RNA_lipid_weight.RData")
#'
#' # Load data
#' rna_data <- load_rna_data()
#' lipidomic_data <- load_lipidomic_data()
#'
#' # Compute correlations
#' cor_matrix <- compute_gene_lipid_correlations_for_structure(
#'   mofa_weights = weight_RNA_and_lip,
#'   rna_data = rna_data,
#'   lipidomic_data = lipidomic_data,
#'   cell_type = "Hepatocytes",
#'   factor_number = 1,
#'   top_percent = 0.5
#' )
#'
#' @export
compute_gene_lipid_correlations_for_structure <- function(mofa_weights,
                                                           rna_data,
                                                           lipidomic_data,
                                                           cell_type = "Hepatocytes",
                                                           factor_number = 1,
                                                           top_percent = 0.5,
                                                           lipid_metabolism_genes = NULL,
                                                           verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Computing Gene-Lipid Correlations for Structural Analysis\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  # Step 1: Filter MOFA weights to selected factor
  # ==============================================
  if (verbose) cat("Step 1: Filtering MOFA weights to factor", factor_number, "...\n")

  factor_weights <- mofa_weights %>%
    dplyr::filter(factor == factor_number)

  if (verbose) {
    cat("  Total features in factor", factor_number, ":", nrow(factor_weights), "\n")
  }

  # Step 2: Select top percent of features by absolute weight
  # ==========================================================
  if (verbose) cat("\nStep 2: Selecting top", top_percent * 100, "% of features by |weight|...\n")

  # Calculate threshold for top percent (both tails)
  weight_threshold <- quantile(abs(factor_weights$weight), 1 - top_percent, na.rm = TRUE)

  top_features <- factor_weights %>%
    dplyr::filter(abs(weight) >= weight_threshold)

  if (verbose) {
    cat("  Weight threshold (|weight| >=):", weight_threshold, "\n")
    cat("  Features selected:", nrow(top_features), "\n")
    cat("  Views represented:\n")
    view_counts <- table(top_features$view)
    for (v in names(view_counts)) {
      cat("    ", v, ":", view_counts[v], "\n")
    }
  }

  # Step 3: Extract gene and lipid feature lists
  # =============================================
  if (verbose) cat("\nStep 3: Extracting gene and lipid feature lists...\n")

  # RNA features (remove view suffix)
  rna_features <- top_features %>%
    dplyr::filter(view == cell_type) %>%
    pull(feature)

  # Remove view suffix to get gene names
  rna_genes <- gsub(paste0("_", cell_type), "", rna_features)

  # Lipid features (remove "_lipidomic" suffix)
  lipid_features <- top_features %>%
    dplyr::filter(view == "lipidomic") %>%
    pull(feature)

  lipid_names <- gsub("_lipidomic", "", lipid_features)

  if (verbose) {
    cat("  RNA genes from", cell_type, ":", length(rna_genes), "\n")
    cat("  Lipids:", length(lipid_names), "\n")
  }

  # Step 4: Filter to lipid metabolism genes
  # =========================================
  if (verbose) cat("\nStep 4: Filtering to lipid metabolism genes...\n")

  if (is.null(lipid_metabolism_genes)) {
    # Load from gene set file
    if (file.exists("support/m5.all.v2023.2.Mm.symbols.gmt")) {
      if (!exists("import_gmt")) {
        source("scripts/support_functions.R")
      }

      # Load GSEABase if needed for import_gmt
      if (!requireNamespace("GSEABase", quietly = TRUE)) {
        warning("GSEABase package not available. Using all genes.")
        lipid_metabolism_genes <- rna_genes
      } else {
        library(GSEABase)
        gene_sets <- import_gmt("support/m5.all.v2023.2.Mm.symbols.gmt")

        if ("GOBP_CELLULAR_LIPID_METABOLIC_PROCESS" %in% names(gene_sets)) {
          lipid_metabolism_genes <- gene_sets[["GOBP_CELLULAR_LIPID_METABOLIC_PROCESS"]]
        } else {
          warning("Lipid metabolism gene set not found. Using all genes.")
          lipid_metabolism_genes <- rna_genes
        }
      }
    } else {
      warning("Gene set file not found. Using all genes.")
      lipid_metabolism_genes <- rna_genes
    }
  }

  # Filter genes
  rna_genes_filtered <- rna_genes[rna_genes %in% lipid_metabolism_genes]

  if (verbose) {
    cat("  Lipid metabolism genes available:", length(lipid_metabolism_genes), "\n")
    cat("  Genes after filtering:", length(rna_genes_filtered), "\n")
  }

  if (length(rna_genes_filtered) == 0) {
    stop("No genes remaining after filtering to lipid metabolism genes")
  }

  # Step 5: Prepare data matrices
  # ==============================
  if (verbose) cat("\nStep 5: Preparing data matrices...\n")

  # Filter RNA data to the specified cell type
  cell_type_pattern <- paste0("_", cell_type)
  rna_celltype_samples <- grep(cell_type_pattern, colnames(rna_data), value = TRUE)

  if (length(rna_celltype_samples) == 0) {
    stop("No samples found for cell type: ", cell_type)
  }

  # Subset RNA to cell type
  rna_celltype <- rna_data[, rna_celltype_samples, drop = FALSE]

  # Remove suffixes from RNA sample names to match lipidomic sample names
  # Remove cell type suffix (e.g., "_Hepatocytes")
  # Remove gonad suffix (e.g., "_Ovaries", "_Testes")
  rna_sample_names_clean <- colnames(rna_celltype)
  rna_sample_names_clean <- gsub("_Testes", "", rna_sample_names_clean)
  rna_sample_names_clean <- gsub("_Ovaries", "", rna_sample_names_clean)
  rna_sample_names_clean <- gsub(cell_type_pattern, "", rna_sample_names_clean)
  colnames(rna_celltype) <- rna_sample_names_clean

  # Find common samples
  common_samples <- intersect(colnames(rna_celltype), colnames(lipidomic_data))

  if (length(common_samples) == 0) {
    stop("No common samples between RNA and lipidomic data after cleaning sample names.\n",
         "RNA samples (first 5): ", paste(head(colnames(rna_celltype), 5), collapse = ", "), "\n",
         "Lipidomic samples (first 5): ", paste(head(colnames(lipidomic_data), 5), collapse = ", "))
  }

  if (verbose) {
    cat("  ", cell_type, "samples:", length(rna_celltype_samples), "\n")
    cat("  Common samples after name cleaning:", length(common_samples), "\n")
  }

  # Subset and filter data
  rna_subset <- rna_celltype[rna_genes_filtered, common_samples, drop = FALSE]
  lipid_subset <- lipidomic_data[lipid_names, common_samples, drop = FALSE]

  # Remove any genes/lipids with missing data
  rna_subset <- rna_subset[complete.cases(rna_subset), , drop = FALSE]
  lipid_subset <- lipid_subset[complete.cases(lipid_subset), , drop = FALSE]

  if (verbose) {
    cat("  Final dimensions:\n")
    cat("    RNA:", nrow(rna_subset), "genes x", ncol(rna_subset), "samples\n")
    cat("    Lipids:", nrow(lipid_subset), "lipids x", ncol(lipid_subset), "samples\n")
  }

  # Step 6: Compute correlations
  # =============================
  if (verbose) cat("\nStep 6: Computing Pearson correlations...\n")

  # Transpose for correlation (features as columns)
  rna_t <- t(rna_subset)
  lipid_t <- t(lipid_subset)

  # Compute correlation matrix (genes x lipids)
  cor_matrix <- cor(rna_t, lipid_t, method = "pearson", use = "pairwise.complete.obs")

  if (verbose) {
    cat("  Correlation matrix:", nrow(cor_matrix), "genes x", ncol(cor_matrix), "lipids\n")
    cat("  Range:", round(min(cor_matrix, na.rm = TRUE), 3), "to",
        round(max(cor_matrix, na.rm = TRUE), 3), "\n")
    cat("  Mean |correlation|:", round(mean(abs(cor_matrix), na.rm = TRUE), 3), "\n")
  }

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Correlation computation complete!\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  return(cor_matrix)
}


#' Load Lipid Structural Features
#'
#' Loads lipid structural data (chain length, double bond position) from CSV files.
#'
#' Expected format: CSV files organized by lipid class in a directory, with columns:
#' - lipid_name: Full lipid name (e.g., "PC(16:0/18:1)")
#' - C: Chain length (total carbon atoms)
#' - DB: Double bond position (number of double bonds)
#' - LogFC: Log fold change from differential analysis
#' - Additional columns (e.g., P.Value, adj.P.Val) are optional
#'
#' @param data_dir Character. Directory containing CSV files with structural data.
#'   Default: "data/metabolomics/for_bubbles_with_names/"
#' @param verbose Logical. Print progress messages?
#'
#' @return Data frame with columns:
#'   - lipid_name: Lipid identifier
#'   - lipid_class: Lipid class extracted from name
#'   - chain_length: Total carbon atoms (C)
#'   - db_position: Number of double bonds (DB)
#'   - logfc: Log fold change
#'   - ... (additional columns from input files)
#'
#' @examples
#' lipid_structure <- load_lipid_structural_features(
#'   data_dir = "data/metabolomics/for_bubbles_with_names/"
#' )
#'
#' @export
load_lipid_structural_features <- function(data_dir = "data/metabolomics/for_bubbles_with_names/",
                                           verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Loading Lipid Structural Features\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  # Find all CSV files in directory
  csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(csv_files) == 0) {
    stop("No CSV files found in ", data_dir)
  }

  if (verbose) {
    cat("Found", length(csv_files), "CSV files\n\n")
  }

  # Read all files and combine
  all_data <- list()

  for (file in csv_files) {
    # Extract lipid class from filename
    filename <- basename(file)
    lipid_class <- gsub("\\.csv$", "", filename)

    # Read data
    data <- read_csv(file, show_col_types = FALSE)

    # Standardize column names (case-insensitive matching)
    colnames(data) <- tolower(colnames(data))

    # Check required columns
    required_cols <- c("c", "db", "logfc")
    missing_cols <- setdiff(required_cols, colnames(data))

    if (length(missing_cols) > 0) {
      warning("File ", filename, " missing columns: ", paste(missing_cols, collapse = ", "),
              ". Skipping.")
      next
    }

    # Add lipid class
    data$lipid_class <- lipid_class

    # Rename columns to standard names
    # Use .data pronoun to avoid confusion with base R functions
    data <- data %>%
      dplyr::mutate(
        chain_length = .data$c,
        db_position = .data$db
      )

    # Try to extract lipid names
    if ("lipid_name" %in% colnames(data)) {
      # Already has lipid names
    } else if ("lipid" %in% colnames(data)) {
      data <- data %>% dplyr::mutate(lipid_name = .data$lipid)
    } else if ("id" %in% colnames(data)) {
      data <- data %>% dplyr::mutate(lipid_name = .data$id)
    } else if ("name" %in% colnames(data)) {
      data <- data %>% dplyr::mutate(lipid_name = .data$name)
    } else {
      # Create placeholder names
      data$lipid_name <- paste0(lipid_class, "_", seq_len(nrow(data)))
    }

    # Select only relevant columns
    data <- data %>%
      dplyr::select(lipid_name, lipid_class, chain_length, db_position, logfc)

    all_data[[lipid_class]] <- data

    if (verbose) {
      cat("  ", lipid_class, ":", nrow(data), "lipids\n")
    }
  }

  # Combine all data
  combined_data <- bind_rows(all_data)

  if (verbose) {
    cat("\n")
    cat("Total lipids loaded:", nrow(combined_data), "\n")
    cat("Lipid classes:", length(unique(combined_data$lipid_class)), "\n")
    cat("  Classes:", paste(unique(combined_data$lipid_class), collapse = ", "), "\n")
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  return(combined_data)
}


#' Run Structural Regression Analysis
#'
#' Tests how lipid structural features (chain length, double bond position)
#' predict gene-lipid correlation strength.
#'
#' Main regression model:
#'   gene_correlation ~ lipid_class + lipid_abundance + chain_length +
#'                      db_position + lipid_abundance:chain_length +
#'                      lipid_abundance:db_position
#'
#' This global model tests effects across all lipid classes simultaneously.
#'
#' @param cor_matrix Matrix. Gene-lipid correlations (genes x lipids).
#' @param lipid_structure Data frame. Lipid structural features from
#'   load_lipid_structural_features().
#' @param min_lipids_per_class Numeric. Minimum lipids required per class (default: 10).
#' @param verbose Logical. Print progress messages?
#'
#' @return Data frame with columns:
#'   - gene: Gene symbol
#'   - term: Regression term (e.g., "chain_length", "db_position")
#'   - estimate: Coefficient estimate
#'   - std_error: Standard error
#'   - t_value: t-statistic
#'   - p_value: P-value
#'   - significant: Logical, TRUE if p < 0.05
#'
#' @examples
#' regression_results <- run_structural_regression(
#'   cor_matrix = gene_lipid_cor,
#'   lipid_structure = lipid_features
#' )
#'
#' @export
run_structural_regression <- function(cor_matrix,
                                      lipid_structure,
                                      min_lipids_per_class = 10,
                                      verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Running Structural Regression Analysis (Global Model)\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  # Prepare data
  # ============
  if (verbose) cat("Preparing data...\n")

  # Get lipid names from correlation matrix
  cor_lipids <- colnames(cor_matrix)
  struct_lipids <- lipid_structure$lipid_name

  # Match lipids (may need fuzzy matching if names don't match exactly)
  common_lipids <- intersect(cor_lipids, struct_lipids)

  if (length(common_lipids) == 0) {
    stop("No matching lipids between correlation matrix and structural data.\n",
         "Check lipid naming conventions.")
  }

  if (verbose) {
    cat("  Lipids in correlation matrix:", length(cor_lipids), "\n")
    cat("  Lipids in structural data:", length(struct_lipids), "\n")
    cat("  Common lipids:", length(common_lipids), "\n\n")
  }

  # Filter to common lipids
  cor_matrix_filtered <- cor_matrix[, common_lipids, drop = FALSE]
  lipid_structure_filtered <- lipid_structure %>%
    dplyr::filter(lipid_name %in% common_lipids)

  # Filter lipid classes with sufficient data
  class_counts <- table(lipid_structure_filtered$lipid_class)
  valid_classes <- names(class_counts)[class_counts >= min_lipids_per_class]

  if (verbose) {
    cat("Lipid class filtering:\n")
    cat("  Classes with >=", min_lipids_per_class, "lipids:", length(valid_classes), "\n")
    cat("  Valid classes:", paste(valid_classes, collapse = ", "), "\n\n")
  }

  lipid_structure_filtered <- lipid_structure_filtered %>%
    dplyr::filter(lipid_class %in% valid_classes)

  # Update common lipids after class filtering
  common_lipids <- lipid_structure_filtered$lipid_name
  cor_matrix_filtered <- cor_matrix_filtered[, common_lipids, drop = FALSE]

  # Run regression for each gene
  # =============================
  if (verbose) {
    cat("Running regression for", nrow(cor_matrix_filtered), "genes...\n")
    cat("Progress: ")
  }

  all_results <- list()
  n_genes <- nrow(cor_matrix_filtered)
  progress_interval <- max(1, floor(n_genes / 20))

  for (i in seq_len(n_genes)) {
    gene <- rownames(cor_matrix_filtered)[i]

    # Show progress
    if (verbose && i %% progress_interval == 0) {
      cat(round(100 * i / n_genes), "% ")
    }

    # Get correlations for this gene
    gene_cors <- cor_matrix_filtered[i, ]

    # Build data frame for regression
    model_df <- data.frame(
      lipid_name = names(gene_cors),
      gene_correlation = as.numeric(gene_cors),
      stringsAsFactors = FALSE
    )

    # Merge with structural data
    model_df <- model_df %>%
      left_join(lipid_structure_filtered, by = "lipid_name") %>%
      dplyr::filter(!is.na(gene_correlation) & !is.na(chain_length) & !is.na(db_position))

    if (nrow(model_df) < 10) {
      next  # Skip genes with insufficient data
    }

    # Run regression
    tryCatch({
      model <- lm(gene_correlation ~ lipid_class + logfc + chain_length + db_position +
                    logfc:chain_length + logfc:db_position + chain_length:db_position,
                  data = model_df)

      # Extract coefficients
      coef_summary <- summary(model)$coefficients

      # Keep only structural terms (not lipid class indicators)
      terms_of_interest <- c("logfc", "chain_length", "db_position",
                             "logfc:chain_length", "logfc:db_position",
                             "chain_length:db_position")

      coef_df <- as.data.frame(coef_summary)
      coef_df$term <- rownames(coef_df)
      coef_df <- coef_df %>%
        dplyr::filter(term %in% terms_of_interest) %>%
        dplyr::mutate(
          gene = gene,
          estimate = Estimate,
          std_error = `Std. Error`,
          t_value = `t value`,
          p_value = `Pr(>|t|)`,
          significant = p_value < 0.05
        ) %>%
        dplyr::select(gene, term, estimate, std_error, t_value, p_value, significant)

      all_results[[gene]] <- coef_df

    }, error = function(e) {
      # Skip genes with model failures
    })
  }

  if (verbose) cat("\n\n")

  # Combine results
  combined_results <- bind_rows(all_results)

  if (verbose) {
    cat("Results summary:\n")
    cat("  Total gene-term pairs:", nrow(combined_results), "\n")

    if (nrow(combined_results) > 0 && "gene" %in% colnames(combined_results)) {
      cat("  Unique genes:", length(unique(combined_results$gene)), "\n")
      if ("significant" %in% colnames(combined_results)) {
        cat("  Significant associations (p < 0.05):",
            sum(combined_results$significant, na.rm = TRUE), "\n")
        cat("\n")
        cat("Significant associations by term:\n")
        sig_by_term <- combined_results %>%
          dplyr::filter(significant) %>%
          dplyr::group_by(term) %>%
          dplyr::summarise(n_significant = n(), .groups = "drop")
        print(sig_by_term)
      }
    } else {
      cat("  No results returned (empty dataframe)\n")
    }
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  return(combined_results)
}


#' Run Class-Specific Structural Regression
#'
#' Tests structural associations separately within each lipid class.
#'
#' For each lipid class, runs regression:
#'   gene_correlation ~ lipid_abundance + chain_length + db_position +
#'                      lipid_abundance:chain_length +
#'                      lipid_abundance:db_position + chain_length:db_position
#'
#' This reveals class-specific structural preferences.
#'
#' @param cor_matrix Matrix. Gene-lipid correlations (genes x lipids).
#' @param lipid_structure Data frame. Lipid structural features.
#' @param min_lipids_per_class Numeric. Minimum lipids per class (default: 10).
#' @param verbose Logical. Print progress messages?
#'
#' @return Data frame with columns:
#'   - gene: Gene symbol
#'   - lipid_class: Lipid class analyzed
#'   - term: Regression term
#'   - estimate: Coefficient
#'   - std_error: Standard error
#'   - t_value: t-statistic
#'   - p_value: P-value
#'   - significant: Logical, TRUE if p < 0.05
#'
#' @examples
#' class_results <- run_class_specific_regression(
#'   cor_matrix = gene_lipid_cor,
#'   lipid_structure = lipid_features
#' )
#'
#' @export
run_class_specific_regression <- function(cor_matrix,
                                          lipid_structure,
                                          min_lipids_per_class = 10,
                                          verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Running Class-Specific Structural Regression\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  # Get valid lipid classes
  common_lipids <- intersect(colnames(cor_matrix), lipid_structure$lipid_name)
  lipid_structure_filtered <- lipid_structure %>%
    dplyr::filter(lipid_name %in% common_lipids)

  class_counts <- table(lipid_structure_filtered$lipid_class)
  valid_classes <- names(class_counts)[class_counts >= min_lipids_per_class]

  if (verbose) {
    cat("Valid lipid classes (n >=", min_lipids_per_class, "):\n")
    for (cls in valid_classes) {
      cat("  ", cls, ":", class_counts[cls], "lipids\n")
    }
    cat("\n")
  }

  # Run analysis for each class
  all_class_results <- list()

  for (cls in valid_classes) {
    if (verbose) cat("Processing", cls, "...\n")

    # Filter to this class
    class_lipids <- lipid_structure_filtered %>%
      dplyr::filter(lipid_class == cls) %>%
      pull(lipid_name)

    class_cor_matrix <- cor_matrix[, class_lipids, drop = FALSE]
    class_structure <- lipid_structure_filtered %>%
      dplyr::filter(lipid_class == cls)

    # Run regression for each gene
    class_results <- list()

    for (i in seq_len(nrow(class_cor_matrix))) {
      gene <- rownames(class_cor_matrix)[i]

      # Get correlations
      gene_cors <- class_cor_matrix[i, ]

      # Build model data
      model_df <- data.frame(
        lipid_name = names(gene_cors),
        gene_correlation = as.numeric(gene_cors),
        stringsAsFactors = FALSE
      )

      model_df <- model_df %>%
        left_join(class_structure, by = "lipid_name") %>%
        dplyr::filter(!is.na(gene_correlation) & !is.na(chain_length) & !is.na(db_position))

      if (nrow(model_df) < 5) {
        next  # Need at least 5 points
      }

      # Run regression (without lipid_class term since it's constant)
      tryCatch({
        model <- lm(gene_correlation ~ logfc + chain_length + db_position +
                      logfc:chain_length + logfc:db_position + chain_length:db_position,
                    data = model_df)

        coef_summary <- summary(model)$coefficients

        terms_of_interest <- c("logfc", "chain_length", "db_position",
                               "logfc:chain_length", "logfc:db_position",
                               "chain_length:db_position")

        coef_df <- as.data.frame(coef_summary)
        coef_df$term <- rownames(coef_df)
        coef_df <- coef_df %>%
          dplyr::filter(term %in% terms_of_interest) %>%
          dplyr::mutate(
            gene = gene,
            lipid_class = cls,
            estimate = Estimate,
            std_error = `Std. Error`,
            t_value = `t value`,
            p_value = `Pr(>|t|)`,
            significant = p_value < 0.05
          ) %>%
          dplyr::select(gene, lipid_class, term, estimate, std_error, t_value, p_value, significant)

        class_results[[gene]] <- coef_df

      }, error = function(e) {
        # Skip on error
      })
    }

    if (length(class_results) > 0) {
      all_class_results[[cls]] <- bind_rows(class_results)
    }
  }

  # Combine all results
  combined_results <- bind_rows(all_class_results)

  if (verbose) {
    cat("\n")
    cat("Class-specific results summary:\n")
    cat("  Total gene-class-term combinations:", nrow(combined_results), "\n")

    if (nrow(combined_results) > 0 && "significant" %in% colnames(combined_results)) {
      cat("  Significant associations (p < 0.05):",
          sum(combined_results$significant, na.rm = TRUE), "\n")
      cat("\n")
      cat("Significant associations by class and term:\n")
      sig_summary <- combined_results %>%
        dplyr::filter(significant) %>%
        dplyr::group_by(lipid_class, term) %>%
        dplyr::summarise(n_significant = n(), .groups = "drop") %>%
        dplyr::arrange(lipid_class, term)
      print(sig_summary)
    } else {
      cat("  No results with 'significant' column\n")
    }
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  return(combined_results)
}


#' Extract Significant Structural Associations
#'
#' Filters regression results to significant associations and summarizes by
#' structural feature type.
#'
#' @param regression_results Data frame from run_structural_regression() or
#'   run_class_specific_regression().
#' @param p_cutoff Numeric. P-value cutoff (default: 0.05).
#'
#' @return List with elements:
#'   - all_significant: All significant associations
#'   - chain_length_genes: Genes with chain length associations
#'   - db_position_genes: Genes with double bond associations
#'   - interaction_genes: Genes with interaction effects
#'   - multi_feature_genes: Genes with multiple feature associations
#'
#' @examples
#' sig_results <- extract_significant_associations(
#'   regression_results = regression_output,
#'   p_cutoff = 0.05
#' )
#'
#' @export
extract_significant_associations <- function(regression_results,
                                             p_cutoff = 0.05) {

  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Extracting Significant Structural Associations\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  # Check if results are empty or invalid
  if (nrow(regression_results) == 0) {
    cat("WARNING: Empty regression results. Returning empty dataframe.\n\n")
    return(data.frame(
      gene = character(),
      chain_length = integer(),
      db_position = integer(),
      interactions = integer(),
      total_significant = integer()
    ))
  }

  if (!"p_value" %in% colnames(regression_results)) {
    cat("WARNING: No p_value column in results. Returning empty dataframe.\n\n")
    return(data.frame(
      gene = character(),
      chain_length = integer(),
      db_position = integer(),
      interactions = integer(),
      total_significant = integer()
    ))
  }

  # Filter to significant
  sig_results <- regression_results %>%
    dplyr::filter(p_value < p_cutoff)

  cat("Total significant associations (p <", p_cutoff, "):", nrow(sig_results), "\n\n")

  # Categorize by term type
  chain_length_terms <- c("chain_length", "logfc:chain_length")
  db_position_terms <- c("db_position", "logfc:db_position")
  interaction_terms <- c("logfc:chain_length", "logfc:db_position", "chain_length:db_position")

  # Genes with chain length associations
  chain_length_genes <- sig_results %>%
    dplyr::filter(term %in% chain_length_terms) %>%
    dplyr::select(gene, term, estimate, p_value) %>%
    dplyr::arrange(p_value)

  # Genes with DB position associations
  db_position_genes <- sig_results %>%
    dplyr::filter(term %in% db_position_terms) %>%
    dplyr::select(gene, term, estimate, p_value) %>%
    dplyr::arrange(p_value)

  # Genes with interaction effects
  interaction_genes <- sig_results %>%
    dplyr::filter(term %in% interaction_terms) %>%
    dplyr::select(gene, term, estimate, p_value) %>%
    dplyr::arrange(p_value)

  # Genes with multiple structural associations
  gene_feature_counts <- sig_results %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_associations = n(),
      terms = paste(unique(term), collapse = ", "),
      min_p_value = min(p_value),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_associations > 1) %>%
    dplyr::arrange(desc(n_associations), min_p_value)

  # Print summaries
  cat("Genes with chain length associations:", nrow(chain_length_genes), "\n")
  cat("Genes with double bond associations:", nrow(db_position_genes), "\n")
  cat("Genes with interaction effects:", nrow(interaction_genes), "\n")
  cat("Genes with multiple feature associations:", nrow(gene_feature_counts), "\n\n")

  if (nrow(gene_feature_counts) > 0) {
    cat("Top genes with multiple associations:\n")
    print(head(gene_feature_counts, 10))
    cat("\n")
  }

  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  return(list(
    all_significant = sig_results,
    chain_length_genes = chain_length_genes,
    db_position_genes = db_position_genes,
    interaction_genes = interaction_genes,
    multi_feature_genes = gene_feature_counts
  ))
}


#' Save Structural Analysis Results
#'
#' Saves all structural analysis results to organized output files.
#'
#' @param cor_matrix Matrix. Gene-lipid correlation matrix.
#' @param lipid_structure Data frame. Lipid structural features.
#' @param regression_results Data frame. Global regression results.
#' @param class_specific_results Data frame. Class-specific regression results.
#' @param significant_results List. Output from extract_significant_associations().
#' @param output_dir Character. Output directory (default: "results/multi-omic/structural_model/").
#' @param verbose Logical. Print progress messages?
#'
#' @return Invisibly returns list of file paths created.
#'
#' @examples
#' save_structural_analysis_results(
#'   cor_matrix = gene_lipid_cor,
#'   lipid_structure = lipid_features,
#'   regression_results = global_results,
#'   class_specific_results = class_results,
#'   significant_results = sig_results,
#'   output_dir = "results/multi-omic/structural_model/"
#' )
#'
#' @export
save_structural_analysis_results <- function(cor_matrix,
                                             lipid_structure,
                                             regression_results,
                                             class_specific_results,
                                             significant_results,
                                             output_dir = "results/multi-omic/structural_model/",
                                             verbose = TRUE) {

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("Saving Structural Analysis Results\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  files_created <- c()

  # Save correlation matrix
  cor_file <- file.path(output_dir, "hepatocyte_cor_mat_RNA_lip.csv")
  write_csv(as.data.frame(cor_matrix) %>% tibble::rownames_to_column("gene"),
            cor_file)
  files_created <- c(files_created, cor_file)
  if (verbose) cat("  Saved:", cor_file, "\n")

  # Save lipid structure
  struct_file <- file.path(output_dir, "lipids_with_structure.csv")
  write_csv(lipid_structure, struct_file)
  files_created <- c(files_created, struct_file)
  if (verbose) cat("  Saved:", struct_file, "\n")

  # Save global regression results
  global_file <- file.path(output_dir, "gene_lipid_regression_coefficients.csv")
  write_csv(regression_results, global_file)
  files_created <- c(files_created, global_file)
  if (verbose) cat("  Saved:", global_file, "\n")

  # Save class-specific results
  class_file <- file.path(output_dir, "gene_lipid_class_specific_all_terms.csv")
  write_csv(class_specific_results, class_file)
  files_created <- c(files_created, class_file)
  if (verbose) cat("  Saved:", class_file, "\n")

  # Save significant associations
  sig_file <- file.path(output_dir, "significant_gene_lipid_class_associations.csv")
  write_csv(significant_results$all_significant, sig_file)
  files_created <- c(files_created, sig_file)
  if (verbose) cat("  Saved:", sig_file, "\n")

  # Save chain length genes
  chain_file <- file.path(output_dir, "genes_associated_with_chain_length.csv")
  write_csv(significant_results$chain_length_genes, chain_file)
  files_created <- c(files_created, chain_file)
  if (verbose) cat("  Saved:", chain_file, "\n")

  # Save DB position genes
  db_file <- file.path(output_dir, "genes_associated_with_db_position.csv")
  write_csv(significant_results$db_position_genes, db_file)
  files_created <- c(files_created, db_file)
  if (verbose) cat("  Saved:", db_file, "\n")

  # Save interaction genes
  int_file <- file.path(output_dir, "genes_with_significant_interactions.csv")
  write_csv(significant_results$interaction_genes, int_file)
  files_created <- c(files_created, int_file)
  if (verbose) cat("  Saved:", int_file, "\n")

  # Save multi-feature genes
  multi_file <- file.path(output_dir, "genes_multiple_structural_associations.csv")
  write_csv(significant_results$multi_feature_genes, multi_file)
  files_created <- c(files_created, multi_file)
  if (verbose) cat("  Saved:", multi_file, "\n")

  if (verbose) {
    cat("\n")
    cat("All results saved to:", output_dir, "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }

  invisible(files_created)
}
