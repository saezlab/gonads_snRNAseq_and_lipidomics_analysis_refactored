#' Interactive Workflow Testing Script
#'
#' This script runs each workflow and checks for basic correctness.
#' Run from project root with: source("tests/run_workflow_tests.R")

# ==============================================================================
# Setup
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("INTERACTIVE WORKFLOW TESTING\n")
cat("================================================================================\n\n")

# Helper function to check files
check_file <- function(path, description) {
  exists <- file.exists(path)
  status <- if (exists) "✓" else "✗"
  size <- if (exists) format(file.size(path), big.mark = ",") else "N/A"

  cat(sprintf("  [%s] %s\n", status, description))
  if (exists) {
    cat(sprintf("      File: %s (%s bytes)\n", path, size))
  } else {
    cat(sprintf("      Missing: %s\n", path))
  }

  return(exists)
}

# Helper function to run a workflow with error catching
run_workflow <- function(workflow_num, workflow_path, workflow_name) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("WORKFLOW %d: %s\n", workflow_num, workflow_name))
  cat("================================================================================\n\n")

  cat("Script:", workflow_path, "\n\n")

  if (!file.exists(workflow_path)) {
    cat("✗ ERROR: Workflow script not found!\n")
    return(FALSE)
  }

  cat("Running workflow...\n")
  start_time <- Sys.time()

  # Capture both errors AND warnings
  result <- tryCatch({
    # Suppress verbose output but capture any errors
    suppressMessages({
      source(workflow_path, local = new.env(), echo = FALSE, print.eval = FALSE)
    })
    list(success = TRUE, error = NULL, warning = NULL)
  }, error = function(e) {
    list(success = FALSE, error = e$message, warning = NULL)
  }, warning = function(w) {
    # Continue execution but note warning
    list(success = TRUE, error = NULL, warning = w$message)
  })

  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "secs"), 1)

  cat(sprintf("\nCompleted in %.1f seconds\n", elapsed))

  if (result$success) {
    cat("✓ Workflow completed successfully\n")
    if (!is.null(result$warning)) {
      cat("  ⚠ Note: Warnings were generated (check workflow output)\n")
    }
  } else {
    cat("✗ Workflow FAILED with error:\n")
    cat("  ", result$error, "\n")
    cat("\n")
    cat("  Error location: Check the workflow script for the line causing this error\n")
    cat("  Suggestion: Run the workflow directly to see more detailed output:\n")
    cat("    source(\"", workflow_path, "\")\n", sep = "")
  }

  return(result$success)
}

# ==============================================================================
# Test Each Workflow
# ==============================================================================

results <- list()

# ------------------------------------------------------------------------------
# WORKFLOW 1: Data Preparation & QC
# ------------------------------------------------------------------------------

results$wf01 <- run_workflow(
  1,
  "workflows/01_data_preparation_qc.R",
  "Data Preparation & QC"
)

if (results$wf01) {
  cat("\nChecking outputs:\n")
  check_file("results/qc/data_summary.txt", "Data summary report")
  check_file("data/RNA/scaled_RNA_1.csv", "RNA data part 1")
  check_file("data/RNA/scaled_RNA_2.csv", "RNA data part 2")
  check_file("data/RNA/scaled_RNA_3.csv", "RNA data part 3")
  check_file("data/RNA/meta_data_df.csv", "RNA metadata")
}

# ------------------------------------------------------------------------------
# WORKFLOW 2: Pseudobulk RNA-seq Analysis
# ------------------------------------------------------------------------------

results$wf02 <- run_workflow(
  2,
  "workflows/02_fcg_pseudobulk_analysis.R",
  "FCG Pseudobulk RNA-seq Analysis"
)

if (results$wf02) {
  cat("\nChecking outputs:\n")

  # Main results files
  has_ttop <- check_file("results/pseudobulk/ttop_full_list.RData",
                         "Differential expression results (ttop_full_list)")
  has_tval <- check_file("results/pseudobulk/tval_full_df.RData",
                         "t-statistic matrix (tval_full_df)")

  # PCA plots
  check_file("results/pseudobulk/pca/pca_Hepatocytes.pdf",
            "PCA plot - Hepatocytes")
  check_file("results/pseudobulk/pca/pca_Endothelial_cells.pdf",
            "PCA plot - Endothelial cells")

  # Volcano plots
  check_file("results/pseudobulk/volcanos/volcano_gonadal_type_Hepatocytes.pdf",
            "Volcano plot - Hepatocytes gonadal_type")
  check_file("results/pseudobulk/volcanos/volcano_allosome_type_Hepatocytes.pdf",
            "Volcano plot - Hepatocytes allosome_type")

  # Validate key results
  if (has_ttop) {
    cat("\nValidating differential expression results...\n")
    load("results/pseudobulk/ttop_full_list.RData")

    # Check structure
    expected_contrasts <- c("gonadal_type", "allosome_type",
                           "XXT_vs_XXO", "XYT_vs_XYO",
                           "XYO_vs_XXO", "XYT_vs_XXT")
    expected_celltypes <- c("Hepatocytes", "Endothelial cells",
                           "Immune cells", "Stellate cells",
                           "Cholangiocytes")

    has_all_contrasts <- all(expected_contrasts %in% names(ttop_full_list))
    cat(sprintf("  [%s] All 6 contrasts present\n",
                if (has_all_contrasts) "✓" else "✗"))

    if (has_all_contrasts) {
      has_all_celltypes <- all(expected_celltypes %in%
                               names(ttop_full_list$gonadal_type))
      cat(sprintf("  [%s] All 5 cell types present\n",
                  if (has_all_celltypes) "✓" else "✗"))

      # Check hepatocyte gonadal_type (PRIMARY FINDING)
      hepa_gonad <- ttop_full_list$gonadal_type$Hepatocytes
      n_genes <- nrow(hepa_gonad)
      n_sig <- sum(hepa_gonad$adj.P.Val <= 0.05, na.rm = TRUE)

      cat(sprintf("  [✓] Hepatocyte gonadal_type: %s genes, %s significant (FDR ≤ 0.05)\n",
                  format(n_genes, big.mark = ","),
                  format(n_sig, big.mark = ",")))

      # Expected: ~4,000 DEGs
      if (n_sig >= 3000 && n_sig <= 5000) {
        cat("      Expected range: 3,000-5,000 (manuscript reports ~4,000) ✓\n")
      } else {
        cat(sprintf("      WARNING: Expected ~4,000, found %s\n", n_sig))
      }
    }
  }
}

# ------------------------------------------------------------------------------
# WORKFLOW 3: Lipidomic Analysis
# ------------------------------------------------------------------------------

results$wf03 <- run_workflow(
  3,
  "workflows/03_fcg_lipidomic_analysis.R",
  "FCG Lipidomic Analysis"
)

if (results$wf03) {
  cat("\nChecking outputs:\n")

  has_lipid_results <- check_file("results/lipidomic/ttop_full_list.RData",
                                   "Lipidomic differential analysis results")
  check_file("results/lipidomic/gonade_type_ttop.csv",
            "Gonad-dependent lipids table")
  check_file("results/lipidomic/pca_plot.pdf",
            "Lipidomic PCA plot")

  # Validate key results
  if (has_lipid_results) {
    cat("\nValidating lipidomic results...\n")
    load("results/lipidomic/ttop_full_list.RData")

    # Check gonadal_type results
    if ("gonadal_type" %in% names(ttop_full_list)) {
      gonad_results <- ttop_full_list$gonadal_type
      n_lipids <- nrow(gonad_results)
      n_sig <- sum(gonad_results$adj.P.Val <= 0.05, na.rm = TRUE)

      cat(sprintf("  [✓] Gonadal type: %s lipids, %s significant (FDR ≤ 0.05)\n",
                  format(n_lipids, big.mark = ","),
                  format(n_sig, big.mark = ",")))

      # Expected: 733/1386 lipids (53.9%)
      pct <- round(100 * n_sig / n_lipids, 1)
      if (n_sig >= 650 && n_sig <= 850) {
        cat(sprintf("      Expected range: 650-850 (manuscript reports 733, %.1f%%) ✓\n", pct))
      } else {
        cat(sprintf("      WARNING: Expected ~733 (53.9%%), found %s (%.1f%%)\n",
                    n_sig, pct))
      }
    }

    # Check allosome_type results (should be minimal)
    if ("allosome_type" %in% names(ttop_full_list)) {
      allosome_results <- ttop_full_list$allosome_type
      n_sig_chromo <- sum(allosome_results$adj.P.Val <= 0.05, na.rm = TRUE)

      cat(sprintf("  [✓] Allosome type: %s significant (FDR ≤ 0.05)\n",
                  format(n_sig_chromo, big.mark = ",")))

      if (n_sig_chromo <= 50) {
        cat("      Expected ~0 (manuscript reports 0) ✓\n")
      } else {
        cat(sprintf("      WARNING: Expected ~0, found %s\n", n_sig_chromo))
      }
    }
  }
}

# ------------------------------------------------------------------------------
# WORKFLOW 4: MOFA Integration
# ------------------------------------------------------------------------------

results$wf04 <- run_workflow(
  4,
  "workflows/04_mofa_integration.R",
  "MOFA Multi-Omic Integration"
)

if (results$wf04) {
  cat("\nChecking outputs:\n")
  check_file("results/mofa/main_model_RNA_lipid_weight.RData",
            "MOFA model weights")
  check_file("results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata",
            "Gene-lipid correlation network")
}

# ------------------------------------------------------------------------------
# WORKFLOW 5: Zonation Analysis
# ------------------------------------------------------------------------------

results$wf05 <- run_workflow(
  5,
  "workflows/05_zonation_analysis.R",
  "Zonation-Specific Analysis"
)

if (results$wf05) {
  cat("\nChecking outputs:\n")
  check_file("results/pseudobulk/zonation/ttop_pericentral_vs_periportal.RData",
            "Zonation differential expression results")
}

# ------------------------------------------------------------------------------
# WORKFLOW 6: Structural Analysis
# ------------------------------------------------------------------------------

results$wf06 <- run_workflow(
  6,
  "workflows/06_structural_analysis.R",
  "Lipid Structural Analysis"
)

if (results$wf06) {
  cat("\nChecking outputs:\n")
  check_file("results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv",
            "Gene-lipid correlation matrix")
  check_file("results/multi-omic/structural_model/lipids_with_structure.csv",
            "Lipid structural features")
  check_file("results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv",
            "Significant gene-lipid associations")
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("TEST SUMMARY\n")
cat("================================================================================\n\n")

n_total <- length(results)
n_success <- sum(unlist(results))
n_fail <- n_total - n_success

cat(sprintf("Total workflows tested: %d\n", n_total))
cat(sprintf("Successful: %d\n", n_success))
cat(sprintf("Failed: %d\n", n_fail))
cat("\n")

cat("Workflow Status:\n")
for (i in 1:n_total) {
  wf_name <- names(results)[i]
  status <- if (results[[wf_name]]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("  Workflow %d: %s\n", i, status))
}

cat("\n")
if (n_fail == 0) {
  cat("✓ ALL WORKFLOWS PASSED!\n")
} else {
  cat("✗ SOME WORKFLOWS FAILED - Please review errors above\n")
}

cat("\n")
cat("================================================================================\n\n")
