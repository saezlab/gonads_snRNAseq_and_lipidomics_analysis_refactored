#' Comprehensive Test Suite for All Workflows
#'
#' This script systematically tests all refactored workflows to ensure they:
#' 1. Run without errors
#' 2. Produce expected output files
#' 3. Generate results consistent with manuscript findings
#'
#' Run this script from project root:
#'   source("tests/test_all_workflows.R")
#'
#' Expected runtime: ~30-60 minutes for all workflows

library(testthat)

# ==============================================================================
# Test Configuration
# ==============================================================================

# Set to TRUE to run full workflows (slow), FALSE for quick checks only
RUN_FULL_WORKFLOWS <- FALSE

# Set to TRUE to clean up test outputs after completion
CLEANUP_AFTER_TESTS <- FALSE

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("COMPREHENSIVE WORKFLOW TEST SUITE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Configuration:\n")
cat("  Run full workflows:", RUN_FULL_WORKFLOWS, "\n")
cat("  Cleanup after tests:", CLEANUP_AFTER_TESTS, "\n")
cat("\n")

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Test that a file exists and is not empty
test_file_exists <- function(filepath, description = NULL) {
  desc <- if (is.null(description)) filepath else description

  test_that(paste("File exists:", desc), {
    expect_true(file.exists(filepath),
                info = paste("File not found:", filepath))
  })

  if (file.exists(filepath)) {
    test_that(paste("File not empty:", desc), {
      expect_gt(file.size(filepath), 0,
                info = paste("File is empty:", filepath))
    })
  }
}

#' Test that an RData file contains expected objects
test_rdata_contents <- function(filepath, expected_objects, description = NULL) {
  desc <- if (is.null(description)) filepath else description

  test_that(paste("RData contains expected objects:", desc), {
    env <- new.env()
    load(filepath, envir = env)

    for (obj_name in expected_objects) {
      expect_true(exists(obj_name, envir = env),
                  info = paste("Object not found:", obj_name, "in", filepath))
    }
  })
}

#' Check that a differential expression result has expected structure
check_de_result_structure <- function(ttop_df, min_genes = 100) {
  # Check required columns
  required_cols <- c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
  expect_true(all(required_cols %in% names(ttop_df)),
              info = paste("Missing required columns:",
                          paste(setdiff(required_cols, names(ttop_df)), collapse = ", ")))

  # Check minimum number of genes
  expect_gte(nrow(ttop_df), min_genes,
             info = paste("Too few genes:", nrow(ttop_df)))

  # Check no NAs in critical columns
  expect_false(any(is.na(ttop_df$logFC)),
               info = "NA values found in logFC")
  expect_false(any(is.na(ttop_df$adj.P.Val)),
               info = "NA values found in adj.P.Val")
}

#' Run a workflow script and capture any errors
run_workflow_safely <- function(workflow_path, workflow_name) {
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("Testing:", workflow_name, "\n")
  cat(paste(rep("-", 80), collapse = ""), "\n\n")

  result <- tryCatch({
    source(workflow_path, local = TRUE)
    list(success = TRUE, error = NULL)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })

  test_that(paste("Workflow runs without error:", workflow_name), {
    expect_true(result$success,
                info = paste("Workflow failed with error:", result$error))
  })

  return(result$success)
}

# ==============================================================================
# Test 1: Data Preparation Workflow
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("TEST 1: Data Preparation & QC Workflow\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

test_that("Workflow 01: Data Preparation script exists", {
  expect_true(file.exists("workflows/01_data_preparation_qc.R"))
})

if (RUN_FULL_WORKFLOWS) {
  workflow_01_success <- run_workflow_safely(
    "workflows/01_data_preparation_qc.R",
    "01_data_preparation_qc"
  )
} else {
  cat("Skipping full workflow run (set RUN_FULL_WORKFLOWS = TRUE to run)\n")
  cat("Testing data loading functions only...\n\n")

  source("scripts/core/data_loading.R")

  # Test RNA data loading
  test_that("RNA data loads correctly", {
    scaled_RNA <- load_rna_data(verbose = FALSE)
    expect_equal(nrow(scaled_RNA), 20123,
                 info = "Expected 20,123 genes")
    expect_equal(ncol(scaled_RNA), 215,
                 info = "Expected 215 pseudobulk samples")
  })

  # Test metadata loading
  test_that("Metadata loads correctly", {
    metadata <- load_metadata(
      rna_metadata_file = "data/RNA/meta_data_df.csv",
      verbose = FALSE
    )
    expect_equal(nrow(metadata), 215,
                 info = "Expected 215 samples")
    expect_true("gonadal_type" %in% names(metadata),
                info = "gonadal_type column missing")
    expect_true("allosome_type" %in% names(metadata),
                info = "allosome_type column missing")
    expect_true("group" %in% names(metadata),
                info = "group column missing")
  })
}

# Check expected output files
cat("\nChecking expected output files...\n")
test_file_exists("results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc/data_summary.txt", "Data summary report")

# ==============================================================================
# Test 2: Pseudobulk RNA-seq Analysis
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("TEST 2: FCG Pseudobulk RNA-seq Analysis\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

test_that("Workflow 02: Pseudobulk analysis script exists", {
  expect_true(file.exists("workflows/02_fcg_pseudobulk_analysis.R"))
})

if (RUN_FULL_WORKFLOWS) {
  workflow_02_success <- run_workflow_safely(
    "workflows/02_fcg_pseudobulk_analysis.R",
    "02_fcg_pseudobulk_analysis"
  )

  if (workflow_02_success) {
    cat("\nValidating differential expression results...\n")

    # Check output files exist
    test_file_exists("results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/ttop_full_list.RData",
                    "DE results (ttop_full_list)")
    test_file_exists("results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/tval_full_df.RData",
                    "t-statistic matrix")

    # Load and validate DE results
    test_that("DE results have correct structure", {
      load("results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/ttop_full_list.RData")

      # Expected contrasts
      expected_contrasts <- c("gonadal_type", "allosome_type",
                             "XXT_vs_XXO", "XYT_vs_XYO",
                             "XYO_vs_XXO", "XYT_vs_XXT")

      expect_true(all(expected_contrasts %in% names(ttop_full_list)),
                  info = "Missing expected contrasts")

      # Expected cell types
      expected_cell_types <- c("Hepatocytes", "Endothelial cells",
                              "Immune cells", "Stellate cells",
                              "Cholangiocytes")

      for (contrast in expected_contrasts) {
        expect_true(all(expected_cell_types %in% names(ttop_full_list[[contrast]])),
                    info = paste("Missing cell types in", contrast))
      }

      # Check hepatocyte gonadal_type results (PRIMARY FINDING)
      hepa_gonad <- ttop_full_list$gonadal_type$Hepatocytes
      check_de_result_structure(hepa_gonad)

      # Should have ~4,000 DEGs at FDR <= 0.05
      n_sig <- sum(hepa_gonad$adj.P.Val <= 0.05)
      expect_gte(n_sig, 3000,
                 info = paste("Expected ~4,000 gonad-dependent genes, found:", n_sig))
      expect_lte(n_sig, 5000,
                 info = paste("Too many DEGs (expected ~4,000), found:", n_sig))
    })

    # Check PCA plots generated
    test_file_exists("results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/pca/pca_Hepatocytes.pdf",
                    "Hepatocyte PCA plot")

    # Check volcano plots generated
    test_file_exists("results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/volcanos/volcano_gonadal_type_Hepatocytes.pdf",
                    "Hepatocyte gonadal volcano plot")
  }
} else {
  cat("Skipping full workflow run (set RUN_FULL_WORKFLOWS = TRUE to run)\n")
  cat("Testing core differential analysis functions only...\n\n")

  source("scripts/core/data_loading.R")
  source("scripts/core/differential_analysis.R")
  source("scripts/config/contrasts.R")

  # Quick test: Run DE on hepatocytes only
  test_that("Differential analysis functions work", {
    scaled_RNA <- load_rna_data(verbose = FALSE)
    metadata <- load_metadata(
      rna_metadata_file = "data/RNA/meta_data_df.csv",
      verbose = FALSE
    )

    # Test gonadal_type contrast on hepatocytes
    hepa_data <- prepare_pseudobulk(
      rna_data = scaled_RNA,
      metadata = metadata,
      cell_type = "Hepatocytes",
      verbose = FALSE
    )

    comparisons <- build_contrast_list("gonadal_type")

    result <- run_differential_analysis(
      data = hepa_data$data,
      targets = hepa_data$targets,
      group_variable = "gonadal_type",
      comparisons = comparisons,
      verbose = FALSE
    )

    expect_true("Testes vs Ovaries" %in% names(result))
    ttop <- result[["Testes vs Ovaries"]]
    check_de_result_structure(ttop)

    # Check for expected number of DEGs
    n_sig <- sum(ttop$adj.P.Val <= 0.05)
    expect_gte(n_sig, 3000,
               info = paste("Expected ~4,000 gonad-dependent genes, found:", n_sig))
  })
}

# ==============================================================================
# Test 3: Lipidomic Analysis
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("TEST 3: FCG Lipidomic Analysis\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

test_that("Workflow 03: Lipidomic analysis script exists", {
  expect_true(file.exists("workflows/03_fcg_lipidomic_analysis.R"))
})

if (RUN_FULL_WORKFLOWS) {
  workflow_03_success <- run_workflow_safely(
    "workflows/03_fcg_lipidomic_analysis.R",
    "03_fcg_lipidomic_analysis"
  )

  if (workflow_03_success) {
    cat("\nValidating lipidomic results...\n")

    # Check output files
    test_file_exists("results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData",
                    "Lipidomic DE results")
    test_file_exists("results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gonade_type_ttop.csv",
                    "Gonad-dependent lipids")

    # Validate results
    test_that("Lipidomic results are correct", {
      load("results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData")

      # Check gonadal_type results
      gonad_results <- ttop_full_list$gonadal_type

      # Should have 733 significant lipids at FDR <= 0.05 (53.9% of 1,386)
      n_sig <- sum(gonad_results$adj.P.Val <= 0.05)
      expect_gte(n_sig, 650,
                 info = paste("Expected ~733 gonad-dependent lipids, found:", n_sig))
      expect_lte(n_sig, 850,
                 info = paste("Too many significant lipids, found:", n_sig))

      # Check allosome_type results (should be minimal)
      allosome_results <- ttop_full_list$allosome_type
      n_sig_chromo <- sum(allosome_results$adj.P.Val <= 0.05)
      expect_lte(n_sig_chromo, 50,
                 info = paste("Expected ~0 chromosome-dependent lipids, found:", n_sig_chromo))
    })
  }
} else {
  cat("Skipping full workflow run\n")
}

# ==============================================================================
# Test 4: MOFA Integration
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("TEST 4: MOFA Multi-Omic Integration\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

test_that("Workflow 04: MOFA integration script exists", {
  expect_true(file.exists("workflows/04_mofa_integration.R"))
})

if (RUN_FULL_WORKFLOWS) {
  workflow_04_success <- run_workflow_safely(
    "workflows/04_mofa_integration.R",
    "04_mofa_integration"
  )

  if (workflow_04_success) {
    # Check MOFA outputs
    test_file_exists("results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/main_model_RNA_lipid_weight.RData",
                    "MOFA model weights")
    test_file_exists("results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata",
                    "Gene-lipid correlation network")
  }
} else {
  cat("Skipping full workflow run (MOFA is computationally intensive)\n")
}

# ==============================================================================
# Test 5: Zonation Analysis
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("TEST 5: Zonation-Specific Analysis\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

test_that("Workflow 05: Zonation analysis script exists", {
  expect_true(file.exists("workflows/05_zonation_analysis.R"))
})

if (RUN_FULL_WORKFLOWS) {
  workflow_05_success <- run_workflow_safely(
    "workflows/05_zonation_analysis.R",
    "05_zonation_analysis"
  )
} else {
  cat("Skipping full workflow run\n")
}

# ==============================================================================
# Test 6: Structural Analysis
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("TEST 6: Lipid Structural Analysis\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

test_that("Workflow 06: Structural analysis script exists", {
  expect_true(file.exists("workflows/06_structural_analysis.R"))
})

if (RUN_FULL_WORKFLOWS) {
  workflow_06_success <- run_workflow_safely(
    "workflows/06_structural_analysis.R",
    "06_structural_analysis"
  )

  if (workflow_06_success) {
    # Check structural analysis outputs
    test_file_exists("results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv",
                    "Gene-lipid correlation matrix")
    test_file_exists("results/multi-omic/structural_model/lipids_with_structure.csv",
                    "Lipid structural features")
  }
} else {
  cat("Skipping full workflow run\n")
}

# ==============================================================================
# Test Summary
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("TEST SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

test_results <- test_results()
cat("Total tests:", test_results$n_tests, "\n")
cat("Passed:", test_results$n_pass, "\n")
cat("Failed:", test_results$n_fail, "\n")
cat("Warnings:", test_results$n_warn, "\n")
cat("Skipped:", test_results$n_skip, "\n")

if (test_results$n_fail > 0) {
  cat("\n❌ TESTS FAILED - Please review errors above\n")
} else {
  cat("\n✓ ALL TESTS PASSED\n")
}

# Cleanup if requested
if (CLEANUP_AFTER_TESTS && RUN_FULL_WORKFLOWS) {
  cat("\nCleaning up test outputs...\n")
  # Add cleanup code here if needed
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
