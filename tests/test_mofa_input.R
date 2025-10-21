#' Test MOFA Input Preparation - CRITICAL TEST
#'
#' ⚠️ **MOST IMPORTANT TEST IN THE ENTIRE REFACTORING** ⚠️
#'
#' This test verifies that the new prepare_mofa_input() function produces
#' EXACTLY IDENTICAL output to the original MOFA_lipidcell_RNA.R code.
#'
#' If this test fails, MOFA models need to be retrained (hours of computation).
#'
#' The test checks:
#' - Exact number of rows
#' - Exact column names and order
#' - Exact sample names
#' - Exact feature names (including view suffixes)
#' - Exact numeric values (within floating point tolerance)
#' - Exact view assignments

library(readr)
library(reshape2)
library(dplyr)

# Source modules
source("scripts/core/data_loading.R")
source("scripts/core/multiomics.R")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("⚠️  CRITICAL TEST: MOFA Input Preparation ⚠️\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("This test verifies EXACT numerical identity between:\n")
cat("  OLD: Original MOFA_lipidcell_RNA.R code\n")
cat("  NEW: prepare_mofa_input() function\n\n")

cat("If this test fails, MOFA models need retraining!\n\n")

# ==============================================================================
# Load Data
# ==============================================================================
cat("Loading data...\n")

# Load RNA data (using new function)
scaled_RNA <- load_rna_data(verbose = FALSE)

# Load lipidomic data
load("data/metabolomics/lipidomic.RData")

cat("  RNA:", nrow(scaled_RNA), "genes x", ncol(scaled_RNA), "samples\n")
cat("  Lipidomic:", nrow(lipidomic), "lipids x", ncol(lipidomic), "samples\n\n")

# ==============================================================================
# OLD METHOD (Exact copy from MOFA_lipidcell_RNA.R lines 29-60)
# ==============================================================================
cat("Running OLD method (original code)...\n")

# Create targets (not actually used in MOFA input, but part of original code)
targets_old <- data.frame(sample = names(scaled_RNA),
                          condition = names(scaled_RNA))
targets_old$condition <- gsub(".*_", "", targets_old$condition)
targets_old$sample <- gsub("_[A-Z].*", "", targets_old$sample)

# Clean sample names
scaled_RNA_old <- scaled_RNA
names(scaled_RNA_old) <- gsub("_Testes", "", names(scaled_RNA_old))
names(scaled_RNA_old) <- gsub("_Ovaries", "", names(scaled_RNA_old))

# Prepare RNA for MOFA
scaled_RNA_old$feature <- row.names(scaled_RNA_old)
RNA_MOFA_old <- melt(scaled_RNA_old)
RNA_MOFA_old$variable <- as.character(RNA_MOFA_old$variable)
RNA_MOFA_old$view <- gsub(".*_", "", RNA_MOFA_old$variable)
RNA_MOFA_old$variable <- gsub("_[A-Z].*", "", RNA_MOFA_old$variable)
names(RNA_MOFA_old) <- c("feature", "sample", "value", "view")
RNA_MOFA_old <- RNA_MOFA_old[, c(2, 1, 4, 3)]
RNA_MOFA_old <- RNA_MOFA_old[RNA_MOFA_old$sample %in% names(lipidomic), ]

# Prepare lipidomic for MOFA
lipidomic_MOFA_old <- lipidomic
lipidomic_MOFA_old$feature <- row.names(lipidomic_MOFA_old)
lipidomic_MOFA_old <- melt(lipidomic_MOFA_old)
lipidomic_MOFA_old$variable <- as.character(lipidomic_MOFA_old$variable)
lipidomic_MOFA_old$view <- "lipidomic"
names(lipidomic_MOFA_old) <- c("feature", "sample", "value", "view")

# Combine
MOFA_data_old <- as.data.frame(rbind(RNA_MOFA_old, lipidomic_MOFA_old))
MOFA_data_old$feature <- paste(MOFA_data_old$feature, MOFA_data_old$view, sep = "_")

cat("  OLD method complete:\n")
cat("    Rows:", nrow(MOFA_data_old), "\n")
cat("    Samples:", length(unique(MOFA_data_old$sample)), "\n")
cat("    Features:", length(unique(MOFA_data_old$feature)), "\n")
cat("    Views:", paste(unique(MOFA_data_old$view), collapse = ", "), "\n\n")

# ==============================================================================
# NEW METHOD (Using prepare_mofa_input())
# ==============================================================================
cat("Running NEW method (prepare_mofa_input)...\n")

MOFA_data_new <- prepare_mofa_input(
  rna_data = scaled_RNA,
  lipidomic_data = lipidomic,
  remove_gonad_suffix = TRUE,
  verbose = FALSE
)

cat("  NEW method complete:\n")
cat("    Rows:", nrow(MOFA_data_new), "\n")
cat("    Samples:", length(unique(MOFA_data_new$sample)), "\n")
cat("    Features:", length(unique(MOFA_data_new$feature)), "\n")
cat("    Views:", paste(unique(MOFA_data_new$view), collapse = ", "), "\n\n")

# ==============================================================================
# RIGOROUS COMPARISON
# ==============================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("RIGOROUS COMPARISON\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

all_tests_pass <- TRUE

# Test 1: Dimensions
# ==================
cat("Test 1: Dimensions\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

identical_rows <- nrow(MOFA_data_old) == nrow(MOFA_data_new)
identical_cols <- ncol(MOFA_data_old) == ncol(MOFA_data_new)

cat("  Rows match:    ", identical_rows, "(", nrow(MOFA_data_old), "vs", nrow(MOFA_data_new), ")\n")
cat("  Columns match: ", identical_cols, "(", ncol(MOFA_data_old), "vs", ncol(MOFA_data_new), ")\n")

test1_pass <- identical_rows && identical_cols
if (test1_pass) {
  cat("  ✓ Test 1 PASSED\n\n")
} else {
  cat("  ✗ Test 1 FAILED\n\n")
  all_tests_pass <- FALSE
}

# Test 2: Column Names and Order
# ===============================
cat("Test 2: Column Names and Order\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

identical_colnames <- identical(names(MOFA_data_old), names(MOFA_data_new))
cat("  Column names match: ", identical_colnames, "\n")
if (!identical_colnames) {
  cat("  OLD:", paste(names(MOFA_data_old), collapse = ", "), "\n")
  cat("  NEW:", paste(names(MOFA_data_new), collapse = ", "), "\n")
}

test2_pass <- identical_colnames
if (test2_pass) {
  cat("  ✓ Test 2 PASSED\n\n")
} else {
  cat("  ✗ Test 2 FAILED\n\n")
  all_tests_pass <- FALSE
}

# Test 3: Unique Samples
# =======================
cat("Test 3: Unique Samples\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

samples_old <- sort(unique(MOFA_data_old$sample))
samples_new <- sort(unique(MOFA_data_new$sample))

identical_samples <- identical(samples_old, samples_new)
cat("  Samples match: ", identical_samples, "\n")
cat("  OLD samples:", length(samples_old), "\n")
cat("  NEW samples:", length(samples_new), "\n")

if (!identical_samples) {
  cat("  Samples in OLD but not NEW:", setdiff(samples_old, samples_new), "\n")
  cat("  Samples in NEW but not OLD:", setdiff(samples_new, samples_old), "\n")
}

test3_pass <- identical_samples
if (test3_pass) {
  cat("  ✓ Test 3 PASSED\n\n")
} else {
  cat("  ✗ Test 3 FAILED\n\n")
  all_tests_pass <- FALSE
}

# Test 4: Unique Features
# ========================
cat("Test 4: Unique Features\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

features_old <- sort(unique(MOFA_data_old$feature))
features_new <- sort(unique(MOFA_data_new$feature))

identical_features <- identical(features_old, features_new)
cat("  Features match: ", identical_features, "\n")
cat("  OLD features:", length(features_old), "\n")
cat("  NEW features:", length(features_new), "\n")

if (!identical_features) {
  n_diff <- length(setdiff(features_old, features_new))
  if (n_diff > 0 && n_diff < 20) {
    cat("  Features in OLD but not NEW:\n")
    print(head(setdiff(features_old, features_new), 20))
  }
  n_diff_new <- length(setdiff(features_new, features_old))
  if (n_diff_new > 0 && n_diff_new < 20) {
    cat("  Features in NEW but not OLD:\n")
    print(head(setdiff(features_new, features_old), 20))
  }
}

test4_pass <- identical_features
if (test4_pass) {
  cat("  ✓ Test 4 PASSED\n\n")
} else {
  cat("  ✗ Test 4 FAILED\n\n")
  all_tests_pass <- FALSE
}

# Test 5: View Assignments
# =========================
cat("Test 5: View Assignments\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

views_old <- sort(unique(MOFA_data_old$view))
views_new <- sort(unique(MOFA_data_new$view))

identical_views <- identical(views_old, views_new)
cat("  Views match: ", identical_views, "\n")
cat("  OLD views:", paste(views_old, collapse = ", "), "\n")
cat("  NEW views:", paste(views_new, collapse = ", "), "\n")

test5_pass <- identical_views
if (test5_pass) {
  cat("  ✓ Test 5 PASSED\n\n")
} else {
  cat("  ✗ Test 5 FAILED\n\n")
  all_tests_pass <- FALSE
}

# Test 6: Numeric Values (MOST CRITICAL)
# ========================================
cat("Test 6: Numeric Values (Row-by-Row Comparison)\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

cat("  Sorting both datasets for comparison...\n")

# Sort both datasets by sample, feature, view for row-by-row comparison
MOFA_data_old_sorted <- MOFA_data_old[order(MOFA_data_old$sample,
                                             MOFA_data_old$feature,
                                             MOFA_data_old$view), ]
MOFA_data_new_sorted <- MOFA_data_new[order(MOFA_data_new$sample,
                                             MOFA_data_new$feature,
                                             MOFA_data_new$view), ]

# Reset row names for comparison
rownames(MOFA_data_old_sorted) <- NULL
rownames(MOFA_data_new_sorted) <- NULL

# Compare values with tolerance
values_comparison <- all.equal(
  MOFA_data_old_sorted$value,
  MOFA_data_new_sorted$value,
  tolerance = 1e-10
)

cat("  Values comparison:\n")
if (isTRUE(values_comparison)) {
  cat("    ✓ ALL VALUES IDENTICAL (within tolerance 1e-10)\n")
  test6_pass <- TRUE
} else {
  cat("    ✗ VALUES DIFFER\n")
  cat("    Difference:", values_comparison, "\n")

  # Find first differing row
  value_diffs <- which(abs(MOFA_data_old_sorted$value - MOFA_data_new_sorted$value) > 1e-10)
  if (length(value_diffs) > 0) {
    cat("    Number of differing rows:", length(value_diffs), "\n")
    cat("    First 5 differing rows:\n")
    for (i in head(value_diffs, 5)) {
      cat("      Row", i, ":\n")
      cat("        OLD:", MOFA_data_old_sorted$sample[i],
          MOFA_data_old_sorted$feature[i],
          MOFA_data_old_sorted$value[i], "\n")
      cat("        NEW:", MOFA_data_new_sorted$sample[i],
          MOFA_data_new_sorted$feature[i],
          MOFA_data_new_sorted$value[i], "\n")
    }
  }
  test6_pass <- FALSE
  all_tests_pass <- FALSE
}

if (test6_pass) {
  cat("  ✓ Test 6 PASSED\n\n")
} else {
  cat("  ✗ Test 6 FAILED\n\n")
}

# Test 7: Complete Data Frame Identity
# =====================================
cat("Test 7: Complete Data Frame Identity\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

complete_identity <- all.equal(MOFA_data_old_sorted, MOFA_data_new_sorted, tolerance = 1e-10)

if (isTRUE(complete_identity)) {
  cat("  ✓ COMPLETE IDENTITY CONFIRMED\n")
  test7_pass <- TRUE
} else {
  cat("  ✗ Differences detected:\n")
  cat("    ", complete_identity, "\n")
  test7_pass <- FALSE
  all_tests_pass <- FALSE
}

if (test7_pass) {
  cat("  ✓ Test 7 PASSED\n\n")
} else {
  cat("  ✗ Test 7 FAILED\n\n")
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("FINAL SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("  Test 1 (Dimensions):             ", ifelse(test1_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 2 (Column names):           ", ifelse(test2_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 3 (Samples):                ", ifelse(test3_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 4 (Features):               ", ifelse(test4_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 5 (Views):                  ", ifelse(test5_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 6 (Values):                 ", ifelse(test6_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 7 (Complete identity):      ", ifelse(test7_pass, "✓ PASS", "✗ FAIL"), "\n")

cat("\n")

if (all_tests_pass) {
  cat("✓✓✓ ALL TESTS PASSED ✓✓✓\n")
  cat("\n")
  cat("🎉 MOFA INPUT IS IDENTICAL! 🎉\n")
  cat("\n")
  cat("The new prepare_mofa_input() function produces EXACTLY the same output\n")
  cat("as the original code. It is SAFE to use in production.\n")
  cat("\n")
  cat("Pre-trained MOFA models can be used without retraining.\n")
} else {
  cat("✗✗✗ SOME TESTS FAILED ✗✗✗\n")
  cat("\n")
  cat("⚠️  WARNING: MOFA INPUT DIFFERS FROM ORIGINAL ⚠️\n")
  cat("\n")
  cat("DO NOT use the new function until differences are resolved.\n")
  cat("Using different MOFA input would require RETRAINING ALL MODELS (hours).\n")
}

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n")
