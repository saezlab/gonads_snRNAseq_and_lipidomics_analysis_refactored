#' Test Differential Analysis Functions
#'
#' This script validates that the new differential analysis functions in
#' core/differential_analysis.R produce identical results to the original code.

library(readr)
library(readxl)
library(dplyr)
library(limma)

# Source modules
source("scripts/core/data_loading.R")
source("scripts/core/differential_analysis.R")
source("scripts/config/contrasts.R")
source("scripts/support_functions.R")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Testing Differential Analysis Functions\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Load data
cat("Loading data...\n")
scaled_RNA <- load_rna_data(verbose = FALSE)
metadata <- load_metadata(rna_metadata_file = "data/RNA/meta_data_df.csv", verbose = FALSE)
cat("  Data loaded: ", nrow(scaled_RNA), "genes x", ncol(scaled_RNA), "samples\n\n")

# Test on Hepatocytes only (most important cell type)
cell_type <- "Hepatocytes"
cat("Testing on", cell_type, "...\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# Prepare data
targets <- metadata[metadata$cell_type == cell_type, ]
scaled_RNA_subset <- scaled_RNA[, targets$sample]

# Test 1: Gonadal Type Comparison
# ================================
cat("Test 1: Gonadal Type (Testes vs Ovaries)\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# OLD METHOD (from pseudobulk_DA.R)
cat("  Running OLD method...\n")
group_variable <- "gonadal_type"
targets_old <- targets[, c("sample", group_variable)]
names(targets_old)[2] <- "condition"
comparisons_old <- list("testes_vs_Ova" = c(-1, 2))

limmaRes_old <- runLimma(scaled_RNA_subset, targets_old, comparisons = comparisons_old)
ttop_old <- limma_res_to_ttop_list(
  limma_res = limmaRes_old,
  comp_names = names(comparisons_old),
  number = nrow(scaled_RNA_subset),
  adjust.method = "fdr"
)[[1]]

cat("    Old method: ", nrow(ttop_old), "features,",
    sum(ttop_old$adj.P.Val < 0.05, na.rm = TRUE), "significant (FDR < 0.05)\n")

# NEW METHOD (using our new functions)
cat("  Running NEW method...\n")
targets_new <- targets[, c("sample", group_variable)]
names(targets_new)[2] <- "condition"
comparisons_new <- build_contrast_list("gonadal_type")

ttop_new <- run_differential_analysis(
  data = scaled_RNA_subset,
  targets = targets_new,
  comparisons = comparisons_new,
  verbose = FALSE
)[[1]]

cat("    New method: ", nrow(ttop_new), "features,",
    sum(ttop_new$adj.P.Val < 0.05, na.rm = TRUE), "significant (FDR < 0.05)\n")

# Compare results
cat("\n  Comparing results...\n")

# Match rows by ID
ttop_old_sorted <- ttop_old[order(ttop_old$ID), ]
ttop_new_sorted <- ttop_new[order(ttop_new$ID), ]

# Check key columns
identical_ids <- identical(ttop_old_sorted$ID, ttop_new_sorted$ID)
identical_logfc <- all.equal(ttop_old_sorted$logFC, ttop_new_sorted$logFC, tolerance = 1e-6)
identical_t <- all.equal(ttop_old_sorted$t, ttop_new_sorted$t, tolerance = 1e-6)
identical_pval <- all.equal(ttop_old_sorted$P.Value, ttop_new_sorted$P.Value, tolerance = 1e-6)
identical_adjpval <- all.equal(ttop_old_sorted$adj.P.Val, ttop_new_sorted$adj.P.Val, tolerance = 1e-6)

cat("    IDs match:        ", identical_ids, "\n")
cat("    logFC match:      ", isTRUE(identical_logfc), "\n")
cat("    t-statistics match:", isTRUE(identical_t), "\n")
cat("    P-values match:   ", isTRUE(identical_pval), "\n")
cat("    Adj P-values match:", isTRUE(identical_adjpval), "\n")

test1_pass <- identical_ids && isTRUE(identical_logfc) && isTRUE(identical_t) &&
              isTRUE(identical_pval) && isTRUE(identical_adjpval)

if (test1_pass) {
  cat("  ✓ Test 1 PASSED: Gonadal type analysis produces identical results\n")
} else {
  cat("  ✗ Test 1 FAILED: Results differ!\n")
  if (!isTRUE(identical_logfc)) cat("    logFC diff:", identical_logfc, "\n")
  if (!isTRUE(identical_t)) cat("    t diff:", identical_t, "\n")
  if (!isTRUE(identical_pval)) cat("    P.Value diff:", identical_pval, "\n")
  if (!isTRUE(identical_adjpval)) cat("    adj.P.Val diff:", identical_adjpval, "\n")
}


# Test 2: Multi-Cell-Type Analysis
# =================================
cat("\n")
cat("Test 2: Multi-Cell-Type Analysis\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

cat("  Testing on 2 cell types (Hepatocytes + Endothelial cells)...\n")

# NEW METHOD (multi-cell-type function)
test_cell_types <- c("Hepatocytes", "Endothelial cells")
comparisons_multi <- build_contrast_list("gonadal_type")

results_multi <- run_differential_analysis_multi_celltype(
  data = scaled_RNA,
  metadata = metadata,
  group_variable = "gonadal_type",
  comparisons = comparisons_multi,
  cell_types = test_cell_types,
  verbose = FALSE
)

# Verify structure
cat("    Results structure:\n")
cat("      Cell types:", paste(names(results_multi), collapse = ", "), "\n")
for (ct in names(results_multi)) {
  cat("      ", ct, "- comparisons:", paste(names(results_multi[[ct]]), collapse = ", "), "\n")
  n_sig <- sum(results_multi[[ct]][[1]]$adj.P.Val < 0.05, na.rm = TRUE)
  cat("         Significant features:", n_sig, "\n")
}

# Check that Hepatocytes result matches Test 1
hepa_result_multi <- results_multi$Hepatocytes[[1]]
hepa_result_multi_sorted <- hepa_result_multi[order(hepa_result_multi$ID), ]

# Remove cell-type-specific t column name for comparison
names(hepa_result_multi_sorted)[4] <- "t"

identical_to_test1 <- all.equal(
  ttop_new_sorted[, c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
  hepa_result_multi_sorted[, c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
  tolerance = 1e-6
)

test2_pass <- isTRUE(identical_to_test1) &&
              length(results_multi) == 2 &&
              all(test_cell_types %in% names(results_multi))

if (test2_pass) {
  cat("  ✓ Test 2 PASSED: Multi-cell-type function works correctly\n")
} else {
  cat("  ✗ Test 2 FAILED\n")
  if (!isTRUE(identical_to_test1)) {
    cat("    Hepatocytes result differs from single-cell-type analysis\n")
    cat("    Difference:", identical_to_test1, "\n")
  }
}


# Test 3: Helper Functions
# =========================
cat("\n")
cat("Test 3: Helper Functions\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Test filter_significant_features
cat("  Testing filter_significant_features()...\n")
sig_features <- filter_significant_features(ttop_new, fdr_cutoff = 0.05)
n_sig_manual <- sum(ttop_new$adj.P.Val < 0.05, na.rm = TRUE)
n_sig_function <- nrow(sig_features)
cat("    Manual count:", n_sig_manual, "\n")
cat("    Function count:", n_sig_function, "\n")
test3a_pass <- n_sig_manual == n_sig_function

# Test ttop_list_to_stat_matrix
cat("  Testing ttop_list_to_stat_matrix()...\n")
ttop_list_test <- list("Comparison1" = ttop_new)
stat_matrix <- ttop_list_to_stat_matrix(ttop_list_test, stat_column = "t")
cat("    Matrix dimensions:", nrow(stat_matrix), "x", ncol(stat_matrix), "\n")
cat("    Expected:", nrow(ttop_new), "x 1\n")
test3b_pass <- nrow(stat_matrix) == nrow(ttop_new) && ncol(stat_matrix) == 1

test3_pass <- test3a_pass && test3b_pass

if (test3_pass) {
  cat("  ✓ Test 3 PASSED: Helper functions work correctly\n")
} else {
  cat("  ✗ Test 3 FAILED\n")
}


# Test 4: Contrast Configuration
# ===============================
cat("\n")
cat("Test 4: Contrast Configuration\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Test get_contrast
cat("  Testing get_contrast()...\n")
gonadal_contrast <- get_contrast("gonadal_type")
cat("    Name:", gonadal_contrast$name, "\n")
cat("    Vector:", paste(gonadal_contrast$vector, collapse = ", "), "\n")
test4a_pass <- gonadal_contrast$name == "Testes vs Ovaries" &&
               identical(gonadal_contrast$vector, c(-1, 2))

# Test build_contrast_list
cat("  Testing build_contrast_list()...\n")
comparisons_built <- build_contrast_list(c("gonadal_type", "allosome_type"))
cat("    Contrasts built:", paste(names(comparisons_built), collapse = ", "), "\n")
test4b_pass <- length(comparisons_built) == 2 &&
               "Testes vs Ovaries" %in% names(comparisons_built) &&
               "XY vs XX" %in% names(comparisons_built)

# Test get_contrast_colors
cat("  Testing get_contrast_colors()...\n")
colors <- get_contrast_colors("gonadal_type")
cat("    Colors:", paste(names(colors), "=", colors, collapse = ", "), "\n")
test4c_pass <- !is.null(colors) && "Ovaries" %in% names(colors) && "Testes" %in% names(colors)

test4_pass <- test4a_pass && test4b_pass && test4c_pass

if (test4_pass) {
  cat("  ✓ Test 4 PASSED: Contrast configuration works correctly\n")
} else {
  cat("  ✗ Test 4 FAILED\n")
}


# Summary
# =======
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Summary of Test Results\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

all_tests_pass <- test1_pass && test2_pass && test3_pass && test4_pass

cat("  Test 1 (Gonadal type DA):             ", ifelse(test1_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 2 (Multi-cell-type function):    ", ifelse(test2_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 3 (Helper functions):            ", ifelse(test3_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 4 (Contrast configuration):      ", ifelse(test4_pass, "✓ PASS", "✗ FAIL"), "\n")

cat("\n")
if (all_tests_pass) {
  cat("✓✓✓ ALL TESTS PASSED ✓✓✓\n")
  cat("The new differential analysis functions produce identical results to the original code.\n")
  cat("It is safe to use these functions in refactored workflows.\n")
} else {
  cat("✗✗✗ SOME TESTS FAILED ✗✗✗\n")
  cat("Review the differences above before using the new functions.\n")
}

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n")
