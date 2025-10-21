#' Test Data Loading Functions
#'
#' This script validates that the new modular data loading functions in
#' core/data_loading.R produce identical results to the original scripts.

library(readr)
library(readxl)
library(dplyr)

# Source the new modules
source("scripts/core/data_loading.R")
source("scripts/config/datasets.R")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Testing Data Loading Functions\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Test 1: RNA Data Loading
# =========================
cat("Test 1: RNA Data Loading (FCG Main Dataset)\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Old way (from pseudobulk_DA.R)
cat("  Loading RNA data using OLD method...\n")
scaled_RNA_1_old <- as.data.frame(read_csv("data/RNA/scaled_RNA_1.csv", show_col_types = FALSE))
scaled_RNA_2_old <- as.data.frame(read_csv("data/RNA/scaled_RNA_2.csv", show_col_types = FALSE))
scaled_RNA_3_old <- as.data.frame(read_csv("data/RNA/scaled_RNA_3.csv", show_col_types = FALSE))

scaled_RNA_old <- merge(scaled_RNA_1_old, scaled_RNA_2_old, by = "gene")
scaled_RNA_old <- merge(scaled_RNA_old, scaled_RNA_3_old, by = "gene")
row.names(scaled_RNA_old) <- scaled_RNA_old$gene
scaled_RNA_old <- scaled_RNA_old[, -1]

rm(scaled_RNA_1_old, scaled_RNA_2_old, scaled_RNA_3_old)

cat("    Old method: ", nrow(scaled_RNA_old), "genes x", ncol(scaled_RNA_old), "samples\n")

# New way (using load_rna_data)
cat("  Loading RNA data using NEW method...\n")
scaled_RNA_new <- load_rna_data(verbose = FALSE)
cat("    New method: ", nrow(scaled_RNA_new), "genes x", ncol(scaled_RNA_new), "samples\n")

# Compare
cat("\n  Comparing results...\n")
identical_dims <- identical(dim(scaled_RNA_old), dim(scaled_RNA_new))
identical_rownames <- identical(row.names(scaled_RNA_old), row.names(scaled_RNA_new))
identical_colnames <- identical(colnames(scaled_RNA_old), colnames(scaled_RNA_new))
identical_values <- all.equal(scaled_RNA_old, scaled_RNA_new)

cat("    Dimensions match:  ", identical_dims, "\n")
cat("    Row names match:   ", identical_rownames, "\n")
cat("    Column names match:", identical_colnames, "\n")
cat("    Values match:      ", isTRUE(identical_values), "\n")

if (!isTRUE(identical_values)) {
  cat("    Difference details:", identical_values, "\n")
}

test1_pass <- identical_dims && identical_rownames && identical_colnames && isTRUE(identical_values)
if (test1_pass) {
  cat("  ✓ Test 1 PASSED: RNA data loading produces identical results\n")
} else {
  cat("  ✗ Test 1 FAILED: Results differ!\n")
}

# Clean up
rm(scaled_RNA_old, scaled_RNA_new)


# Test 2: Lipidomic Data Loading
# ===============================
cat("\n")
cat("Test 2: Lipidomic Data Loading\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Old way (from lipidomic_DA.R)
cat("  Loading lipidomic data using OLD method...\n")
lipidomic_old <- as.data.frame(read_csv("data/metabolomics/Lipid_data_file_JPC_cleaned.csv", show_col_types = FALSE))
row.names(lipidomic_old) <- lipidomic_old$`Unnamed: 0`
lipidomic_old <- as.data.frame(t(lipidomic_old[, -1]))
lipidomic_old <- log2(lipidomic_old)
names(lipidomic_old) <- gsub("-0", "-", names(lipidomic_old))
names(lipidomic_old) <- gsub("-", "_", names(lipidomic_old))

cat("    Old method: ", nrow(lipidomic_old), "lipids x", ncol(lipidomic_old), "samples\n")

# New way (using load_lipidomic_data)
cat("  Loading lipidomic data using NEW method...\n")
lipidomic_new <- load_lipidomic_data(verbose = FALSE)
cat("    New method: ", nrow(lipidomic_new), "lipids x", ncol(lipidomic_new), "samples\n")

# Compare
cat("\n  Comparing results...\n")
identical_dims <- identical(dim(lipidomic_old), dim(lipidomic_new))
identical_rownames <- identical(row.names(lipidomic_old), row.names(lipidomic_new))
identical_colnames <- identical(sort(colnames(lipidomic_old)), sort(colnames(lipidomic_new)))
# Sort columns for comparison (order may differ)
lipidomic_old_sorted <- lipidomic_old[, sort(colnames(lipidomic_old))]
lipidomic_new_sorted <- lipidomic_new[, sort(colnames(lipidomic_new))]
identical_values <- all.equal(lipidomic_old_sorted, lipidomic_new_sorted)

cat("    Dimensions match:  ", identical_dims, "\n")
cat("    Row names match:   ", identical_rownames, "\n")
cat("    Column names match:", identical_colnames, "\n")
cat("    Values match:      ", isTRUE(identical_values), "\n")

if (!isTRUE(identical_values)) {
  cat("    Difference details:", identical_values, "\n")
}

test2_pass <- identical_dims && identical_rownames && identical_colnames && isTRUE(identical_values)
if (test2_pass) {
  cat("  ✓ Test 2 PASSED: Lipidomic data loading produces identical results\n")
} else {
  cat("  ✗ Test 2 FAILED: Results differ!\n")
}

# Clean up
rm(lipidomic_old, lipidomic_new, lipidomic_old_sorted, lipidomic_new_sorted)


# Test 3: Metadata Loading (for pseudobulk RNA analysis)
# =======================================================
cat("\n")
cat("Test 3: Metadata Loading (with RNA metadata)\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Old way (from pseudobulk_DA.R)
cat("  Loading metadata using OLD method...\n")
Suppl_Table_1_OverviewData <- as.data.frame(read_excel("data/Suppl_Table_1_OverviewData.xlsx"))
Suppl_Table_1_OverviewData <- Suppl_Table_1_OverviewData[, -4]
meta_data_df_old <- as.data.frame(read_csv("data/RNA/meta_data_df.csv", show_col_types = FALSE))
meta_data_df_old <- merge(meta_data_df_old, Suppl_Table_1_OverviewData,
                          by.x = "sample_name", by.y = "sequencing_name")

cat("    Old method: ", nrow(meta_data_df_old), "samples\n")

# New way (using load_metadata)
cat("  Loading metadata using NEW method...\n")
meta_data_df_new <- load_metadata(rna_metadata_file = "data/RNA/meta_data_df.csv",
                                   verbose = FALSE)
cat("    New method: ", nrow(meta_data_df_new), "samples\n")

# Compare
cat("\n  Comparing results...\n")
identical_dims <- identical(dim(meta_data_df_old), dim(meta_data_df_new))
identical_samples <- identical(sort(meta_data_df_old$sample), sort(meta_data_df_new$sample))

# Column names might differ in order, so check sets
cols_old <- sort(colnames(meta_data_df_old))
cols_new <- sort(colnames(meta_data_df_new))
identical_colnames <- identical(cols_old, cols_new)

cat("    Dimensions match:  ", identical_dims, "\n")
cat("    Samples match:     ", identical_samples, "\n")
cat("    Column names match:", identical_colnames, "\n")

test3_pass <- identical_dims && identical_samples && identical_colnames
if (test3_pass) {
  cat("  ✓ Test 3 PASSED: Metadata loading produces identical results\n")
} else {
  cat("  ✗ Test 3 FAILED: Results differ!\n")
}

# Clean up
rm(meta_data_df_old, meta_data_df_new, Suppl_Table_1_OverviewData)


# Test 4: Metadata Loading (for lipidomic analysis)
# ==================================================
cat("\n")
cat("Test 4: Metadata Loading (overview only, for lipidomics)\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Old way (from lipidomic_DA.R)
cat("  Loading metadata using OLD method...\n")
Suppl_Table_1_OverviewData <- as.data.frame(read_excel("data/Suppl_Table_1_OverviewData.xlsx"))
targets_all_old <- Suppl_Table_1_OverviewData[, c(2, 3, 4, 5)]
names(targets_all_old)[1] <- "sample"

cat("    Old method: ", nrow(targets_all_old), "samples\n")

# New way (using load_metadata without RNA metadata)
cat("  Loading metadata using NEW method...\n")
targets_all_new <- load_metadata(rna_metadata_file = NULL, verbose = FALSE)
cat("    New method: ", nrow(targets_all_new), "samples\n")

# Compare
cat("\n  Comparing results...\n")
identical_dims <- identical(dim(targets_all_old), dim(targets_all_new))
identical_samples <- identical(targets_all_old$sample, targets_all_new$sample)
identical_colnames <- identical(colnames(targets_all_old), colnames(targets_all_new))
identical_values <- all.equal(targets_all_old, targets_all_new)

cat("    Dimensions match:  ", identical_dims, "\n")
cat("    Samples match:     ", identical_samples, "\n")
cat("    Column names match:", identical_colnames, "\n")
cat("    Values match:      ", isTRUE(identical_values), "\n")

test4_pass <- identical_dims && identical_samples && identical_colnames && isTRUE(identical_values)
if (test4_pass) {
  cat("  ✓ Test 4 PASSED: Lipidomic metadata loading produces identical results\n")
} else {
  cat("  ✗ Test 4 FAILED: Results differ!\n")
}


# Summary
# =======
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Summary of Test Results\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

all_tests_pass <- test1_pass && test2_pass && test3_pass && test4_pass

cat("  Test 1 (RNA data loading):           ", ifelse(test1_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 2 (Lipidomic data loading):     ", ifelse(test2_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 3 (Metadata with RNA):          ", ifelse(test3_pass, "✓ PASS", "✗ FAIL"), "\n")
cat("  Test 4 (Metadata for lipidomics):    ", ifelse(test4_pass, "✓ PASS", "✗ FAIL"), "\n")

cat("\n")
if (all_tests_pass) {
  cat("✓✓✓ ALL TESTS PASSED ✓✓✓\n")
  cat("The new data loading functions produce identical results to the original code.\n")
  cat("It is safe to use these functions in refactored workflows.\n")
} else {
  cat("✗✗✗ SOME TESTS FAILED ✗✗✗\n")
  cat("Review the differences above before using the new functions.\n")
}

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n")
