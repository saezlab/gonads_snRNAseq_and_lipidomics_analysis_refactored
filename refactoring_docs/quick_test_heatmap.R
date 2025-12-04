# QUICK TEST: Verify pheatmap namespace fix
# This creates a simple test heatmap to verify the fix works

cat("\n=== QUICK HEATMAP TEST ===\n\n")

# Load both packages (this creates the conflict)
library(pheatmap)
library(ComplexHeatmap)

cat("Both packages loaded (this causes the namespace conflict)\n\n")

# Create test data
set.seed(123)
test_matrix <- matrix(rnorm(200, mean = 0, sd = 5), nrow = 20, ncol = 10)
rownames(test_matrix) <- paste0("Pathway_", 1:20)
colnames(test_matrix) <- paste0("Contrast_", 1:10)

cat("Created test matrix:", nrow(test_matrix), "x", ncol(test_matrix), "\n\n")

# Test 1: Using unqualified pheatmap (BROKEN - uses ComplexHeatmap version)
cat("Test 1: Unqualified pheatmap() [BROKEN]\n")
test_file1 <- "test_unqualified.pdf"
if (file.exists(test_file1)) unlink(test_file1)

tryCatch({
  pheatmap(test_matrix, filename = test_file1, width = 8, height = 6)

  if (file.exists(test_file1)) {
    size1 <- file.info(test_file1)$size
    cat("  File created, size:", size1, "bytes\n")
    if (size1 < 1000) {
      cat("  ✗ CORRUPTED (size too small)\n")
    } else {
      cat("  ? May or may not work\n")
    }
  } else {
    cat("  ✗ File not created\n")
  }
}, error = function(e) {
  cat("  ✗ ERROR:", conditionMessage(e), "\n")
})

cat("\n")

# Test 2: Using explicit namespace (FIXED - uses original pheatmap)
cat("Test 2: pheatmap::pheatmap() [FIXED]\n")
test_file2 <- "test_qualified.pdf"
if (file.exists(test_file2)) unlink(test_file2)

tryCatch({
  pheatmap::pheatmap(test_matrix, filename = test_file2, width = 8, height = 6)

  if (file.exists(test_file2)) {
    size2 <- file.info(test_file2)$size
    cat("  File created, size:", size2, "bytes\n")
    if (size2 > 5000) {
      cat("  ✓ SUCCESS (valid PDF size)\n")
    } else if (size2 > 1000) {
      cat("  ? Possibly OK\n")
    } else {
      cat("  ✗ CORRUPTED (size too small)\n")
    }
  } else {
    cat("  ✗ File not created\n")
  }
}, error = function(e) {
  cat("  ✗ ERROR:", conditionMessage(e), "\n")
})

cat("\n")
cat("=== TEST COMPLETE ===\n")
cat("\nCompare file sizes:\n")
if (file.exists(test_file1)) {
  cat("  ", test_file1, ":", file.info(test_file1)$size, "bytes\n")
}
if (file.exists(test_file2)) {
  cat("  ", test_file2, ":", file.info(test_file2)$size, "bytes\n")
}

cat("\nTry opening both PDFs to see which one works!\n")
cat("The second one (test_qualified.pdf) should open correctly.\n\n")
