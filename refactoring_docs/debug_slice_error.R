# Debug script for the slice() error

library(dplyr)
library(tidyr)
library(tibble)

cat("Loading pathway results...\n")
load("results/pseudobulk/pathways/pathway_scores_df_all.RData")

cat("Checking data structure:\n")
cat("  Class:", class(pathway_results), "\n")
cat("  Dimensions:", nrow(pathway_results), "x", ncol(pathway_results), "\n")
cat("  Column names:", paste(colnames(pathway_results), collapse = ", "), "\n")
cat("  Column classes:\n")
for (col in colnames(pathway_results)) {
  cat("    ", col, ":", class(pathway_results[[col]]), "\n")
}

cat("\nFirst few rows:\n")
print(head(pathway_results))

# Filter to gonadal contrasts
cat("\n\nFiltering to gonadal contrasts...\n")
pathway_gonadal <- pathway_results %>%
  filter(grepl("Testes_vs_Ovaries_tval", condition))

cat("After filter:\n")
cat("  Class:", class(pathway_gonadal), "\n")
cat("  Dimensions:", nrow(pathway_gonadal), "x", ncol(pathway_gonadal), "\n")

# Test the problematic pattern
cat("\n\nTesting pattern matching...\n")
pattern <- "RIBOSOM"
cat("Pattern:", pattern, "\n")

matches <- pathway_gonadal %>%
  filter(grepl(pattern, source, ignore.case = TRUE))

cat("Matches found:", nrow(matches), "\n")
cat("Matches class:", class(matches), "\n")

if (nrow(matches) > 0) {
  cat("First match:\n")
  print(head(matches, 1))

  cat("\n\nTrying arrange...\n")
  tryCatch({
    arranged <- matches %>%
      arrange(desc(abs(score)))
    cat("Arrange successful!\n")
    cat("Arranged class:", class(arranged), "\n")

    cat("\n\nTrying slice...\n")
    tryCatch({
      sliced <- arranged %>% slice(1)
      cat("Slice successful!\n")
      print(sliced)
    }, error = function(e) {
      cat("ERROR in slice:\n")
      print(e)

      cat("\n\nTrying alternative: convert to data.frame first\n")
      arranged_df <- as.data.frame(arranged)
      sliced_df <- arranged_df[1, , drop = FALSE]
      cat("Alternative successful!\n")
      print(sliced_df)
    })

  }, error = function(e) {
    cat("ERROR in arrange:\n")
    print(e)
  })
}

cat("\n\n=== TESTING FIX ===\n")
cat("Method 1: Using as.data.frame() before slice\n")
tryCatch({
  best_match <- matches %>%
    as.data.frame() %>%
    arrange(desc(abs(score))) %>%
    slice(1)
  cat("SUCCESS with as.data.frame()!\n")
  print(best_match)
}, error = function(e) {
  cat("FAILED:\n")
  print(e)
})

cat("\nMethod 2: Using indexing instead of slice\n")
tryCatch({
  arranged <- matches %>%
    arrange(desc(abs(score)))
  best_match <- arranged[1, , drop = FALSE]
  cat("SUCCESS with indexing!\n")
  print(best_match)
}, error = function(e) {
  cat("FAILED:\n")
  print(e)
})

cat("\nMethod 3: Using head(1)\n")
tryCatch({
  best_match <- matches %>%
    arrange(desc(abs(score))) %>%
    head(1)
  cat("SUCCESS with head()!\n")
  print(best_match)
}, error = function(e) {
  cat("FAILED:\n")
  print(e)
})

cat("\nMethod 4: Using dplyr::slice explicitly\n")
tryCatch({
  best_match <- matches %>%
    arrange(desc(abs(score))) %>%
    dplyr::slice(1)
  cat("SUCCESS with dplyr::slice()!\n")
  print(best_match)
}, error = function(e) {
  cat("FAILED:\n")
  print(e)
})

cat("\n=== DEBUG COMPLETE ===\n")
