# Verify key numbers for figure legend
library(dplyr)
source("scripts/core/data_loading.R")

# 1. Load lipidomic data using the same function as the workflow
lipidomic <- load_lipidomic_data(
  file_path = "data/metabolomics/Lipid_data_file_JPC_cleaned.csv",
  log_transform = TRUE,
  verbose = TRUE
)
cat("Lipidomic after loading:", nrow(lipidomic), "lipids x", ncol(lipidomic), "samples\n")
cat("Sample names (first 5):", paste(head(colnames(lipidomic), 5), collapse=", "), "\n")

# Match to RNA samples
rna_meta <- read.csv("data/RNA/meta_data_df.csv")
rna_mice <- unique(gsub("-", "_", rna_meta$sample_name))
cat("RNA mice (first 5):", paste(head(rna_mice, 5), collapse=", "), "\n")

# Try both formats
keep1 <- intersect(colnames(lipidomic), rna_mice)
cat("Direct match:", length(keep1), "\n")

# Try with dash
rna_mice2 <- unique(rna_meta$sample_name)
keep2 <- intersect(colnames(lipidomic), rna_mice2)
cat("Dash match:", length(keep2), "\n")

# Use whichever works
keep <- if (length(keep1) > length(keep2)) keep1 else keep2
cat("Using:", length(keep), "samples\n")

if (length(keep) > 0) {
  lip_filtered <- lipidomic[, keep]
} else {
  # Just use all samples
  lip_filtered <- lipidomic
  cat("Using all", ncol(lipidomic), "samples for PCA check\n")
}

# PCA
pca <- prcomp(t(lip_filtered), center = TRUE)
var_explained <- summary(pca)$importance[2, 1:3] * 100
cat("\nPCA variance explained:\n")
cat("  PC1:", round(var_explained[1], 2), "%\n")
cat("  PC2:", round(var_explained[2], 2), "%\n")
cat("  PC3:", round(var_explained[3], 2), "%\n")

# 2. DE results
load("results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData")
ttop_gonadal <- ttop_full_list[["Testes vs Ovaries"]][["Testes vs Ovaries"]]
cat("\nLipids in DE analysis:", nrow(ttop_gonadal), "\n")
n_sig <- sum(ttop_gonadal$adj.P.Val <= 0.05)
cat("Significant (FDR <= 0.05):", n_sig, "\n")
cat("Significant fraction:", round(n_sig / nrow(ttop_gonadal) * 100, 1), "%\n")

# 3. Enrichment
enrich <- read.csv("results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv")
cat("\nLipid classes tested:", nrow(enrich), "\n")
sig_classes <- enrich[enrich$p_value < 0.05, ]
cat("Significant classes:", nrow(sig_classes), "\n")
cat("  Male-biased:", sum(sig_classes$enrichment_score > 0), "\n")
cat("  Female-biased:", sum(sig_classes$enrichment_score < 0), "\n")
