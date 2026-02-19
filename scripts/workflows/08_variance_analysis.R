################################################################################
# Workflow 8: Variance Analysis - Comparing Lipid Variability Between Groups
################################################################################
#
# This workflow compares the variance/variability of lipid abundances between
# mice with ovaries versus mice with testes.
#
# Question: Are lipid abundances more variable in male or female mice?
#
# Statistical approaches:
# 1. Levene's Test (Brown-Forsythe variant) - tests equality of variances
#    - Uses median (robust to outliers and non-normality)
#    - Applied per-lipid to identify specific lipids with differential variance
#
# 2. Coefficient of Variation (CV) - easy-to-interpret measure of relative variability
#    - CV = SD / |Mean| (expressed as percentage)
#    - Allows comparison even when group means differ
#    - Global paired test to assess systematic differences
#
# Output:
# - Per-lipid variance statistics (CV, Levene's p-value)
# - Global summary: Do lipids systematically have higher variance in one group?
# - Visualizations: CV comparison, volcano-style plot of differential variance
#
################################################################################

# Clear workspace
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)  # For labeling points
  library(car)  # For Levene's test
  library(readr)
})

################################################################################
# 1. CONFIGURATION
################################################################################

cat(strrep("=", 80), "\n")
cat("WORKFLOW 8: VARIANCE ANALYSIS\n")
cat("Comparing Lipid Variability: Ovaries vs Testes\n")
cat(strrep("=", 80), "\n\n")

# Load configuration and core functions
source("scripts/config/analysis_config.R")
source("scripts/config/plotting_themes.R")
source("scripts/core/data_loading.R")

# Define output directory
variance_dir <- "results/lipidomic/variance_analysis"
if (!dir.exists(variance_dir)) {
  dir.create(variance_dir, recursive = TRUE)
}

################################################################################
# 2. LOAD DATA
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("2. LOADING DATA\n")
cat(strrep("-", 80), "\n\n")

# Load lipidomic data (log2 transformed)
cat("Loading lipidomic data...\n")
lipidomic <- load_lipidomic_data(
  file_path = "data/metabolomics/Lipid_data_file_JPC_cleaned.csv",
  log_transform = TRUE,
  verbose = TRUE
)

# Load metadata
cat("\nLoading metadata...\n")
metadata <- load_metadata(verbose = TRUE)

# Filter metadata to samples with lipidomic data
metadata_lipid <- metadata %>%
  dplyr::filter(sample %in% colnames(lipidomic))

cat("\nData summary:\n")
cat("  Lipids:", nrow(lipidomic), "\n")
cat("  Samples:", ncol(lipidomic), "\n")
cat("  Gonadal types:", paste(unique(metadata_lipid$gonadal_type), collapse = ", "), "\n")

# Count samples per group
group_counts <- table(metadata_lipid$gonadal_type)
cat("\nSamples per group:\n")
for (g in names(group_counts)) {
  cat("  ", g, ":", group_counts[g], "\n")
}

################################################################################
# 3. SEPARATE SAMPLES BY GONADAL TYPE
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("3. SEPARATING SAMPLES BY GONADAL TYPE\n")
cat(strrep("-", 80), "\n\n")

# Get sample names for each group
ovaries_samples <- metadata_lipid %>%
  dplyr::filter(gonadal_type == "Ovaries") %>%
  dplyr::pull(sample)

testes_samples <- metadata_lipid %>%
  dplyr::filter(gonadal_type == "Testes") %>%
  dplyr::pull(sample)

cat("Ovaries samples:", length(ovaries_samples), "\n")
cat("Testes samples:", length(testes_samples), "\n")

# Subset lipidomic data
lipidomic_ovaries <- lipidomic[, ovaries_samples]
lipidomic_testes <- lipidomic[, testes_samples]

################################################################################
# 4. CALCULATE COEFFICIENT OF VARIATION (CV) FOR EACH LIPID
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("4. CALCULATING COEFFICIENT OF VARIATION (CV)\n")
cat(strrep("-", 80), "\n\n")

cat("CV = SD / |Mean| × 100 (expressed as percentage)\n")
cat("Higher CV indicates more relative variability.\n\n")

# Function to calculate CV (handles negative means from log-transformed data)
calc_cv <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA)
  mean_val <- mean(x)
  sd_val <- sd(x)
  if (abs(mean_val) < 1e-10) return(NA)  # Avoid division by near-zero
  return((sd_val / abs(mean_val)) * 100)
}

# Calculate CV for each lipid in each group
cv_results <- data.frame(
  lipid = rownames(lipidomic),
  cv_ovaries = apply(lipidomic_ovaries, 1, calc_cv),
  cv_testes = apply(lipidomic_testes, 1, calc_cv),
  mean_ovaries = rowMeans(lipidomic_ovaries, na.rm = TRUE),
  mean_testes = rowMeans(lipidomic_testes, na.rm = TRUE),
  sd_ovaries = apply(lipidomic_ovaries, 1, sd, na.rm = TRUE),
  sd_testes = apply(lipidomic_testes, 1, sd, na.rm = TRUE),
  var_ovaries = apply(lipidomic_ovaries, 1, var, na.rm = TRUE),
  var_testes = apply(lipidomic_testes, 1, var, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Calculate CV difference and ratio
cv_results <- cv_results %>%
  dplyr::mutate(
    cv_diff = cv_testes - cv_ovaries,  # Positive = more variable in testes
    cv_ratio = cv_testes / cv_ovaries,  # >1 = more variable in testes
    cv_log2_ratio = log2(cv_testes / cv_ovaries),  # Log2 ratio for symmetry
    var_ratio = var_testes / var_ovaries,
    var_log2_ratio = log2(var_testes / var_ovaries)
  )

cat("CV calculation complete.\n\n")

# Summary statistics
cat("CV Summary Statistics:\n")
cat("  Ovaries - Mean CV:", round(mean(cv_results$cv_ovaries, na.rm = TRUE), 2), "%\n")
cat("  Ovaries - Median CV:", round(median(cv_results$cv_ovaries, na.rm = TRUE), 2), "%\n")
cat("  Testes - Mean CV:", round(mean(cv_results$cv_testes, na.rm = TRUE), 2), "%\n")
cat("  Testes - Median CV:", round(median(cv_results$cv_testes, na.rm = TRUE), 2), "%\n")
cat("\n")
cat("  Lipids with higher CV in Testes:", sum(cv_results$cv_diff > 0, na.rm = TRUE), "\n")
cat("  Lipids with higher CV in Ovaries:", sum(cv_results$cv_diff < 0, na.rm = TRUE), "\n")

################################################################################
# 5. RUN LEVENE'S TEST (BROWN-FORSYTHE VARIANT) FOR EACH LIPID
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("5. RUNNING LEVENE'S TEST (BROWN-FORSYTHE VARIANT)\n")
cat(strrep("-", 80), "\n\n")

cat("Testing equality of variances using absolute deviations from group MEDIAN.\n")
cat("This is robust to non-normality and outliers.\n\n")

# Function to run Levene's test for one lipid
run_levene_test <- function(lipid_name, lipidomic_data, metadata) {
  # Get lipid values
  values <- as.numeric(lipidomic_data[lipid_name, ])

  # Create data frame for test
  test_df <- data.frame(
    value = values,
    group = factor(metadata$gonadal_type[match(colnames(lipidomic_data), metadata$sample)])
  )

  # Remove NA values

test_df <- test_df[!is.na(test_df$value) & !is.na(test_df$group), ]

  if (nrow(test_df) < 4) return(list(p_value = NA, statistic = NA))

  # Run Levene's test (center = "median" for Brown-Forsythe)
  tryCatch({
    levene_result <- car::leveneTest(value ~ group, data = test_df, center = "median")
    return(list(
      p_value = levene_result$`Pr(>F)`[1],
      statistic = levene_result$`F value`[1]
    ))
  }, error = function(e) {
    return(list(p_value = NA, statistic = NA))
  })
}

# Run Levene's test for all lipids
cat("Running Levene's test for", nrow(lipidomic), "lipids...\n")

levene_results <- lapply(rownames(lipidomic), function(lipid) {
  run_levene_test(lipid, lipidomic, metadata_lipid)
})

# Extract results
cv_results$levene_p <- sapply(levene_results, function(x) x$p_value)
cv_results$levene_F <- sapply(levene_results, function(x) x$statistic)

# Adjust p-values for multiple testing
cv_results$levene_padj <- p.adjust(cv_results$levene_p, method = "BH")

cat("Levene's test complete.\n\n")

# Count significant results
n_sig_nominal <- sum(cv_results$levene_p < 0.05, na.rm = TRUE)
n_sig_adj <- sum(cv_results$levene_padj < 0.05, na.rm = TRUE)

cat("Significant differential variance (Levene's test):\n")
cat("  Nominal p < 0.05:", n_sig_nominal, "lipids\n")
cat("  FDR-adjusted p < 0.05:", n_sig_adj, "lipids\n")

################################################################################
# 6. GLOBAL PAIRED TEST: SYSTEMATIC DIFFERENCE IN CV?
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("6. GLOBAL TEST: SYSTEMATIC DIFFERENCE IN VARIABILITY?\n")
cat(strrep("-", 80), "\n\n")

cat("Testing whether lipids systematically have higher CV in one group.\n")
cat("Using paired Wilcoxon signed-rank test (each lipid is a paired observation).\n\n")

# Paired Wilcoxon test on CVs
wilcox_cv <- wilcox.test(cv_results$cv_testes, cv_results$cv_ovaries,
                          paired = TRUE, alternative = "two.sided")

cat("Paired Wilcoxon signed-rank test (CV: Testes vs Ovaries):\n")
cat("  V statistic:", wilcox_cv$statistic, "\n")
cat("  p-value:", format.pval(wilcox_cv$p.value, digits = 4), "\n")

# Also test with paired t-test for comparison
ttest_cv <- t.test(cv_results$cv_testes, cv_results$cv_ovaries,
                   paired = TRUE, alternative = "two.sided")

cat("\nPaired t-test (CV: Testes vs Ovaries):\n")
cat("  t statistic:", round(ttest_cv$statistic, 3), "\n")
cat("  p-value:", format.pval(ttest_cv$p.value, digits = 4), "\n")
cat("  Mean CV difference (Testes - Ovaries):", round(ttest_cv$estimate, 3), "%\n")

# Interpretation
cat("\n")
cat(strrep("-", 40), "\n")
cat("INTERPRETATION:\n")
cat(strrep("-", 40), "\n")

if (ttest_cv$estimate > 0) {
  cat("Lipids tend to be MORE VARIABLE in TESTES (male) mice.\n")
} else {
  cat("Lipids tend to be MORE VARIABLE in OVARIES (female) mice.\n")
}

if (wilcox_cv$p.value < 0.05) {
  cat("This difference is STATISTICALLY SIGNIFICANT (p < 0.05).\n")
} else {
  cat("This difference is NOT statistically significant (p >= 0.05).\n")
}

################################################################################
# 7. IDENTIFY LIPIDS WITH SIGNIFICANT DIFFERENTIAL VARIANCE
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("7. IDENTIFYING LIPIDS WITH DIFFERENTIAL VARIANCE\n")
cat(strrep("-", 80), "\n\n")

# Add significance flags
cv_results <- cv_results %>%
  dplyr::mutate(
    sig_levene = levene_padj < 0.05,
    more_var_in = case_when(
      sig_levene & cv_diff > 0 ~ "Testes",
      sig_levene & cv_diff < 0 ~ "Ovaries",
      TRUE ~ "NS"
    )
  )

# Summary of differential variance
cat("Lipids with significantly different variance (FDR < 0.05):\n")
cat("  More variable in Testes:", sum(cv_results$more_var_in == "Testes", na.rm = TRUE), "\n")
cat("  More variable in Ovaries:", sum(cv_results$more_var_in == "Ovaries", na.rm = TRUE), "\n")
cat("  Not significant:", sum(cv_results$more_var_in == "NS", na.rm = TRUE), "\n")

# Top lipids with higher variance in Testes
cat("\nTop 10 lipids with higher variance in TESTES:\n")
top_testes_var <- cv_results %>%
  dplyr::filter(more_var_in == "Testes") %>%
  dplyr::arrange(levene_padj) %>%
  dplyr::slice_head(n = 10)

if (nrow(top_testes_var) > 0) {
  print(top_testes_var %>%
          dplyr::select(lipid, cv_ovaries, cv_testes, cv_diff, levene_padj) %>%
          dplyr::mutate(across(where(is.numeric), ~round(., 3))))
} else {
  cat("  None\n")
}

# Top lipids with higher variance in Ovaries
cat("\nTop 10 lipids with higher variance in OVARIES:\n")
top_ovaries_var <- cv_results %>%
  dplyr::filter(more_var_in == "Ovaries") %>%
  dplyr::arrange(levene_padj) %>%
  dplyr::slice_head(n = 10)

if (nrow(top_ovaries_var) > 0) {
  print(top_ovaries_var %>%
          dplyr::select(lipid, cv_ovaries, cv_testes, cv_diff, levene_padj) %>%
          dplyr::mutate(across(where(is.numeric), ~round(., 3))))
} else {
  cat("  None\n")
}

################################################################################
# 8. VISUALIZATIONS
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("8. CREATING VISUALIZATIONS\n")
cat(strrep("-", 80), "\n\n")

# Color palette
ovaries_color <- "#e54141"
testes_color <- "#4173d8"

#-------------------------------------------------------------------------------
# Lipid-level Volcano Plot
#-------------------------------------------------------------------------------
cat("Creating CV difference volcano plot...\n")

# Add label column for significant lipids
cv_results <- cv_results %>%
  dplyr::mutate(label = ifelse(sig_levene, lipid, NA_character_))

# Find the nominal p-value threshold corresponding to FDR = 0.05
# This is the largest nominal p-value among those with FDR < 0.05
fdr_threshold_pval <- cv_results %>%
  dplyr::filter(levene_padj < 0.05) %>%
  dplyr::pull(levene_p) %>%
  max()

# If no significant results, use a very small value (line won't show in plot range)
if (length(fdr_threshold_pval) == 0 || is.infinite(fdr_threshold_pval)) {
  fdr_threshold_pval <- 1e-10
}

p2 <- ggplot(cv_results, aes(x = cv_diff, y = -log10(levene_p))) +
  geom_hline(yintercept = -log10(fdr_threshold_pval), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = more_var_in), alpha = 0.6, size = 2) +
  geom_text_repel(
    aes(label = label),
    size = 3.5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 20,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c("Testes" = testes_color, "Ovaries" = ovaries_color, "NS" = "grey70"),
    name = "Higher variance in"
  ) +
  labs(
    title = "Differential Variance: Volcano Plot",
    subtitle = paste0("Levene's test (Brown-Forsythe); ", n_sig_adj, " lipids FDR < 0.05"),
    x = "CV Difference (Testes - Ovaries, %)",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(variance_dir, "cv_difference_volcano.pdf"), p2,
       width = 10, height = 7)

cat("\nLipid-level visualization saved to:", variance_dir, "\n")

################################################################################
# 9. LIPID CLASS-LEVEL VARIANCE ANALYSIS
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("9. LIPID CLASS-LEVEL VARIANCE ANALYSIS\n")
cat(strrep("-", 80), "\n\n")

cat("Pooling lipids into classes and repeating the variance analysis.\n")
cat("This reveals whether specific lipid classes show differential variability.\n\n")

# Extract lipid class from lipid names
# Pattern: everything before the first parenthesis
extract_lipid_class <- function(lipid_name) {
  gsub("\\(.*", "", lipid_name)
}

# Add lipid class to results
cv_results$lipid_class <- sapply(cv_results$lipid, extract_lipid_class)

cat("Lipid classes identified:", length(unique(cv_results$lipid_class)), "\n")

# Calculate class-level statistics by averaging across lipids within each class
# For each sample, calculate the mean abundance of all lipids in each class
cat("\nCalculating class-level abundances (mean across lipids per class)...\n")

# Get unique lipid classes
lipid_classes <- unique(cv_results$lipid_class)

# Create class-level abundance matrices
class_abundance_list <- lapply(lipid_classes, function(lclass) {
  # Get lipids in this class
  class_lipids <- cv_results$lipid[cv_results$lipid_class == lclass]

  if (length(class_lipids) < 2) return(NULL)

  # Calculate mean abundance across lipids for each sample
  class_data <- lipidomic[class_lipids, , drop = FALSE]
  class_mean <- colMeans(class_data, na.rm = TRUE)

  return(class_mean)
})
names(class_abundance_list) <- lipid_classes

# Remove NULL entries (classes with only 1 lipid)
class_abundance_list <- class_abundance_list[!sapply(class_abundance_list, is.null)]

# Convert to matrix
class_abundance <- do.call(rbind, class_abundance_list)

cat("  Classes with >= 2 lipids:", nrow(class_abundance), "\n")

# Separate by gonadal type
class_ovaries <- class_abundance[, ovaries_samples]
class_testes <- class_abundance[, testes_samples]

# Calculate CV for each class in each group
cv_class_results <- data.frame(
  lipid_class = rownames(class_abundance),
  n_lipids = sapply(rownames(class_abundance), function(lc) {
    sum(cv_results$lipid_class == lc)
  }),
  cv_ovaries = apply(class_ovaries, 1, calc_cv),
  cv_testes = apply(class_testes, 1, calc_cv),
  mean_ovaries = rowMeans(class_ovaries, na.rm = TRUE),
  mean_testes = rowMeans(class_testes, na.rm = TRUE),
  sd_ovaries = apply(class_ovaries, 1, sd, na.rm = TRUE),
  sd_testes = apply(class_testes, 1, sd, na.rm = TRUE),
  var_ovaries = apply(class_ovaries, 1, var, na.rm = TRUE),
  var_testes = apply(class_testes, 1, var, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Calculate differences
cv_class_results <- cv_class_results %>%
  dplyr::mutate(
    cv_diff = cv_testes - cv_ovaries,
    cv_ratio = cv_testes / cv_ovaries
  )

# Run Levene's test for each class
cat("Running Levene's test for each lipid class...\n")

run_levene_class <- function(class_name, class_abundance, metadata) {
  values <- as.numeric(class_abundance[class_name, ])

  test_df <- data.frame(
    value = values,
    group = factor(metadata$gonadal_type[match(colnames(class_abundance), metadata$sample)])
  )

  test_df <- test_df[!is.na(test_df$value) & !is.na(test_df$group), ]

  if (nrow(test_df) < 4) return(list(p_value = NA, statistic = NA))

  tryCatch({
    levene_result <- car::leveneTest(value ~ group, data = test_df, center = "median")
    return(list(
      p_value = levene_result$`Pr(>F)`[1],
      statistic = levene_result$`F value`[1]
    ))
  }, error = function(e) {
    return(list(p_value = NA, statistic = NA))
  })
}

levene_class_results <- lapply(rownames(class_abundance), function(lc) {
  run_levene_class(lc, class_abundance, metadata_lipid)
})

cv_class_results$levene_p <- sapply(levene_class_results, function(x) x$p_value)
cv_class_results$levene_F <- sapply(levene_class_results, function(x) x$statistic)
cv_class_results$levene_padj <- p.adjust(cv_class_results$levene_p, method = "BH")

# Add significance flags
# Use FDR-adjusted p-values for coloring (conservative approach)
cv_class_results <- cv_class_results %>%
  dplyr::mutate(
    sig_levene = levene_padj < 0.05,
    sig_nominal = levene_p < 0.05,
    more_var_in = case_when(
      sig_levene & cv_diff > 0 ~ "Testes",
      sig_levene & cv_diff < 0 ~ "Ovaries",
      TRUE ~ "NS"
    ),
    label = ifelse(sig_levene, lipid_class, NA_character_)
  )

# Sort by CV difference
cv_class_results <- cv_class_results %>%
  dplyr::arrange(cv_diff)

cat("\nLipid Class Variance Summary:\n")
cat(strrep("-", 60), "\n")

# Print results table
print_table <- cv_class_results %>%
  dplyr::select(lipid_class, n_lipids, cv_ovaries, cv_testes, cv_diff, levene_p, levene_padj) %>%
  dplyr::mutate(
    cv_ovaries = round(cv_ovaries, 2),
    cv_testes = round(cv_testes, 2),
    cv_diff = round(cv_diff, 2),
    levene_p = round(levene_p, 4),
    levene_padj = round(levene_padj, 4)
  )

print(as.data.frame(print_table), row.names = FALSE)

# Global test for class-level
cat("\n")
cat("Global test (class-level):\n")
wilcox_class <- wilcox.test(cv_class_results$cv_testes, cv_class_results$cv_ovaries,
                             paired = TRUE, alternative = "two.sided")
ttest_class <- t.test(cv_class_results$cv_testes, cv_class_results$cv_ovaries,
                      paired = TRUE, alternative = "two.sided")

cat("  Paired Wilcoxon p-value:", format.pval(wilcox_class$p.value, digits = 4), "\n")
cat("  Paired t-test p-value:", format.pval(ttest_class$p.value, digits = 4), "\n")
cat("  Mean CV difference:", round(ttest_class$estimate, 3), "%\n")

# Count
n_class_sig <- sum(cv_class_results$sig_levene, na.rm = TRUE)
n_class_nominal <- sum(cv_class_results$sig_nominal, na.rm = TRUE)

cat("\nClasses with differential variance:\n")
cat("  Nominal p < 0.05:", n_class_nominal, "\n")
cat("  FDR-adjusted p < 0.05:", n_class_sig, "\n")

#-------------------------------------------------------------------------------
# Class-level Visualizations
#-------------------------------------------------------------------------------
cat("\nCreating class-level visualizations...\n")

# Find the nominal p-value threshold corresponding to FDR = 0.05 for classes
fdr_threshold_pval_class <- cv_class_results %>%
  dplyr::filter(levene_padj < 0.05) %>%
  dplyr::pull(levene_p) %>%
  max()

# If no significant results, use a very small value
if (length(fdr_threshold_pval_class) == 0 || is.infinite(fdr_threshold_pval_class)) {
  fdr_threshold_pval_class <- 1e-10
}

# Volcano plot for classes
p_class_volcano <- ggplot(cv_class_results, aes(x = cv_diff, y = -log10(levene_p))) +
  geom_hline(yintercept = -log10(fdr_threshold_pval_class), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = more_var_in, size = n_lipids), alpha = 0.7) +
  geom_text_repel(
    aes(label = label),
    size = 3.5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 30,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c("Testes" = testes_color, "Ovaries" = ovaries_color, "NS" = "grey70"),
    name = "Higher variance in\n(FDR < 0.05)"
  ) +
  scale_size_continuous(name = "# Lipids", range = c(2, 8)) +
  labs(
    title = "Differential Variance by Lipid Class",
    subtitle = paste0("Levene's test; ", n_class_sig, " classes FDR < 0.05"),
    x = "CV Difference (Testes - Ovaries, %)",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(variance_dir, "cv_difference_volcano_by_class.pdf"), p_class_volcano,
       width = 12, height = 8)

cat("Class-level visualization saved.\n")

################################################################################
# 10. SAVE RESULTS
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("10. SAVING RESULTS\n")
cat(strrep("-", 80), "\n\n")

# Save full results table (lipid-level)
write_csv(cv_results, file.path(variance_dir, "variance_analysis_results.csv"))
cat("Saved: variance_analysis_results.csv\n")

# Save significant differential variance lipids
sig_var_lipids <- cv_results %>%
  dplyr::filter(sig_levene) %>%
  dplyr::arrange(levene_padj)

write_csv(sig_var_lipids, file.path(variance_dir, "significant_differential_variance_lipids.csv"))
cat("Saved: significant_differential_variance_lipids.csv\n")

# Save class-level results
write_csv(cv_class_results, file.path(variance_dir, "variance_analysis_by_class.csv"))
cat("Saved: variance_analysis_by_class.csv\n")

# Save summary statistics
summary_stats <- data.frame(
  Metric = c(
    "--- LIPID-LEVEL ANALYSIS ---",
    "Total lipids analyzed",
    "Ovaries samples (n)",
    "Testes samples (n)",
    "Mean CV - Ovaries (%)",
    "Mean CV - Testes (%)",
    "Median CV - Ovaries (%)",
    "Median CV - Testes (%)",
    "Mean CV difference (Testes - Ovaries)",
    "Median CV difference (Testes - Ovaries)",
    "Paired Wilcoxon p-value",
    "Paired t-test p-value",
    "Lipids with sig. diff. variance (FDR < 0.05)",
    "  - More variable in Testes",
    "  - More variable in Ovaries",
    "--- CLASS-LEVEL ANALYSIS ---",
    "Total lipid classes analyzed",
    "Mean CV - Ovaries (class-level, %)",
    "Mean CV - Testes (class-level, %)",
    "Paired Wilcoxon p-value (class-level)",
    "Paired t-test p-value (class-level)",
    "Classes with sig. diff. variance (FDR < 0.05)",
    "Classes with nominal sig. (p < 0.05)"
  ),
  Value = c(
    "",
    nrow(cv_results),
    length(ovaries_samples),
    length(testes_samples),
    round(mean(cv_results$cv_ovaries, na.rm = TRUE), 3),
    round(mean(cv_results$cv_testes, na.rm = TRUE), 3),
    round(median(cv_results$cv_ovaries, na.rm = TRUE), 3),
    round(median(cv_results$cv_testes, na.rm = TRUE), 3),
    round(ttest_cv$estimate, 3),
    round(median(cv_results$cv_diff, na.rm = TRUE), 3),
    format.pval(wilcox_cv$p.value, digits = 4),
    format.pval(ttest_cv$p.value, digits = 4),
    n_sig_adj,
    sum(cv_results$more_var_in == "Testes", na.rm = TRUE),
    sum(cv_results$more_var_in == "Ovaries", na.rm = TRUE),
    "",
    nrow(cv_class_results),
    round(mean(cv_class_results$cv_ovaries, na.rm = TRUE), 3),
    round(mean(cv_class_results$cv_testes, na.rm = TRUE), 3),
    format.pval(wilcox_class$p.value, digits = 4),
    format.pval(ttest_class$p.value, digits = 4),
    n_class_sig,
    n_class_nominal
  )
)

write_csv(summary_stats, file.path(variance_dir, "variance_analysis_summary.csv"))
cat("Saved: variance_analysis_summary.csv\n")

################################################################################
# 11. FINAL SUMMARY
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("VARIANCE ANALYSIS COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("KEY FINDINGS:\n")
cat(strrep("-", 40), "\n")

cat("LIPID-LEVEL ANALYSIS:\n")
mean_cv_diff <- ttest_cv$estimate
if (mean_cv_diff > 0) {
  cat("1. Lipids are systematically MORE VARIABLE in TESTES (male) mice\n")
  cat("   Mean CV difference: +", round(mean_cv_diff, 2), "% (Testes > Ovaries)\n\n", sep = "")
} else {
  cat("1. Lipids are systematically MORE VARIABLE in OVARIES (female) mice\n")
  cat("   Mean CV difference: ", round(mean_cv_diff, 2), "% (Ovaries > Testes)\n\n", sep = "")
}

cat("2. Statistical significance (lipid-level):\n")
cat("   - Paired Wilcoxon test p-value:", format.pval(wilcox_cv$p.value, digits = 4), "\n")
cat("   - Paired t-test p-value:", format.pval(ttest_cv$p.value, digits = 4), "\n\n")

cat("3. Per-lipid differential variance (Levene's test):\n")
cat("   - ", n_sig_adj, " lipids show significantly different variance (FDR < 0.05)\n", sep = "")
cat("   - ", sum(cv_results$more_var_in == "Testes", na.rm = TRUE),
    " more variable in Testes\n", sep = "")
cat("   - ", sum(cv_results$more_var_in == "Ovaries", na.rm = TRUE),
    " more variable in Ovaries\n\n", sep = "")

cat("CLASS-LEVEL ANALYSIS:\n")
mean_cv_diff_class <- ttest_class$estimate
if (mean_cv_diff_class > 0) {
  cat("4. At class level: lipid classes tend to be MORE VARIABLE in TESTES\n")
  cat("   Mean CV difference: +", round(mean_cv_diff_class, 2), "%\n\n", sep = "")
} else {
  cat("4. At class level: lipid classes tend to be MORE VARIABLE in OVARIES\n")
  cat("   Mean CV difference: ", round(mean_cv_diff_class, 2), "%\n\n", sep = "")
}

cat("5. Statistical significance (class-level):\n")
cat("   - Paired Wilcoxon test p-value:", format.pval(wilcox_class$p.value, digits = 4), "\n")
cat("   - Paired t-test p-value:", format.pval(ttest_class$p.value, digits = 4), "\n\n")

cat("6. Per-class differential variance (Levene's test):\n")
cat("   - ", n_class_nominal, " classes show nominally different variance (p < 0.05)\n", sep = "")
cat("   - ", n_class_sig, " classes significant after FDR correction\n\n", sep = "")

cat("OUTPUT FILES:\n")
cat("  Data files:\n")
cat("    - ", file.path(variance_dir, "variance_analysis_results.csv"), "\n", sep = "")
cat("    - ", file.path(variance_dir, "significant_differential_variance_lipids.csv"), "\n", sep = "")
cat("    - ", file.path(variance_dir, "variance_analysis_by_class.csv"), "\n", sep = "")
cat("    - ", file.path(variance_dir, "variance_analysis_summary.csv"), "\n", sep = "")
cat("  Plots:\n")
cat("    - ", file.path(variance_dir, "cv_difference_volcano.pdf"), "\n", sep = "")
cat("    - ", file.path(variance_dir, "cv_difference_volcano_by_class.pdf"), "\n", sep = "")

cat("\n")
cat(strrep("=", 80), "\n\n")
