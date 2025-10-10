library(ggplot2)

all_lipid_class_files <- list.files("data/metabolomics/for_bubbles_with_names/")

all_lipid_class_files <- paste("data/metabolomics/for_bubbles_with_names/",all_lipid_class_files,sep = "")

df_list <- list()

for(lipid_class_file in all_lipid_class_files)
{
  df <- read.csv(lipid_class_file)
  
  lipid_class <- gsub("_with_name","",gsub("data/.*/","",gsub(".csv$","",lipid_class_file)))
  
  df$lipid_class <- lipid_class
  
  df <- df[,c(4,6,1,2,3)]
  
  df_list[[lipid_class]] <- df
}

full_df <- do.call(rbind, df_list)

######

#get the correlation network and check which correlations are associated with structural features
load("results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata")
sum(row.names(hepatocyte_cor_mat_RNA_lip) %in% full_df$lipid)

full_df <- full_df[full_df$lipid %in% row.names(hepatocyte_cor_mat_RNA_lip),]

hepatocyte_cor_mat_RNA_lip <- hepatocyte_cor_mat_RNA_lip[full_df$lipid,]

row.names(full_df) <- full_df$lipid
full_df <- full_df[row.names(hepatocyte_cor_mat_RNA_lip),]

write_csv(hepatocyte_cor_mat_RNA_lip, file = "results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv")
write_csv(full_df, file = "results/multi-omic/structural_model/lipids_with_structure.csv")

############

# Load required libraries
library(tidyverse)
library(broom)

# Read data
lipid_structure <- read_csv("results/multi-omic/structural_model/lipids_with_structure.csv")
gene_lipid_cor <- read_csv("results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv")

# Ensure row order matches
stopifnot(nrow(lipid_structure) == nrow(gene_lipid_cor))

# Prepare lipid data
lipid_data <- lipid_structure %>%
  select(lipid, lipid_class, C, DB, LogFC) %>%
  rename(chain_length = C,
         db_position = DB,
         lipid_abundance = LogFC)

# Get gene names
gene_names <- colnames(gene_lipid_cor)

# Function to run regression with lipid class and test class-specific effects
run_regression_with_class <- function(gene_name, lipid_df, cor_matrix) {
  # Get correlation values for this gene
  gene_cors <- cor_matrix[[gene_name]]
  
  # Create modeling dataframe
  model_df <- lipid_df %>%
    mutate(gene_correlation = gene_cors) %>%
    filter(!is.na(gene_correlation))
  
  # Check if we have enough data
  if (nrow(model_df) < 10) {
    return(NULL)
  }
  
  # Fit main model with lipid class
  model <- tryCatch({
    lm(gene_correlation ~ lipid_class + lipid_abundance + chain_length + 
         db_position + lipid_abundance:chain_length + 
         lipid_abundance:db_position,
       data = model_df)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(model)) {
    return(NULL)
  }
  
  # Extract main model results
  results <- tidy(model) %>%
    mutate(gene = gene_name)
  
  model_stats <- glance(model) %>%
    mutate(gene = gene_name)
  
  # Now test class-specific effects for ALL parameters
  # For each lipid class, run a separate model
  class_specific_results <- list()
  
  lipid_classes <- unique(model_df$lipid_class)
  
  for (lclass in lipid_classes) {
    # Subset to this class
    class_df <- model_df %>% filter(lipid_class == lclass)
    
    # Need at least 10 observations
    if (nrow(class_df) < 10) next
    
    # Fit model for this class only
    class_model <- tryCatch({
      lm(gene_correlation ~ lipid_abundance + chain_length + db_position + 
           lipid_abundance:chain_length + lipid_abundance:db_position,
         data = class_df)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(class_model)) {
      class_results <- tidy(class_model) %>%
        filter(term != "(Intercept)") %>%  # Exclude intercept
        mutate(
          gene = gene_name,
          lipid_class = lclass,
          n_lipids = nrow(class_df)
        )
      
      class_specific_results[[lclass]] <- class_results
    }
  }
  
  # Combine class-specific results
  if (length(class_specific_results) > 0) {
    class_specific_df <- bind_rows(class_specific_results)
  } else {
    class_specific_df <- NULL
  }
  
  return(list(
    coefficients = results, 
    stats = model_stats,
    class_specific = class_specific_df
  ))
}

# Run analysis for all genes
cat("Running regression analysis for", length(gene_names), "genes...\n")

results_with_class <- map(gene_names,
                          ~run_regression_with_class(.x, lipid_data, gene_lipid_cor),
                          .progress = TRUE)

# Filter out NULL results
results_with_class <- compact(results_with_class)

# Extract results
coefficients_with_class <- map_dfr(results_with_class, ~.x$coefficients)
stats_with_class <- map_dfr(results_with_class, ~.x$stats)
class_specific_all <- map_dfr(results_with_class, ~.x$class_specific)

# Save results
write_csv(coefficients_with_class, "results/multi-omic/structural_model/gene_lipid_regression_coefficients.csv")
write_csv(stats_with_class, "results/multi-omic/structural_model/gene_lipid_regression_stats.csv")
write_csv(class_specific_all, "results/multi-omic/structural_model/gene_lipid_class_specific_all_terms.csv")

cat("\nAnalysis complete!\n")
cat("Successfully analyzed", length(results_with_class), "genes\n")
cat("Skipped", length(gene_names) - length(results_with_class), "genes due to insufficient data\n")
cat("Results saved to:\n")
cat("  - gene_lipid_regression_coefficients.csv (main model)\n")
cat("  - gene_lipid_regression_stats.csv (model fit statistics)\n")
cat("  - gene_lipid_class_specific_all_terms.csv (ALL class-specific effects)\n")

# Create summary table of significant associations by class
cat("\n=== SIGNIFICANT ASSOCIATIONS SUMMARY ===\n")

sig_associations <- class_specific_all %>%
  filter(p.value < 0.05) %>%
  select(gene, lipid_class, term, estimate, std.error, statistic, p.value, n_lipids) %>%
  arrange(term, p.value)

write_csv(sig_associations, "results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv")

cat("\nFound", nrow(sig_associations), "significant class-specific associations\n")
cat("Saved to: significant_gene_lipid_class_associations.csv\n")

# Breakdown by term type
term_summary <- sig_associations %>%
  mutate(term_type = case_when(
    term == "chain_length" ~ "Chain Length",
    term == "db_position" ~ "Double Bond Position",
    term == "lipid_abundance" ~ "Lipid Abundance",
    str_detect(term, ":") ~ "Interaction",
    TRUE ~ "Other"
  )) %>%
  count(term_type, lipid_class) %>%
  pivot_wider(names_from = lipid_class, values_from = n, values_fill = 0)

cat("\nSignificant associations by term type and lipid class:\n")
print(term_summary)

# Genes associated with chain length (by class)
chain_length_genes <- sig_associations %>%
  filter(term == "chain_length") %>%
  select(gene, lipid_class, estimate, p.value) %>%
  arrange(lipid_class, p.value)

cat("\n=== GENES WITH SIGNIFICANT CHAIN LENGTH ASSOCIATIONS ===\n")
cat("Total:", nrow(chain_length_genes), "gene-class pairs\n\n")

chain_summary <- chain_length_genes %>%
  group_by(lipid_class) %>%
  summarise(
    n_genes = n(),
    n_positive = sum(estimate > 0),
    n_negative = sum(estimate < 0),
    .groups = "drop"
  )

print(chain_summary)
write_csv(chain_length_genes, "results/multi-omic/structural_model/genes_associated_with_chain_length.csv")

# Genes associated with DB position (by class)
db_genes <- sig_associations %>%
  filter(term == "db_position") %>%
  select(gene, lipid_class, estimate, p.value) %>%
  arrange(lipid_class, p.value)

cat("\n=== GENES WITH SIGNIFICANT DB POSITION ASSOCIATIONS ===\n")
cat("Total:", nrow(db_genes), "gene-class pairs\n\n")

db_summary <- db_genes %>%
  group_by(lipid_class) %>%
  summarise(
    n_genes = n(),
    n_positive = sum(estimate > 0),
    n_negative = sum(estimate < 0),
    .groups = "drop"
  )

print(db_summary)
write_csv(db_genes, "results/multi-omic/structural_model/genes_associated_with_db_position.csv")

# Genes with significant interactions (by class)
interaction_genes <- sig_associations %>%
  filter(str_detect(term, ":")) %>%
  select(gene, lipid_class, term, estimate, p.value) %>%
  arrange(term, lipid_class, p.value)

cat("\n=== GENES WITH SIGNIFICANT INTERACTION EFFECTS ===\n")
cat("Total:", nrow(interaction_genes), "gene-class pairs\n\n")

interaction_summary <- interaction_genes %>%
  group_by(term, lipid_class) %>%
  summarise(
    n_genes = n(),
    mean_estimate = mean(estimate),
    .groups = "drop"
  ) %>%
  arrange(term, desc(n_genes))

print(interaction_summary)
write_csv(interaction_genes, "results/multi-omic/structural_model/genes_with_significant_interactions.csv")

# Genes appearing in multiple categories
cat("\n=== GENES WITH MULTIPLE STRUCTURAL ASSOCIATIONS ===\n")

multi_structural <- sig_associations %>%
  filter(term %in% c("chain_length", "db_position") | str_detect(term, ":")) %>%
  group_by(gene, lipid_class) %>%
  summarise(
    terms = paste(term, collapse = ", "),
    n_terms = n(),
    .groups = "drop"
  ) %>%
  filter(n_terms > 1) %>%
  arrange(desc(n_terms), lipid_class)

if (nrow(multi_structural) > 0) {
  cat("Found", nrow(multi_structural), "gene-class pairs with multiple structural effects\n")
  print(multi_structural, n = 50)
  write_csv(multi_structural, "results/multi-omic/structural_model/genes_multiple_structural_associations.csv")
} else {
  cat("No genes found with multiple structural associations in same class\n")
}

cat("\n=== OUTPUT FILES CREATED ===\n")
cat("1. significant_gene_lipid_class_associations.csv - All significant associations\n")
cat("2. genes_associated_with_chain_length.csv - Chain length specific\n")
cat("3. genes_associated_with_db_position.csv - DB position specific\n")
cat("4. genes_with_significant_interactions.csv - Interaction effects\n")
cat("5. genes_multiple_structural_associations.csv - Genes with multiple effects\n")