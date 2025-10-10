# Load required libraries
library(tidyverse)
library(patchwork)

# Read data
lipid_structure <- read_csv("results/multi-omic/structural_model/lipids_with_structure.csv")
gene_lipid_cor <- read_csv("results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv")
sig_associations <- read_csv("results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv")

# Prepare lipid data
lipid_data <- lipid_structure %>%
  select(lipid, lipid_class, C, DB, LogFC) %>%
  rename(chain_length = C,
         db_position = DB,
         fold_change = LogFC)

# Function to create scatterplot canvas for a gene
plot_gene_lipid_structure <- function(gene_name, max_classes = 5, p_threshold = 0.05) {
  
  # Get correlation values for this gene
  gene_cors <- gene_lipid_cor[[gene_name]]
  
  # Combine with lipid data
  plot_data <- lipid_data %>%
    mutate(gene_correlation = gene_cors) %>%
    filter(!is.na(gene_correlation))
  
  # Get significant classes for this gene
  sig_classes <- sig_associations %>%
    filter(gene == gene_name) %>%
    distinct(lipid_class) %>%
    pull(lipid_class)
  
  if (length(sig_classes) == 0) {
    cat("No significant associations found for gene:", gene_name, "\n")
    return(NULL)
  }
  
  # Limit to max_classes
  if (length(sig_classes) > max_classes) {
    cat("Gene", gene_name, "has", length(sig_classes), "significant classes.",
        "Showing top", max_classes, "\n")
    sig_classes <- sig_classes[1:max_classes]
  }
  
  # Filter data to significant classes
  plot_data <- plot_data %>%
    filter(lipid_class %in% sig_classes)
  
  # Create list to store plots
  plot_list <- list()
  
  # For each significant class, create 3 plots
  for (lclass in sig_classes) {
    class_data <- plot_data %>% filter(lipid_class == lclass)
    
    # Check significance for this gene-class combination
    sig_data <- sig_associations %>%
      filter(gene == gene_name, lipid_class == lclass)
    
    # Check which terms are significant at the specified threshold
    fc_sig <- any(sig_data$term == "lipid_abundance" & sig_data$p.value < p_threshold)
    chain_sig <- any(sig_data$term == "chain_length" & sig_data$p.value < p_threshold)
    db_sig <- any(sig_data$term == "db_position" & sig_data$p.value < p_threshold)
    fc_chain_sig <- any(sig_data$term == "lipid_abundance:chain_length" & sig_data$p.value < p_threshold)
    fc_db_sig <- any(sig_data$term == "lipid_abundance:db_position" & sig_data$p.value < p_threshold)
    
    # Set colors and y-axis labels based on significance
    fc_color <- if(fc_sig) "#00BA38" else "black"  # green if sig, black otherwise
    chain_color <- if(chain_sig) "#00BA38" else "black"
    db_color <- if(db_sig) "#00BA38" else "black"
    fc_chain_color <- if(fc_chain_sig) "#00BA38" else "black"
    fc_db_color <- if(fc_db_sig) "#00BA38" else "black"
    
    # Plot 1: Correlation vs Fold Change
    p1 <- ggplot(class_data, aes(x = gene_correlation, y = fold_change)) +
      geom_point(alpha = 0.6, size = 2, color = "#2E86AB") +
      geom_smooth(method = "lm", se = TRUE, color = fc_color, linewidth = 1) +
      labs(
        title = paste0(lclass),
        x = "Gene-Lipid Correlation",
        y = "Fold Change (LogFC)"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = fc_color, face = "bold")
      )
    
    # Plot 2: Correlation vs Chain Length
    p2 <- ggplot(class_data, aes(x = gene_correlation, y = chain_length)) +
      geom_point(alpha = 0.6, size = 2, color = "#2E86AB") +
      geom_smooth(method = "lm", se = TRUE, color = chain_color, linewidth = 1) +
      labs(
        x = "Gene-Lipid Correlation",
        y = "Carbon Chain Length"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = chain_color, face = "bold")
      )
    
    # Plot 3: Correlation vs DB Position
    p3 <- ggplot(class_data, aes(x = gene_correlation, y = db_position)) +
      geom_point(alpha = 0.6, size = 2, color = "#2E86AB") +
      geom_smooth(method = "lm", se = TRUE, color = db_color, linewidth = 1) +
      labs(
        x = "Gene-Lipid Correlation",
        y = "Double Bond Position"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = db_color, face = "bold")
      )
    
    # Plot 4: Correlation vs Fold Change * Chain Length (interaction)
    class_data <- class_data %>%
      mutate(fc_chain_interaction = fold_change * chain_length)
    
    p4 <- ggplot(class_data, aes(x = gene_correlation, y = fc_chain_interaction)) +
      geom_point(alpha = 0.6, size = 2, color = "#2E86AB") +
      geom_smooth(method = "lm", se = TRUE, color = fc_chain_color, linewidth = 1) +
      labs(
        x = "Gene-Lipid Correlation",
        y = "FC × Chain Length"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = fc_chain_color, face = "bold")
      )
    
    # Plot 5: Correlation vs Fold Change * DB Position (interaction)
    class_data <- class_data %>%
      mutate(fc_db_interaction = fold_change * db_position)
    
    p5 <- ggplot(class_data, aes(x = gene_correlation, y = fc_db_interaction)) +
      geom_point(alpha = 0.6, size = 2, color = "#2E86AB") +
      geom_smooth(method = "lm", se = TRUE, color = fc_db_color, linewidth = 1) +
      labs(
        x = "Gene-Lipid Correlation",
        y = "FC × DB Position"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = fc_db_color, face = "bold")
      )
    
    # Add plots to list
    plot_list <- c(plot_list, list(p1, p2, p3, p4, p5))
  }
  
  # Combine all plots
  n_classes <- length(sig_classes)
  combined_plot <- wrap_plots(plot_list, ncol = 5, nrow = n_classes) +
    plot_annotation(
      title = paste0("Gene: ", gene_name),
      subtitle = paste0("Showing ", n_classes, " significant lipid classes (p < ", p_threshold, ")"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )
  
  return(combined_plot)
}

# Example usage - replace with your gene of interest
gene_of_interest <- "Cyp2c38"  # Change this to any gene
p_value_threshold <- 0.05    # Adjust significance threshold (e.g., 0.01, 0.05, 0.1)

# Create the plot
p <- plot_gene_lipid_structure(gene_of_interest, max_classes = 5, p_threshold = p_value_threshold)

if (!is.null(p)) {
  print(p)
  
  # Save the plot
  ggsave(
    filename = paste0("results/multi-omic/structural_model/gene_lipid_structure_", gene_of_interest, ".pdf"),
    plot = p,
    width = 18,
    height = 3 * min(5, length(unique(sig_associations %>% 
                                        filter(gene == gene_of_interest) %>% 
                                        pull(lipid_class)))),
    units = "in"
  )
  
  cat("\nPlot saved as:", paste0("results/multi-omic/structural_model/gene_lipid_structure_", gene_of_interest, ".pdf\n"))

  # Save plot data for testing
  if (Sys.getenv("TESTING_MODE") == "TRUE") {
    plot_data_for_gene <- lipid_data %>%
      mutate(gene_correlation = gene_lipid_cor[[gene_of_interest]]) %>%
      filter(!is.na(gene_correlation)) %>%
      filter(lipid_class %in% (sig_associations %>%
                                filter(gene == gene_of_interest) %>%
                                distinct(lipid_class) %>%
                                pull(lipid_class)))

    write_csv(plot_data_for_gene,
              paste0("tests/outputs/lipid_structure_vis_", gene_of_interest, "_data.csv"))
    cat("Test data saved to: tests/outputs/lipid_structure_vis_", gene_of_interest, "_data.csv\n")
  }
}

# Function to generate plots for multiple genes
generate_plots_for_genes <- function(gene_list, max_classes = 5, p_threshold = 0.05, output_dir = "results/multi-omic/structural_model/gene_plots") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  for (gene in gene_list) {
    cat("\nGenerating plot for", gene, "...\n")
    
    p <- plot_gene_lipid_structure(gene, max_classes = max_classes, p_threshold = p_threshold)
    
    if (!is.null(p)) {
      # Determine height based on number of classes
      n_classes <- length(unique(sig_associations %>% 
                                   filter(gene == gene) %>% 
                                   pull(lipid_class)))
      n_classes <- min(n_classes, max_classes)
      
      ggsave(
        filename = file.path(output_dir, paste0(gene, "_lipid_structure.pdf")),
        plot = p,
        width = 18,
        height = 3 * n_classes,
        units = "in"
      )
    }
  }
  
  cat("\nAll plots saved in:", output_dir, "\n")
}

# Example: Generate plots for top genes with structural associations
top_genes <- sig_associations %>%
  filter(term %in% c("chain_length", "db_position")) %>%
  count(gene, sort = TRUE) %>%
  head(10) %>%
  pull(gene)

cat("\nTop genes with structural associations:\n")
print(top_genes)

# Uncomment to generate plots for all top genes
# generate_plots_for_genes(top_genes, max_classes = 5, p_threshold = 0.05)