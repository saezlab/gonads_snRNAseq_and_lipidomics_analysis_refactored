################################################################################
# Workflow 7: Create Summary Figure Panel for Structural Analysis
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(patchwork)
library(grid)
library(gridExtra)
library(RColorBrewer)

cat("Creating summary figure panel for structural analysis...\n\n")

# Load data
class_results <- read.csv("results/multi-omic/structural_model/class_specific_regression_results.csv")
lipid_structure <- read.csv("results/multi-omic/structural_model/lipids_with_structure.csv")
cor_matrix <- read.csv("results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv",
                      row.names = 1, check.names = FALSE)

# Key validated genes
key_genes <- c("Cyp2c38", "Cyp2b9", "Cyp2d9", "Cyp2b13")

#===============================================================================
# PANEL A: Heatmap of Chain Length Structural Specificity
#===============================================================================

cat("Creating Panel A: Structural specificity heatmap...\n")

# Create matrix for chain_length significance
chain_data <- class_results %>%
  filter(term == "chain_length", gene %in% key_genes, p_value < 0.05) %>%
  mutate(sig_value = estimate) %>%
  select(gene, lipid_class, sig_value)

# Clean lipid class names
chain_data$lipid_class <- gsub("_with_name", "", chain_data$lipid_class)

# Pivot to wide format
chain_wide <- chain_data %>%
  pivot_wider(names_from = lipid_class, values_from = sig_value, values_fill = NA)

chain_matrix <- as.matrix(chain_wide[, -1])
rownames(chain_matrix) <- chain_wide$gene

# Order columns by median association
col_order <- order(colMeans(chain_matrix, na.rm = TRUE), decreasing = TRUE)
chain_matrix <- chain_matrix[, col_order]

# Create annotation for gene direction
gene_annot <- data.frame(
  Direction = c("Long-chain", "Long-chain", "Short-chain", "Long-chain"),
  row.names = rownames(chain_matrix)
)

# Plot heatmap
pdf("results/multi-omic/structural_model/summary_figure_panelA.pdf", width = 10, height = 4)
pheatmap(chain_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-0.1, 0.1, length.out = 101),
         na_col = "grey90",
         main = "Chain Length Structural Specificity",
         fontsize = 10,
         fontsize_row = 11,
         fontsize_col = 10,
         angle_col = 45,
         annotation_row = gene_annot,
         annotation_colors = list(
           Direction = c("Long-chain" = "#B2182B", "Short-chain" = "#2166AC")
         ))
dev.off()

#===============================================================================
# PANEL B: Cyp2c38-BMP Detailed Example (Raw vs Partial)
#===============================================================================

cat("Creating Panel B: Cyp2c38-BMP example...\n")

# Get BMP data
bmp_lipids <- lipid_structure[lipid_structure$lipid_class == "BMP_with_name", ]
common_lipids <- intersect(bmp_lipids$lipid_name, colnames(cor_matrix))
bmp_data <- bmp_lipids[bmp_lipids$lipid_name %in% common_lipids, ]
bmp_data$gene_correlation <- as.numeric(cor_matrix["Cyp2c38", bmp_data$lipid_name])

# Get p-values for coloring
cyp_bmp_pvals <- class_results %>%
  filter(gene == "Cyp2c38", lipid_class == "BMP_with_name") %>%
  select(term, p_value)

fc_color <- ifelse(cyp_bmp_pvals$p_value[cyp_bmp_pvals$term == "logfc"] < 0.05,
                   "#00BA38", "black")
chain_color <- ifelse(cyp_bmp_pvals$p_value[cyp_bmp_pvals$term == "chain_length"] < 0.05,
                      "#00BA38", "black")
db_color <- ifelse(cyp_bmp_pvals$p_value[cyp_bmp_pvals$term == "db_position"] < 0.05,
                   "#00BA38", "black")

# Compute partial residuals
model_without_fc <- lm(gene_correlation ~ chain_length + db_position, data = bmp_data)
model_fc_on_others <- lm(logfc ~ chain_length + db_position, data = bmp_data)
partial_gene_cor_fc <- residuals(model_without_fc)
partial_fc <- residuals(model_fc_on_others)

model_without_chain <- lm(gene_correlation ~ logfc + db_position, data = bmp_data)
model_chain_on_others <- lm(chain_length ~ logfc + db_position, data = bmp_data)
partial_gene_cor_chain <- residuals(model_without_chain)
partial_chain <- residuals(model_chain_on_others)

model_without_db <- lm(gene_correlation ~ logfc + chain_length, data = bmp_data)
model_db_on_others <- lm(db_position ~ logfc + chain_length, data = bmp_data)
partial_gene_cor_db <- residuals(model_without_db)
partial_db <- residuals(model_db_on_others)

partial_data <- data.frame(
  partial_gene_cor_fc, partial_fc,
  partial_gene_cor_chain, partial_chain,
  partial_gene_cor_db, partial_db
)

# Raw plots
p1_raw <- ggplot(bmp_data, aes(x = gene_correlation, y = logfc)) +
  geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB") +
  geom_smooth(method = "lm", se = TRUE, color = fc_color, linewidth = 1.2) +
  labs(title = "RAW: Fold Change", x = "Gene-Lipid Correlation", y = "LogFC") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_text(color = fc_color, face = "bold"))

p2_raw <- ggplot(bmp_data, aes(x = gene_correlation, y = chain_length)) +
  geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB") +
  geom_smooth(method = "lm", se = TRUE, color = chain_color, linewidth = 1.2) +
  labs(title = "RAW: Chain Length", x = "Gene-Lipid Correlation", y = "Chain Length (C)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_text(color = chain_color, face = "bold"))

p3_raw <- ggplot(bmp_data, aes(x = gene_correlation, y = db_position)) +
  geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB") +
  geom_smooth(method = "lm", se = TRUE, color = db_color, linewidth = 1.2) +
  labs(title = "RAW: DB Position", x = "Gene-Lipid Correlation", y = "Double Bonds") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_text(color = db_color, face = "bold"))

# Partial plots
p1_partial <- ggplot(partial_data, aes(x = partial_gene_cor_fc, y = partial_fc)) +
  geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB") +
  geom_smooth(method = "lm", se = TRUE, color = fc_color, linewidth = 1.2) +
  labs(title = "PARTIAL: Fold Change", x = "Gene Correlation | others", y = "LogFC | others") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_text(color = fc_color, face = "bold"))

p2_partial <- ggplot(partial_data, aes(x = partial_gene_cor_chain, y = partial_chain)) +
  geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB") +
  geom_smooth(method = "lm", se = TRUE, color = chain_color, linewidth = 1.2) +
  labs(title = "PARTIAL: Chain Length", x = "Gene Correlation | others", y = "Chain Length | others") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_text(color = chain_color, face = "bold"))

p3_partial <- ggplot(partial_data, aes(x = partial_gene_cor_db, y = partial_db)) +
  geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB") +
  geom_smooth(method = "lm", se = TRUE, color = db_color, linewidth = 1.2) +
  labs(title = "PARTIAL: DB Position", x = "Gene Correlation | others", y = "DB Position | others") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_text(color = db_color, face = "bold"))

# Combine
combined <- (p1_raw | p2_raw | p3_raw) / (p1_partial | p2_partial | p3_partial) +
  plot_annotation(
    title = "Cyp2c38 - BMP Lipid Structural Specificity",
    subtitle = "Raw bivariate relationships (top) vs. Partial effects controlling for other predictors (bottom)\nGreen = p < 0.05, Black = not significant",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 11))
  )

ggsave("results/multi-omic/structural_model/summary_figure_panelB.pdf",
       combined, width = 15, height = 8)

#===============================================================================
# PANEL C: Summary Statistics
#===============================================================================

cat("Creating Panel C: Summary statistics...\n")

# Count associations by lipid class
summary_by_class <- class_results %>%
  filter(term == "chain_length", p_value < 0.05) %>%
  mutate(lipid_class = gsub("_with_name", "", lipid_class)) %>%
  group_by(lipid_class) %>%
  summarise(
    n_genes = length(unique(gene)),
    mean_effect = mean(abs(estimate))
  ) %>%
  arrange(desc(n_genes)) %>%
  head(15)

p_summary <- ggplot(summary_by_class, aes(x = reorder(lipid_class, n_genes), y = n_genes)) +
  geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Genes with Chain Length Structural Specificity",
    subtitle = "Number of genes per lipid class (p < 0.05)",
    x = "Lipid Class",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

ggsave("results/multi-omic/structural_model/summary_figure_panelC.pdf",
       p_summary, width = 8, height = 6)

cat("\n=== Summary figure panels created ===\n")
cat("  - Panel A: results/multi-omic/structural_model/summary_figure_panelA.pdf\n")
cat("  - Panel B: results/multi-omic/structural_model/summary_figure_panelB.pdf\n")
cat("  - Panel C: results/multi-omic/structural_model/summary_figure_panelC.pdf\n\n")

cat("Recommended figure layout for publication:\n")
cat("----------------------------------------\n")
cat("|  Panel A        |  Panel C          |\n")
cat("|  (Heatmap)      |  (Summary bars)   |\n")
cat("|                 |                   |\n")
cat("|---------------------------------------|\n")
cat("|                                       |\n")
cat("|  Panel B: Cyp2c38-BMP Example        |\n")
cat("|  (6-panel: Raw top, Partial bottom)  |\n")
cat("|                                       |\n")
cat("-----------------------------------------\n\n")
