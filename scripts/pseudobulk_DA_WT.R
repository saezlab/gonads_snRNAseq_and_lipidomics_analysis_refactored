# library(Seurat)
library(limma)
library(readxl)
library(ggplot2)
library(ggrepel)
library(readr)
library(cowplot)
library(decoupleR)
library(ggridges)
library(dplyr)
library(viridis)
library(reshape2)
library(gridExtra)
library(grid)

source("scripts/support_functions.R")


Suppl_Table_1_OverviewData <- as.data.frame(read_excel("data/Suppl_Table_1_OverviewData.xlsx"))
Suppl_Table_1_OverviewData <- Suppl_Table_1_OverviewData[,-4]

#I sploit the data in three parts so that the file size doesn't exceed github limit
scaled_RNA_1 <- as.data.frame(read_csv("data/RNA/scaled_RNA_withWT_1.csv"))
scaled_RNA_2 <- as.data.frame(read_csv("data/RNA/scaled_withWT_RNA_2.csv"))
scaled_RNA_3 <- as.data.frame(read_csv("data/RNA/scaled_withWT_RNA_3.csv"))

scaled_RNA <- merge(scaled_RNA_1, scaled_RNA_2, by = "gene")
scaled_RNA <- merge(scaled_RNA, scaled_RNA_3, by = "gene")

row.names(scaled_RNA) <- scaled_RNA$gene
scaled_RNA <- scaled_RNA[,-1]

rm(scaled_RNA_1)
rm(scaled_RNA_2)
rm(scaled_RNA_3)

meta_data_df <- as.data.frame(read_csv("data/RNA/meta_data_df_withWT.csv"))

meta_data_df <- merge(meta_data_df, Suppl_Table_1_OverviewData[,-3], by.x = "sample_name",by.y = "sequencing_name")

color_mapping <- c("Hepatocytes" = "#ca3659",
                   "Endothelial cells" = "#8f7b2f",
                   "Immune cells" = "#0d934d",
                   "Stellate cells" = "#1892b1",
                   "Cholangiocytes" = "#834f9f")

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(2,3)], no_label = T, color_mapping = color_mapping))

print(unique(meta_data_df$cell_type))
print(unique(meta_data_df$gonadal_type))
print(unique(meta_data_df$group))

#############
target <-  meta_data_df
# Loop through each cell type to subset the data and perform PCA
cell_types <- unique(target$cell_type)

# Prepare a list to store PCA results for each cell type
pca_results <- list()
factor_enrichment_list <- list()

for (cell_type in cell_types) {
  print(cell_type)
  # Subset the scaled_RNA for the current cell type
  samples_in_celltype <- target %>%
    filter(cell_type == !!cell_type) %>%
    pull(sample)  # Get sample names for the current cell type
  
  # Subset the scaled_RNA data
  scaled_RNA_subset <- scaled_RNA[, colnames(scaled_RNA) %in% samples_in_celltype]
  
  # Perform PCA
  pca <- prcomp(t(scaled_RNA_subset), center = TRUE)  # Transpose the data to have samples as rows
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)
  
  
  # Extract the PCA scores for the first principal component (PC1)
  pca_scores <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1]  # Taking the first principal component
  )
  
  # Merge with the target metadata (cell type, gonadal type)
  pca_merged <- pca_scores %>%
    inner_join(target, by = "sample") %>%
    mutate(cell_type = cell_type)  # Add the cell type information
  
  # Store the result in the list
  pca_results[[cell_type]] <- pca_merged
  
  which(explained < 0.025)[1]
  sample_coordinates <- pca$x[,-c((which(explained < 0.025)[1]):length(scaled_RNA_subset[1,]))]
  
  ### factor enrichment
  to_enrich <- melt(target[,c(2,4,6,7)], id.vars = "sample")
  to_enrich <- to_enrich[,-2]
  names(to_enrich) <- c("target","source")
  to_enrich$mor <- 1
  
  factor_enrichment <- decoupleR::run_ulm(sample_coordinates, to_enrich)
  factor_enrichment <- dcast(factor_enrichment, source~condition, value.var = "score")
  row.names(factor_enrichment) <- factor_enrichment[,1]
  factor_enrichment <- factor_enrichment[,-1]
  
  #rescale with variance epxlained
  factor_enrichment <- factor_enrichment * explained[1:(which(explained < 0.025)[1])-1]
  factor_enrichment <- factor_enrichment[,order(as.numeric(gsub("PC","",names(factor_enrichment))))]
  names(factor_enrichment) <- paste(names(factor_enrichment), cell_type, sep = "_")
  factor_enrichment$condition <- row.names(factor_enrichment)
  
  factor_enrichment_list[[cell_type]] <- factor_enrichment
  
}

# Combine PCA results into a single data frame
pca_combined <- bind_rows(pca_results)
pca_combined <- pca_combined[rev(seq_len(nrow(pca_combined))), ]
pca_combined$cell_type <- factor(pca_combined$cell_type, levels = unique(pca_combined$cell_type))

# Create the joy plot (ridgeline plot) using ggridges
ggplot(pca_combined, aes(x = PC1, y = cell_type, fill = allosome_type, color = NULL)) +
  geom_density_ridges(scale = 2, alpha = 0.7) +
  labs(x = "PC1 Scores",
       y = "Density (per cell type)") +
  theme_ridges() +
  scale_fill_manual(values = c("brown", "cyan")) +
  theme(legend.position = "right")

#############


cell_type <- "Hepatocytes"
targets <- meta_data_df[meta_data_df$cell_type == cell_type,]

scaled_RNA_subset <- scaled_RNA[,targets$sample]
pca_groups <- targets[,c(2,4)]
names(pca_groups) <- c("sample","condition")

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA_subset, pca_groups, no_label = T))

pcares <- prcomp(t(scaled_RNA_subset))
sample_coordinates <- pcares$x[,-48]


pdf_size <- 9

groups_ttop_list <- list()
for(cell_type in unique(meta_data_df$cell_type))
{
  print(cell_type)
  group_variable <- "group" #gonadal_type or group
  targets <- meta_data_df[meta_data_df$cell_type == cell_type,]
  
  scaled_RNA_subset <- scaled_RNA[,targets$sample]
  
  targets <- targets[,c("sample",group_variable)]
  names(targets)[2] <- "condition"
  print(unique(targets$condition))
  
  comparisons <- list("XYT-WT_vs_XXO" = c(5,-1),
                      "XYO_vs_XXO" = c(3,-1))
  
  limmaRes <- runLimma(scaled_RNA_subset, targets, comparisons = comparisons)
  
  ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                      comp_names = names(comparisons),
                                      number = length(scaled_RNA_subset[,1]),
                                      adjust.method = "fdr")
  
  names(ttop_list[[1]])[4] <- paste("t",cell_type, sep = "_")
  
  groups_ttop_list[[cell_type]] <- ttop_list
  # volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40) + ggtitle(paste(cell_type, names(comparisons), sep = " "))
  # name_of_file <- paste("results/pseudobulk/volcanos/",paste(paste(cell_type, names(comparisons), sep = "_"),".pdf",sep = ""),sep = "")
  # ggsave(filename = name_of_file, plot = volca_plt, device = "pdf", width = pdf_size, height = pdf_size)
}

ttop_WT <- groups_ttop_list$Hepatocytes$`XYT-WT_vs_XXO`

ttop_Ynosry <- groups_ttop_list$Hepatocytes$XYO_vs_XXO

load("C:/Users/dugourd/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis/results/pseudobulk/ttop_full_list.RData")

ttop_gonad <-  ttop_full_list$testes_vs_ova$Hepatocytes

t_table <- merge(ttop_WT[,c(1,4)],ttop_Ynosry[,c(1,4)], by = "ID")
t_table <- merge(t_table,ttop_gonad[,c(1,4)], by = "ID")

names(t_table) <- c("ID","WT","Ynosry","gonad")

plot(t_table$WT, t_table$Ynosry)
plot(t_table$WT, t_table$gonad)

t_table$label_Ynosry <- ifelse(t_table$WT * t_table$Ynosry > 50,t_table$ID,NA)
t_table$label_gonad <- ifelse(t_table$WT * t_table$gonad > 50 | abs(t_table$WT) > 25,t_table$ID,NA)

ggplot(t_table, aes(x = WT, y = Ynosry,label = label_Ynosry)) + 
  geom_point() +
  geom_label_repel() +
  theme_minimal() +
  xlab("XYT-WT_vs_XXO") + ylab("XYO_vs_XXO")

ggplot(t_table, aes(x = WT, y = gonad, label = label_gonad)) + 
  geom_point() +
  geom_label_repel() +
  theme_minimal() +
  xlab("XYT-WT_vs_XXO") + ylab("testes_vs_ova")

#what % of the genes in the XYT−WT_vs_XXO DEG list (i.e. sex-dependent genes) are also in the testes_vs_ova DEG list (i.e. gonad-dependent genes)?
sum(ttop_WT[ttop_WT$adj.P.Val < 0.05,1] %in% ttop_gonad[ttop_gonad$adj.P.Val < 0.05,1]) / length(ttop_WT[ttop_WT$adj.P.Val < 0.05,1])

# library(plotly)
# 
# plot_ly(
#   data = t_table,
#   x = ~WT,
#   y = ~Ynosry,
#   z = ~gonad,
#   text = ~ID,            # hover text
#   type = "scatter3d",
#   mode = "markers",
#   marker = list(size = 4)
# ) %>%
#   layout(
#     scene = list(
#       xaxis = list(title = "WT"),
#       yaxis = list(title = "Ynosry"),
#       zaxis = list(title = "Gonad")
#     )
#   )
# 
# 
