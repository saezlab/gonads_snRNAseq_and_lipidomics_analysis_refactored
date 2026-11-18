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
library(patchwork)

source("scripts/support_functions.R")


Suppl_Table_1_OverviewData <- as.data.frame(read_excel("data/Suppl_Table_1_OverviewData.xlsx"))
Suppl_Table_1_OverviewData <- Suppl_Table_1_OverviewData[,-4]

#I sploit the data in three parts so that the file size doesn't exceed github limit
scaled_RNA_1 <- as.data.frame(read_csv("data/RNA/scaled_RNA_1.csv"))
scaled_RNA_2 <- as.data.frame(read_csv("data/RNA/scaled_RNA_2.csv"))
scaled_RNA_3 <- as.data.frame(read_csv("data/RNA/scaled_RNA_3.csv"))

scaled_RNA <- merge(scaled_RNA_1, scaled_RNA_2, by = "gene")
scaled_RNA <- merge(scaled_RNA, scaled_RNA_3, by = "gene")

row.names(scaled_RNA) <- scaled_RNA$gene
scaled_RNA <- scaled_RNA[,-1]

rm(scaled_RNA_1)
rm(scaled_RNA_2)
rm(scaled_RNA_3)

meta_data_df <- as.data.frame(read_csv("data/RNA/meta_data_df.csv"))
meta_data_df$sample_name <- gsub("-","_",meta_data_df$sample_name)

meta_data_df <- merge(meta_data_df, Suppl_Table_1_OverviewData, by.x = "sample_name",by.y = "sequencing_name")

color_mapping <- c("Hepatocytes" = "#ca3659",
                   "Endothelial cells" = "#8f7b2f",
                   "Immune cells" = "#0d934d",
                   "Stellate cells" = "#1892b1",
                   "Cholangiocytes" = "#834f9f")

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(2,4)], no_label = T, color_mapping = color_mapping))
grid.draw(nicePCA(scaled_RNA_pca, meta_data_df_pca[,c(2,4)], no_label = T, color_mapping = color_mapping))

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
ggplot(pca_combined, aes(x = PC1, y = cell_type, fill = gonadal_type, color = NULL)) +
  geom_density_ridges(scale = 2, alpha = 0.7) +
  labs(x = "PC1 Scores",
       y = "Density (per cell type)") +
  theme_ridges() +
  scale_fill_manual(values = c("red", "royalblue")) +
  theme(legend.position = "right")

# Create the joy plot (ridgeline plot) using ggridges
ggplot(pca_combined, aes(x = PC1, y = cell_type, fill = allosome_type, color = NULL)) +
  geom_density_ridges(scale = 2, alpha = 0.7) +
  labs(x = "PC1 Scores",
       y = "Density (per cell type)") +
  theme_ridges() +
  scale_fill_manual(values = c("red", "royalblue")) +
  # scale_fill_manual(values = c("grey30", "grey70")) +
  # scale_fill_manual(values = c("brown", "cyan")) +
  theme(legend.position = "right")


sub_pca_combined <- pca_combined[pca_combined$group %in% c("XXO","XYO"),]
# Create the joy plot (ridgeline plot) using ggridges
ggplot(sub_pca_combined, aes(x = PC1, y = cell_type, fill = group, color = NULL)) +
  geom_density_ridges(scale = 2, alpha = 0.7) +
  labs(x = "PC1 Scores",
       y = "Density (per cell type)") +
  theme_ridges() +
  scale_fill_manual(values = c("brown", "cyan")) +
  theme(legend.position = "right")

factor_enrichment_df <- Reduce(function(x, y) merge(x, y, all=TRUE), factor_enrichment_list)
row.names(factor_enrichment_df) <- factor_enrichment_df$condition

factor_enrichment_df_annotation <- as.data.frame(names(factor_enrichment_df))[-1,,drop = F]
names(factor_enrichment_df_annotation) <- "cell_type"
# row.names(factor_enrichment_df_annotation) <- gsub("_.*","",factor_enrichment_df_annotation$cell_type)
row.names(factor_enrichment_df_annotation) <- factor_enrichment_df_annotation$cell_type
factor_enrichment_df_annotation$cell_type <- gsub(".*_","",factor_enrichment_df_annotation$cell_type)

# names(factor_enrichment_df) <- gsub("_.*","",names(factor_enrichment_df))

t <- as.vector(t(factor_enrichment_df[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(factor_enrichment_df[,-1], display_numbers = F,
                   show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315, annotation_col = factor_enrichment_df_annotation)

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

to_write <- as.data.frame(sample_coordinates)
to_write$samples <- row.names(to_write)

write_csv(to_write, file = "results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv")
# 
# 
# to_enrich <- melt(targets[,c(2,4,6,7)], id.vars = "sample")
# to_enrich <- to_enrich[,-2]
# names(to_enrich) <- c("target","source")
# to_enrich$mor <- 1
# 
# factor_enrichment <- decoupleR::run_ulm(sample_coordinates, to_enrich)
# factor_enrichment <- dcast(factor_enrichment, source~condition, value.var = "score")
# row.names(factor_enrichment) <- factor_enrichment[,1] 
# factor_enrichment <- factor_enrichment[,-1]
# 
# #rescale with variance epxlained
# factor_enrichment <- factor_enrichment * c((pcares$sdev)^2 / sum(pcares$sdev^2))[-48]
# factor_enrichment <- factor_enrichment[,order(as.numeric(gsub("PC","",names(factor_enrichment))))]
# 
# t <- as.vector(t(factor_enrichment))
# palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
# palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
# palette <- c(palette1, palette2)
# pheatmap::pheatmap(factor_enrichment, display_numbers = F,
#                    show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315)
# 
# factor_enrichment_long <- factor_enrichment
# factor_enrichment_long$condition <- row.names(factor_enrichment_long)
# 
# factor_enrichment_long <- melt(factor_enrichment_long)
# 
# ggplot(factor_enrichment_long, aes(x = variable, y = value, group = condition, color = condition)) +
#   geom_line(size = 2, alpha = 0.5) + 
#   geom_point(size = 3) +
#   theme_minimal()
# 
# write_csv(factor_enrichment_long, file = "results/pseudobulk/factor_enrichment_long.csv")
# 

pdf_size <- 9

gonade_type_ttop_list <- list()
for(cell_type in unique(meta_data_df$cell_type))
{
  print(cell_type)
  group_variable <- "gonadal_type" #gonadal_type or group or allosome_type
  targets <- meta_data_df[meta_data_df$cell_type == cell_type,]
  
  scaled_RNA_subset <- scaled_RNA[,targets$sample]
  
  targets <- targets[,c("sample",group_variable)]
  names(targets)[2] <- "condition"
  print(unique(targets$condition))
  
  #just to make some presentation plots
  # gene_expr <- as.data.frame(t(scaled_RNA_subset["Dbp",]))
  # gene_expr <- merge(gene_expr, targets, by.x = "row.names", by.y = "sample")
  # 
  # ggplot(gene_expr, aes(x = condition, y = Dbp, group = condition, fill = condition, color = condition)) + 
  #   geom_jitter() + theme_minimal()
  # 
  # 
  # gene_expr <- as.data.frame(t(scaled_RNA_subset["Esr1",]))
  # gene_expr <- merge(gene_expr, targets, by.x = "row.names", by.y = "sample")
  # 
  # ggplot(gene_expr, aes(x = condition, y = Esr1, group = condition, fill = condition, color = condition)) + 
  #   geom_jitter() + theme_minimal()
  
  comparisons <- list("testes_vs_Ova" = c(-1,2))
  
  limmaRes <- runLimma(scaled_RNA_subset, targets, comparisons = comparisons)
  
  ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                      comp_names = names(comparisons),
                                      number = length(scaled_RNA_subset[,1]),
                                      adjust.method = "fdr")
  
  names(ttop_list[[1]])[4] <- paste("t",cell_type, sep = "_")
  
  gonade_type_ttop_list[[cell_type]] <- ttop_list[[1]]
  volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40) + ggtitle(paste(cell_type, names(comparisons), sep = " "))
  name_of_file <- paste("results/pseudobulk/volcanos/",paste(paste(cell_type, names(comparisons), sep = "_"),".pdf",sep = ""),sep = "")
  ggsave(filename = name_of_file, plot = volca_plt, device = "pdf", width = pdf_size, height = pdf_size)
}

allosome_type_ttop_list <- list()
for(cell_type in unique(meta_data_df$cell_type))
{
  print(cell_type)
  group_variable <- "allosome_type" #gonadal_type or group
  targets <- meta_data_df[meta_data_df$cell_type == cell_type,]
  
  scaled_RNA_subset <- scaled_RNA[,targets$sample]
  
  targets <- targets[,c("sample",group_variable)]
  names(targets)[2] <- "condition"
  print(unique(targets$condition))
  
  comparisons <- list("XY_vs_XX" = c(-1,2))
  
  limmaRes <- runLimma(scaled_RNA_subset, targets, comparisons = comparisons)
  
  ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                      comp_names = names(comparisons),
                                      number = length(scaled_RNA_subset[,1]),
                                      adjust.method = "fdr")
  
  names(ttop_list[[1]])[4] <- paste("t",cell_type, sep = "_")
  
  allosome_type_ttop_list[[cell_type]] <- ttop_list[[1]]
  volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40) + ggtitle(paste(cell_type, names(comparisons), sep = " "))
  name_of_file <- paste("results/pseudobulk/volcanos/",paste(paste(cell_type, names(comparisons), sep = "_"),".pdf",sep = ""),sep = "")
  ggsave(filename = name_of_file, plot = volca_plt, device = "pdf", width = pdf_size, height = pdf_size)
}

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
  
  comparisons <- list("XXT_vs_XXO" = c(1,-4))
  
  limmaRes <- runLimma(scaled_RNA_subset, targets, comparisons = comparisons)
  
  ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                      comp_names = names(comparisons),
                                      number = length(scaled_RNA_subset[,1]),
                                      adjust.method = "fdr")
  
  names(ttop_list[[1]])[4] <- paste("t",cell_type, sep = "_")
  
  groups_ttop_list[[cell_type]] <- ttop_list[[1]]
  volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40) + ggtitle(paste(cell_type, names(comparisons), sep = " "))
  name_of_file <- paste("results/pseudobulk/volcanos/",paste(paste(cell_type, names(comparisons), sep = "_"),".pdf",sep = ""),sep = "")
  ggsave(filename = name_of_file, plot = volca_plt, device = "pdf", width = pdf_size, height = pdf_size)
}

ttop_full_list <- list("testes_vs_ova" = gonade_type_ttop_list, "XY_vs_XX" = allosome_type_ttop_list, "XXT_vs_XXO" = groups_ttop_list)
save(ttop_full_list, file = "results/pseudobulk/ttop_full_list.RData")

X9WZDEU_mMulti_100g_gene_list <- as.data.frame(read_csv("support/FCGL_Xenium_Panel.csv"))
X9WZDEU_mMulti_100g_gene_list <- X9WZDEU_mMulti_100g_gene_list[X9WZDEU_mMulti_100g_gene_list$Reason != "cell type marker",]

ttop <- ttop_full_list$testes_vs_ova$Hepatocytes

volcano_nice(ttop, IDIndex = 1, FCIndex = 2, pValIndex = 5, label = T, nlabels = 30)  +ylab("-log10(p-value)") +
  xlim(-4,4) + ylim(0,60)


ttop <- ttop[ttop$ID %in% X9WZDEU_mMulti_100g_gene_list$Gene,]

volcano_nice(ttop, IDIndex = 1, FCIndex = 2, pValIndex = 5, label = T, nlabels = 50) + ylab("-log10(p-value)")

ttop <- ttop_full_list$XY_vs_XX$Hepatocytes

# volcano_nice(ttop, IDIndex = 1, FCIndex = 2, pValIndex = 5, label = T, nlabels = 30, col_high_neut_low = c("brown","grey30","cyan")) +ylab("-log10(p-value)")
volcano_nice(ttop, IDIndex = 1, FCIndex = 2, pValIndex = 5, label = T, nlabels = 30, col_high_neut_low = c("grey30","black","grey70")) +ylab("-log10(p-value)") + 
  xlim(-4,4) + ylim(0,60)
volcano_nice(ttop, IDIndex = 1, FCIndex = 2, pValIndex = 5, label = T, nlabels = 30) +ylab("-log10(p-value)") +
  xlim(-4,4) + ylim(0,60)

ttop <- ttop[ttop$ID %in% X9WZDEU_mMulti_100g_gene_list$Gene,]

volcano_nice(ttop, IDIndex = 1, FCIndex = 2, pValIndex = 5, label = T, nlabels = 50) + ylab("-log10(p-value)")

ttop <- ttop_full_list$XXT_vs_XXO$Hepatocytes

volcano_nice(ttop, IDIndex = 1, FCIndex = 2, pValIndex = 5, label = T, nlabels = 30)

ttop <- ttop[ttop$ID %in% X9WZDEU_mMulti_100g_gene_list$Gene,]

volcano_nice(ttop, IDIndex = 1, FCIndex = 2, pValIndex = 5, label = T, nlabels = 50) + ylab("-log10(p-value)")

  
sum(ttop$adj.P.Val <= 0.05)

row.names(ttop) <- ttop$ID

Sry_Target_GeneCard <- as.data.frame(read_csv("support/Sry_Target_GeneCard.csv"))

Sry_Target_GeneCard$target <- capitalize_first(Sry_Target_GeneCard$target)
Sry_Target_GeneCard <- unique(Sry_Target_GeneCard)

Sry_act <- decoupleR::run_ulm(mat = ttop[,4,drop = F], network = Sry_Target_GeneCard)

length(which(ttop$adj.P.Val <= 0.05 & ttop$t_Hepatocytes > 0 & ttop$ID %in% Sry_Target_GeneCard$target))

##TF analysis
load("support/TF_regulons.RData")
TF_regulons <- TF_regulons[complete.cases(TF_regulons),]

Gonadal_Hormone_Genes <- as.data.frame(read_csv("support/Gonadal_Hormone_Genes.csv"))
Gonadal_Hormone_Genes$source <- gsub("_spec","",Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("ESR1","Esr1",Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("AR","Ar",Gonadal_Hormone_Genes$source)

TF_regulons <- TF_regulons[!TF_regulons$source %in% c("Ar", "Esr1"),]

TF_regulons <- rbind(TF_regulons, Gonadal_Hormone_Genes)

TF_act <- decoupleR::run_ulm(mat = ttop[,4,drop = F],network = TF_regulons)

to_plot <- merge(ttop[,c(1,4)], TF_regulons[TF_regulons$source %in% c("Esr1","Ar"),],by.x = "ID",by.y = "target", all.x = T)
to_plot$mor <- ifelse(is.na(to_plot$source),"background",to_plot$source)

ggplot(to_plot, aes(x = mor, y = t_Hepatocytes)) + geom_boxplot() + theme_minimal()
