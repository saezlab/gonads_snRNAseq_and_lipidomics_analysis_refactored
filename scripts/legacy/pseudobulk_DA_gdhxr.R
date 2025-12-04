# library(Seurat)
library(limma)
library(readxl)
library(ggplot2)
library(ggrepel)
library(readr)
library(cowplot)
library(gridExtra)
library(grid)
library(decoupleR)
library(ggridges)
library(dplyr)
library(viridis)
library(reshape2)
library(pheatmap)


source("scripts/support_functions.R")

#I sploit the data in three parts so that the file size doesn't exceed github limit
scaled_RNA_1 <- as.data.frame(read_csv("data/RNA/scaled_RNA_gdxhr_1.csv"))
scaled_RNA_2 <- as.data.frame(read_csv("data/RNA/scaled_RNA_gdxhr_2.csv"))
scaled_RNA_3 <- as.data.frame(read_csv("data/RNA/scaled_RNA_gdxhr_3.csv"))

scaled_RNA <- merge(scaled_RNA_1, scaled_RNA_2, by = "gene")
scaled_RNA <- merge(scaled_RNA, scaled_RNA_3, by = "gene")

row.names(scaled_RNA) <- scaled_RNA$gene
scaled_RNA <- scaled_RNA[,-1]

rm(scaled_RNA_1)
rm(scaled_RNA_2)
rm(scaled_RNA_3)

meta_data_df <- as.data.frame(read_csv("data/RNA/meta_data_df_gdxhr.csv"))

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,3)], no_label = T))

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
  to_enrich <- melt(target[,c(1,4)], id.vars = "sample")
  to_enrich <- to_enrich[,-2]
  names(to_enrich) <- c("target","source")
  to_enrich$mor <- 1
  
  factor_enrichment <- decoupleR::run_ulm(sample_coordinates, to_enrich, minsize = 2)
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

# Create the joy plot (ridgeline plot) using ggridges
ggplot(pca_combined, aes(x = PC1, y = cell_type, fill = group, color = NULL)) +
  geom_density_ridges(scale = 2, alpha = 0.7) +
  labs(x = "PC1 Scores",
       y = "Density (per cell type)") +
  theme_ridges() +
  # scale_fill_manual(values = c("red", "royalblue")) +
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
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 100)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 100)
palette <- c(palette1, palette2)
pheatmap::pheatmap(factor_enrichment_df[,-1], display_numbers = F,
                   show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315, annotation_col = factor_enrichment_df_annotation)

#############


cell_type <- "Hepatocytes"
targets <- meta_data_df[meta_data_df$cell_type == cell_type,]
#
targets$mfsg <- gsub("-[etv]$","",targets$group)
targets <- targets[!grepl("-[et]$",targets$group),]

scaled_RNA_subset <- scaled_RNA[,targets$sample]

pca_groups <- targets[,c(1,5)]
names(pca_groups) <- c("sample","condition")



grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA_subset, pca_groups, no_label = T))

pcares <- prcomp(t(scaled_RNA_subset))
sample_coordinates <- pcares$x[,-38]

to_write <- as.data.frame(sample_coordinates)
to_write$samples <- row.names(to_write)

write_csv(to_write, file = "results/pseudobulk/gdxrh/PCA_sample_coord_RNA_gdxhr_HEPA.csv")

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
  
  if(cell_type == "Hepatocytes")
  {
    condensed_targets <- targets
    # condensed_targets$condition <- gsub("-[etv]$","",condensed_targets$condition)
    # condensed_targets <- condensed_targets[condensed_targets$condition %in% c("m-g-v","m-s-v","f-g-v","f-s-v"),]
    grid.newpage()  # Start a new page
    grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
    grid.draw(nicePCA(scaled_RNA_subset, targets = condensed_targets, no_label = T))
  }
  
  comparisons <- list("testes_gonadosectomy" = c(5,-1),
                      "testes_oestro" = c(8,-5),
                      "testes_testo" = c(7,-5),
                      "ova_gonadosectomy" = c(6,-3),
                      "ova_oestro" = c(4,-6),
                      "ova_testo" = c(2,-6))
  
  limmaRes <- runLimma(scaled_RNA_subset, targets, comparisons = comparisons)
  
  ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                      comp_names = names(comparisons),
                                      number = length(scaled_RNA_subset[,1]),
                                      adjust.method = "fdr")
  
  groups_ttop_list[[cell_type]] <- ttop_list
}

save(groups_ttop_list, file = "results/pseudobulk/gdhxr_groups_ttop_list.RData")

hepatocyte_ttop <- groups_ttop_list$Hepatocytes

tval_df <- merge(hepatocyte_ttop[[1]][,c(1,4)],hepatocyte_ttop[[2]][,c(1,4)], by = "ID")
for(i in 3:length(hepatocyte_ttop))
{
  tval_df <- merge(tval_df,hepatocyte_ttop[[i]][,c(1,4)], by = "ID")
}
names(tval_df)[-1] <- names(hepatocyte_ttop)

write_csv(tval_df, file = "results/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv")

pheatmap(cor(tval_df[,-1]), display_numbers = T)

