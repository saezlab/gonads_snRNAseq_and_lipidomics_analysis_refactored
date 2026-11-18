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
library(tidyr)


source("scripts/support_functions.R")


#I sploit the data in three parts so that the file size doesn't exceed github limit
scaled_RNA <- as.data.frame(read_csv("data/RNA/gdxhr_scaled_RNA_zonation.csv"))


row.names(scaled_RNA) <- scaled_RNA$gene
scaled_RNA <- scaled_RNA[,-1]

meta_data_df <- as.data.frame(read_csv("data/RNA/gdxhr_meta_data_df_zonation.csv"))

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,4)], no_label = T))

pcares <- prcomp(t(scaled_RNA))

View(pcares$x)
# 34, 36 and 41 are looking like strong outliers on PC1 and 2. Need to remove them.

scaled_RNA <- scaled_RNA[,!grepl("_34_", names(scaled_RNA))]
scaled_RNA <- scaled_RNA[,!grepl("_36_", names(scaled_RNA))]
scaled_RNA <- scaled_RNA[,!grepl("_41_", names(scaled_RNA))]

meta_data_df <- meta_data_df[meta_data_df$sample %in% names(scaled_RNA),]

scaled_RNA <- scaled_RNA[,meta_data_df$sample]

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,4)], no_label = T))

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,3)], no_label = T))

# Separate the column into three new columns
meta_data_df <- meta_data_df %>%
  separate(group, into = c("col1", "col2", "col3"), sep = "-", remove = F)
meta_data_df <- as.data.frame(meta_data_df)

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,5)], no_label = T))

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,6)], no_label = T))

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,7)], no_label = T))

meta_data_df$no_fm <- paste(meta_data_df$col2, meta_data_df$col3, sep = "-")

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,8)], no_label = T))

meta_data_df$no_vte <- paste(meta_data_df$col1, meta_data_df$col2, sep = "-")
meta_data_df$no_sg <- paste(meta_data_df$col1, meta_data_df$col3, sep = "-")

pcares <- prcomp(t(scaled_RNA))

Z_matrix <- as.data.frame(pcares$x)[,c(1:15)]

meta_data_df_long <- melt(meta_data_df[,-2], id.vars = "sample") 
meta_data_df_long <- meta_data_df_long[,-2]
meta_data_df_long$mor <- 1
names(meta_data_df_long) <- c("target","source","mor")
# meta_data_df_long <- unique(meta_data_df_long)

meta_data_enrichment <- decoupleR::run_ulm(mat = Z_matrix, network = meta_data_df_long)
meta_data_enrichment <- dcast(meta_data_enrichment, source~condition, value.var = "score")

meta_data_enrichment <- meta_data_enrichment[,c("source",paste("PC",1:15,sep = ""))]
row.names(meta_data_enrichment) <- meta_data_enrichment$source

explained <- c((pcares$sdev)^2 / sum(pcares$sdev^2))[1:15]
names(explained) <- names(meta_data_enrichment)[-1]
explained <- as.data.frame(explained)

t <- as.vector(t(meta_data_enrichment[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(meta_data_enrichment[,-1], display_numbers = T,
                   show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315, annotation_col = explained)

#e vs t
#e-t vs v

meta_data_df$hormonal_zonation <- paste(meta_data_df$zonation_type, meta_data_df$col3, sep = "_") 
#############

pdf_size <- 9

group_variable <- "zonation_type" #gonadal_type or group or allosome_type
targets <- meta_data_df

targets <- targets[,c("sample",group_variable)]
names(targets)[2] <- "condition"
print(unique(targets$condition))

comparisons <- list("Central_vs_portal" = c(1,-2))

limmaRes <- runLimma(scaled_RNA, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(scaled_RNA[,1]),
                                    adjust.method = "fdr")



volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40)
volca_plt + ggtitle("Central_vs_portal")

###################
group_variable <- "hormonal_zonation" #gonadal_type or group or allosome_type
targets <- meta_data_df

targets <- targets[,c("sample",group_variable)]
names(targets)[2] <- "condition"
print(unique(targets$condition))

comparisons <- list("Pericentral_t_vs_Pericentral_e" = c(3,-5),
                    "Periportal_t_vs_Periportal_e" = c(4,-6),
                    "Central_vs_portal" = c(1,3,5,-2,-4,-6),
                    "Pericentral_te_vs_Pericentral_v" = c(3,5,-1,-1),
                    "Periportal_te_vs_Periportal_v" = c(4,6,-2,-2),
                    "Pericentral_t_vs_Pericentral_v" = c(3,-1),
                    "Pericentral_e_vs_Pericentral_v" = c(5,-1),
                    "Periportal_t_vs_Periportal_v" = c(4,-2),
                    "Periportal_e_vs_Periportal_v" = c(6,-2))

limmaRes <- runLimma(scaled_RNA, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(scaled_RNA[,1]),
                                    adjust.method = "fdr")

volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40)

t_table <- ttop_list_to_t_table(ttop_list)
names(t_table) <- c("gene",names(comparisons))

write_csv(t_table, file = "results/pseudobulk/RNA_zonation/t_table_zonation.csv")
