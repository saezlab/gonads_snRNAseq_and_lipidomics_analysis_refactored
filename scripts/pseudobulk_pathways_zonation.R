library(decoupleR)
library(readr)
library(GSEABase)
library(reshape2)
library(pheatmap)

source("scripts/support_functions.R")

t_table_zonation <- as.data.frame(
  read_csv("results/pseudobulk/rna_zonation/t_table_zonation.csv"))

pathways <- import_gmt("support/m5.all.v2023.2.Mm.symbols.gmt")
names(pathways) <- c("target","source")
pathways$mor <- 1

tval_full_matrix <- t_table_zonation
row.names(tval_full_matrix) <- tval_full_matrix$gene
tval_full_matrix <- tval_full_matrix[,-1]

pathway_scores <- decoupleR::run_ulm(tval_full_matrix, network = pathways)
pathway_scores_df <- dcast(pathway_scores, source~condition, value.var = "score")
row.names(pathway_scores_df)<- tolower(gsub("_"," ",gsub("GO[A-Z]{2}_","",pathway_scores_df$source)))

to_heatmap <- pathway_scores_df[,c(2,3,5,4,8)]

# to_heatmap <- to_heatmap[apply(to_heatmap,1,function(x){max(abs(x)) > 4}),]
to_heatmap <- to_heatmap[to_heatmap[,5] > 7 | to_heatmap[,5] < -7,]
to_heatmap[to_heatmap > 30] <- 30
to_heatmap[to_heatmap < -30] <- -30



t <- as.vector(t(to_heatmap))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(to_heatmap, show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, gaps_col = c(2,4))
pheatmap::pheatmap(to_heatmap, show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, gaps_col = c(2,4), 
                   filename = "results/pseudobulk/rna_zonation/pathway_heatmap.pdf", width = 7, height = 7)
dev.off()

optimality_tasks <- as.data.frame(
  read_csv("support/Cleaned_Full_Task_Names_and_Single_Standardized_Gene_Table.csv"))[,-1]

names(optimality_tasks) <- c("source","target")
optimality_tasks$mor <- 1

tasks_scores <- decoupleR::run_ulm(tval_full_matrix, network = optimality_tasks)
tasks_scores_df <- dcast(tasks_scores, source~condition, value.var = "score")

to_heatmap <- tasks_scores_df[,c(1,3,8,5,4,7)]
row.names(to_heatmap) <- to_heatmap$source
to_heatmap <- to_heatmap[,-1]

t <- as.vector(t(to_heatmap))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(to_heatmap, show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, gaps_col = c(2,4))
pheatmap::pheatmap(to_heatmap, show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, gaps_col = c(2,4), filename = "results/pseudobulk/rna_zonation/pathway_heatmap.pdf", width = 14, height = 6)
dev.off()

pathway_scores <- decoupleR::run_ulm(tval_full_matrix, network = pathways)
pathway_scores_df <- dcast(pathway_scores, source~condition, value.var = "score")

write_csv(pathway_scores_df, file = "results/pseudobulk/rna_zonation/pathway_scores_df.csv")

to_heatmap <- pathway_scores_df[,c(1,3,8,5,4,7)]
row.names(to_heatmap) <- to_heatmap$source
to_heatmap <- to_heatmap[,-1]

# to_heatmap <- to_heatmap[apply(to_heatmap,1,function(x){max(abs(x)) > 4}),]
to_heatmap <- to_heatmap[to_heatmap[,5] > 4 | to_heatmap[,5] < -9.5,]

t <- as.vector(t(to_heatmap))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(to_heatmap, show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, gaps_col = c(2,4))
pheatmap::pheatmap(to_heatmap, show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, gaps_col = c(2,4), filename = "results/pseudobulk/rna_zonation/pathway_heatmap.pdf", width = 14, height = 6)
dev.off()

#rest is just sanity check

plot(t_table_zonation$portaltestes_vs_centraltestes, t_table_zonation$portalova_vs_centralova)

t_table_diff <- data.frame("testes_zonation_effect_p_c" = t_table_zonation$portaltestes_vs_centraltestes - t_table_zonation$portalova_vs_centralova)
t_table_diff$gene <- t_table_zonation$gene
row.names(t_table_diff) <- t_table_diff$gene

mean_diff <- rowMeans(abs(t_table_zonation[,c("portaltestes_vs_centraltestes","portalova_vs_centralova")]))
plot(density(log2(mean_diff)))
abline(v = -0.25)

t_table_diff <- t_table_diff[log2(mean_diff) > -0.25,]
mean_diff <- mean_diff[log2(mean_diff) > -0.25]

plot(t_table_diff$testes_zonation_effect_p_c, mean_diff)

names(mean_diff) <- t_table_diff$gene

t_table_diff$testes_zonation_effect_p_c <- t_table_diff$testes_zonation_effect_p_c / mean_diff



pathway_scores_diff <- decoupleR::run_ulm(t_table_diff[,-2,drop = F], network = pathways)
pathway_scores_diff_df <- dcast(pathway_scores_diff, source~condition, value.var = "score")

write_csv(pathway_scores_diff_df, file = "results/pseudobulk/rna_zonation/pathway_scores_diff.csv")

write_csv(t_table_diff, file = "results/pseudobulk/rna_zonation/t_table_diff.csv")



temp <- merge(pathway_scores_df, pathway_scores_diff_df, by = "source")
plot(temp$testes_effect_zonation, temp$testes_zonation_effect_p_c)