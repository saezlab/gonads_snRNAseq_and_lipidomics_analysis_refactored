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

View(pathway_scores_df[,c(1,2,3,5,7,9)])
View(t_table_zonation[,c(1,4,7,8,9,10)])

plot(t_table_zonation$Periportal_t_vs_Periportal_v, t_table_zonation$Pericentral_t_vs_Pericentral_v)

t_table_diff <- data.frame("testo_zonation_effect_p_c" = t_table_zonation$Periportal_t_vs_Periportal_v - t_table_zonation$Pericentral_t_vs_Pericentral_v,
                           "oestro_zonation_effect_p_c" = t_table_zonation$Periportal_e_vs_Periportal_v - t_table_zonation$Pericentral_e_vs_Pericentral_v)
t_table_diff$gene <- t_table_zonation$gene
row.names(t_table_diff) <- t_table_diff$gene

pathway_scores_diff <- decoupleR::run_ulm(t_table_diff[,-3], network = pathways)
pathway_scores_diff_df <- dcast(pathway_scores_diff, source~condition, value.var = "score")

write_csv(pathway_scores_diff_df, file = "results/pseudobulk/rna_zonation/pathway_scores_diff.csv")

write_csv(t_table_diff, file = "results/pseudobulk/rna_zonation/t_table_diff.csv")
