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
scaled_RNA <- as.data.frame(read_csv("data/RNA/scaled_RNA_zonation.csv"))


row.names(scaled_RNA) <- scaled_RNA$gene
scaled_RNA <- scaled_RNA[,-1]

meta_data_df <- as.data.frame(read_csv("data/RNA/meta_data_df_zonation.csv"))

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,3)], no_label = T))

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,4)], no_label = T))

meta_data_df$condition <- paste(meta_data_df$zonation_type, meta_data_df$gonadal_type, sep = "_")

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(scaled_RNA, meta_data_df[,c(1,5)], no_label = T))

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

#############

pdf_size <- 9

group_variable <- "condition" #gonadal_type or group or allosome_type
targets <- meta_data_df

targets <- targets[,c("sample",group_variable)]
names(targets)[2] <- "condition"
print(unique(targets$condition))

comparisons <- list("portal_vs_Central" = c(-1,-3,2,4),
                    "testes_vs_ova" = c(-1,-2,3,4),
                    "portalova_vs_centralova" = c(-1,2),
                    "portaltestes_vs_centraltestes" = c(-3,4),
                    "centraltestes_vs_centralova" = c(-1,3),
                    "portaltestes_vs_portalova" = c(-2,4),
                    "testes_effect_zonation" = c(-3,4,1,-2))

limmaRes <- runLimma(scaled_RNA, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(scaled_RNA[,1]),
                                    adjust.method = "fdr")



k <- 3
volca_plt <- volcano_nice(ttop_list[[k]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40)
plot(volca_plt) + ggtitle(names(comparisons[k]))

t_table <- ttop_list_to_t_table(ttop_list)
names(t_table) <- c("gene",names(comparisons))

write_csv(t_table, file = "results/pseudobulk/rna_zonation/t_table_zonation.csv")

##TF analysis
load("support/TF_regulons.RData")
TF_regulons <- TF_regulons[complete.cases(TF_regulons),]

Gonadal_Hormone_Genes <- as.data.frame(read_csv("support/Gonadal_Hormone_Genes.csv"))
Gonadal_Hormone_Genes$source <- gsub("_spec","",Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("ESR1","Esr1",Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("AR","Ar",Gonadal_Hormone_Genes$source)

TF_regulons <- TF_regulons[!TF_regulons$source %in% c("Ar", "Esr1"),]

TF_regulons <- rbind(TF_regulons, Gonadal_Hormone_Genes)

TF_scores <- decoupleR::run_ulm(mat = scaled_RNA, network = TF_regulons[,], minsize = 20)
TF_scores_df <- dcast(TF_scores, source~condition, value.var = "score")
row.names(TF_scores_df) <- TF_scores_df$source
TF_scores_df <- TF_scores_df[,-1]

decoupler_input_t_table <- t_table[,-1]
row.names(decoupler_input_t_table) <- t_table$gene

TF_scores_contrasts <- decoupleR::run_ulm(mat = decoupler_input_t_table, network = TF_regulons[,], minsize = 20)
TF_scores_contrasts_df <- dcast(TF_scores_contrasts, source~condition, value.var = "score")
row.names(TF_scores_contrasts_df) <- TF_scores_contrasts_df$source
# TF_scores_contrasts_df <- TF_scores_contrasts_df[,-1]

load("data/metabolomics/lipidomic.RData")


testes_zonation_baias_genes <- t_table[abs(t_table$testes_effect_zonation) > 5 & abs(t_table$portalova_vs_centralova) < 1,]
ovaries_zonation_baias_genes <- t_table[abs(t_table$testes_effect_zonation) > 5 & abs(t_table$portaltestes_vs_centraltestes) < 1,]
zonation_baias_genes <- t_table[abs(t_table$testes_effect_zonation) < 1 & abs(t_table$portal_vs_Central) > 5,]

testes_pericentral_genes <- testes_zonation_baias_genes[testes_zonation_baias_genes$portaltestes_vs_centraltestes < 0,]
testes_periportal_genes <- testes_zonation_baias_genes[testes_zonation_baias_genes$portaltestes_vs_centraltestes > 0,]

ovaries_pericentral_genes <- ovaries_zonation_baias_genes[ovaries_zonation_baias_genes$portalova_vs_centralova < 0,]
ovaries_periportal_genes <- ovaries_zonation_baias_genes[ovaries_zonation_baias_genes$portalova_vs_centralova > 0,]

top_genes <- as.data.frame(do.call(rbind,list(testes_pericentral_genes,testes_periportal_genes,ovaries_pericentral_genes,ovaries_periportal_genes)))
row.names(top_genes) <- top_genes$gene
top_genes <- top_genes[,c(4,5)]

plot(top_genes$portalova_vs_centralova, top_genes$portaltestes_vs_centraltestes)

t <- as.vector(t(top_genes))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10, my_colors = c("lightgreen","white","orange"))
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10, my_colors = c("lightgreen","white","orange"))
palette <- c(palette1, palette2)
pheatmap::pheatmap(t(top_genes), display_numbers = T,
                   show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 0)


########################################################################################################################
########################################################################################################################
########################################################################################################################

#1#######################################################################################################################

#for paper = 0.65
threshold_cor <- 0.65

#testes_portal
scaled_RNA_testes_portal <- scaled_RNA[testes_periportal_genes$gene,grepl("Periportal_Testes",names(scaled_RNA))]
TF_scores_df_testes_portal <- TF_scores_df[,grepl("Periportal_Testes",names(TF_scores_df))]
lipidomic_testes <- lipidomic[,gsub("_Peri.*_Testes","",names(scaled_RNA_testes_portal))]
names(scaled_RNA_testes_portal) <- gsub("_Peri.*_Testes","",names(scaled_RNA_testes_portal))
names(TF_scores_df_testes_portal) <- gsub("_Peri.*_Testes","",names(TF_scores_df_testes_portal))


cormat_testes_portal <- as.data.frame(rbind(make_cor_network(lipidomic_testes, scaled_RNA_testes_portal,theshold_cor = threshold_cor),
                                            make_cor_network(scaled_RNA_testes_portal,TF_scores_df_testes_portal, theshold_cor = threshold_cor, TF = T)))
names(cormat_testes_portal) <- c("node1","node2","weight")
cormat_testes_portal_ATT <- expand_attributes(cormat_testes_portal,t_table,TF_scores_df)
cormat_testes_portal_ATT$type <- ifelse(cormat_testes_portal_ATT$gene %in% row.names(TF_scores_df),"TF",
                                         ifelse(cormat_testes_portal_ATT$gene %in% t_table$gene,"gene","lipid"))


write_csv(cormat_testes_portal, file = "results/multi-omic/cormat_testes_portal_SIF.csv")
write_csv(cormat_testes_portal_ATT, file = "results/multi-omic/cormat_testes_portal_ATT.csv")

cormat_testes_portal$node2 <- gsub("[(].*","",cormat_testes_portal$node2)
cormat_testes_portal$sign <- sign(cormat_testes_portal$weight)
cormat_testes_portal <- aggregate(weight ~ node1+node2+sign, data = cormat_testes_portal, FUN = mean)
cormat_testes_portal$abs_weight <- abs(cormat_testes_portal$weight)

cormat_testes_portal_ATT$gene <- gsub("[(].*","",cormat_testes_portal_ATT$gene)
cormat_testes_portal_ATT <- aggregate(cbind(portal_vs_Central,testes_vs_ova,portalova_vs_centralova,portaltestes_vs_centraltestes,centraltestes_vs_centralova,portaltestes_vs_portalova,testes_effect_zonation) ~ gene+type, data = cormat_testes_portal_ATT, FUN = mean)

write_csv(cormat_testes_portal, file = "results/multi-omic/compressed_cormat_testes_portal_SIF.csv")
write_csv(cormat_testes_portal_ATT, file = "results/multi-omic/compressed_cormat_testes_portal_ATT.csv")

#2#######################################################################################################################

#testes_central
scaled_RNA_testes_central <- scaled_RNA[testes_pericentral_genes$gene,grepl("Pericentral_Testes",names(scaled_RNA))]
TF_scores_df_testes_central <- TF_scores_df[,grepl("Pericentral_Testes",names(TF_scores_df))]
lipidomic_testes <- lipidomic[,gsub("_Peri.*_Testes","",names(scaled_RNA_testes_central))]
names(scaled_RNA_testes_central) <- gsub("_Peri.*_Testes","",names(scaled_RNA_testes_central))
names(TF_scores_df_testes_central) <- gsub("_Peri.*_Testes","",names(TF_scores_df_testes_central))


cormat_testes_central <- as.data.frame(rbind(make_cor_network(lipidomic_testes,scaled_RNA_testes_central,theshold_cor = threshold_cor),
                                            make_cor_network(scaled_RNA_testes_central,TF_scores_df_testes_central, theshold_cor = threshold_cor, TF = T)))
names(cormat_testes_central) <- c("node1","node2","weight")
cormat_testes_central_ATT <- expand_attributes(cormat_testes_central,t_table,TF_scores_df)
cormat_testes_central_ATT$type <- ifelse(cormat_testes_central_ATT$gene %in% row.names(TF_scores_df),"TF",
                                         ifelse(cormat_testes_central_ATT$gene %in% t_table$gene,"gene","lipid"))

write_csv(cormat_testes_central, file = "results/multi-omic/cormat_testes_central_SIF.csv")
write_csv(cormat_testes_central_ATT, file = "results/multi-omic/cormat_testes_central_ATT.csv")

cormat_testes_central$node2 <- gsub("[(].*","",cormat_testes_central$node2)
cormat_testes_central$sign <- sign(cormat_testes_central$weight)
cormat_testes_central <- aggregate(weight ~ node1+node2+sign, data = cormat_testes_central, FUN = mean)
cormat_testes_central$abs_weight <- abs(cormat_testes_central$weight)

cormat_testes_central_ATT$gene <- gsub("[(].*","",cormat_testes_central_ATT$gene)
cormat_testes_central_ATT <- aggregate(cbind(portal_vs_Central,testes_vs_ova,portalova_vs_centralova,portaltestes_vs_centraltestes,centraltestes_vs_centralova,portaltestes_vs_portalova,testes_effect_zonation) ~ gene+type, data = cormat_testes_central_ATT, FUN = mean)

write_csv(cormat_testes_central, file = "results/multi-omic/compressed_cormat_testes_central_SIF.csv")
write_csv(cormat_testes_central_ATT, file = "results/multi-omic/compressed_cormat_testes_central_ATT.csv")

#3#######################################################################################################################

#ovaries_portal
scaled_RNA_ovaries_portal <- scaled_RNA[ovaries_periportal_genes$gene,grepl("Periportal_Ovaries",names(scaled_RNA))]
TF_scores_df_Ovaries_portal <- TF_scores_df[,grepl("Periportal_Ovaries",names(TF_scores_df))]
lipidomic_Ovaries <- lipidomic[,gsub("_Peri.*_Ovaries","",names(scaled_RNA_ovaries_portal))]
names(scaled_RNA_ovaries_portal) <- gsub("_Peri.*_Ovaries","",names(scaled_RNA_ovaries_portal))
names(TF_scores_df_Ovaries_portal) <- gsub("_Peri.*_Ovaries","",names(TF_scores_df_Ovaries_portal))


cormat_Ovaries_portal <- as.data.frame(rbind(make_cor_network(lipidomic_Ovaries,scaled_RNA_ovaries_portal,theshold_cor = threshold_cor),
                                            make_cor_network(scaled_RNA_ovaries_portal,TF_scores_df_Ovaries_portal, theshold_cor = threshold_cor, TF = T)))

names(cormat_Ovaries_portal) <- c("node1","node2","weight")
cormat_Ovaries_portal_ATT <- expand_attributes(cormat_Ovaries_portal,t_table,TF_scores_df)
cormat_Ovaries_portal_ATT$type <- ifelse(cormat_Ovaries_portal_ATT$gene %in% row.names(TF_scores_df),"TF",
                                      ifelse(cormat_Ovaries_portal_ATT$gene %in% t_table$gene,"gene","lipid"))

write_csv(cormat_Ovaries_portal, file = "results/multi-omic/cormat_Ovaries_portal_SIF.csv")
write_csv(cormat_Ovaries_portal_ATT, file = "results/multi-omic/cormat_Ovaries_portal_ATT.csv")

cormat_Ovaries_portal$node2 <- gsub("[(].*","",cormat_Ovaries_portal$node2)
cormat_Ovaries_portal$sign <- sign(cormat_Ovaries_portal$weight)
cormat_Ovaries_portal <- aggregate(weight ~ node1+node2+sign, data = cormat_Ovaries_portal, FUN = mean)
cormat_Ovaries_portal$abs_weight <- abs(cormat_Ovaries_portal$weight)

cormat_Ovaries_portal_ATT$gene <- gsub("[(].*","",cormat_Ovaries_portal_ATT$gene)
cormat_Ovaries_portal_ATT <- aggregate(cbind(portal_vs_Central,testes_vs_ova,portalova_vs_centralova,portaltestes_vs_centraltestes,centraltestes_vs_centralova,portaltestes_vs_portalova,testes_effect_zonation) ~ gene+type, data = cormat_Ovaries_portal_ATT, FUN = mean)

write_csv(cormat_Ovaries_portal, file = "results/multi-omic/compressed_cormat_Ovaries_portal_SIF.csv")
write_csv(cormat_Ovaries_portal_ATT, file = "results/multi-omic/compressed_cormat_Ovaries_portal_ATT.csv")

#4#######################################################################################################################

#ova_central
scaled_RNA_ova_central <- scaled_RNA[ovaries_pericentral_genes$gene,grepl("Pericentral_Ovaries",names(scaled_RNA))]
TF_scores_df_ova_central <- TF_scores_df[,grepl("Pericentral_Ovaries",names(TF_scores_df))]
lipidomic_Ovaries <- lipidomic[,gsub("_Peri.*_Ovaries","",names(scaled_RNA_ova_central))]
names(scaled_RNA_ova_central) <- gsub("_Peri.*_Ovaries","",names(scaled_RNA_ova_central))
names(TF_scores_df_ova_central) <- gsub("_Peri.*_Ovaries","",names(TF_scores_df_ova_central))


cormat_ova_central <- as.data.frame(rbind(make_cor_network(lipidomic_Ovaries, scaled_RNA_ova_central,theshold_cor = threshold_cor),
                                             make_cor_network(scaled_RNA_ova_central,TF_scores_df_ova_central, theshold_cor = threshold_cor, TF = T)))
names(cormat_ova_central) <- c("node1","node2","weight")
cormat_ova_central_ATT <- expand_attributes(cormat_ova_central,t_table,TF_scores_df)
cormat_ova_central_ATT$type <- ifelse(cormat_ova_central_ATT$gene %in% row.names(TF_scores_df),"TF",
                                      ifelse(cormat_ova_central_ATT$gene %in% t_table$gene,"gene","lipid"))

write_csv(cormat_ova_central, file = "results/multi-omic/cormat_ova_central_SIF.csv")
write_csv(cormat_ova_central_ATT, file = "results/multi-omic/cormat_ova_central_ATT.csv")

cormat_ova_central$node2 <- gsub("[(].*","",cormat_ova_central$node2)
cormat_ova_central$sign <- sign(cormat_ova_central$weight)
cormat_ova_central <- aggregate(weight ~ node1+node2+sign, data = cormat_ova_central, FUN = mean)
cormat_ova_central$abs_weight <- abs(cormat_ova_central$weight)

cormat_ova_central_ATT$gene <- gsub("[(].*","",cormat_ova_central_ATT$gene)
cormat_ova_central_ATT <- aggregate(cbind(portal_vs_Central,testes_vs_ova,portalova_vs_centralova,portaltestes_vs_centraltestes,centraltestes_vs_centralova,portaltestes_vs_portalova,testes_effect_zonation) ~ gene+type, data = cormat_ova_central_ATT, FUN = mean)

write_csv(cormat_ova_central, file = "results/multi-omic/compressed_cormat_ova_central_SIF.csv")
write_csv(cormat_ova_central_ATT, file = "results/multi-omic/compressed_cormat_ova_central_ATT.csv")
