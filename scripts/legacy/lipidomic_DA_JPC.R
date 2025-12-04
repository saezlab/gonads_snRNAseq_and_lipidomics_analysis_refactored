library(readr)
library(readxl)
library(pheatmap)
library(limma)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(cowplot)
library(gridExtra)
library(decoupleR)
library(grid)

setwd("~/Dropbox/postdoc/Manuscripts/MS_FourCoreGenotypesLiver/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis")
source("scripts/support_functions.R")

Suppl_Table_1_OverviewData <- as.data.frame(
  read_excel("data/Suppl_Table_1_OverviewData.xlsx"))

#from https://github.com/saezlab/liver_sex_differences
lipidomic <- as.data.frame(
  read_csv("data/metabolomics/Lipid_data_file_JPC_cleaned.csv"))
row.names(lipidomic) <- lipidomic$`Unnamed: 0`
lipidomic <- as.data.frame(t(lipidomic[,-1]))

pheatmap(cor(log2(lipidomic)))
pheatmap(log2(lipidomic), show_rownames = F)

targets_all <- Suppl_Table_1_OverviewData[,c(2,3,4,5)]
names(targets_all)[1] <- "sample"

lipidomic <- log2(lipidomic)
names(lipidomic) <- gsub("-0","-",names(lipidomic))
names(lipidomic) <- gsub("-","_",names(lipidomic))

lipidomic <- lipidomic[,targets_all$sample]

targets_all_noWT <- targets_all[!grepl("-WT",targets_all$group),]

color_mapping <- c("Ovaries" = "#ED2024",
                   "Testes" = "#4D69B1")

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(lipidomic, targets_all_noWT[,c(1,3)], no_label = T, color_mapping = color_mapping))

pcares <- prcomp(t(lipidomic))
sample_coordinates <- pcares$x[,-48]

to_write <- as.data.frame(sample_coordinates)
to_write$samples <- row.names(to_write)

write_csv(to_write, file = "results/lipidomic/PCA_sample_coord_lipid.csv")

to_enrich <- melt(targets_all, id.vars = "sample")
to_enrich <- to_enrich[,-2]
names(to_enrich) <- c("target","source")
to_enrich$mor <- 1

factor_enrichment <- decoupleR::run_ulm(sample_coordinates, to_enrich)
factor_enrichment <- dcast(factor_enrichment, source~condition, value.var = "score")
row.names(factor_enrichment) <- factor_enrichment[,1] 
factor_enrichment <- factor_enrichment[,-1]

#rescale with variance epxlained
factor_enrichment <- factor_enrichment * c((pcares$sdev)^2 / sum(pcares$sdev^2))[-48]
factor_enrichment <- factor_enrichment[,order(as.numeric(gsub("PC","",names(factor_enrichment))))]

t <- as.vector(t(factor_enrichment))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(factor_enrichment, display_numbers = F,
                   show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315)

factor_enrichment_long <- factor_enrichment
factor_enrichment_long$condition <- row.names(factor_enrichment_long)

factor_enrichment_long <- melt(factor_enrichment_long)

ggplot(factor_enrichment_long, aes(x = variable, y = value, group = condition, color = condition)) +
  geom_line(size = 2, alpha = 0.5) + 
  geom_point(size = 3) +
  theme_minimal()

write_csv(factor_enrichment_long, file = "results/lipidomic/factor_enrichment_long.csv")

save(lipidomic, file = "data/metabolomics/lipidomic.RData")

group_variable <- "gonadal_type" #gonadal_type or group or allosome_type

targets <- targets_all[,c("sample",group_variable)]
names(targets)[2] <- "condition"
print(unique(targets$condition))

# comparisons <- list("testes_vs_ova" = c(1,-2))
comparisons <- list("testes_vs_ova" = c(-1,2))

limmaRes <- runLimma(lipidomic, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(lipidomic[,1]),
                                    adjust.method = "fdr")

gonade_type_ttop <- ttop_list[[1]]

write_csv(gonade_type_ttop, file = "results/lipidomic/gonade_type_ttop.csv")

volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40)
volca_plt

name_of_file <- paste("results/lipidomic/volcanos/",paste(names(comparisons),".pdf",sep = ""),sep = "")
pdf_size <- 9
ggsave(filename = name_of_file, plot = volca_plt, device = "pdf", width = pdf_size, height = pdf_size)
gonade_type_ttop <- gonade_type_ttop[order(gonade_type_ttop$P.Value, decreasing = F),]
gonade_type_ttop_20 <- as.data.frame(gonade_type_ttop[1:20,])
write_csv(gonade_type_ttop_20, file = "results/lipidomic/gonade_type_ttop_20.csv")
###

group_variable <- "allosome_type" #gonadal_type or group or allosome_type

targets <- targets_all[,c("sample",group_variable)]
names(targets)[2] <- "condition"
print(unique(targets$condition))

comparisons <- list("XX_vs_XY" = c(1,-2))

limmaRes <- runLimma(lipidomic, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(lipidomic[,1]),
                                    adjust.method = "fdr")

allosome_type_ttop <- ttop_list[[1]]
volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40)
volca_plt

name_of_file <- paste("results/lipidomic/volcanos/",paste(names(comparisons),".pdf",sep = ""),sep = "")
pdf_size <- 9
ggsave(filename = name_of_file, plot = volca_plt, device = "pdf", width = pdf_size, height = pdf_size)


###
group_variable <- "group" #gonadal_type or group or allosome_type

targets <- targets_all[,c("sample",group_variable)]
names(targets)[2] <- "condition"
print(unique(targets$condition))

comparisons <- list("XXO_vs_XXT" = c(1,-4),
                    "XYT_vs_XYT-WT" = c(2,-5),
                    "ovaries_vs_XYT" = c(2,2,-1,-3),
                    "ovaries_vs_XYT-WT" = c(5,5,-1,-3))

limmaRes <- runLimma(lipidomic, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(lipidomic[,1]),
                                    adjust.method = "fdr")

group_ttop <- ttop_list[[1]]
volca_plt <- volcano_nice(ttop_list[[1]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40)
volca_plt

name_of_file <- paste("results/lipidomic/volcanos/",paste(names(comparisons)[1],".pdf",sep = ""),sep = "")
pdf_size <- 9
ggsave(filename = name_of_file, plot = volca_plt, device = "pdf", width = pdf_size, height = pdf_size)

group_ttop <- ttop_list[[2]]
volca_plt <- volcano_nice(ttop_list[[2]], FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 40)
volca_plt

name_of_file <- paste("results/lipidomic/volcanos/",paste(names(comparisons)[2],".pdf",sep = ""),sep = "")
pdf_size <- 9
ggsave(filename = name_of_file, plot = volca_plt, device = "pdf", width = pdf_size, height = pdf_size)

t_val_fcgl_vs_WT <- merge(ttop_list[[3]][,c(1,4)],ttop_list[[4]][,c(1,4)], by = "ID")
names(t_val_fcgl_vs_WT) <- c("ID",names(comparisons)[c(3,4)])
t_val_fcgl_vs_WT$diff <- abs(t_val_fcgl_vs_WT[,2] - t_val_fcgl_vs_WT[,3])
t_val_fcgl_vs_WT$label <- ifelse(t_val_fcgl_vs_WT$diff > 6, t_val_fcgl_vs_WT$ID, NA)

ggplot(t_val_fcgl_vs_WT, aes(x = `ovaries_vs_XYT-WT`, y = ovaries_vs_XYT, label = label)) + 
  geom_point() +
  theme_minimal() +
  stat_summary(fun.data=mean_cl_normal) +
  geom_smooth(method='lm', formula= y~x) + 
  geom_label_repel() + ggtitle("Pearson corelation coeficient = 0.88")

cor.test(t_val_fcgl_vs_WT$ovaries_vs_XYT, t_val_fcgl_vs_WT$`ovaries_vs_XYT-WT`)

ttop_full_list <- list("testes_vs_ova" = gonade_type_ttop, "XX_vs_XY" = allosome_type_ttop, "XXO_vs_XXT" = group_ttop)
save(ttop_full_list, file = "results/lipidomic/ttop_full_list.RData")

sum(ttop_full_list$testes_vs_ova$adj.P.Val <= 0.05)
sum(ttop_full_list$XX_vs_XY$adj.P.Val <= 0.05)


sum(ttop_full_list$testes_vs_ova$adj.P.Val <= 0.05 & ttop_full_list$testes_vs_ova$t > 0)
sum(ttop_full_list$testes_vs_ova$adj.P.Val <= 0.05 & ttop_full_list$testes_vs_ova$t < 0)
### DO THE LIPID CLASSE ENRICHMENT
#CONNECT LIPID CLASSE WITH ENZYMES

lipid_classes <- as.data.frame(cbind(ttop_full_list$testes_vs_ova$ID,gsub("[(].*","",ttop_full_list$testes_vs_ova$ID)))
names(lipid_classes) <- c("target","source")
barplot(table(lipid_classes))
lipid_classes$mor <- 1

t_table <- merge(ttop_full_list$testes_vs_ova[,c(1,4)], ttop_full_list$XX_vs_XY[,c(1,4)], by = "ID")
t_table <- merge(t_table, ttop_full_list$XXO_vs_XXT[,c(1,4)])
names(t_table) <- c("ID",names(ttop_full_list))
row.names(t_table) <- t_table$ID

lipid_classe_enrichment <- decoupleR::run_ulm(mat = t_table[,-1], network = lipid_classes, minsize = 10)
lipid_classe_enrichment <- dcast(lipid_classe_enrichment, source~condition, value.var = "score")
row.names(lipid_classe_enrichment) <- lipid_classe_enrichment$source

lipid_counts_in_classes <- table(lipid_classes$source)
lipid_counts_in_classes <- lipid_counts_in_classes[lipid_counts_in_classes >= 10]

row.names(lipid_classe_enrichment) <- paste(row.names(lipid_classe_enrichment), as.numeric(lipid_counts_in_classes), sep = " n=")

t <- as.vector(t(lipid_classe_enrichment[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(lipid_classe_enrichment[,-1], display_numbers = T,
                   show_rownames = T, cluster_cols = T, cluster_rows = T,color = palette, angle_col = 315)

lipid_classe_enrichment$source

lipid_classification <- as.data.frame(read_csv("support/lipid_classification.csv"))

lipid_enrichment_gonade <- lipid_classe_enrichment[,c(1,2)]
lipid_enrichment_gonade$n <- row.names(lipid_enrichment_gonade)

lipid_enrichment_gonade <- merge(lipid_enrichment_gonade, lipid_classification, by.x = "source", by.y = "Lipid")
lipid_enrichment_gonade <- lipid_enrichment_gonade[order(abs(lipid_enrichment_gonade$testes_vs_ova), decreasing = F),]
lipid_enrichment_gonade$n <- factor(lipid_enrichment_gonade$n, levels = lipid_enrichment_gonade$n)

ggplot(lipid_enrichment_gonade, aes(x = n, y = testes_vs_ova, color = testes_vs_ova)) + 
  geom_point(size = 4) + 
  geom_segment(aes(x=n, xend=n, y=0, yend=testes_vs_ova), size = 1) + 
  coord_flip() + 
  theme_minimal() +
  scale_color_gradient2(
    midpoint = 0,    # zero is the midpoint
    low = "red",     # negative values
    mid = "grey",   # zero
    high = "royalblue"    # positive values
  )

ggplot(lipid_enrichment_gonade, aes(x = n, y = testes_vs_ova, color = testes_vs_ova > 0)) + 
  geom_point(size = 4) + 
  geom_segment(aes(x = n, xend = n, y = 0, yend = testes_vs_ova), size = 1) + 
  scale_color_manual(values = c("TRUE" = "royalblue", "FALSE" = "red")) +
  coord_flip() + 
  theme_minimal()

ggplot(lipid_enrichment_gonade, aes(x = n, y = testes_vs_ova, color = testes_vs_ova > 0)) + 
  geom_point(size = 4) + 
  geom_segment(aes(x = n, xend = n, y = 0, yend = testes_vs_ova), size = 1) + 
  scale_color_manual(values = c("TRUE" = "royalblue", "FALSE" = "red")) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal()


#from https://github.com/saezlab/liver_sex_differences/blob/main/linex_export.ipynb (without the filtering step)
linex2_adjacency_matrix_OT <- as.data.frame(read_csv("support/linex2_adjacency_matrix_OT.csv"))
row.names(linex2_adjacency_matrix_OT) <- linex2_adjacency_matrix_OT$...1
linex2_adjacency_matrix_OT <- linex2_adjacency_matrix_OT[,-1]

linex2_SIF <- melt(linex2_adjacency_matrix_OT)
linex2_SIF$lipid_B <- rep(names(linex2_adjacency_matrix_OT), length(linex2_adjacency_matrix_OT[1,]))
names(linex2_SIF)[1] <- "lipid_A"
linex2_SIF <- linex2_SIF[linex2_SIF$value != 0,]
linex2_SIF <- linex2_SIF[,c(1,3,2)]
names(linex2_SIF)[3] <- "sign"

significant_lipids <- gonade_type_ttop[gonade_type_ttop$P.Value < 10e-7,]
significant_lipids_network <- linex2_SIF[linex2_SIF$lipid_A %in% significant_lipids$ID &
                                           linex2_SIF$lipid_B %in% significant_lipids$ID,]

significant_lipids_network$to_keep <- TRUE
for(i in 1:length(significant_lipids_network[,1]))
{
  edge_ID <- paste(significant_lipids_network[i,1],significant_lipids_network[i,2],sep = "_")
  print(edge_ID)
  edges_ID_to_compare <- paste(significant_lipids_network[1:i,2],significant_lipids_network[1:i,1], sep = "_")
  significant_lipids_network[i,"to_keep"] <- !(edge_ID %in% edges_ID_to_compare)
}
significant_lipids_network <- significant_lipids_network[significant_lipids_network$to_keep, -4]

write_csv(significant_lipids_network, file = "results/lipidomic/significant_lipids_network_SIF.csv")
write_csv(significant_lipids, file = "results/lipidomic/significant_lipids_network_ATT.csv")
