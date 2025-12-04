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

source("scripts/support_functions.R")

lipidomic <- as.data.frame(
  read_csv("data/metabolomics/lipidomics_20250325transfections.csv"))

lipidomic <- as.data.frame(clean_lipid_data(lipidomic))

targets <- as.data.frame(cbind(names(lipidomic),names(lipidomic)))[-1,]
targets[,2] <- gsub("_.*","",targets[,2])
names(targets) <- c("sample","condition")

row.names(lipidomic) <- lipidomic$Lipid
lipidomic <- lipidomic[,-1]

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(lipidomic, targets, no_label = T))

unique(targets$condition)

comparisons <- list("sult3a1" = c(2,-1),
                    "abcb1a" = c(3,-1),
                    "cyp2a4" = c(4,-1),
                    "mtmr11" = c(5,-1),
                    "mindy3" = c(6,-1),
                    "esco2" = c(7,-1),
                    "kitl" = c(8,-1))

limmaRes <- runLimma(lipidomic, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(lipidomic[,1]),
                                    adjust.method = "fdr")

t_table <- ttop_list_to_t_table(ttop_list)
row.names(t_table) <- t_table$ID

t <- as.vector(t(t_table[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(t_table[,-1], show_rownames = F, cluster_cols = T, cluster_rows = T,color = palette, angle_col = 315, display_numbers = F)
dev.off()

testes_portal <- as.data.frame(read_csv("results/multi-omic/cormat_testes_portal_SIF.csv")) 
testes_portal <- testes_portal[grepl("[()]",testes_portal$node2),]
testes_portal$node1 <- paste(testes_portal$node1,"testes_portal",sep = "_")

testes_central <- as.data.frame(read_csv("results/multi-omic/cormat_testes_central_SIF.csv")) 
testes_central <- testes_central[grepl("[()]",testes_central$node2),]
testes_central$node1 <- paste(testes_central$node1,"testes_central",sep = "_")

ova_portal <- as.data.frame(read_csv("results/multi-omic/cormat_Ovaries_portal_SIF.csv")) 
ova_portal <- ova_portal[grepl("[()]",ova_portal$node2),]
ova_portal$node1 <- paste(ova_portal$node1,"ova_portal",sep = "_")

ova_central <- as.data.frame(read_csv("results/multi-omic/cormat_ova_central_SIF.csv")) 
ova_central <- ova_central[grepl("[()]",ova_central$node2),]
ova_central$node1 <- paste(ova_central$node1,"ova_central",sep = "_")

sszb_network <- as.data.frame(do.call(rbind,list(testes_portal,testes_central,ova_portal,ova_central)))
names(sszb_network) <- c("source","target","mor")

sszb_validation_scores <- decoupleR::run_ulm(t_table[,-1], network = sszb_network, minsize = 2)
sszb_validation_scores$gonade_zonation <- sub("^[^_]*_", "", sszb_validation_scores$source)
sszb_validation_scores$source <- tolower(gsub("_.*","",sszb_validation_scores$source))

validation_set_only <- sszb_validation_scores[sszb_validation_scores$source == sszb_validation_scores$condition,]

write_csv(validation_set_only, file = "results/in_vitro_validation/validation_scores_df_1_sszb.csv")


KO_index <- 8
zzsb <- as.character(validation_set_only[validation_set_only$source == names(t_table)[KO_index],"gonade_zonation"])
sszb_network_sub <- sszb_network[which(grepl(zzsb, sszb_network$source)),]
lipid_scoring_Scatter <- merge(t_table[,c(1,KO_index)], sszb_network_sub, by.x = "ID", by.y = "target")
lipid_scoring_Scatter$source <- gsub("_.*","",lipid_scoring_Scatter$source)
lipid_scoring_Scatter$mor <- ifelse(tolower(lipid_scoring_Scatter$source) == names(lipid_scoring_Scatter)[2],lipid_scoring_Scatter$mor,0)

names(lipid_scoring_Scatter) <- c("ID","KO_tval","KO_target","predicted_association_weight")
lipid_scoring_Scatter$ID <- ifelse(lipid_scoring_Scatter$predicted_association_weight != 0, lipid_scoring_Scatter$ID, NA)

ggplot(lipid_scoring_Scatter, aes(x = predicted_association_weight, y = KO_tval, label = ID)) + 
  geom_point() + 
  geom_label_repel() +
  ggtitle(paste(names(t_table)[KO_index],zzsb,sep = "_")) +
  theme_minimal() + 
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x)









