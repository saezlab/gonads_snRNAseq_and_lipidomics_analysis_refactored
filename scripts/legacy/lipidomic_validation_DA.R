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
source("scripts/support_functions.R")

Suppl_Table_1_OverviewData <- as.data.frame(
  read_excel("data/Suppl_Table_1_OverviewData.xlsx"))

#cleaned according to from https://github.com/saezlab/liver_sex_differences
lipidomic <- read_csv("data/metabolomics/lipidomics_expt2.csv")

lipidomic <- lipidomic[,c(4,43:63)]

lipidomic <- as.data.frame(clean_lipid_data(lipidomic))

row.names(lipidomic) <- lipidomic$Lipid
lipidomic <- as.data.frame(lipidomic[,-1])

targets <- Suppl_Table_1_OverviewData[,c(2,3,4,5)]
names(targets)[1] <- "sample"
targets <- targets[targets$sample %in% names(lipidomic),]
lipidomic <- lipidomic[,targets$sample]
lipidomic <- log2(lipidomic)

unique(targets$gonadal_type)

plot(nicePCA(lipidomic, targets = targets[,c(1,3)]))

targets <- targets[,c(1,3)]
names(targets) <- c("sample","condition")

comparisons <- list("testes_vs_ova" = c(2,-1))

limmaRes <- runLimma(lipidomic, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(lipidomic[,1]),
                                    adjust.method = "fdr")

tval_df_original <- ttop_list[[1]][,c(1,4)]

##############################

Suppl_Table_1_OverviewData <- as.data.frame(
  read_excel("data/Suppl_Table_1_OverviewData.xlsx"))

#from https://github.com/saezlab/liver_sex_differences
lipidomic <- as.data.frame(
  read_csv("data/metabolomics/Lipid_data_file_JPC_cleaned.csv"))
row.names(lipidomic) <- lipidomic$`Unnamed: 0`
lipidomic <- as.data.frame(t(lipidomic[,-1]))

lipidomic <- log2(lipidomic)
names(lipidomic) <- gsub("-0","-",names(lipidomic))
names(lipidomic) <- gsub("-","_",names(lipidomic))

lipidomic <- lipidomic[,targets$sample]


limmaRes <- runLimma(lipidomic, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(lipidomic[,1]),
                                    adjust.method = "fdr")

tval_df_old <- ttop_list[[1]][,c(1,4)]


t_comparison <- merge(tval_df_original,tval_df_old, by = "ID")

plot(t_comparison$t.x, t_comparison$t.y)

t_comparison_TG <- t_comparison[grepl("TG[(]",t_comparison$ID),]

plot(t_comparison_TG$t.x, t_comparison_TG$t.y)
cor.test(t_comparison_TG$t.x, t_comparison_TG$t.y)

# t_comparison_TG$residuals <- summary(lm(t.x~t.y, data = t_comparison_TG))$residuals

t_comparison_TG$bad <- abs(t_comparison_TG$t.y) > 4 & abs(t_comparison_TG$t.x) < 2

hepatocyte_TF_spec_gene_lipid_network_mofa_ATT <- as.data.frame(
  read_csv("results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_ATT_full.csv"))

sum(t_comparison_TG[t_comparison_TG$bad,"ID"] %in% hepatocyte_TF_spec_gene_lipid_network_mofa_ATT$node)
sum(t_comparison_TG[!t_comparison_TG$bad,"ID"] %in% hepatocyte_TF_spec_gene_lipid_network_mofa_ATT$node)


t_comparison_DG <- t_comparison[grepl("[TD]G[(]",t_comparison$ID),]

plot(t_comparison_DG$t.x, t_comparison_DG$t.y)
cor.test(t_comparison_DG$t.x, t_comparison_DG$t.y)

names(t_comparison_DG) <- c("ID","validation","fcgl")

t_comparison_DG$labels <- ifelse(
  abs(t_comparison_DG$validation) > 4 & 
    abs(t_comparison_DG$fcgl) > 2.5,
  t_comparison_DG$ID,
  NA)


ggplot(t_comparison_DG, aes(x = fcgl, y = validation, label = labels)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = T, color = "royalblue") +
  geom_label_repel() +
  theme_minimal()
