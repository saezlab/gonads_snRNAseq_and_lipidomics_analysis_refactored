library(readr)
library(vsn)
library(limma)
library(ggplot2)
library(ggrepel)

source("scripts/support_functions.R")

proteo <- as.data.frame(
  read.delim("data/proteomics/P3225_fractionated_protein.tsv"))

proteo <- proteo[-which(grepl("contam",proteo[,1])),]
proteo <- proteo[,c(4,24:39)]

dubs <- proteo[duplicated(proteo$Gene),1]

proteo <- proteo[-which(proteo$Gene %in% dubs),]

row.names(proteo) <- proteo$Gene

proteo <- proteo[,-1]

proteo[proteo == 0] <- NA

fit <- vsnMatrix(as.matrix(proteo))

# Create a plot of the mean versus the standard deviation of the VSN-transformed data
meanSdPlot(fit)

# Apply the VSN transformation to the original data frame
proteo <- as.data.frame(vsn::predict(fit,as.matrix(proteo)))


metadata <- as.data.frame(
  read_csv("data/proteomics/metadata.csv"))

targets <- metadata[,c(1,2)]
names(targets)[1] <- "sample"
unique(targets$condition)

proteo <- proteo[,targets$sample]

comparisons <- list("testes_vs_Ova" = c(1,-2))

limmaRes <- runLimma(proteo, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(proteo[,1]),
                                    adjust.method = "fdr")

ttop <- ttop_list$testes_vs_Ova

tval_full_df <- as.data.frame(read_csv("results/pseudobulk/tval_full_df.csv"))

t_val_compararison <- merge(ttop[,c(1,4)],tval_full_df[,c(1,2)])

plot(t_val_compararison$t, t_val_compararison$testes_vs_ova_tval_Hepatocytes)

names(t_val_compararison) <- c("ID","proteo","RNA")

t_val_compararison$labels <- ifelse(
  abs(t_val_compararison$proteo) > 20,
  t_val_compararison$ID,
  NA)

ggplot(t_val_compararison, aes(x = RNA, y = proteo, label = labels)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = T, color = "royalblue") +
  geom_label_repel() +
  theme_minimal()

cor.test(t_val_compararison$proteo, t_val_compararison$RNA)

load("results/pseudobulk/ttop_full_list.RData")

ttop_RNA <- ttop_full_list$testes_vs_ova$Hepatocytes

pval_merged <- merge(ttop[,c(1,6)], ttop_RNA[,c(1,6)], by = "ID")

sum(pval_merged$adj.P.Val.x <= 0.05 & pval_merged$adj.P.Val.y <= 0.05) / sum(pval_merged$adj.P.Val.y <= 0.05)

common_sig_genes <- pval_merged[pval_merged$adj.P.Val.x <= 0.05 & pval_merged$adj.P.Val.y <= 0.05,1]

t_val_compararison$common_sig_genes <- t_val_compararison$ID %in% common_sig_genes

t_val_compararison$common_sig_genes <- ifelse(sign(t_val_compararison$proteo) == sign(t_val_compararison$RNA),
                                              t_val_compararison$common_sig_genes,FALSE)

sum(t_val_compararison$common_sig_genes) / sum(pval_merged$adj.P.Val.y <= 0.05)

ggplot(t_val_compararison, aes(x = RNA, y = proteo, label = labels, color = common_sig_genes)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = T, color = "royalblue") +
  geom_label_repel() +
  theme_minimal()
