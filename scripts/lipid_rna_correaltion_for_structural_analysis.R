library(readr)
library(GSEABase)
library(reshape2)
library(decoupleR)
library(cosmosR)
library(dplyr)

source("scripts/support_functions.R")

load("results/mofa/main_model_RNA_lipid_weight.RData")
RNA_lipid_weights <- weights

rm(weights)

# RNA_lipid_weights <- lapply(RNA_lipid_weights, function(factor_weights){
#   apply(factor_weights,2,function(x){
#     x <- x[order(abs(x), decreasing = T)]
#     top_n <- round(length(x) * 0.05)
#     x[top_n:length(x)] <- NA
#     return(x)
#   }, simplify = F)
# })

RNA_lipid_weights <- lapply(RNA_lipid_weights, function(factor_weights){
  apply(factor_weights,2,function(x){
    x <- x[order(x, decreasing = T)]
    top_n <- round(length(x) * 0.5)
    print(top_n)
    x[top_n:(length(x) - top_n)] <- NA
    return(x)
  }, simplify = F)
})

factor_n <- 1

RNA_lipid_factor_n <- as.data.frame(c(
  RNA_lipid_weights$Hepatocytes[[factor_n]],
  RNA_lipid_weights$`Immune cells`[[factor_n]],
  RNA_lipid_weights$lipidomic[[factor_n]],
  RNA_lipid_weights$Hepatocytes_TFs[[factor_n]]))

RNA_lipid_factor_n <- RNA_lipid_factor_n[complete.cases(RNA_lipid_factor_n),,drop = F]
names(RNA_lipid_factor_n) <- "weight"
RNA_lipid_factor_n$feature <- gsub("_[A-Za-z].*","",row.names(RNA_lipid_factor_n))

#### corelation RNA_lipids
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

names(scaled_RNA) <- gsub("_Testes","",names(scaled_RNA))
names(scaled_RNA) <- gsub("_Ovaries","",names(scaled_RNA))

load("data/metabolomics/lipidomic.RData")

## Lipid hepatocyte RNA
cell_type <- "_Hepatocytes"

RNA_celltype <- scaled_RNA[,grepl(cell_type,names(scaled_RNA))]
names(RNA_celltype) <- gsub(cell_type,"",names(RNA_celltype))
RNA_celltype <- RNA_celltype[row.names(RNA_celltype) %in% RNA_lipid_factor_n[grepl(cell_type,row.names(RNA_lipid_factor_n)),"feature"],]

lipidomic_RNA <- lipidomic[row.names(lipidomic) %in% RNA_lipid_factor_n$feature,]

lipidomic_RNA <- lipidomic_RNA[,names(RNA_celltype)]

hepatocyte_cor_mat_RNA_lip <- as.data.frame(cor(t(rbind(lipidomic_RNA, RNA_celltype))))
hepatocyte_cor_mat_RNA_lip <- hepatocyte_cor_mat_RNA_lip[row.names(lipidomic_RNA), row.names(RNA_celltype)]

##need to extract the top correlated genes for each lipid 
pathways <- import_gmt("support/m5.all.v2023.2.Mm.symbols.gmt")
names(pathways) <- c("target","source")

lipid_genes <- pathways[pathways$source == "GOBP_CELLULAR_LIPID_METABOLIC_PROCESS","target"]

hepatocyte_cor_mat_RNA_lip <- hepatocyte_cor_mat_RNA_lip[,names(hepatocyte_cor_mat_RNA_lip) %in% lipid_genes]

save(hepatocyte_cor_mat_RNA_lip, file = "results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata")
