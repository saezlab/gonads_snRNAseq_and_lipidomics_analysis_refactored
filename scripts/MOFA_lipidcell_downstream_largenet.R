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

n_upstream <- 500
top_percentile <- 0.33

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
    top_n <- round(length(x) * top_percentile)
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

hepatocyte_cor_mat_RNA_lip_sparse <- as.data.frame(t(apply(hepatocyte_cor_mat_RNA_lip, 1, function(x){
  x_top <- x[order(abs(x), decreasing = T)]
  x_top <- x_top[n_upstream] #for top 3
  # x_top <- x_top[20]
  x[abs(x) < abs(x_top)] <- NA
  return(x)
})))


hepatocyte_cor_mat_RNA_lip_network <- melt(cbind(hepatocyte_cor_mat_RNA_lip_sparse, row.names(hepatocyte_cor_mat_RNA_lip_sparse)))
hepatocyte_cor_mat_RNA_lip_network <- hepatocyte_cor_mat_RNA_lip_network[complete.cases(hepatocyte_cor_mat_RNA_lip_network),]
names(hepatocyte_cor_mat_RNA_lip_network) <- c("node1","node2","mor")

##TF analysis
load("support/TF_regulons.RData")
TF_regulons <- TF_regulons[complete.cases(TF_regulons),]

Gonadal_Hormone_Genes <- as.data.frame(read_csv("support/Gonadal_Hormone_Genes.csv"))
Gonadal_Hormone_Genes$source <- gsub("_spec","",Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("ESR1","Esr1",Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("AR","Ar",Gonadal_Hormone_Genes$source)

TF_regulons <- TF_regulons[!TF_regulons$source %in% c("Ar", "Esr1"),]

TF_regulons <- rbind(TF_regulons, Gonadal_Hormone_Genes)

TF_lipid_factor_n <- data.frame(RNA_lipid_weights$Hepatocytes_TFs[[factor_n]])
TF_lipid_factor_n <- TF_lipid_factor_n[complete.cases(TF_lipid_factor_n),,drop = F]
#hepatocytes
TF_regulons_hepatocytes <- TF_regulons[TF_regulons$source %in% gsub("_Hepatocytes","",row.names(TF_lipid_factor_n)),]
TF_regulons_hepatocytes <- TF_regulons_hepatocytes[TF_regulons_hepatocytes$target %in% hepatocyte_cor_mat_RNA_lip_network$enzyme,]

TF_scores <- decoupleR::run_ulm(mat = scaled_RNA, network = TF_regulons[,], minsize = 20)
TF_scores_samples <- dcast(TF_scores, source~condition, value.var = "score")


TF_scores_samples <- TF_scores_samples[TF_scores_samples$source %in% gsub("_Hepatocytes_TFs","",row.names(TF_lipid_factor_n)),]
row.names(TF_scores_samples) <- TF_scores_samples$source
TF_scores_samples <- TF_scores_samples[,-1]


TF_scores_samples_celltype <- TF_scores_samples[,grepl(cell_type,names(TF_scores_samples))]
names(TF_scores_samples_celltype) <- gsub(cell_type,"",names(TF_scores_samples_celltype))

hepatocyte_cor_mat_RNA_TF <- as.data.frame(cor(t(rbind(TF_scores_samples_celltype, RNA_celltype))))
hepatocyte_cor_mat_RNA_TF <- hepatocyte_cor_mat_RNA_TF[row.names(TF_scores_samples_celltype), row.names(RNA_celltype)]
hepatocyte_cor_mat_RNA_TF <- hepatocyte_cor_mat_RNA_TF[,names(hepatocyte_cor_mat_RNA_TF) %in% hepatocyte_cor_mat_RNA_lip_network$node2]

hepatocyte_cor_mat_RNA_TF_sparse <- as.data.frame(apply(hepatocyte_cor_mat_RNA_TF, 2, function(x){
  x_top <- x[order(abs(x), decreasing = T)]
  x_top <- x_top[n_upstream] #for top 3
  x[abs(x) < abs(x_top)] <- NA
  return(x)
}))

hepatocyte_cor_mat_RNA_TF_network <- melt(cbind(hepatocyte_cor_mat_RNA_TF_sparse, row.names(hepatocyte_cor_mat_RNA_TF_sparse)))
hepatocyte_cor_mat_RNA_TF_network <- hepatocyte_cor_mat_RNA_TF_network[complete.cases(hepatocyte_cor_mat_RNA_TF_network),]
names(hepatocyte_cor_mat_RNA_TF_network) <- c("node1","node2","mor")

hepa_TF_gene_lipid_network <- as.data.frame(rbind(hepatocyte_cor_mat_RNA_lip_network, hepatocyte_cor_mat_RNA_TF_network))

nodes_lipnetwork <- unique(c(as.character(hepa_TF_gene_lipid_network$node1), as.character(hepa_TF_gene_lipid_network$node2)))
nodes_lipnetwork <- data.frame("node" = nodes_lipnetwork)
nodes_lipnetwork$type <- ifelse(grepl(":",nodes_lipnetwork$node), "lipid", 
                                ifelse(nodes_lipnetwork$node %in% TF_regulons$source, "TF", "enzyme"))

TF_lipid_factor_n$feature <- row.names(TF_lipid_factor_n)
names(TF_lipid_factor_n) <- c("weight","feature")
factor_weights_unique <- rbind(TF_lipid_factor_n, RNA_lipid_factor_n)
factor_weights_unique <- factor_weights_unique[!duplicated(factor_weights_unique$feature),]
nodes_lipnetwork <- merge(nodes_lipnetwork, factor_weights_unique, all.x = T, by.x = "node", by.y = "feature")

write_csv(hepa_TF_gene_lipid_network, file = "results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full_largenet.csv")
write_csv(nodes_lipnetwork, file = "results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_ATT_full_largenet.csv")
