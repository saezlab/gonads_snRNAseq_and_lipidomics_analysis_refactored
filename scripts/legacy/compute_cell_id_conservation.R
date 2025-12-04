library(readr)
# library(Seurat)
library(limma)
library(readxl)
library(ggplot2)
library(ggrepel)
library(readr)
library(cowplot)
library(decoupleR)
library(ggridges)
library(dplyr)
library(viridis)
library(reshape2)
library(gridExtra)
library(grid)

source("scripts/support_functions.R")


Suppl_Table_1_OverviewData <- as.data.frame(read_excel("data/Suppl_Table_1_OverviewData.xlsx"))
Suppl_Table_1_OverviewData <- Suppl_Table_1_OverviewData[,-4]

#I sploit the data in three parts so that the file size doesn't exceed github limit
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

meta_data_df <- as.data.frame(read_csv("data/RNA/meta_data_df.csv"))

meta_data_df <- merge(meta_data_df, Suppl_Table_1_OverviewData, by.x = "sample_name",by.y = "sequencing_name")

scaled_RNA <- scaled_RNA[,which(grepl("Hepatocyte",names(scaled_RNA)))]

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74767
SC3seq_Ms_ProcessedData <- as.data.frame(
  read_delim("data/RNA/SC3seq_Ms_ProcessedData.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE))
mouse_ESC <- SC3seq_Ms_ProcessedData[,c(2,which(grepl("2iLESC",names(SC3seq_Ms_ProcessedData))))]

validation_RNA <- as.data.frame(
  read_csv("data/RNA/merged_counts_matrix_250508_VH00211_433_AACVFM3HV_with_gene_symbols_and_groups.csv"))
#only control samples
validation_RNA <- validation_RNA[,c(1:5)]
names(validation_RNA) <- c("gene","in_vitro_hepa_neg_1","in_vitro_hepa_neg_2","in_vitro_hepa_neg_3","in_vitro_hepa_neg_4")

full_count <- merge(scaled_RNA, validation_RNA, by.x = "row.names", by.y = "gene")
names(full_count)[1] <- "gene"
full_count <- merge(full_count, mouse_ESC, by.x = "gene", by.y = "mm10_gene_symbol")

row.names(full_count) <- full_count$gene


# 1. Pull out the three groups of columns
hep_cols   <- grep("^sx42_.*_Hepatocytes_",   colnames(full_count), value=TRUE)
vitro_cols <- grep("^in_vitro_hepa_neg_",     colnames(full_count), value=TRUE)
stem_cols  <- grep("^X2iLESC_",                colnames(full_count), value=TRUE)

# 2. Subset to common genes (if you have any other references, intersect here)
expr_mat <- as.matrix(full_count[ , c(hep_cols, vitro_cols, stem_cols)])
expr_mat <- apply(expr_mat,2,as.numeric)
rownames(expr_mat) <- full_count$gene

# 3. Normalize across all samples to make them comparable
#    Here we use between-arrays quantile normalization from limma
expr_norm <- normalizeBetweenArrays(expr_mat, method="quantile")
# expr_norm <- apply(expr_mat,2,rank)

# 4. Build the two “centroids” (mean profile) for your anchors
centroid_hep  <- rowMeans(expr_norm[, hep_cols])
centroid_stem <- rowMeans(expr_norm[, stem_cols])

# 5. Compute Pearson correlations of each sample to each centroid
cor_to_hep  <- apply(expr_norm, 2, function(s) cor(s, centroid_hep,  method="pearson"))
cor_to_stem <- apply(expr_norm, 2, function(s) cor(s, centroid_stem, method="pearson"))

# 6. Convert into a 0–1 score:
#    score = cor_to_hep / (cor_to_hep + cor_to_stem)
identity_score <- cor_to_hep / (cor_to_hep + cor_to_stem)

# 7. Inspect your in vitro samples
cor_to_hep[vitro_cols]
cor_to_stem[vitro_cols]
identity_score[vitro_cols]




