library(Seurat)
library(readr)
library(ggplot2)

fcgl_pseudobulk_gonad <- readRDS("data/RNA/fcgl_pseudobulks_gonad_withxytwt.rds")

scaled_RNA <- as.data.frame(fcgl_pseudobulk_gonad@assays$RNA@layers$scale.data)

row.names(scaled_RNA) <- Features(fcgl_pseudobulk_gonad[["RNA"]])

meta_data_df <- as.data.frame(fcgl_pseudobulk_gonad@meta.data)
names(meta_data_df)[1] <- "sample"
meta_data_df$sample <- gsub("-","_",meta_data_df$sample)
meta_data_df$sample_name <- gsub("-","_",meta_data_df$sample_name)

names(scaled_RNA) <- meta_data_df$sample
scaled_RNA$gene <- row.names(scaled_RNA)
scaled_RNA <- scaled_RNA[,c(241,1:240)]

write_csv(scaled_RNA[,c(1:80)], file = "data/RNA/scaled_RNA_withWT_1.csv")
write_csv(scaled_RNA[,c(1,81:160)], file = "data/RNA/scaled_withWT_RNA_2.csv")
write_csv(scaled_RNA[,c(1,161:241)], file = "data/RNA/scaled_withWT_RNA_3.csv")

write_csv(meta_data_df, file = "data/RNA/meta_data_df_withWT.csv")