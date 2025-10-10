library(Seurat)
library(readr)

gdxhr_heps_zonationsubclustered <- readRDS("data/RNA/gdxhr_heps_zonationsubclustered_pseudobulk_group_20241028.rds")

scaled_RNA <- as.data.frame(gdxhr_heps_zonationsubclustered@assays$RNA@layers$scale.data)

row.names(scaled_RNA) <- Features(gdxhr_heps_zonationsubclustered[["RNA"]])

meta_data_df <- as.data.frame(gdxhr_heps_zonationsubclustered@meta.data)
names(meta_data_df)[1] <- "sample"
meta_data_df$sample <- gsub("-","_",meta_data_df$sample)
meta_data_df$sample_name <- gsub("-","_",meta_data_df$sample_name)

names(scaled_RNA) <- meta_data_df$sample
scaled_RNA$gene <- row.names(scaled_RNA)
scaled_RNA <- scaled_RNA[,c(77,1:76)]

write_csv(scaled_RNA, file = "data/RNA/scaled_RNA_zonation.csv")

write_csv(meta_data_df, file = "data/RNA/meta_data_df_zonation.csv")
