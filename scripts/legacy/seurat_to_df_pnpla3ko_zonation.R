library(Seurat)
library(readr)

fcgl_pseudobulk_gonad <- readRDS("data/RNA/pnpla3ko_sci_heps_zonationsubsetted_pseudobulks.rds")

scaled_RNA <- as.data.frame(fcgl_pseudobulk_gonad@assays$RNA@layers$scale.data)

row.names(scaled_RNA) <- Features(fcgl_pseudobulk_gonad[["RNA"]])

meta_data_df <- as.data.frame(fcgl_pseudobulk_gonad@meta.data)
names(meta_data_df)[1] <- "sample"
meta_data_df$sample <- gsub("-","_",meta_data_df$sample)
meta_data_df$sample_name <- gsub("-","_",meta_data_df$sample_name)

names(scaled_RNA) <- meta_data_df$sample
scaled_RNA$gene <- row.names(scaled_RNA)
scaled_RNA <- scaled_RNA[,c(43,1:42)]

write_csv(scaled_RNA, file = "data/RNA/scaled_RNA_pnpla3ko_zonation.csv")

write_csv(meta_data_df, file = "data/RNA/meta_data_df_pnpla3ko_zonation.csv")
