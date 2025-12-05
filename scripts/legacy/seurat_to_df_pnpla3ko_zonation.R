library(Seurat)
library(readr)

# Load compatibility functions for Seurat v4/v5
source("scripts/seurat_compatibility.R")

fcgl_pseudobulk_gonad <- readRDS("data/RNA/pnpla3ko_sci_heps_zonationsubsetted_pseudobulks.rds")

# Extract scaled data (compatible with both Seurat v4 and v5)
scaled_RNA <- as.data.frame(get_layer_data(fcgl_pseudobulk_gonad, assay = "RNA", layer = "scale.data"))

row.names(scaled_RNA) <- Features(fcgl_pseudobulk_gonad[["RNA"]])

meta_data_df <- as.data.frame(fcgl_pseudobulk_gonad@meta.data)
names(meta_data_df)[1] <- "sample"
meta_data_df$sample <- gsub("-","_",meta_data_df$sample)
meta_data_df$sample_name <- gsub("-","_",meta_data_df$sample_name)

names(scaled_RNA) <- meta_data_df$sample
scaled_RNA$gene <- row.names(scaled_RNA)
scaled_RNA <- scaled_RNA[,c(43,1:42)]

write_csv(scaled_RNA, file = "data/RNA/scaled_RNA_pnpla3ko_zonation.csv")

# Save metadata without row names to preserve column structure
rownames(meta_data_df) <- NULL
write_csv(meta_data_df, file = "data/RNA/meta_data_df_pnpla3ko_zonation.csv")
