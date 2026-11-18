library(Seurat)
library(readr)

# Load compatibility functions for Seurat v4/v5
source("scripts/seurat_compatibility.R")

heps_zonationsubclustered <- readRDS("data/RNA/fcgl_heps_zonationsubclustered_pseudobulk_gonad_20241023.rds")

# Extract scaled data (compatible with both Seurat v4 and v5)
scaled_RNA <- as.data.frame(get_layer_data(heps_zonationsubclustered, assay = "RNA", layer = "scale.data"))

row.names(scaled_RNA) <- Features(heps_zonationsubclustered[["RNA"]])

meta_data_df <- as.data.frame(heps_zonationsubclustered@meta.data)
names(meta_data_df)[1] <- "sample"
meta_data_df$sample <- gsub("-","_",meta_data_df$sample)
meta_data_df$sample_name <- gsub("-","_",meta_data_df$sample_name)

names(scaled_RNA) <- meta_data_df$sample
scaled_RNA$gene <- row.names(scaled_RNA)
scaled_RNA <- scaled_RNA[,c(87,1:86)]

write_csv(scaled_RNA, file = "data/RNA/scaled_RNA_zonation.csv")

# Save metadata without row names to preserve column structure
rownames(meta_data_df) <- NULL
write_csv(meta_data_df, file = "data/RNA/meta_data_df_zonation.csv")
