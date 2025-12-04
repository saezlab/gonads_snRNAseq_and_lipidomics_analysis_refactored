library(Seurat)
library(readr)
library(ggplot2)

fcgl_pseudobulk_gonad <- readRDS("data/RNA/fcgl_pseudobulk_gonad_20241023.rds")

scaled_RNA <- as.data.frame(fcgl_pseudobulk_gonad@assays$RNA@layers$scale.data)

row.names(scaled_RNA) <- Features(fcgl_pseudobulk_gonad[["RNA"]])

meta_data_df <- as.data.frame(fcgl_pseudobulk_gonad@meta.data)
names(meta_data_df)[1] <- "sample"
meta_data_df$sample <- gsub("-","_",meta_data_df$sample)
meta_data_df$sample_name <- gsub("-","_",meta_data_df$sample_name)

names(scaled_RNA) <- meta_data_df$sample
scaled_RNA$gene <- row.names(scaled_RNA)
scaled_RNA <- scaled_RNA[,c(216,1:215)]

write_csv(scaled_RNA[,c(1:80)], file = "data/RNA/scaled_RNA_1.csv")
write_csv(scaled_RNA[,c(1,81:160)], file = "data/RNA/scaled_RNA_2.csv")
write_csv(scaled_RNA[,c(1,161:216)], file = "data/RNA/scaled_RNA_3.csv")

write_csv(meta_data_df, file = "data/RNA/meta_data_df.csv")

colnames(fcgl_pseudobulk_gonad@meta.data)
#"sample_name"

fcgl_pseudobulk_gonad_umap <- FindVariableFeatures(fcgl_pseudobulk_gonad)

fcgl_pseudobulk_gonad_umap <- RunPCA(fcgl_pseudobulk_gonad_umap)

fcgl_pseudobulk_gonad_umap <- RunUMAP(fcgl_pseudobulk_gonad_umap, dims = 1:20)

# Plot UMAP colored by mouse ID
DimPlot(fcgl_pseudobulk_gonad_umap, group.by = "sample_name", label = FALSE, reduction = "umap") +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(unique(fcgl_pseudobulk_gonad_umap$sample_name)))) +
  theme_minimal()
  