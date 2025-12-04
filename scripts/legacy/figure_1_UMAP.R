library(Seurat)
library(readr)
library(ggplot2)

seurat_obj <- readRDS("../original_datasets/snRNAseq_seurat_objects/fcgl_annotated.rds")
# Make sure you've done PCA and UMAP first
# seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# View available metadata columns
colnames(seurat_obj@meta.data)

# Suppose the column storing mouse IDs is called "mouse_id"
DimPlot(seurat_obj, group.by = "sample_name", label = FALSE) +
  ggtitle("UMAP colored by mouse ID (466412 cells)") + 
  labs(x = "UMAP Axis 1 (don't overthink it)", y = "UMAP Axis 2 (still arbitrary)")
