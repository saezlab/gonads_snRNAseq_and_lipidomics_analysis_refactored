### 1.0 Setup ###

# Load packages and color palette
library(Seurat)
library(dplyr)
# suppressPackageStartupMessages(library(SCpubr))
# library(gridExtra)
# library(ggplotify)
# library(scir)
library(glmGamPoi)
# library(scDblFinder)
# library(Nebulosa)
# library(EnhancedVolcano)
library(ggplot2)
# library(tidyr)
# library(patchwork)
# library(ggpubr)
# library(svglite)
# library(future.apply)
# library(purrr)
# library(stringr)

source("/omics/odcf/analysis/OE0538_projects/DO-0012/james_misc/misc_functions.R")
source("/omics/odcf/analysis/OE0538_projects/DO-0012/james_misc/theme_color_palette.R")



### 2.0 FCG pre-processing ###

setwd("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci")

## Import data ##

metadata <- readxl::read_xlsx(
  path = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/sci_metadata.xlsx",
  sheet = "metadata"
) %>%
  dplyr::filter(!is.na(sequencing_name)) %>%
  dplyr::mutate(sample_name = sequencing_name) %>%
  dplyr::arrange(sample_name)

fcg <- import_samples_seurat_raw(
  folder = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_data/20240524_36047_sx42b/36047/alignment/",
  samples = unique(metadata$sample_name),
  gtf = "/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf",
  metadata = metadata %>% dplyr::distinct(sample_name, group, gonadal_type, mouse_label, allosome_type)
)
saveRDS(fcg, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/imported.rds")

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/imported.rds")

## Remove UMI outliers ##

lower_umi_threshold = 500
fcg <- subset(fcg, subset = nCount_RNA > lower_umi_threshold)
upper_umi_threshold <- base::mean(fcg$nCount_RNA) + (3 * stats::sd(fcg$nCount_RNA))
fcg <- subset(fcg, subset = nCount_RNA < upper_umi_threshold)

mean(fcg$nCount_RNA)
mean(fcg$nFeature_RNA)

saveRDS(fcg, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/filtered.rds")

## Clustering ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/filtered.rds")

fcg <- NormalizeData(fcg)
fcg <- FindVariableFeatures(fcg)
fcg <- ScaleData(fcg)
saveRDS(fcg,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/scaled.rds")

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/scaled.rds")

fcg <- RunPCA(fcg)
ElbowPlot(fcg, ndims = 50)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/elbowplot.pdf", width = 10, height = 8)
fcg <- FindNeighbors(fcg, dims = 1:20, reduction = "pca")
fcg <- FindClusters(fcg, resolution = 0.03, cluster.name = "unintegrated_clusters")
fcg <- RunUMAP(fcg, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
saveRDS(fcg,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/clustered.rds")

SCpubr::do_DimPlot(fcg, pt.size = 0.03, legend.position = "right")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/clustered_dimplot.pdf", width = 8, height = 8)
SCpubr::do_DimPlot(fcg, pt.size = 0.03, legend.position = "right", group.by = "group")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/clustered_dimplot_group.pdf", width = 9, height = 8)

## Doublet removal ##

sce <- scDblFinder(GetAssayData(fcg, slot="counts"),
                   dbr = .02,
                   clusters=Idents(fcg))
fcg$scDblFinder.score <- sce$scDblFinder.score
fcg$scDblFinder.class <- sce$scDblFinder.class
rm(sce)

knitr::kable(table(fcg$scDblFinder.class), align = "c")
doublets <- rownames(fcg[,fcg@meta.data$scDblFinder.class == "doublet"]@meta.data)

cells_prior <- base::ncol(fcg)
fcg <- subset(fcg, subset = scDblFinder.class == "singlet")
cells_post <- base::ncol(fcg)
cells_lost <- cells_prior - cells_post

knitr::kable(
  base::data.frame(
    cells_prior = cells_prior,
    cells_post = cells_post,
    cells_lost = cells_lost,
    perc_lost = round(cells_lost / cells_prior * 100, 1)
  )
)

SCpubr::do_DimPlot(fcg, pt.size = 0.03, legend.position = "right", label = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_dimplot.pdf", width = 8, height = 8)

# known cell type markers
cell_type_markers <- list("Hepatocytes" = c("Mlxipl","Apoa2","Hnf4a"),
                          "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
                          "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
                          "Stellate cells" = c("Lrat","Des","Reln"),
                          "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
SCpubr::do_DotPlot(sample = fcg, features = cell_type_markers, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_celltypemarkersdotplot.pdf", width = 10, height = 8)

# unbiased cluster marker detection
fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")
cluster_markers <- FindAllMarkers(fcg, only.pos = TRUE)
cluster_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup() -> top5
saveRDS(cluster_markers, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_clustermarkers.rds")
fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")
cluster_markers <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_clustermarkers.rds")
cluster_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup() -> top5
SCpubr::do_DotPlot(sample = fcg, features = top5$gene, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_clustermarkersdotplot.pdf", width = 3, height = 15)

mean(fcg$nCount_RNA)
mean(fcg$nFeature_RNA)

saveRDS(fcg, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")

## Broad cell type annotation and pseudobulking ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")

new.cluster.ids <- c("Hepatocytes","Hepatocytes","Endothelial cells","Immune cells","Stellate cells","Immune cells","Immune cells","Hepatocytes","Cholangiocytes")
names(new.cluster.ids) <- levels(fcg)
fcg <- RenameIdents(fcg, new.cluster.ids)

# dim plot
SCpubr::do_DimPlot(fcg, pt.size = 0.03, legend.position = "right", label = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated_dimplot.pdf", width = 8, height = 8)

# biased marker dotplot
cell_type_markers <- list("Hepatocytes" = c("Mlxipl","Apoa2","Hnf4a"),
                          "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
                          "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
                          "Stellate cells" = c("Lrat","Des","Reln"),
                          "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
SCpubr::do_DotPlot(sample = fcg, features = cell_type_markers, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated_celltypemarkersdotplot.pdf", width = 10, height = 8)

fcg$cell_type <- Seurat::Idents(fcg)
saveRDS(fcg, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated.rds")

# pseudobulk
fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated.rds")
fcg <- subset(x = fcg, group != "XYT-WT")
fcg$cell_type <- Seurat::Idents(fcg)
fcg_pseudobulks <- AggregateExpression(fcg, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "allosome_type"))
saveRDS(fcg_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/pseudobulks_scc.rds")

## Hepatocyte zonation annotation and pseudobulking ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated.rds")

fcg$cell_type <- Seurat::Idents(fcg)
fcg_heps <- subset(x = fcg, subset = cell_type == "Hepatocytes")
fcg_heps
rm(fcg)
saveRDS(fcg_heps, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps.rds")

fcg_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps.rds")
fcg_heps_pc <- subset(x = fcg_heps, subset = Cyp2e1 > 1 & Cyp2f2 < 1)
fcg_heps_pp <- subset(x = fcg_heps, subset = Cyp2e1 < 1 & Cyp2f2 > 1)
fcg_heps <- merge(fcg_heps_pc, y = fcg_heps_pp)
rm(fcg_heps_pc)
fcg_heps <- JoinLayers(fcg_heps, assay = "RNA")
fcg_heps$zonation_type <- ifelse(colnames(fcg_heps) %in% colnames(fcg_heps_pp),"Periportal","Pericentral")
saveRDS(fcg_heps,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted.rds")

# pseudobulk
fcg_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted.rds")
fcg_heps_zonationsubsetted_pseudobulks <- AggregateExpression(fcg_heps, assays = "RNA", return.seurat = T, group.by = c("sample_name", "zonation_type", "gonadal_type"))
saveRDS(fcg_heps_zonationsubsetted_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted_pseudobulks_gonad.rds")



### 3.0 GDX pre-processing ###

setwd("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/")

## Import data ##

metadata <- readxl::read_xlsx(
  path = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/sci_metadata.xlsx",
  sheet = "metadata"
) %>%
  dplyr::filter(!is.na(sequencing_name)) %>%
  dplyr::mutate(sample_name = sequencing_name) %>%
  dplyr::arrange(sample_name)

gdxhr <- import_samples_seurat_raw(
  folder = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_data/41064_gdxhr/41064/alignment",
  samples = unique(metadata$sample_name),
  gtf = "/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf",
  metadata = metadata %>% dplyr::distinct(sample_name, group, pump_status))
saveRDS(gdxhr, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/imported_20250410.rds")

gdxhr_old <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/imported.rds")
gdxhr_new <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/imported_20250410.rds")

mean(gdxhr_old$nCount_RNA)
mean(gdxhr_old$nFeature_RNA)
mean(gdxhr_new$nCount_RNA)
mean(gdxhr_new$nFeature_RNA)

## Remove UMI outliers ##

lower_umi_threshold = 500
gdxhr <- subset(gdxhr, subset = nCount_RNA > lower_umi_threshold)
upper_umi_threshold <- base::mean(gdxhr$nCount_RNA) + (3 * stats::sd(gdxhr$nCount_RNA))
gdxhr <- subset(gdxhr, subset = nCount_RNA < upper_umi_threshold)

mean(gdxhr$nCount_RNA)
mean(gdxhr$nFeature_RNA)

# Save rds
saveRDS(gdxhr, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/filtered.rds")

## Clustering ##

gdxhr <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/filtered.rds")

gdxhr <- NormalizeData(gdxhr)
gdxhr <- FindVariableFeatures(gdxhr)
gdxhr <- ScaleData(gdxhr)
saveRDS(gdxhr,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/scaled.rds")

gdxhr <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/scaled.rds")

gdxhr <- RunPCA(gdxhr)
ElbowPlot(gdxhr, ndims = 50)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/elbowplot.pdf", width = 10, height = 8)
gdxhr <- FindNeighbors(gdxhr, dims = 1:20, reduction = "pca")
gdxhr <- FindClusters(gdxhr, resolution = 0.03, cluster.name = "unintegrated_clusters")
gdxhr <- RunUMAP(gdxhr, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
saveRDS(gdxhr,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/clustered.rds")

SCpubr::do_DimPlot(gdxhr, pt.size = 0.03, legend.position = "right")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/clustered_dimplot.pdf", width = 8, height = 8)
SCpubr::do_DimPlot(gdxhr, pt.size = 0.03, legend.position = "right", group.by = "group")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/clustered_dimplot_group.pdf", width = 9, height = 8)

## Doublet removal ##

sce <- scDblFinder(GetAssayData(gdxhr, slot="counts"),
                   dbr = .02,
                   clusters=Idents(gdxhr))
gdxhr$scDblFinder.score <- sce$scDblFinder.score
gdxhr$scDblFinder.class <- sce$scDblFinder.class
rm(sce)

knitr::kable(table(gdxhr$scDblFinder.class), align = "c")
doublets <- rownames(gdxhr[,gdxhr@meta.data$scDblFinder.class == "doublet"]@meta.data)

cells_prior <- base::ncol(gdxhr)
gdxhr <- subset(gdxhr, subset = scDblFinder.class == "singlet")
cells_post <- base::ncol(gdxhr)
cells_lost <- cells_prior - cells_post

knitr::kable(
  base::data.frame(
    cells_prior = cells_prior,
    cells_post = cells_post,
    cells_lost = cells_lost,
    perc_lost = round(cells_lost / cells_prior * 100, 1)
  )
)

SCpubr::do_DimPlot(gdxhr, pt.size = 0.03, legend.position = "right", label = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/doubletless_dimplot.pdf", width = 8, height = 8)

# known cell type markers
cell_type_markers <- list("Hepatocytes" = c("Mlxipl","Apoa2","Hnf4a"),
                          "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
                          "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
                          "Stellate cells" = c("Lrat","Des","Reln"),
                          "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
SCpubr::do_DotPlot(sample = gdxhr, features = cell_type_markers, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/doubletless_celltypemarkersdotplot.pdf", width = 10, height = 8)

# unbiased cluster marker detection
cluster_markers <- FindAllMarkers(gdxhr, only.pos = TRUE)
cluster_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup() -> top5
SCpubr::do_DotPlot(sample = gdxhr, features = top5$gene, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/doubletless_clustermarkersdotplot.pdf", width = 10, height = 15)

mean(gdxhr$nCount_RNA)
mean(gdxhr$nFeature_RNA)

saveRDS(gdxhr, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/doubletless.rds")

## Broad cell type annotation and pseudobulking ##

gdxhr <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/doubletless.rds")

new.cluster.ids <- c("Hepatocytes","Endothelial cells","Immune cells","Stellate cells","Hepatocytes", "Cholangiocytes")
names(new.cluster.ids) <- levels(gdxhr)
gdxhr <- RenameIdents(gdxhr, new.cluster.ids)

SCpubr::do_DimPlot(gdxhr, pt.size = 0.03, legend.position = "right", label = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/annotated_dimplot.pdf", width = 8, height = 8)

# biased marker dotplot
cell_type_markers <- list("Hepatocytes" = c("Mlxipl","Apoa2","Hnf4a"),
                          "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
                          "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
                          "Stellate cells" = c("Lrat","Des","Reln"),
                          "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
SCpubr::do_DotPlot(sample = gdxhr, features = cell_type_markers, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/annotated_celltypemarkersdotplot.pdf", width = 10, height = 8)

saveRDS(gdxhr, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/annotated.rds")

# pseudobulk

gdxhr <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/annotated.rds")
gdxhr$cell_type <- Seurat::Idents(gdxhr)
gdxhr_pseudobulks <- AggregateExpression(gdxhr, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "group"))
saveRDS(gdxhr_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/pseudobulks_withoutpumpstatussubsetting.rds")

gdxhr <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/annotated.rds")
gdxhr <- subset(x = gdxhr, pump_status == "yes")
gdxhr$cell_type <- Seurat::Idents(gdxhr)
gdxhr_pseudobulks <- AggregateExpression(gdxhr, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "group"))
saveRDS(gdxhr_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/pseudobulks_withpumpstatussubsetting.rds")

## Hepatocyte zonation annotation and pseudobulking ##

gdxhr <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/gdxhr_annotated.rds")

gdxhr$cell_type <- Seurat::Idents(gdxhr)
gdxhr_heps <- subset(x = gdxhr, subset = cell_type == "Hepatocytes")
gdxhr_heps
rm(gdxhr)
saveRDS(gdxhr_heps, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/heps.rds")

gdxhr_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/heps.rds")
gdxhr_heps_pc <- subset(x = gdxhr_heps, subset = Cyp2e1 > 1 & Cyp2f2 < 1)
gdxhr_heps_pp <- subset(x = gdxhr_heps, subset = Cyp2e1 < 1 & Cyp2f2 > 1)
gdxhr_heps <- merge(gdxhr_heps_pc, y = gdxhr_heps_pp)
rm(gdxhr_heps_pc)
gdxhr_heps <- JoinLayers(gdxhr_heps, assay = "RNA")
gdxhr_heps$zonation_type <- ifelse(colnames(gdxhr_heps) %in% colnames(gdxhr_heps_pp),"Periportal","Pericentral")
saveRDS(gdxhr_heps,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/heps_zonationsubsetted.rds")

# pseudobulk

gdxhr_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/heps_zonationsubsetted.rds")
gdxhr_heps_zonationsubsetted_pseudobulks <- AggregateExpression(gdxhr_heps, assays = "RNA", return.seurat = T, group.by = c("sample_name", "zonation_type", "group"))
saveRDS(gdxhr_heps_zonationsubsetted_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/heps_zonationsubsetted_pseudobulks_withoutpumpstatussubsetting.rds")

gdxhr_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/heps_zonationsubsetted.rds")
gdxhr_heps <- subset(x = gdxhr_heps, pump_status == "yes")
gdxhr_heps_zonationsubsetted_pseudobulks <- AggregateExpression(gdxhr_heps, assays = "RNA", return.seurat = T, group.by = c("sample_name", "zonation_type", "group"))
saveRDS(gdxhr_heps_zonationsubsetted_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/heps_zonationsubsetted_pseudobulks_withpumpstatussubsetting.rds")



### 4.0 P14 pre-processing ###

setwd("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci")

## Import data ##

metadata <- readxl::read_xlsx(
  path = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/metadata.xlsx",
  sheet = "metadata"
) %>%
  dplyr::filter(!is.na(sequencing_name)) %>%
  dplyr::mutate(sample_name = sequencing_name) %>%
  dplyr::arrange(sample_name)

csb303a <- import_samples_seurat_raw(
  folder = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_sci-seq_data/20251006_44766_csb303/44766/alignment",
  samples = metadata %>%
    dplyr::pull(sample_name) %>%
    unique(),
  gtf = "/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf",
  metadata = metadata %>% dplyr::distinct(sample_name, group)
)
mean(csb303a$nCount_RNA)
mean(csb303a$nFeature_RNA)
saveRDS(csb303a, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/csb303a.rds")

csb303b <- import_samples_seurat_raw(
  folder = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_sci-seq_data/20251015_44847_csb303b_and_immortheps/44847/alignment",
  samples = metadata %>%
    dplyr::pull(sample_name) %>%
    unique(),
  gtf = "/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf",
  metadata = metadata %>% dplyr::distinct(sample_name, group)
)
mean(csb303b$nCount_RNA)
mean(csb303b$nFeature_RNA)
saveRDS(csb303b, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/csb303b.rds")

csb303a <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/csb303a.rds")
csb303b <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/csb303b.rds")

p14 <- merge(csb303a, y = c(csb303b))
p14 <- JoinLayers(p14, assay = "RNA")
p14
mean(p14$nCount_RNA)
mean(p14$nFeature_RNA)
saveRDS(p14, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/imported.rds")

## Remove UMI outliers ##

p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/imported.rds")

lower_umi_threshold = 500
p14 <- subset(p14, subset = nCount_RNA > lower_umi_threshold)
upper_umi_threshold <- base::mean(p14$nCount_RNA) + (3 * stats::sd(p14$nCount_RNA))
p14 <- subset(p14, subset = nCount_RNA < upper_umi_threshold)

mean(p14$nCount_RNA)
mean(p14$nFeature_RNA)

saveRDS(p14, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/filtered.rds")

## Clustering ##

p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/filtered.rds")

p14 <- NormalizeData(p14)
p14 <- FindVariableFeatures(p14)
p14 <- ScaleData(p14)
saveRDS(p14,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/scaled.rds")

p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/scaled.rds")

p14 <- RunPCA(p14)
ElbowPlot(p14, ndims = 50)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/elbowplot.pdf", width = 10, height = 8)
p14 <- FindNeighbors(p14, dims = 1:20, reduction = "pca")
p14 <- FindClusters(p14, resolution = 0.06, cluster.name = "unintegrated_clusters")
p14 <- RunUMAP(p14, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
saveRDS(p14,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/clustered.rds")

p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/clustered.rds")

SCpubr::do_DimPlot(p14, pt.size = 0.03, legend.position = "right")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/clustered_dimplot.pdf", width = 8, height = 8)
SCpubr::do_DimPlot(p14, pt.size = 0.03, legend.position = "right", group.by = "group")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/clustered_dimplot_group.pdf", width = 9, height = 8)

## Doublet removal ##

p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/clustered.rds")

sce <- scDblFinder(GetAssayData(p14, slot="counts"),
                   dbr = .02,
                   clusters=Idents(p14))
p14$scDblFinder.score <- sce$scDblFinder.score
p14$scDblFinder.class <- sce$scDblFinder.class
rm(sce)

knitr::kable(table(p14$scDblFinder.class), align = "c")
doublets <- rownames(p14[,p14@meta.data$scDblFinder.class == "doublet"]@meta.data)

cells_prior <- base::ncol(p14)
p14 <- subset(p14, subset = scDblFinder.class == "singlet")
cells_post <- base::ncol(p14)
cells_lost <- cells_prior - cells_post

knitr::kable(
  base::data.frame(
    cells_prior = cells_prior,
    cells_post = cells_post,
    cells_lost = cells_lost,
    perc_lost = round(cells_lost / cells_prior * 100, 1)
  )
)

SCpubr::do_DimPlot(p14, pt.size = 0.03, legend.position = "right", label = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless_dimplot.pdf", width = 8, height = 8)

# known cell type markers
cell_type_markers <- list("Hepatocytes" = c("Mlxipl","Apoa2","Hnf4a"),
                          "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
                          "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
                          "Stellate cells" = c("Lrat","Des","Reln"),
                          "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
SCpubr::do_DotPlot(sample = p14, features = cell_type_markers, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless_celltypemarkersdotplot.pdf", width = 10, height = 8)

# unbiased cluster marker detection
p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless.rds")
cluster_markers <- FindAllMarkers(p14, only.pos = TRUE)
cluster_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup() -> top5
saveRDS(cluster_markers, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless_clustermarkers.rds")
p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless.rds")
cluster_markers <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless_clustermarkers.rds")
cluster_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup() -> top5
SCpubr::do_DotPlot(sample = p14, features = top5$gene, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless_clustermarkersdotplot.pdf", width = 3, height = 15)

mean(p14$nCount_RNA)
mean(p14$nFeature_RNA)

saveRDS(p14, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless.rds")

## Broad cell type annotation and pseudobulking ##

p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/doubletless.rds")

new.cluster.ids <- c("Hepatocytes","Endothelial cells","Stellate cells","Immune cells","Immune cells","Hepatocytes","Hepatocytes","Cholangiocytes","Hepatocytes")
names(new.cluster.ids) <- levels(p14)
p14 <- RenameIdents(p14, new.cluster.ids)

# dim plot
SCpubr::do_DimPlot(p14, pt.size = 0.03, legend.position = "right", label = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/annotated_dimplot.pdf", width = 8, height = 8)

# biased marker dotplot
cell_type_markers <- list("Hepatocytes" = c("Mlxipl","Apoa2","Hnf4a"),
                          "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
                          "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
                          "Stellate cells" = c("Lrat","Des","Reln"),
                          "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
SCpubr::do_DotPlot(sample = p14, features = cell_type_markers, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/annotated_celltypemarkersdotplot.pdf", width = 10, height = 8)

p14$cell_type <- Seurat::Idents(p14)
saveRDS(p14, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/annotated.rds")

# pseudobulk
p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/annotated.rds")
p14$cell_type <- Seurat::Idents(p14)
p14_pseudobulks <- AggregateExpression(p14, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "group"))
saveRDS(p14_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/pseudobulks.rds")

## Hepatocyte zonation annotation and pseudobulking ##

p14 <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/annotated.rds")

p14$cell_type <- Seurat::Idents(p14)
p14_heps <- subset(x = p14, subset = cell_type == "Hepatocytes")
p14_heps
rm(p14)
saveRDS(p14_heps, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/heps.rds")

p14_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/heps.rds")
p14_heps_pc <- subset(x = p14_heps, subset = Cyp2e1 > 1 & Cyp2f2 < 1)
p14_heps_pp <- subset(x = p14_heps, subset = Cyp2e1 < 1 & Cyp2f2 > 1)
p14_heps <- merge(p14_heps_pc, y = p14_heps_pp)
rm(p14_heps_pc)
p14_heps <- JoinLayers(p14_heps, assay = "RNA")
p14_heps$zonation_type <- ifelse(colnames(p14_heps) %in% colnames(p14_heps_pp),"Periportal","Pericentral")
saveRDS(p14_heps,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/heps_zonationsubsetted.rds")

# pseudobulk
p14_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/heps_zonationsubsetted.rds")
p14_heps_zonationsubsetted_pseudobulks <- AggregateExpression(p14_heps, assays = "RNA", return.seurat = T, group.by = c("sample_name", "zonation_type", "group"))
saveRDS(p14_heps_zonationsubsetted_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/p14_sci/heps_zonationsubsetted_pseudobulks.rds")



### 5.0 Human snRNAseq pre-processing ###

setwd("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett")

## Import data ##

library(future.apply)  # for parallel loading

# Set your directory
data_dir <- "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/raw_data"
setwd(data_dir)

# Use multiple cores (adjust workers to your node)
plan(multicore, workers = 8)  # on Windows use: plan(multisession, workers = 8)

# 1. Identify all matrix files → 1 per sample
mtx_files <- list.files(pattern = "_matrix\\.mtx\\.gz$")

# 2. Remove suffix to get sample prefixes (e.g. GSM6964935_hea01)
sample_prefixes <- sub("_matrix\\.mtx\\.gz$", "", mtx_files)

# 3. Load each sample in parallel
seurat_list <- future_lapply(sample_prefixes, function(prefix) {
  
  mtx      <- paste0(prefix, "_matrix.mtx.gz")
  features <- paste0(prefix, "_features.tsv.gz")
  barcodes <- paste0(prefix, "_barcodes.tsv.gz")

  # Read matrix (gene symbols in column 2)
  mat <- ReadMtx(
    mtx      = mtx,
    features = features,
    cells    = barcodes,
    feature.column = 2
  )
  
  # Sample ID like "hea01" (last part after underscore)
  sample_id <- sub(".*_", "", prefix)
  
  obj <- CreateSeuratObject(mat, project = sample_id)
  obj$sample <- sample_id
  obj
})

# Name the list by sample IDs (hea01, hea02, ...)
names(seurat_list) <- sub(".*_", "", sample_prefixes)

combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

combined

human <- combined

mean(human$nCount_RNA)
mean(human$nFeature_RNA)

saveRDS(human, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/imported.rds")

## Remove UMI outliers ##

human <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/imported.rds")

lower_umi_threshold = 500
human <- subset(human, subset = nCount_RNA > lower_umi_threshold)
upper_umi_threshold <- base::mean(human$nCount_RNA) + (3 * stats::sd(human$nCount_RNA))
human <- subset(human, subset = nCount_RNA < upper_umi_threshold)

mean(human$nCount_RNA)
mean(human$nFeature_RNA)

human <- JoinLayers(human, assay = "RNA")

saveRDS(human, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/filtered_500.rds")

## Clustering ##

human <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/filtered_500.rds")

human <- SplitObject(human, split.by = "sample")

human <- lapply(human, function(x) {
  SCTransform(
    x,
    assay = "RNA",
    variable.features.n = 3000,
    verbose = TRUE
  )
})

features <- SelectIntegrationFeatures(human, nfeatures = 3000)

human <- PrepSCTIntegration(human, anchor.features = features)

anchors <- FindIntegrationAnchors(
  human,
  normalization.method = "SCT",
  anchor.features = features
)

human <- IntegrateData(
  anchors,
  normalization.method = "SCT"
)

DefaultAssay(human) <- "integrated"

human <- RunPCA(human, verbose = FALSE)

saveRDS(human, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/clustered_beforefindneighbours.rds")
human <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/clustered_beforefindneighbours.rds")

human <- FindNeighbors(human, dims = 1:30, verbose = FALSE)
human <- FindClusters(human, resolution = 0.1, verbose = FALSE)
human <- RunUMAP(human, dims = 1:30, verbose = FALSE)

saveRDS(human, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/clustered.rds")
human <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/clustered.rds")

DimPlot(human, reduction = "umap", label = TRUE, raster = FALSE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/dimplot_clusters.pdf", width = 8, height = 8)

DimPlot(human, reduction = "umap", group.by = "sample", shuffle = TRUE, raster = FALSE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/dimplot_sample.pdf", width = 8, height = 8)

## Broad cell type annotation and pseudobulking ##

human <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/clustered.rds")

# biased marker dotplot
cell_type_markers <- list("Hepatocytes" = c("ALB","APOA2","CYP2E1"),
                          "Immune cells" = c("PTPRC","MARCO"),
                          "Endothelial cells" = c("BTNL9","PECAM1"),
                          "Stellate cells" = c("RELN"))
SCpubr::do_DotPlot(sample = human, features = cell_type_markers, flip = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/annotated_celltypemarkersdotplot.pdf", width = 10, height = 8)

new.cluster.ids <- c("Hepatocytes","Hepatocytes","Hepatocytes","Other","Other","Other","Other","Other","Other","Other")
names(new.cluster.ids) <- levels(human)
human <- RenameIdents(human, new.cluster.ids)

# dim plot
DimPlot(human, reduction = "umap", shuffle = TRUE, raster = FALSE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/dimplot_annotated.pdf", width = 8, height = 8)

human$cell_type <- Seurat::Idents(human)
saveRDS(human, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/annotated.rds")

## Hepatocyte zonation annotation and pseudobulking ##

human <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/annotated.rds")

human_heps <- subset(x = human, subset = cell_type == "Hepatocytes")
human_heps
rm(human)
saveRDS(human_heps, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/heps.rds")

human_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/heps.rds")

human_heps$zonation_type <- NA
human_heps$zonation_type[FetchData(human_heps, "CYP3A4") > 1 &
                  FetchData(human_heps, "HAL") < 1] <- "Pericentral"
human_heps$zonation_type[FetchData(human_heps, "CYP3A4") < 1 &
                  FetchData(human_heps, "HAL") > 1] <- "Periportal"
human_heps <- subset(human_heps, subset = !is.na(zonation_type))
saveRDS(human_heps,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/heps_zonationsubsetted.rds")

human_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/heps_zonationsubsetted.rds")
VlnPlot(human_heps, group.by = "sample", features = "XIST", assay = "SCT")
sex_by_sample <- c(
  "hea01" = "male",
  "hea02" = "male",
  "hea03" = "male",
  "hea04" = "male",
  "hea05" = "male",
  "hea06" = "male",
  "hea07" = "male",
  "hea08" = "female",
  "hea09" = "female"
)
human_heps$sex <- unname(sex_by_sample[as.character(human_heps$sample)])
human_heps_zonationsubsetted_pseudobulks <- AggregateExpression(
  human_heps,
  assays        = "RNA",
  slot          = "counts",  # <- makes intent clear
  group.by      = c("sample", "zonation_type", "sex"),
  return.seurat = TRUE)
saveRDS(human_heps_zonationsubsetted_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/heps_zonationsubsetted_pseudobulks.rds")

# check expression of key mouse zonation-dependent genes
pseudobulks <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/heps_zonationsubsetted_pseudobulks.rds")
pseudobulks_male = subset(pseudobulks, sex == "male")
VlnPlot(pseudobulks_male, group.by = "zonation_type",   pt.size = 2, features = "LCP1")



### 6.0 Figure panels ###

## Fig 1B ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated.rds")
Idents(fcg) <- "sample_name"
color_palette <- c("sx42_1" = "black", "sx42_14" = "black", "sx42_27" = "black", "sx42_13" = "black")

# Create individual plots
p1 <- SCpubr::do_DimPlot(sample = fcg,
                         idents.keep = c("sx42_1"),
                         legend.position = "none", label = FALSE,
                         pt.size = 0.04,
                         colors.use = color_palette,
                         na.value = "grey90")

p2 <- SCpubr::do_DimPlot(sample = fcg,
                         idents.keep = c("sx42_14"),
                         legend.position = "none", label = FALSE,
                         pt.size = 0.04,
                         colors.use = color_palette,
                         na.value = "grey90")

p3 <- SCpubr::do_DimPlot(sample = fcg,
                         idents.keep = c("sx42_27"),
                         legend.position = "none", label = FALSE,
                         pt.size = 0.04,
                         colors.use = color_palette,
                         na.value = "grey90")

p4 <- SCpubr::do_DimPlot(sample = fcg,
                         idents.keep = c("sx42_13"),
                         legend.position = "none", label = FALSE,
                         pt.size = 0.04,
                         colors.use = color_palette,
                         na.value = "grey90")

# build the grid and keep a handle to it
grid_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# pass that grob explicitly to ggsave
ggsave(
  filename = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated_dimplot_indivres.pdf",
  plot     = grid_plot,      # <-- the crucial bit
  width    = 8,
  height   = 8,
  dpi      = 300             # optional: hi-res export
)

## Fig. 1C ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated.rds")
Idents(fcg) <- "cell_type"
SCpubr::do_DimPlot(fcg, pt.size = 0.04, legend.position = "none", label = FALSE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated_dimplot_celltype.pdf", width = 8, height = 8)

## Fig. 1D ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated.rds")
cluster_markers <- FindAllMarkers(fcg, only.pos = TRUE)
saveRDS(cluster_markers, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated_clustermarkers.rds")
fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated.rds")
cluster_markers <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated_clustermarkers.rds")
cluster_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup() -> top5
feature_list <- split(top5$gene, top5$cluster)
SCpubr::do_DotPlot(
  sample             = fcg,
  features           = feature_list,
  flip               = TRUE,
  dot.scale          = 6,
  legend.length      = 10,
  use_viridis        = FALSE,
  sequential.palette = "Greys"
)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated_clustermarkersdotplot_wide.pdf", width = 2, height = 8)

## Fig. 1F ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated.rds")

fcg$cell_type <- Seurat::Idents(fcg)

p_gonad <- SCpubr::do_BarPlot(
  fcg,
  group.by        = "cell_type",
  split.by        = "gonadal_type",
  legend.position = "none",
  flip            = FALSE,
  position        = "fill"
)

p_allosome <- SCpubr::do_BarPlot(
  fcg,
  group.by        = "cell_type",
  split.by        = "allosome_type",   # <- key change
  legend.position = "none",
  flip            = FALSE,
  position        = "fill"
)

combined <- p_gonad | p_allosome   # combine horizontally (use / for vertical)

ggsave(
  "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/annotated_proportions_gonad_vs_allosome.pdf",
  plot   = combined,
  width  = 8,      # doubled width
  height = 6
)

## Fig. 1 Supp - XX vs XY genotyping ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_xytwtexcluded.rds")

# Define chrY genes.
genes_Y <- fcg[["RNA"]][[]] %>% dplyr::filter(gene_chr == "chrY", !gene_name %in% c("Sry", "Erdr1y")) %>% dplyr::pull(gene_name)
genes_Y <- fcg[["RNA"]][[]] %>% dplyr::filter(gene_chr == "chrY", !gene_name %in% c("Erdr1y")) %>% dplyr::pull(gene_name)
write.csv(genes_Y, file = "genes_Y.csv", row.names = FALSE)

# Summarize sample-wise expression.
sample_exprs <- tibble::tibble(
  sample = fcg$sample_name,
  group = fcg$group,
  "Σ(chrY)" = rowSums(Seurat::FetchData(fcg, vars = genes_Y, layer = "counts")),
  "Xist" = Seurat::FetchData(fcg, vars = "Xist", layer = "counts")[,1],
  "Sry" = Seurat::FetchData(fcg, vars = "Sry", layer = "counts")[,1],
) %>%
  dplyr::group_by(sample, group) %>%
  dplyr::summarize_all(.funs = sum) %>%
  dplyr::ungroup()

sample_exprs <- sample_exprs %>% filter(group != "XYT-WT")

# Define color palette
color_palette <- c(
  "XXO" = "#ED6058",
  "XYO" = "#ED6058",
  "XYT" = "#6D98CB",
  "XXT" = "#6D98CB"
)

x_vals <- sample_exprs$Xist
y_vals <- sample_exprs$`Σ(chrY)`

x_range <- range(x_vals)
y_range <- range(y_vals)

gap_x <- 0.025 * diff(x_range)
gap_y <- 0.025 * diff(y_range)

ggplot(sample_exprs,
       aes(x     = Xist,
           y     = `Σ(chrY)`,
           fill  = group,
           shape = group,    # map shape to group
           label = group)) +
  geom_point(size = 6, color = "black") +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = c(
    "XXO" = 21,   # circle
    "XYO" = 21,   # circle
    "XXT" = 24,   # filled triangle point-up
    "XYT" = 24    # filled triangle point-up
  )) +
  labs(x     = "Xist expression",
       y     = "Y-linked gene expression",
       fill  = NULL,
       shape = NULL) +   # remove the “shape” legend title
  coord_cartesian(
    xlim = c(x_range[1] - gap_x, 110000),
    ylim = c(y_range[1] - gap_y, y_range[2] + gap_y)
  ) +
  theme_classic() +
  scale_x_continuous(breaks = c(0, 50000, 100000),labels = c("0", "50000", "100000")) +
  theme(
    legend.position = c(0.9, 0.8),
    legend.text     = element_text(size = 14),
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 16)
  )

ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_xist_y.pdf", width = 5, height = 5)

## Fig 5B ##

fcg_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted.rds")
zonation_genes <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0012/james_misc/gene_lists/Zonation_Genes.csv")
pericentral_genes <- subset(zonation_genes, zonation_genes$Zonation.status == "Pericentral")
pericentral_genes <- pericentral_genes[,1]
periportal_genes <- subset(zonation_genes, zonation_genes$Zonation.status == "Periportal")
periportal_genes <- periportal_genes[,1]
zonation_genes_grouped <- list("Pericentral" = pericentral_genes,
                               "Periportal" = periportal_genes)

#---------------------------
# 2. Average expression by zonation_type
#    (genes x zonation_type matrix)
#---------------------------
all_zonation_genes <- unique(unlist(zonation_genes_grouped))

avg_exp <- AverageExpression(
  fcg_heps,
  features = all_zonation_genes,
  group.by = "zonation_type",
  slot = "data"                 # log-normalized expression
)[[1]]

# Make sure columns are in the expected order
# and that both Pericentral and Periportal exist
stopifnot(all(c("Pericentral", "Periportal") %in% colnames(avg_exp)))

#---------------------------
# 3. Pick top genes by zonation specificity
#    Using log2 fold-change between Pericentral and Periportal
#---------------------------
eps <- 1e-6  # small number to avoid division by zero

top10_zonation_specific <- map(
  names(zonation_genes_grouped),
  function(zone) {
    genes <- zonation_genes_grouped[[zone]]
    genes <- intersect(genes, rownames(avg_exp))   # keep only genes present

    if (length(genes) == 0) return(character(0))

    # Extract expression in both zones for these genes
    expr_pc <- avg_exp[genes, "Pericentral", drop = TRUE]
    expr_pp <- avg_exp[genes, "Periportal",  drop = TRUE]

    # Define a zone-specific log2 fold-change
    specificity <- if (zone == "Pericentral") {
      log2((expr_pc + eps) / (expr_pp + eps))  # PC-specific
    } else {
      log2((expr_pp + eps) / (expr_pc + eps))  # PP-specific
    }

    # Order by specificity (most specific first)
    specificity <- sort(specificity, decreasing = TRUE)

    n_keep <- min(length(specificity), 10)
    names(specificity)[seq_len(n_keep)]
  }
)

names(top10_zonation_specific) <- names(zonation_genes_grouped)

# Optional: inspect selected genes
print(top10_zonation_specific)

#---------------------------
# 4. Dot plot using top specificity markers
#---------------------------
p <- SCpubr::do_DotPlot(
  sample             = fcg_heps,
  features           = top10_zonation_specific,   # top specific markers
  flip               = TRUE,
  dot.scale          = 6,
  legend.length      = 10,
  use_viridis        = FALSE,
  group.by           = "zonation_type",
  sequential.palette = "Greys",
  dot.min            = 0
)

ggplot2::ggsave(
  filename = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted_zonationmarkersdotplot.pdf",
  width    = 2,
  height   = 9
)

## Fig. 5C ##

fcg_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted.rds")

fcg_heps$zonation_type <- factor(
  fcg_heps$zonation_type,
  levels = c("Pericentral","Periportal")   # <-- your preferred order
)
SCpubr::do_BarPlot(
  fcg_heps,
  group.by = "zonation_type",
  split.by = "gonadal_type",
  legend.position = "right",
  flip = FALSE,
  position = "fill",
  colors = c("Pericentral" = "#3BC900", "Periportal" = "#FF9100")
)
ggplot2::ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted_proportions.pdf", width = 4, height = 5)

means_by_sex <- fcg_heps@meta.data %>%
  group_by(gonadal_type, zonation_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(gonadal_type) %>%
  mutate(prop = n / sum(n))
means_by_sex
means_by_zone <- fcg_heps@meta.data %>%
  group_by(zonation_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop = n / sum(n))
means_by_zone

## Fig. 5D ##

# Load object
fcg_heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted.rds")

# Get counts per sample & zonation, then join gonadal_type
heps_counts <- fcg_heps@meta.data %>%
  filter(zonation_type %in% c("Pericentral", "Periportal")) %>%
  count(zonation_type, sample_name) %>%
  left_join(
    fcg_heps@meta.data[, c("sample_name", "gonadal_type")] %>% distinct(),
    by = "sample_name"
  )

pc_counts <- heps_counts %>%
  filter(zonation_type == "Pericentral")

p_pc <- ggplot(pc_counts, aes(
    x = gonadal_type,
    y = n
  )) +
  geom_boxplot(fill = "#2C9400", color = "#3BC900", outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8, color = "#3BC900") +
  ylab("Pericentral hepatocyte count") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 16)
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Ovaries", "Testes")),
    method      = "wilcox.test",
    label       = "p.format",
    hide.ns     = FALSE
  )

p_pc

pp_counts <- heps_counts %>%
  filter(zonation_type == "Periportal")

p_pp <- ggplot(pp_counts, aes(
    x = gonadal_type,
    y = n
  )) +
  geom_boxplot(fill = "orange", color = "#FF9100", outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8, color = "#FF9100") +
  ylab("Periportal hepatocyte count") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 16)
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Ovaries", "Testes")),
    method      = "wilcox.test",
    label       = "p.format",
    hide.ns     = FALSE
  )

p_pp

p_pc+p_pp
ggplot2::ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/heps_zonationsubsetted_abscounts.pdf", width = 5, height = 5)

## Fig 1. Supp - QC metrics ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")

# Per sample - No. UMI
plot_umi <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = reorder(sample_name, -nCount_RNA, FUN = median), y = nCount_RNA, fill = group)) +
  ggplot2::geom_boxplot(outlier.size = .1, outlier.shape = 21) +
  ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
  ggplot2::labs(x = "Samples", y = "No. of UMI") +
  ggplot2::guides(fill = "none") +
  scir::theme_ggplot() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1))

# Per sample - No. genes
plot_genes <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = reorder(sample_name, -nFeature_RNA, FUN = median), y = nFeature_RNA, fill = group)) +
  ggplot2::geom_boxplot(outlier.size = .1, outlier.shape = 21) +
  ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
  ggplot2::labs(x = "Samples", y = "No. of genes") +
  ggplot2::guides(fill = "none") +
  scir::theme_ggplot() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1))

# Per sample - No. cells.
plot_cells <- fcg@meta.data %>% dplyr::add_count(sample_name) %>% dplyr::distinct(sample_name, group, n) %>%
  ggplot2::ggplot(., ggplot2::aes(x = reorder(sample_name, -n), y = n, fill = group, label = n)) +
  ggplot2::geom_bar(stat = 'identity', width = .8, color = 'black') +
  ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, 20000)) +
  ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
  ggplot2::geom_text(size = 2, angle = 90, nudge_y = -1750, color = 'white') +
  ggplot2::labs(x = "Samples", y = "No. of cells") +
  scir::theme_ggplot() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank()
  )

# Per cell - No. UMI vs. genes.
plot_umi_genes <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = sample_name)) +
  ggplot2::geom_point(size = .15, alpha = .1) +
  ggplot2::geom_smooth(lwd = .8) +
  ggplot2::guides(color = "none", fill = "none") +
  ggplot2::labs(x = "No. of UMI", y = "No. of genes") +
  scir::theme_ggplot()

# All samples - No. of UMI + threshold.
plot_umi_cohort <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = 'Cohort', y = nCount_RNA)) +
  ggplot2::geom_violin(fill = '#E38C94', alpha = .66) +
  ggplot2::labs(x = NULL, y = "No. of UMI") +
  scir::theme_ggplot()

plot_umi_cohort <- add_threshold(plot_umi_cohort)

# All samples - No. of genes + threshold.
plot_genes_cohort <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = 'Cohort', y = nFeature_RNA)) +
  ggplot2::geom_violin(fill = '#63ACC0', alpha = .66) +
  ggplot2::labs(x = NULL, y = "No. of genes") +
  scir::theme_ggplot()

# All samples - MT content.
genes_mt <- base::rownames(fcg[fcg[['RNA']][[]]$gene_chr == 'chrM'])
fcg[["percent.mt"]] <- Seurat::PercentageFeatureSet(fcg, features = genes_mt) / 100

plot_MT_cohort <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = 'Cohort', y = percent.mt)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(), breaks = scales::pretty_breaks(), trans = scales::sqrt_trans()) +
  ggplot2::geom_violin(fill = '#FEF2DD', alpha = .66) +
  ggplot2::labs(x = NULL, y = "Rel. MT expression<br><sub>(√-scaled)</sub>") +
  scir::theme_ggplot()

design = "
ACE
BCF
DDG
"
plot_umi + plot_genes + plot_umi_genes + plot_cells + plot_umi_cohort + plot_genes_cohort + plot_MT_cohort + patchwork::plot_layout(design = design, widths = c(2, 1.25, .5)) + patchwork::plot_annotation(tag_levels = 'a')
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_qc1.pdf", width = 12, height = 8)

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")
fcg <- subset(x = fcg, group != "XYT-WT")

saveRDS(fcg, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_xytwtexcluded.rds")

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")
Idents(fcg) <- "sample_name"

ident_levels <- levels(Idents(fcg))
p <- VlnPlot(fcg, features = "nCount_RNA", pt.size = 0) + NoLegend() + ggtitle(NULL)
means <- fcg@meta.data %>%
  group_by(ident = factor(Idents(fcg), levels = ident_levels)) %>%
  summarize(mean_umis = mean(nCount_RNA), .groups = "drop") %>%
  mutate(x = as.numeric(ident))  # assign numeric positions
p +
  geom_segment(data = means,
               aes(x = x - 0.4, xend = x + 0.4,
                   y = mean_umis, yend = mean_umis),
               color = "black", size = 1) +
  geom_text(data = means,
            aes(x = x, y = mean_umis, label = round(mean_umis, 1)),
            vjust = -0.5, size = 4, color = "black")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/umis_per_nuc.pdf", width = 4, height = 5)

ident_levels <- levels(Idents(fcg))
p <- VlnPlot(fcg, features = "nFeature_RNA", pt.size = 0) + NoLegend() + ggtitle(NULL)
means <- fcg@meta.data %>%
  group_by(ident = factor(Idents(fcg), levels = ident_levels)) %>%
  summarize(mean_genes = mean(nFeature_RNA), .groups = "drop") %>%
  mutate(x = as.numeric(ident))  # assign numeric positions
p +
  geom_segment(data = means,
               aes(x = x - 0.4, xend = x + 0.4,
                   y = mean_genes, yend = mean_genes),
               color = "black", size = 1) +
  geom_text(data = means,
            aes(x = x, y = mean_genes, label = round(mean_genes, 1)),
            vjust = -0.5, size = 4, color = "black")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/genes_per_nuc.pdf", width = 4, height = 5)

ident_levels <- levels(Idents(fcg))
genes_mt <- base::rownames(fcg[fcg[['RNA']][[]]$gene_chr == 'chrM'])
fcg[["percent.mt"]] <- Seurat::PercentageFeatureSet(fcg, features = genes_mt) / 100
p <- VlnPlot(fcg, features = "percent.mt", pt.size = 0) + NoLegend() + ggtitle(NULL)
means <- fcg@meta.data %>%
  group_by(ident = factor(Idents(fcg), levels = ident_levels)) %>%
  summarize(mean_mt = mean(percent.mt), .groups = "drop") %>%
  mutate(x = as.numeric(ident))  # assign numeric positions
p +
  geom_segment(data = means,
               aes(x = x - 0.4, xend = x + 0.4,
                   y = mean_mt, yend = mean_mt),
               color = "black", size = 1) +
  geom_text(data = means,
            aes(x = x, y = mean_mt, label = round(mean_mt, 1)),
            vjust = -0.5, size = 4, color = "black")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/mt_per_nuc.pdf", width = 4, height = 5)## Dataset QC ##

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")

# Per sample - No. UMI
plot_umi <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = reorder(sample_name, -nCount_RNA, FUN = median), y = nCount_RNA, fill = group)) +
  ggplot2::geom_boxplot(outlier.size = .1, outlier.shape = 21) +
  ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
  ggplot2::labs(x = "Samples", y = "No. of UMI") +
  ggplot2::guides(fill = "none") +
  scir::theme_ggplot() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1))

# Per sample - No. genes
plot_genes <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = reorder(sample_name, -nFeature_RNA, FUN = median), y = nFeature_RNA, fill = group)) +
  ggplot2::geom_boxplot(outlier.size = .1, outlier.shape = 21) +
  ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
  ggplot2::labs(x = "Samples", y = "No. of genes") +
  ggplot2::guides(fill = "none") +
  scir::theme_ggplot() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1))

# Per sample - No. cells.
plot_cells <- fcg@meta.data %>% dplyr::add_count(sample_name) %>% dplyr::distinct(sample_name, group, n) %>%
  ggplot2::ggplot(., ggplot2::aes(x = reorder(sample_name, -n), y = n, fill = group, label = n)) +
  ggplot2::geom_bar(stat = 'identity', width = .8, color = 'black') +
  ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, 20000)) +
  ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
  ggplot2::geom_text(size = 2, angle = 90, nudge_y = -1750, color = 'white') +
  ggplot2::labs(x = "Samples", y = "No. of cells") +
  scir::theme_ggplot() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank()
  )

# Per cell - No. UMI vs. genes.
plot_umi_genes <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = sample_name)) +
  ggplot2::geom_point(size = .15, alpha = .1) +
  ggplot2::geom_smooth(lwd = .8) +
  ggplot2::guides(color = "none", fill = "none") +
  ggplot2::labs(x = "No. of UMI", y = "No. of genes") +
  scir::theme_ggplot()

# All samples - No. of UMI + threshold.
plot_umi_cohort <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = 'Cohort', y = nCount_RNA)) +
  ggplot2::geom_violin(fill = '#E38C94', alpha = .66) +
  ggplot2::labs(x = NULL, y = "No. of UMI") +
  scir::theme_ggplot()

plot_umi_cohort <- add_threshold(plot_umi_cohort)

# All samples - No. of genes + threshold.
plot_genes_cohort <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = 'Cohort', y = nFeature_RNA)) +
  ggplot2::geom_violin(fill = '#63ACC0', alpha = .66) +
  ggplot2::labs(x = NULL, y = "No. of genes") +
  scir::theme_ggplot()

# All samples - MT content.
genes_mt <- base::rownames(fcg[fcg[['RNA']][[]]$gene_chr == 'chrM'])
fcg[["percent.mt"]] <- Seurat::PercentageFeatureSet(fcg, features = genes_mt) / 100

plot_MT_cohort <- ggplot2::ggplot(fcg@meta.data, ggplot2::aes(x = 'Cohort', y = percent.mt)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(), breaks = scales::pretty_breaks(), trans = scales::sqrt_trans()) +
  ggplot2::geom_violin(fill = '#FEF2DD', alpha = .66) +
  ggplot2::labs(x = NULL, y = "Rel. MT expression<br><sub>(√-scaled)</sub>") +
  scir::theme_ggplot()

design = "
ACE
BCF
DDG
"
plot_umi + plot_genes + plot_umi_genes + plot_cells + plot_umi_cohort + plot_genes_cohort + plot_MT_cohort + patchwork::plot_layout(design = design, widths = c(2, 1.25, .5)) + patchwork::plot_annotation(tag_levels = 'a')
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_qc1.pdf", width = 12, height = 8)

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")
fcg <- subset(x = fcg, group != "XYT-WT")

saveRDS(fcg, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless_xytwtexcluded.rds")

fcg <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/doubletless.rds")
Idents(fcg) <- "sample_name"

ident_levels <- levels(Idents(fcg))
p <- VlnPlot(fcg, features = "nCount_RNA", pt.size = 0) + NoLegend() + ggtitle(NULL)
means <- fcg@meta.data %>%
  group_by(ident = factor(Idents(fcg), levels = ident_levels)) %>%
  summarize(mean_umis = mean(nCount_RNA), .groups = "drop") %>%
  mutate(x = as.numeric(ident))  # assign numeric positions
p +
  geom_segment(data = means,
               aes(x = x - 0.4, xend = x + 0.4,
                   y = mean_umis, yend = mean_umis),
               color = "black", size = 1) +
  geom_text(data = means,
            aes(x = x, y = mean_umis, label = round(mean_umis, 1)),
            vjust = -0.5, size = 4, color = "black")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/umis_per_nuc.pdf", width = 4, height = 5)

ident_levels <- levels(Idents(fcg))
p <- VlnPlot(fcg, features = "nFeature_RNA", pt.size = 0) + NoLegend() + ggtitle(NULL)
means <- fcg@meta.data %>%
  group_by(ident = factor(Idents(fcg), levels = ident_levels)) %>%
  summarize(mean_genes = mean(nFeature_RNA), .groups = "drop") %>%
  mutate(x = as.numeric(ident))  # assign numeric positions
p +
  geom_segment(data = means,
               aes(x = x - 0.4, xend = x + 0.4,
                   y = mean_genes, yend = mean_genes),
               color = "black", size = 1) +
  geom_text(data = means,
            aes(x = x, y = mean_genes, label = round(mean_genes, 1)),
            vjust = -0.5, size = 4, color = "black")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/genes_per_nuc.pdf", width = 4, height = 5)

ident_levels <- levels(Idents(fcg))
genes_mt <- base::rownames(fcg[fcg[['RNA']][[]]$gene_chr == 'chrM'])
fcg[["percent.mt"]] <- Seurat::PercentageFeatureSet(fcg, features = genes_mt) / 100
p <- VlnPlot(fcg, features = "percent.mt", pt.size = 0) + NoLegend() + ggtitle(NULL)
means <- fcg@meta.data %>%
  group_by(ident = factor(Idents(fcg), levels = ident_levels)) %>%
  summarize(mean_mt = mean(percent.mt), .groups = "drop") %>%
  mutate(x = as.numeric(ident))  # assign numeric positions
p +
  geom_segment(data = means,
               aes(x = x - 0.4, xend = x + 0.4,
                   y = mean_mt, yend = mean_mt),
               color = "black", size = 1) +
  geom_text(data = means,
            aes(x = x, y = mean_mt, label = round(mean_mt, 1)),
            vjust = -0.5, size = 4, color = "black")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/mt_per_nuc.pdf", width = 4, height = 5)

## Fig 5 Supp - Human zonation marker dot plot ##

heps <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/heps_zonationsubsetted.rds")
zonation_genes <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0012/james_misc/gene_lists/Human_Zonation_Genes.csv")
pericentral_genes <- subset(zonation_genes, zonation_genes$Zonation.status == "Pericentral")
pericentral_genes <- pericentral_genes[,1]
periportal_genes <- subset(zonation_genes, zonation_genes$Zonation.status == "Periportal")
periportal_genes <- periportal_genes[,1]
zonation_genes_grouped <- list("Pericentral" = pericentral_genes,
                               "Periportal" = periportal_genes)

DefaultAssay(heps) <- "SCT"

p <- SCpubr::do_DotPlot(
  sample             = heps,
  features           = zonation_genes_grouped,   # your grouped Pericentral / Periportal genes
  flip               = TRUE,
  dot.scale          = 6,
  legend.length      = 10,
  use_viridis        = FALSE,
  group.by           = "zonation_type",          # assumes this column exists in meta.data
  sequential.palette = "Greys"
)
p
ggplot2::ggsave(
  filename = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/human_matchett/heps_zonationsubsetted_zonationmarkersdotplot.pdf",
  plot     = p,
  width    = 2.5,
  height   = 5
)

