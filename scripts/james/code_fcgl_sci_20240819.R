### Setup

# Load packages and color palette
library(dplyr)
library(Seurat)
suppressPackageStartupMessages(library(SCpubr))
library(glmGamPoi)
library(scir)
library(scDblFinder)
library(Nebulosa)
library(EnhancedVolcano)
# library(clusterProfiler)
# library(qtl2)
# library(decoupleR)
# library(dorothea)
# library(tibble)
# library(tidyr)
# library(pheatmap)
# source("/omics/groups/OE0538/internal/users/j083d/sci-fcg/R/theme_color_palette.R")
# 
# setwd("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci")
# 


### Data import and UMI filtering

# Import sample sheet
metadata <- readxl::read_xlsx(
  path = "/omics/groups/OE0538/internal/projects/liver_fcg/tables/Suppl_Table_1_OverviewData.xlsx",
  sheet = "Sample overview"
) %>%
  dplyr::filter(!is.na(sequencing_name)) %>%
  dplyr::mutate(sample_name = sequencing_name) %>%
  dplyr::arrange(sample_name)

# Import raw sci-rocket output
# fcgl <- import_samples_seurat_raw(
#   folder = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_data/20240524_36047_sx42b/36047/alignment/",
#   samples = unique(metadata$sample_name),
#   gtf = "/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf",
#   metadata = metadata %>% dplyr::distinct(sample_name, group, gonadal_type, mouse_label, allosome_type)
# )

# Save rds
# saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_imported.rds")

# Load rds
fcgl <- readRDS("data/RNA/fcgl_pseudobulk_gonad_20240819.rds")

# Remove UMI outliers.
lower_umi_threshold = 500
fcgl <- subset(fcgl, subset = nCount_RNA > lower_umi_threshold)
#AJFD:could use median instead of mean
upper_umi_threshold <- base::mean(fcgl$nCount_RNA) + (3 * stats::sd(fcgl$nCount_RNA))
fcgl <- subset(fcgl, subset = nCount_RNA < upper_umi_threshold)

# Save rds
saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_filtered.rds")



### Pre-processing

fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_filtered.rds")

fcgl <- NormalizeData(fcgl)
fcgl <- FindVariableFeatures(fcgl, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fcgl)
fcgl <- ScaleData(fcgl, features = all.genes)
saveRDS(fcgl,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_scaled.rds")

fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_scaled.rds")

fcgl <- RunPCA(fcgl, features = VariableFeatures(object = fcgl))
ElbowPlot(fcgl, ndims = 50)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_elbowplot.pdf", width = 5, height = 5)

fcgl <- FindNeighbors(fcgl, dims = 1:20)
fcgl <- FindClusters(fcgl, resolution = 0.15)

fcgl <- RunUMAP(fcgl, dims = 1:20, n.neighbors = 30, min.dist = 0.1)

# Save Dim plot
SCpubr::do_DimPlot(fcgl,
                   pt.size = 0.03,
                   legend.position = "right")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_clustered_dimplot.pdf", width = 7, height = 7)

saveRDS(fcgl,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_clustered.rds")



### Doublet removal

# Load rds
#fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_clustered.rds")

# Predict doublets
sce <- scDblFinder(GetAssayData(fcgl, slot="counts"),
                   dbr = .02,
                   clusters=Idents(fcgl))
fcgl$scDblFinder.score <- sce$scDblFinder.score
fcgl$scDblFinder.class <- sce$scDblFinder.class
rm(sce)

# Visualize doublets
knitr::kable(table(fcgl$scDblFinder.class), align = "c")
doublets <- rownames(fcgl[,fcgl@meta.data$scDblFinder.class == "doublet"]@meta.data)
SCpubr::do_DimPlot(fcgl, cells.highlight = doublets, pt.size = .05) +
  ggplot2::labs(x = "UMAP1", y = "UMAP2", title = NULL) +
  scir::theme_ggplot() +
  ggplot2::theme(legend.position = "right")

# Remove doublets
cells_prior <- base::ncol(fcgl)
fcgl <- subset(fcgl, subset = scDblFinder.class == "singlet")
cells_post <- base::ncol(fcgl)
cells_lost <- cells_prior - cells_post

# Print table
knitr::kable(
  base::data.frame(
    cells_prior = cells_prior,
    cells_post = cells_post,
    cells_lost = cells_lost,
    perc_lost = round(cells_lost / cells_prior * 100, 1)
  )
)

# Save rds
saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_doubletless.rds")

# Save Dim plot
SCpubr::do_DimPlot(fcgl,
                   pt.size = 0.03,
                   legend.position = "right")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_doubletless_dimplot.pdf", width = 7, height = 7)

# UMI/gene stats
mean(fcgl$nCount_RNA)
mean(fcgl$nFeature_RNA)



### Re-pre-processing

#fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_doubletless.rds")

fcgl <- NormalizeData(fcgl)
fcgl <- FindVariableFeatures(fcgl, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fcgl)
fcgl <- ScaleData(fcgl, features = all.genes)
fcgl <- RunPCA(fcgl, features = VariableFeatures(object = fcgl))
fcgl <- FindNeighbors(fcgl, dims = 1:20)
fcgl <- FindClusters(fcgl, resolution = 0.15)
fcgl <- RunUMAP(fcgl, dims = 1:20, n.neighbors = 30, min.dist = 0.1)

# Save Dim plot
SCpubr::do_DimPlot(fcgl,
                   pt.size = 0.03,
                   legend.position = "right")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_reclustered_dimplot.pdf", width = 7, height = 7)

saveRDS(fcgl,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_reclustered.rds")



### Cell type annotation

fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_reclustered.rds")

# Dot plot
cell_type_markers <- list("Hepatocytes" = c("Hnf4a", "Cyp2e1", "Cyp2f2"),
                          "Endothelial cells" = c("Btnl9", "Lyve1", "Clec4g"),
                          "Macrophages" = c("Vsig4", "Adgre1", "Clec4f"),
                          "Stellate cells" = c("Lrat", "Des"),
                          "Cholangiocytes" = c("Krt19", "Itgb4"),
                          "T cells" = c("Cd3g","Cd4", "Cd8a"),
                          "B cells" = c("Cd19", "Ms4a1"),
                          "NK cells" = c("Itga2", "Ncam1", "Nkg7"))
SCpubr::do_DotPlot(sample = fcgl, features = cell_type_markers)

# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_reclustered_dotplot.pdf", width = 20, height = 5)

# Assign cell type ID to each cluster

new.cluster.ids <- c("Hepatocytes","Hepatocytes","Hepatocytes","Endothelial cells","Hepatocytes","Hepatocytes","Immune cells","Stellate cells",
                     "Hepatocytes","Immune cells","Immune cells","Immune cells","Immune cells","Hepatocytes","Stellate cells","Cholangiocytes",
                     "Endothelial cells","Cholangiocytes","Endothelial cells")

names(new.cluster.ids) <- levels(fcgl)
# fcgl <- RenameIdents(fcgl, new.cluster.ids)

# Make new Dot plot with cell types instead of cluster numbers, and with cell type labels remove from along the top to avoid redundancy
SCpubr::do_DotPlot(sample = fcgl, features = cell_type_markers)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated_dotplot.pdf", width = 15, height = 3)

Save Dim plot
SCpubr::do_DimPlot(fcgl,
                   pt.size = 0.03,
                   legend.position = "right")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated_dimplot.pdf", width = 8, height = 7)

# Save rds
saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated.rds")



# ### Immune cell subclustering
# 
# # Load rds and pull out non-Heps
# fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated.rds")
# fcgl$cell_type <- Idents(fcgl)
# fcgl <- subset(x = fcgl, subset = cell_type == "Immune cells")
# fcgl <- NormalizeData(fcgl)
# fcgl <- FindVariableFeatures(fcgl, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(fcgl)
# fcgl <- ScaleData(fcgl, features = all.genes)
# fcgl <- RunPCA(fcgl, features = VariableFeatures(object = fcgl))
# fcgl <- FindNeighbors(fcgl, dims = 1:20)
# fcgl <- FindClusters(fcgl, resolution = 0.15)
# fcgl <- RunUMAP(fcgl, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
# 
# # Save Dim plot
# SCpubr::do_DimPlot(fcgl,
#                    pt.size = 0.03,
#                    legend.position = "right")
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_immunecells_dimplot.pdf", width = 7, height = 7)
# 
# saveRDS(fcgl,"/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_immunecells.rds")
# 
# fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_immunecells.rds")
# 
# # Dot plot
# cell_type_markers <- list("Hepatocytes" = c("Hnf4a","Cyp2e1","Cyp2f2"),
#                           "Endothelial cells" = c("Btnl9","Lyve1","Clec4g"),
#                           "Macrophages" = c("Vsig4","Adgre1","Clec4f"),
#                           "Stellate cells" = c("Lrat","Des"),
#                           "Cholangiocytes" = c("Krt19","Itgb4"),
#                           "T cells" = c("Cd3d","Cd3e","Cd3g","Cd4","Cd8a"),
#                           "B cells" = c("Cd79a","Cd79b","Cd19","Ms4a1"),
#                           "NK cells" = c("Itga2","Ncam1","Nkg7"),
#                           "Monocytes" = c("Ccl5","Cd14"))
# SCpubr::do_DotPlot(sample = fcgl, features = cell_type_markers)
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_immunecells_dotplot.pdf", width = 15, height = 6)
# 

# Assign cell type ID to each cluster

# new.cluster.ids <- c("XXX")
#
# names(new.cluster.ids) <- levels(fcgl)
# fcgl <- RenameIdents(fcgl, new.cluster.ids)

# # Make new Dot plot with cell types instead of cluster numbers, and with cell type labels remove from along the top to avoid redundancy
# SCpubr::do_DotPlot(sample = fcgl, features = cell_type_markers)
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated_dotplot.pdf", width = 15, height = 3)

# Save Dim plot
# SCpubr::do_DimPlot(fcgl,
#                    pt.size = 0.03,
#                    legend.position = "right")
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated_dimplot.pdf", width = 8, height = 7)

# # Save rds
# saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated.rds")
#
#
#
# ### Cell type proportions
#
# # Load rds
# fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_reannotatednonhepatocytes.rds")
#
# fcgl$cell_type <- Idents(fcgl)
#
# scc <- SCpubr::do_BarPlot(fcgl,
#                           group.by = "cell_type",
#                           split.by = "allosome_type",
#                           legend.position = "none",
#                           flip = FALSE,
#                           position = "fill")
#
# gonad <- SCpubr::do_BarPlot(fcgl,
#                             group.by = "cell_type",
#                             split.by = "gonadal_type",
#                             legend.position = "none",
#                             flip = FALSE,
#                             position = "fill")
#
# scc + gonad
#
#
#
# ### Differential expression
#
# # Load rds
# fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated.rds")
#
# # Exclude XYT-WT from subsequent analysis
# fcgl <- subset(x = fcgl, group != "XYT-WT")
#
# # Add cell type annotations to metadata
# fcgl$cell_type <- Seurat::Idents(fcgl)
#
# # Pseudobulk
# fcgl_pseudobulk_gonad <- AggregateExpression(fcgl, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "gonadal_type"))
# fcgl_pseudobulk_scc <- AggregateExpression(fcgl, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "allosome_type"))
#
# # Save rds
# saveRDS(fcgl_pseudobulk_gonad, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_gonad.rds")
# saveRDS(fcgl_pseudobulk_scc, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_scc.rds")
#
# # Load rds
# pseudobulk <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_gonad.rds")
# pseudobulk <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_scc.rds")
#
# # Subset on desired first-level cell type
# pseudobulk_hepatocytes <- subset(x = pseudobulk, subset = cell_type == "Hepatocytes")
# pseudobulk_immunecells <- subset(x = pseudobulk, subset = cell_type == "Immune cells")
# pseudobulk_endothelialcells <- subset(x = pseudobulk, subset = cell_type == "Endothelial cells")
# pseudobulk_stellatecells <- subset(x = pseudobulk, subset = cell_type == "Stellate cells")
# pseudobulk_cholangiocytes <- subset(x = pseudobulk, subset = cell_type == "Cholangiocytes")
#
# # Find markers
# Idents(pseudobulk_hepatocytes) <- "gonadal_type"
# de_genes <- FindMarkers(pseudobulk_hepatocytes, ident.1="Ovaries", ident.2="Testes", test.use = "DESeq2")
#
# # # Plot (SCpubr Expression heatmap)
# # de_genes_ovaries <- subset(x = de_genes, subset = p_val_adj < 0.05 & avg_log2FC > 2)
# # de_genes_testes <- subset(x = de_genes, subset = p_val_adj < 0.05 & avg_log2FC < -2)
# # SCpubr::do_ExpressionHeatmap(sample = pseudobulk_hepatocytes,
# #                              features = rownames(de_genes_testes),
# #                              group.by = c("gonadal_type"),
# #                              flip = FALSE,
# #                              viridis.direction = -1)
#
# # Plot (Enhancedvolcano)
# EnhancedVolcano(de_genes,
#                 lab = rownames(de_genes),
#                 x = 'avg_log2FC',
#                 y = 'p_val_adj',
#                 title = NULL,
#                 subtitle = NULL,
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 #xlim = c(-8,8.5),
#                 #ylim = c(0,300),
#                 pointSize = 2,
#                 labSize = 4,
#                 col=c('grey', 'grey', 'grey', 'red3'),
#                 colAlpha = 1,
#                 legendPosition = 'none',
#                 gridlines.major = FALSE,
#                 gridlines.minor = FALSE,
#                 caption = NULL,
#                 titleLabSize = 15,
#                 axisLabSize = 15,
#                 drawConnectors = TRUE,
#                 widthConnectors = 0.3,
#                 max.overlaps = 10,
#                 boxedLabels = TRUE,
#                 #cutoffLineType="blank",
#                 ylab = bquote(-log10padj),
#                 xlab = bquote(log2FC))
#
#

# ## Sci-seq analysis targeted toward Xenium panel
#
# # Pull GDGs out of full Xenium panel list
# xenium_genes <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/FCGL_Xenium_Panel.csv")
# xenium_genes <- subset(xenium_genes, xenium_genes$Reason == "gdg")
# xenium_genes <- xenium_genes[,1]
#
# # Load rds and subset
# pseudobulk <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_allcelltypes_gonad.rds")
# pseudobulk_hepatocytes <- subset(x = pseudobulk, subset = cell_type == "Hepatocytes")
# pseudobulk_hepatocytes_xenium <- subset(x = pseudobulk_hepatocytes, features = xenium_genes)
#
# Idents(pseudobulk_hepatocytes_xenium) <- "gonadal_type"
# de_genes <- FindMarkers(pseudobulk_hepatocytes_xenium, ident.1="Ovaries", ident.2="Testes", test.use = "DESeq2")
#
# # Plot (Enhancedvolcano)
# EnhancedVolcano(de_genes,
#                 lab = rownames(de_genes),
#                 x = 'avg_log2FC',
#                 y = 'p_val_adj',
#                 title = "Xenium panel - Hepatocytes - Ovaries vs Testes",
#                 subtitle = NULL,
#                 pCutoff = 0.05,
#                 FCcutoff = 0,
#                 xlim = c(-8,8),
#                 #ylim = c(0,10),
#                 pointSize = 1.5,
#                 labSize = 4.0,
#                 col=c('grey', 'grey', 'grey', 'red3'),
#                 colAlpha = 1,
#                 legendPosition = 'none',
#                 gridlines.major = FALSE,
#                 gridlines.minor = FALSE,
#                 caption = NULL,
#                 titleLabSize = 15,
#                 axisLabSize = 15,
#                 drawConnectors = TRUE,
#                 widthConnectors = 0.3,
#                 max.overlaps = 10,
#                 #cutoffLineType="blank",
#                 ylab = bquote(-log10padj),
#                 xlab = bquote(log2FC))
#
#
#
# ### TF activity inference (https://saezlab.github.io/decoupleR/articles/tf_sc.html)
#
# # Load rds and subset on hepatocytes
# pseudobulk <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_scipseudobulk_gonad.rds")
# pseudobulk_hepatocytes <- subset(x = pseudobulk, subset = cell_type == "Hepatocytes 1")
#
# # extract the normalised log-transformed counts
# mat <- as.matrix(pseudobulk_hepatocytes@assays$RNA$data)
#
# # get CollecTRI network
# net <- get_collectri(organism='mouse', split_complexes=FALSE)
#
# # Run ulm
# acts <- run_ulm(mat=mat, net=net, .source='source', .target='target', .mor='mor', minsize = 5)
#
# # Extract ulm and store it in tfsulm in Seurat object
# pseudobulk_hepatocytes[['tfsulm']] <- acts %>%
#   pivot_wider(id_cols = 'source', names_from = 'condition',
#               values_from = 'score') %>%
#   column_to_rownames('source') %>%
#   Seurat::CreateAssayObject(.)
#
# # Change assay
# DefaultAssay(object = pseudobulk_hepatocytes) <- "tfsulm"
#
# # Scale the data
# pseudobulk_hepatocytes <- ScaleData(pseudobulk_hepatocytes)
# pseudobulk_hepatocytes@assays$tfsulm@data <- pseudobulk_hepatocytes@assays$tfsulm@scale.data
#
# #set Idents to gonadal_type instead of celltype
# Idents(pseudobulk_hepatocytes) <- "allosome_type"
#
# # Extract activities from object as a long dataframe
# n_tfs <- 25
# df <- t(as.matrix(pseudobulk_hepatocytes@assays$tfsulm@data)) %>%
#   as.data.frame() %>%
#   mutate(cluster = Idents(pseudobulk_hepatocytes)) %>%
#   pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
#   group_by(cluster, source) %>%
#   summarise(mean = mean(score))
#
# # Get top tfs with more variable means across clusters
# tfs <- df %>%
#   group_by(source) %>%
#   summarise(std = sd(mean)) %>%
#   arrange(-abs(std)) %>%
#   head(n_tfs) %>%
#   pull(source)
#
# # Subset long data frame to top tfs and transform to wide matrix
# top_acts_mat <- df %>%
#   filter(source %in% tfs) %>%
#   pivot_wider(id_cols = 'cluster', names_from = 'source',
#               values_from = 'mean') %>%
#   column_to_rownames('cluster') %>%
#   as.matrix()
#
# # Choose color palette
# palette_length = 100
# my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
#
# my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
#                seq(0.05, 3, length.out=floor(palette_length/2)))
#
# # Plot
# tfheatmap_gonad <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
# tfheatmap_scc <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
