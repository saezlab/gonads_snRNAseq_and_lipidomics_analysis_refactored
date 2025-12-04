### Setup

# Load packages and color palette
library(dplyr)
library(Seurat)
suppressPackageStartupMessages(library(SCpubr))
library(glmGamPoi)
library(scir)
library(scDblFinder)
library(Nebulosa)
library(qtl2)
library(decoupleR)
library(dorothea)
library(tibble)
library(tidyr)
library(EnhancedVolcano)
library(pheatmap)
source("/omics/groups/OE0538/internal/users/j083d/sci-fcg/R/theme_color_palette.R")



### Import raw data and do initial UMI filtering

# Import sample sheet
metadata <- readxl::read_xlsx(
  path = "/omics/groups/OE0538/internal/projects/liver_fcg/tables/Suppl_Table_1_OverviewData.xlsx",
  sheet = "Sample overview"
) %>% 
  dplyr::filter(!is.na(sequencing_name)) %>% 
  dplyr::mutate(sample_name = sequencing_name) %>% 
  dplyr::arrange(sample_name)

# Import raw sci-rocket output
fcgl <- import_samples_seurat_raw(
  folder = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_data/20240524_36047_sx42b/36047/alignment/",
  samples = unique(metadata$sample_name),
  gtf = "/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf",
  metadata = metadata %>% dplyr::distinct(sample_name, group, gonadal_type, mouse_label, allosome_type)
)

# Remove UMI outliers.
upper_umi_threshold <- base::mean(fcgl$nCount_RNA) + (3 * stats::sd(fcgl$nCount_RNA))
lower_umi_threshold = 500
fcgl <- subset(fcgl, subset = nCount_RNA < upper_umi_threshold)
fcgl <- subset(fcgl, subset = nCount_RNA > lower_umi_threshold)

# Save rds
saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_imported.rds")



### Generate QC plots

# Load rds
fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_imported.rds")

# UMI/gene stats
mean(fcgl$nCount_RNA)
mean(fcgl$nFeature_RNA)

# Plot stats
SCpubr::do_ViolinPlot(fcgl, features = "nCount_RNA", group.by = "group")

# # Per sample - No. UMI
# plot_umi <- ggplot2::ggplot(fcgl@meta.data, ggplot2::aes(x = reorder(sample_name, -nCount_RNA, FUN = median), y = nCount_RNA, fill = group)) +
#   ggplot2::geom_boxplot(outlier.size = .1, outlier.shape = 21) +
#   ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
#   ggplot2::labs(x = "Samples", y = "No. of UMI") +
#   ggplot2::guides(fill = "none") +
#   scir::theme_ggplot() +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1))
# 
# # Per sample - No. genes
# plot_genes <- ggplot2::ggplot(fcgl@meta.data, ggplot2::aes(x = reorder(sample_name, -nFeature_RNA, FUN = median), y = nFeature_RNA, fill = group)) +
#   ggplot2::geom_boxplot(outlier.size = .1, outlier.shape = 21) +
#   ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
#   ggplot2::labs(x = "Samples", y = "No. of genes") +
#   ggplot2::guides(fill = "none") +
#   scir::theme_ggplot() +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1))
# 
# # Per sample - No. cells.
# plot_cells <- fcgl@meta.data %>% dplyr::add_count(sample_name) %>% dplyr::distinct(sample_name, group, n) %>% 
#   ggplot2::ggplot(., ggplot2::aes(x = reorder(sample_name, -n), y = n, fill = group, label = n)) +
#   ggplot2::geom_bar(stat = 'identity', width = .8, color = 'black') +
#   ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, 20000)) +
#   ggplot2::scale_fill_manual(values = color_palette, name = NULL) +
#   ggplot2::geom_text(size = 2, angle = 90, nudge_y = -1750, color = 'white') +
#   ggplot2::labs(x = "Samples", y = "No. of cells") +
#   scir::theme_ggplot() +
#   ggplot2::theme(
#     axis.text.x = ggplot2::element_text(angle = 90, size = 6, vjust = .5, hjust = 1),
#     panel.grid.major.x = ggplot2::element_blank()
#   )
# 
# # Per cell - No. UMI vs. genes.
# plot_umi_genes <- ggplot2::ggplot(fcgl@meta.data, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = sample_name)) +
#   ggplot2::geom_point(size = .15, alpha = .1) +
#   ggplot2::geom_smooth(lwd = .8) +
#   ggplot2::guides(color = "none", fill = "none") +
#   ggplot2::labs(x = "No. of UMI", y = "No. of genes") +
#   scir::theme_ggplot()
# 
# # All samples - No. of UMI + threshold.
# plot_umi_cohort <- ggplot2::ggplot(fcgl@meta.data, ggplot2::aes(x = 'Cohort', y = nCount_RNA)) +
#   ggplot2::geom_violin(fill = '#E38C94', alpha = .66) +
#   ggplot2::labs(x = NULL, y = "No. of UMI") +
#   scir::theme_ggplot()
# 
# plot_umi_cohort <- add_threshold(plot_umi_cohort)
# 
# # All samples - No. of genes + threshold.
# plot_genes_cohort <- ggplot2::ggplot(fcgl@meta.data, ggplot2::aes(x = 'Cohort', y = nFeature_RNA)) +
#   ggplot2::geom_violin(fill = '#63ACC0', alpha = .66) +
#   ggplot2::labs(x = NULL, y = "No. of genes") +
#   scir::theme_ggplot()
# 
# # All samples - MT content.
# genes_mt <- base::rownames(fcgl[fcgl[['RNA']][[]]$gene_chr == 'chrM'])
# fcgl[["percent.mt"]] <- Seurat::PercentageFeatureSet(fcgl, features = genes_mt) / 100
# 
# plot_MT_cohort <- ggplot2::ggplot(fcgl@meta.data, ggplot2::aes(x = 'Cohort', y = percent.mt)) +
#   ggplot2::scale_y_continuous(labels = scales::percent_format(), breaks = scales::pretty_breaks(), trans = scales::sqrt_trans()) +
#   ggplot2::geom_violin(fill = '#FEF2DD', alpha = .66) +
#   ggplot2::labs(x = NULL, y = "Rel. MT expression<br><sub>(√-scaled)</sub>") +
#   scir::theme_ggplot()
# 
# design = "
# ACE
# BCF
# DDG
# "
# plot_umi + plot_genes + plot_umi_genes + plot_cells + plot_umi_cohort + plot_genes_cohort + plot_MT_cohort + patchwork::plot_layout(design = design, widths = c(2, 1.25, .5)) + patchwork::plot_annotation(tag_levels = 'a')



### Normalization and scaling

# Load rds
fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_imported.rds")

# Normalise and find variable features
fcgl <- fcgl %>% 
  Seurat::NormalizeData(., normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  Seurat::FindVariableFeatures(., selection.method = "vst", verbose = FALSE, nfeatures = 1000)

# Add known cell-marker genes
known_markers <- readxl::read_xlsx("/omics/groups/OE0538/internal/projects/liver_fcg/tables/Suppl_Table_1_OverviewData.xlsx", sheet = "Known markers")

# Retrieve variable scores
data_variability <- SeuratObject::HVFInfo(fcgl) %>% tibble::rownames_to_column(var = 'gene_name')

# Select features with a min. standardized variance or is known marker
data_variability <- data_variability %>% 
  dplyr::mutate(
    Status = dplyr::case_when(
      variance.standardized >= 2 & gene_name %in% known_markers$gene_name ~ 'Variance + Known marker',
      variance.standardized >= 2 & !(gene_name %in% known_markers$gene_name) ~ 'Variance',
      variance.standardized < 2 & gene_name %in% known_markers$gene_name ~ 'Known marker',
      .default = 'Not selected'
    )
  )
selected_features <- data_variability %>% dplyr::filter(Status != 'Not selected') %>% dplyr::pull(gene_name)

# Regress out total no. of UMI (log10) and MT content
fcgl$nCount_RNA_log <- log10(fcgl$nCount_RNA)
genes_mt <- base::rownames(fcgl[fcgl[['RNA']][[]]$gene_chr == 'chrM'])
fcgl$percent_mt <- Seurat::PercentageFeatureSet(fcgl, features = genes_mt) / 100

fcgl <- Seurat::ScaleData(fcgl, features = selected_features, vars.to.regress = c("percent_mt", "nCount_RNA_log")) %>% 
  Seurat::RunPCA(., features = selected_features)

# Save rds
saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_scaled.rds")



### Sanity check FCG group assignment

# To sanity-check the FCG sub-types, we can compare the (summarized) expression of *Xist* vs. the total expression of all genes on chromosome Y; *Xist* signal should be found only within the XXT and XXO and, conversely, Y-signal only in XYT and XYO.
# Chromosome Y has several known paralogs on chromosome X (e.g., *Ddx3y*, *Uty*, *Eif2s3y*, *Erdr1x* and *Kdm5d*) and some leakage (e.g., multi-mappers) are to be expected but this signal should be **significantly** lower in XX-samples compared to XY samples.
# As Sry is transplanted onto chromosome 3 within the FCG, we disregard Sry for calculating overall chrY expression. In addition, *Erdr1x* / *Erdr1y* is masked.

# Load rds
fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_scaled.rds")

# Define chrY genes
genes_Y <- fcgl[["RNA"]][[]] %>% dplyr::filter(gene_chr == "chrY", !gene_name %in% c("Sry", "Erdr1y")) %>% dplyr::pull(gene_name)

# Summarize sample-wise expression
sample_exprs <- tibble::tibble(
  sample = fcgl$sample_name,
  group = fcgl$group,
  "Σ(chrY)" = rowSums(Seurat::FetchData(fcgl, vars = genes_Y, layer = "counts")),
  "Xist" = Seurat::FetchData(fcgl, vars = "Xist", layer = "counts")[,1],
  "Sry" = Seurat::FetchData(fcgl, vars = "Sry", layer = "counts")[,1],
) %>% 
  dplyr::group_by(sample, group) %>%
  dplyr::summarize_all(.funs = sum) %>% 
  dplyr::ungroup()

# Visualize Xist vs. overall chrY expression
ggplot2::ggplot(sample_exprs, ggplot2::aes(x = `Xist`, y = `Σ(chrY)`, fill = group, label = group)) +
  ggplot2::geom_point(shape = 21, color = "black", size = 3) +
  ggforce::geom_mark_ellipse() +
  ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000), limits = c(0, 1000000))+
  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000), limits = c(0, 1000000))+
  ggplot2::labs(x = "Xist<br><sub>Expression</sub>", y = "Σ(chrY)<sup>-Sry/-Erdr1y</sup><br><sub>Expression</sub>") +
  ggplot2::scale_fill_manual(values = color_palette, guide = "none") +
  scir::theme_ggplot()



### Dimension reduction

# Load rds
fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_scaled.rds")

# Plot PCs
plot_screeplot(fcgl, pc_total = 50, pc_chosen = 30, pc_threshold = NULL)

# UMAP [LSF]
fcgl <- Seurat::RunUMAP(fcgl, umap.method = "uwot", dims = 1:30, n.neighbors = 50, metric = "euclidean", min.dist = 0.1, seed.use = 708813, verbose = FALSE) %>% Seurat::FindNeighbors(., dims = 1:30, verbose = FALSE)

# Cluster with a single resolution [LSF]
fcgl <- Seurat::FindClusters(fcgl, resolution = 0.025, verbose = FALSE)

# Save rds
saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_clustered.rds")

DimPlot(fcgl,
        reduction = "umap",
        label=TRUE,
        raster=TRUE) +
  NoLegend()
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/dimplot_beforedoubletremoval.pdf", width = 5, height = 5)



### Detect and remove putative doublets

# Load rds
fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_clustered.rds")

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

DimPlot(fcgl,
        reduction = "umap",
        label=TRUE,
        raster=TRUE) +
  NoLegend()
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/seuratdimplot_afterdoubletremoval.pdf", width = 5, height = 5)

SCpubr::do_DimPlot(fcgl,
                   label = TRUE, 
                   label.color = "black",
                   repel = TRUE,
                   pt.size = 0.05)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/scpubrdimplot_afterdoubletremoval.pdf", width = 4, height = 5)

# UMI/gene stats
mean(fcgl$nCount_RNA)
mean(fcgl$nFeature_RNA)

# Repeat UMAP and clustering
fcgl <- Seurat::RunUMAP(fcgl, umap.method = "uwot", dims = 1:30, n.neighbors = 50, metric = "euclidean", min.dist = 0.1, seed.use = 708813, verbose = FALSE) %>% Seurat::FindNeighbors(., dims = 1:30, verbose = FALSE)
fcgl <- Seurat::FindClusters(fcgl, resolution = 0.025, verbose = FALSE)

# Save rds
saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_reclustered.rds")

SCpubr::do_DimPlot(fcgl,
                   label = TRUE, 
                   label.color = "black",
                   repel = TRUE,
                   pt.size = 0.05)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/dimplot_reclustered.pdf", width = 4, height = 5)



### Cell type annotation

# Load rds
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

# Assign cell type ID to each cluster
new.cluster.ids <- c("Hepatocytes", "Hepatocytes", "Endothelial cells", "Macrophages", "Stellate cells", "B cells", 
                     "NK and T cells", "Hepatocytes", "Hepatocytes", "Hepatocytes", "Cholangiocytes", "Hepatocytes",
                     "Hepatocytes", "Hepatocytes", "Hepatocytes", "Hepatocytes", "Hepatocytes", "Hepatocytes",
                     "Endothelial cells")
names(new.cluster.ids) <- levels(fcgl)
fcgl <- RenameIdents(fcgl, new.cluster.ids)

# Make new Dot plot with cell types instead of cluster numbers, and with cell type labels remove from along the top to avoid redundancy
SCpubr::do_DotPlot(sample = fcgl, features = cell_type_markers)

# Plot
SCpubr::do_DimPlot(fcgl,
                   pt.size = 0.05,
                   plot.axes = TRUE,
                   legend.position = "right")
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/dimplot_annotated.pdf", width = 7, height = 7)

# Save rds
saveRDS(fcgl, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated.rds")



### Cell type proportions

# Load rds
fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated.rds")

fcgl$cell_type <- Idents(fcgl)

scc <- SCpubr::do_BarPlot(fcgl,
                            group.by = "cell_type",
                            split.by = "allosome_type",
                            legend.position = "none",
                            flip = FALSE,
                            position = "fill")

gonad <- SCpubr::do_BarPlot(fcgl,
                   group.by = "cell_type",
                   split.by = "gonadal_type",
                   legend.position = "none",
                   flip = FALSE,
                   position = "fill")

scc + gonad



### Separate subclustering of Heps and non-Heps to try and get more separated non-Heps and less separated Heps for respective figures

# Load rds
fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated.rds")

fcgl$cell_type <- Idents(fcgl)

# Separate Heps and non-Heps
fcgl_heps <- subset(x = fcgl, cell_type == "Hepatocytes")
fcgl_nonheps <- subset(x = fcgl, cell_type != "Hepatocytes")

# Hep subclustering (with different min.dist)
fcgl_heps <- Seurat::RunUMAP(fcgl_heps, umap.method = "uwot", dims = 1:30, n.neighbors = 50, metric = "euclidean", min.dist = 0.3, seed.use = 708813, verbose = FALSE) %>% Seurat::FindNeighbors(., dims = 1:30, verbose = FALSE)
fcgl_heps <- Seurat::FindClusters(fcgl_heps, resolution = 0.025, verbose = FALSE)
saveRDS(fcgl_heps, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_heps_subclustered.rds")
SCpubr::do_DimPlot(fcgl_heps,
                   label = TRUE, 
                   label.color = "black",
                   repel = TRUE,
                   pt.size = 0.05)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/dimplot_heps_subclustered.pdf", width = 4, height = 5)

# Non-hep subclustering
fcgl_nonheps <- Seurat::RunUMAP(fcgl_nonheps, umap.method = "uwot", dims = 1:30, n.neighbors = 50, metric = "euclidean", min.dist = 0.1, seed.use = 708813, verbose = FALSE) %>% Seurat::FindNeighbors(., dims = 1:30, verbose = FALSE)
fcgl_nonheps <- Seurat::FindClusters(fcgl_nonheps, resolution = 0.025, verbose = FALSE)
saveRDS(fcgl_nonheps, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_nonheps_subclustered.rds")
SCpubr::do_DimPlot(fcgl_nonheps,
                   pt.size = 0.05,
                   plot.axes = TRUE)
ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/dimplot_nonheps_subclustered.pdf", width = 4, height = 5)

Seurat::FeaturePlot(fcgl_nonheps, features = c("Vsig4", "Lrat", "Btnl9", "Cd19", "Cd3g"), raster = FALSE)



### Differential expression

# Load rds
fcgl <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_annotated.rds")

# Exclude XYT-WT from subsequent analysis
fcgl <- subset(x = fcgl, group != "XYT-WT")

# Add cell type annotations to metadata
fcgl$cell_type <- Seurat::Idents(fcgl)

# Pseudobulk
fcgl_sci_pseudobulk_gonad <- AggregateExpression(fcgl, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "gonadal_type"))
fcgl_sci_pseudobulk_scc <- AggregateExpression(fcgl, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "allosome_type"))

# Save rds
saveRDS(fcgl_sci_pseudobulk_gonad, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_gonad.rds")
saveRDS(fcgl_sci_pseudobulk_scc, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_scc.rds")

# Load rds
pseudobulk <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_scc.rds")

# Subset on desired first-level cell type
pseudobulk_hepatocytes <- subset(x = pseudobulk, subset = cell_type == "Hepatocytes")
pseudobulk_macrophages <- subset(x = pseudobulk, subset = cell_type == "Macrophages")
pseudobulk_endothelialcells <- subset(x = pseudobulk, subset = cell_type == "Endothelial cells")
pseudobulk_stellatecells <- subset(x = pseudobulk, subset = cell_type == "Stellate cells")
pseudobulk_cholangiocytes <- subset(x = pseudobulk, subset = cell_type == "Cholangiocytes")
pseudobulk_nktcells <- subset(x = pseudobulk, subset = cell_type == "NK and T cells")
pseudobulk_bcells <- subset(x = pseudobulk, subset = cell_type == "B cells")

# Find markers
Idents(pseudobulk_endothelialcells) <- "allosome_type"
de_genes <- FindMarkers(pseudobulk_endothelialcells, ident.1="XX", ident.2="XY", test.use = "DESeq2")

# # Plot (SCpubr Expression heatmap)
# de_genes_ovaries <- subset(x = de_genes, subset = p_val_adj < 0.05 & avg_log2FC > 2)
# de_genes_testes <- subset(x = de_genes, subset = p_val_adj < 0.05 & avg_log2FC < -2)
# SCpubr::do_ExpressionHeatmap(sample = pseudobulk_hepatocytes,
#                              features = rownames(de_genes_testes),
#                              group.by = c("gonadal_type"),
#                              flip = FALSE,
#                              viridis.direction = -1)

# Plot (Enhancedvolcano)
EnhancedVolcano(de_genes,
                lab = rownames(de_genes),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Endothelial cells - XX vs XY",
                subtitle = NULL,
                pCutoff = 0.05,
                FCcutoff = 2,
                xlim = c(-10,10),
                ylim = c(0,305),
                pointSize = 1.5,
                labSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colAlpha = 1,
                legendPosition = 'none',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                caption = NULL,
                titleLabSize = 15,
                axisLabSize = 15,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                max.overlaps = 5,
                #cutoffLineType="blank",
                ylab = bquote(-log10padj),
                xlab = bquote(log2FC)
                )



## Sci-seq analysis targeted toward Xenium panel

# Load Xenium panel genes
xenium_genes <- readxl::read_xlsx("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/XeniumTargetList_Mouse-Panel-1.xlsx")
xenium_genes <- xenium_genes[,2]
xenium_genes <- pull(xenium_genes)

# Load rds and subset
pseudobulk <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_pseudobulk_gonad.rds")
pseudobulk_hepatocytes <- subset(x = pseudobulk, subset = cell_type == "Hepatocytes")
pseudobulk_hepatocytes_xenium <- subset(x = pseudobulk_hepatocytes, features = xenium_genes)

Idents(pseudobulk_hepatocytes_xenium) <- "gonadal_type"
de_genes <- FindMarkers(pseudobulk_hepatocytes_xenium, ident.1="Ovaries", ident.2="Testes", test.use = "DESeq2")

EnhancedVolcano(de_genes,
                lab = rownames(de_genes),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Xenium panel - Hepatocytes - Ovaries vs Testes",
                # ylim = c(0,50),
                # xlim = c(-2,2),
                subtitle = NULL,
                pCutoff = 0.05,
                FCcutoff = 0,
                pointSize = 1.5,
                labSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colAlpha = 1,
                legendPosition = 'none',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                caption = NULL,
                titleLabSize = 15,
                axisLabSize = 15,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                max.overlaps = 15,
                xlab = bquote(log2FC))



### TF activity inference (https://saezlab.github.io/decoupleR/articles/tf_sc.html)

# Load rds and subset on hepatocytes
pseudobulk <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_sci/fcgl_scipseudobulk_gonad.rds")
pseudobulk_hepatocytes <- subset(x = pseudobulk, subset = cell_type == "Hepatocytes 1")

# extract the normalised log-transformed counts
mat <- as.matrix(pseudobulk_hepatocytes@assays$RNA$data)

# get CollecTRI network
net <- get_collectri(organism='mouse', split_complexes=FALSE)

# Run ulm
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target', .mor='mor', minsize = 5)

# Extract ulm and store it in tfsulm in Seurat object
pseudobulk_hepatocytes[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = pseudobulk_hepatocytes) <- "tfsulm"

# Scale the data
pseudobulk_hepatocytes <- ScaleData(pseudobulk_hepatocytes)
pseudobulk_hepatocytes@assays$tfsulm@data <- pseudobulk_hepatocytes@assays$tfsulm@scale.data

#set Idents to gonadal_type instead of celltype
Idents(pseudobulk_hepatocytes) <- "allosome_type"

# Extract activities from object as a long dataframe
n_tfs <- 25
df <- t(as.matrix(pseudobulk_hepatocytes@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(pseudobulk_hepatocytes)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
tfheatmap_gonad <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
tfheatmap_scc <- pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
