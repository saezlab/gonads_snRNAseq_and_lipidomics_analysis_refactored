################################################################################
# Workflow 0: snRNA-seq Pre-processing Pipeline
################################################################################
#
# This workflow documents the pre-processing steps for raw sci-RNA-seq data
# that produce the pseudobulk Seurat objects used in downstream analyses.
#
# Author: James Omichinski / Refactored codebase
# Date: 2025
#
# Replaces legacy scripts:
#   - james_first_submission_code.R (sections 2.0 and 3.0)
#
# Datasets processed:
#   1. FCG (Four Core Genotypes) - Main study cohort
#   2. GDX (Gonadectomy Hormone Replacement) - Validation cohort
#
# Key steps for each dataset:
#   1. Import raw sci-RNA-seq data from alignment output
#   2. Quality control: UMI filtering (remove outliers)
#   3. Standard Seurat workflow (normalize, scale, PCA, cluster, UMAP)
#   4. Doublet detection and removal (scDblFinder)
#   5. Cell type annotation using known markers
#   6. Pseudobulk aggregation by sample and cell type
#   7. Hepatocyte zonation annotation (periportal vs pericentral)
#
# NOTE: This workflow requires raw alignment data and is designed for
#       execution on the original HPC environment. For users starting
#       with the pre-processed data, begin with workflow 01.
#
# Output files:
#   - FCG: fcgl_pseudobulk_gonad_*.rds
#   - GDX: gdxhr_pseudobulks_*.rds
#   - Zonation: heps_zonationsubsetted_pseudobulks_*.rds
#
################################################################################

# Clear workspace
rm(list = ls())

################################################################################
# 0. CONFIGURATION
################################################################################

cat(strrep("=", 80), "\n")
cat("WORKFLOW 0: snRNA-seq PRE-PROCESSING PIPELINE\n")
cat(strrep("=", 80), "\n\n")

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
library(glmGamPoi)  # For faster SCTransform/scaling

# Optional: Load scDblFinder for doublet detection
# This requires BiocManager::install("scDblFinder")
has_scDblFinder <- requireNamespace("scDblFinder", quietly = TRUE)
if (has_scDblFinder) {
  library(scDblFinder)
  cat("scDblFinder loaded for doublet detection\n")
} else {
  cat("WARNING: scDblFinder not installed. Doublet detection will be skipped.\n")
  cat("Install with: BiocManager::install('scDblFinder')\n")
}

# Optional: Load SCpubr for nice visualizations
has_SCpubr <- requireNamespace("SCpubr", quietly = TRUE)
if (has_SCpubr) {
  library(SCpubr)
}

################################################################################
# CONFIGURATION: Adjust these paths for your environment
################################################################################

# Base directories (HPC paths - adjust for your system)
# These are the original paths from the DKFZ HPC environment
config <- list(
  # FCG dataset paths
  fcg = list(
    base_dir = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci",
    alignment_dir = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_data/20240524_36047_sx42b/36047/alignment/",
    metadata_file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/sci_metadata.xlsx",
    gtf_file = "/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf"
  ),

  # GDX dataset paths
  gdx = list(
    base_dir = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci",
    alignment_dir = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_data/41064_gdxhr/41064/alignment",
    metadata_file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/sci_metadata.xlsx"
  )
)

# Cell type markers for annotation
cell_type_markers <- list(
  "Hepatocytes" = c("Mlxipl", "Apoa2", "Hnf4a"),
  "Immune cells" = c("Ptprc", "Vsig4", "Adgre1"),
  "Endothelial cells" = c("Btnl9", "Clec4g", "Pecam1"),
  "Stellate cells" = c("Lrat", "Des", "Reln"),
  "Cholangiocytes" = c("Krt19", "Ankrd1", "Epcam")
)

# Zonation markers
zonation_markers <- list(
  "Pericentral" = "Cyp2e1",
  "Periportal" = "Cyp2f2"
)

# QC parameters
qc_params <- list(
  lower_umi_threshold = 500,      # Minimum UMI count
  upper_umi_sd_multiplier = 3,    # Upper threshold = mean + (n * sd)
  doublet_rate = 0.02,            # Expected doublet rate for scDblFinder
  pca_dims = 20,                  # Number of PCA dimensions for clustering
  cluster_resolution = 0.03       # Clustering resolution
)

cat("\nConfiguration loaded.\n")
cat("To run this workflow, ensure paths are correct for your environment.\n\n")


################################################################################
# HELPER FUNCTIONS
################################################################################

#' Import raw sci-RNA-seq samples into Seurat
#'
#' This function is a placeholder for the custom import function used in the
#' original analysis. Users should implement their own import logic based on
#' their alignment pipeline output format.
#'
#' @param folder Path to alignment output folder
#' @param samples Vector of sample names to import
#' @param gtf Path to GTF annotation file
#' @param metadata Metadata dataframe with sample information
#' @return Seurat object with merged samples
import_samples_seurat_raw <- function(folder, samples, gtf, metadata) {
  stop("This is a placeholder function. Implement based on your alignment pipeline output format.

       The original function:
       1. Read 10X-style matrix files from each sample folder
       2. Created individual Seurat objects
       3. Added metadata (group, gonadal_type, etc.)
       4. Merged all samples into a single Seurat object

       Alternative: Use Seurat's Read10X() or ReadMtx() functions for standard formats.")
}


#' Remove UMI outliers
#'
#' @param seurat_obj Seurat object
#' @param lower_threshold Minimum UMI count (default: 500)
#' @param upper_sd_mult Upper threshold multiplier (mean + n*sd)
#' @return Filtered Seurat object
filter_umi_outliers <- function(seurat_obj,
                                 lower_threshold = 500,
                                 upper_sd_mult = 3) {

  cat("Filtering UMI outliers...\n")
  cells_before <- ncol(seurat_obj)

  # Lower threshold
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA > lower_threshold)

  # Upper threshold (mean + n*sd)
  upper_threshold <- mean(seurat_obj$nCount_RNA) +
                     (upper_sd_mult * sd(seurat_obj$nCount_RNA))
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA < upper_threshold)

  cells_after <- ncol(seurat_obj)

  cat("  Cells before:", cells_before, "\n")
  cat("  Cells after:", cells_after, "\n")
  cat("  Cells removed:", cells_before - cells_after,
      "(", round(100 * (cells_before - cells_after) / cells_before, 1), "%)\n")
  cat("  UMI range:", lower_threshold, "-", round(upper_threshold), "\n")
  cat("  Mean UMI:", round(mean(seurat_obj$nCount_RNA)), "\n")
  cat("  Mean features:", round(mean(seurat_obj$nFeature_RNA)), "\n\n")

  return(seurat_obj)
}


#' Standard Seurat clustering workflow
#'
#' @param seurat_obj Seurat object (normalized and scaled)
#' @param pca_dims Number of PCA dimensions to use
#' @param resolution Clustering resolution
#' @return Seurat object with clusters and UMAP
run_clustering_workflow <- function(seurat_obj,
                                     pca_dims = 20,
                                     resolution = 0.03) {

  cat("Running clustering workflow...\n")

  # Normalize
  cat("  Normalizing data...\n")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

  # Find variable features
  cat("  Finding variable features...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)

  # Scale
  cat("  Scaling data...\n")
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

  # PCA
  cat("  Running PCA...\n")
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

  # Find neighbors and clusters
  cat("  Finding neighbors and clusters (resolution:", resolution, ")...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pca_dims,
                               reduction = "pca", verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution,
                              cluster.name = "unintegrated_clusters",
                              verbose = FALSE)

  # UMAP
  cat("  Running UMAP...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:pca_dims,
                         reduction = "pca",
                         reduction.name = "umap.unintegrated",
                         verbose = FALSE)

  cat("  Clustering complete.\n")
  cat("  Number of clusters:", length(unique(Idents(seurat_obj))), "\n\n")

  return(seurat_obj)
}


#' Detect and remove doublets using scDblFinder
#'
#' @param seurat_obj Seurat object
#' @param doublet_rate Expected doublet rate (default: 0.02)
#' @return Seurat object with doublets removed
remove_doublets <- function(seurat_obj, doublet_rate = 0.02) {

  if (!has_scDblFinder) {
    cat("scDblFinder not available. Skipping doublet detection.\n\n")
    return(seurat_obj)
  }

  cat("Detecting doublets with scDblFinder...\n")

  cells_before <- ncol(seurat_obj)

  # Run scDblFinder
  sce <- scDblFinder::scDblFinder(
    GetAssayData(seurat_obj, slot = "counts"),
    dbr = doublet_rate,
    clusters = Idents(seurat_obj)
  )

  # Add scores to Seurat object
  seurat_obj$scDblFinder.score <- sce$scDblFinder.score
  seurat_obj$scDblFinder.class <- sce$scDblFinder.class
  rm(sce)

  # Report doublet counts
  doublet_table <- table(seurat_obj$scDblFinder.class)
  cat("  Singlets:", doublet_table["singlet"], "\n")
  cat("  Doublets:", doublet_table["doublet"], "\n")

  # Remove doublets
  seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

  cells_after <- ncol(seurat_obj)
  cat("  Cells removed:", cells_before - cells_after,
      "(", round(100 * (cells_before - cells_after) / cells_before, 1), "%)\n\n")

  return(seurat_obj)
}


#' Annotate cell types based on cluster identities
#'
#' @param seurat_obj Seurat object with clusters
#' @param cluster_to_celltype Named vector mapping cluster IDs to cell types
#' @return Seurat object with cell_type annotation
annotate_cell_types <- function(seurat_obj, cluster_to_celltype) {

  cat("Annotating cell types...\n")

  # Ensure cluster IDs match
  cluster_ids <- levels(Idents(seurat_obj))
  if (length(cluster_to_celltype) != length(cluster_ids)) {
    stop("Number of cell type labels (", length(cluster_to_celltype),
         ") does not match number of clusters (", length(cluster_ids), ")")
  }

  names(cluster_to_celltype) <- cluster_ids
  seurat_obj <- RenameIdents(seurat_obj, cluster_to_celltype)
  seurat_obj$cell_type <- Idents(seurat_obj)

  # Report cell type distribution
  cat("  Cell type distribution:\n")
  print(table(seurat_obj$cell_type))
  cat("\n")

  return(seurat_obj)
}


#' Create pseudobulk expression profiles
#'
#' @param seurat_obj Seurat object with cell_type annotation
#' @param group_by Columns to group by for aggregation
#' @return Pseudobulk Seurat object
create_pseudobulk <- function(seurat_obj, group_by) {

  cat("Creating pseudobulk profiles...\n")
  cat("  Grouping by:", paste(group_by, collapse = ", "), "\n")

  seurat_obj$cell_type <- Idents(seurat_obj)

  pseudobulk <- AggregateExpression(
    seurat_obj,
    assays = "RNA",
    return.seurat = TRUE,
    group.by = group_by
  )

  cat("  Pseudobulk samples created:", ncol(pseudobulk), "\n\n")

  return(pseudobulk)
}


#' Annotate hepatocyte zonation (periportal vs pericentral)
#'
#' Uses expression of zonation markers:
#' - Pericentral: Cyp2e1 > 1 AND Cyp2f2 < 1
#' - Periportal: Cyp2e1 < 1 AND Cyp2f2 > 1
#'
#' @param seurat_obj Seurat object (hepatocytes only)
#' @param pc_marker Pericentral marker (default: Cyp2e1)
#' @param pp_marker Periportal marker (default: Cyp2f2)
#' @return Seurat object with zonation annotation
annotate_zonation <- function(seurat_obj,
                               pc_marker = "Cyp2e1",
                               pp_marker = "Cyp2f2") {

  cat("Annotating hepatocyte zonation...\n")
  cells_before <- ncol(seurat_obj)

  # Subset pericentral cells
  heps_pc <- subset(seurat_obj,
                     subset = get(pc_marker) > 1 & get(pp_marker) < 1)
  cat("  Pericentral cells:", ncol(heps_pc), "\n")

  # Subset periportal cells
  heps_pp <- subset(seurat_obj,
                     subset = get(pc_marker) < 1 & get(pp_marker) > 1)
  cat("  Periportal cells:", ncol(heps_pp), "\n")

  # Merge
  seurat_obj <- merge(heps_pc, y = heps_pp)
  rm(heps_pc)
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

  # Annotate zonation type
  seurat_obj$zonation_type <- ifelse(
    colnames(seurat_obj) %in% colnames(heps_pp),
    "Periportal",
    "Pericentral"
  )
  rm(heps_pp)

  cells_after <- ncol(seurat_obj)
  cat("  Cells with clear zonation:", cells_after,
      "(", round(100 * cells_after / cells_before, 1), "% of hepatocytes)\n\n")

  return(seurat_obj)
}


################################################################################
# 1. FCG (FOUR CORE GENOTYPES) PRE-PROCESSING
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("1. FCG (FOUR CORE GENOTYPES) PRE-PROCESSING\n")
cat(strrep("=", 80), "\n\n")

cat("This section processes the main FCG study cohort.\n")
cat("43 FCG samples + 5 wild-type male controls.\n\n")

# NOTE: The following code is commented out because it requires access to
# the original raw data on the HPC system. Uncomment and run if you have
# access to the raw alignment data.

if (FALSE) {  # Set to TRUE to run FCG processing

  #---------------------------------------------------------------------------
  # 1.1 Import raw data
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("1.1 Importing raw data\n")
  cat(strrep("-", 60), "\n\n")

  setwd(config$fcg$base_dir)

  # Load metadata
  metadata <- readxl::read_xlsx(path = config$fcg$metadata_file, sheet = "metadata") %>%
    dplyr::filter(!is.na(sequencing_name)) %>%
    dplyr::mutate(sample_name = sequencing_name) %>%
    dplyr::arrange(sample_name)

  cat("Metadata loaded:", nrow(metadata), "samples\n")
  print(table(metadata$group))

  # Import samples
  fcg <- import_samples_seurat_raw(
    folder = config$fcg$alignment_dir,
    samples = unique(metadata$sample_name),
    gtf = config$fcg$gtf_file,
    metadata = metadata %>%
      dplyr::distinct(sample_name, group, gonadal_type, mouse_label, allosome_type)
  )

  saveRDS(fcg, file.path(config$fcg$base_dir, "imported.rds"))
  cat("Saved: imported.rds\n\n")

  #---------------------------------------------------------------------------
  # 1.2 Quality control: UMI filtering
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("1.2 UMI filtering\n")
  cat(strrep("-", 60), "\n\n")

  fcg <- filter_umi_outliers(
    fcg,
    lower_threshold = qc_params$lower_umi_threshold,
    upper_sd_mult = qc_params$upper_umi_sd_multiplier
  )

  saveRDS(fcg, file.path(config$fcg$base_dir, "filtered.rds"))
  cat("Saved: filtered.rds\n\n")

  #---------------------------------------------------------------------------
  # 1.3 Clustering
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("1.3 Clustering\n")
  cat(strrep("-", 60), "\n\n")

  fcg <- run_clustering_workflow(
    fcg,
    pca_dims = qc_params$pca_dims,
    resolution = qc_params$cluster_resolution
  )

  # Generate diagnostic plots
  if (has_SCpubr) {
    SCpubr::do_DimPlot(fcg, pt.size = 0.03, legend.position = "right")
    ggsave(file.path(config$fcg$base_dir, "clustered_dimplot.pdf"),
           width = 8, height = 8)

    SCpubr::do_DimPlot(fcg, pt.size = 0.03, legend.position = "right",
                       group.by = "group")
    ggsave(file.path(config$fcg$base_dir, "clustered_dimplot_group.pdf"),
           width = 9, height = 8)
  }

  saveRDS(fcg, file.path(config$fcg$base_dir, "clustered.rds"))
  cat("Saved: clustered.rds\n\n")

  #---------------------------------------------------------------------------
  # 1.4 Doublet detection and removal
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("1.4 Doublet removal\n")
  cat(strrep("-", 60), "\n\n")

  fcg <- remove_doublets(fcg, doublet_rate = qc_params$doublet_rate)

  # Visualize cell type markers
  if (has_SCpubr) {
    SCpubr::do_DotPlot(sample = fcg, features = cell_type_markers, flip = TRUE)
    ggsave(file.path(config$fcg$base_dir, "doubletless_celltypemarkersdotplot.pdf"),
           width = 10, height = 8)
  }

  # Find cluster markers (unbiased)
  cluster_markers <- FindAllMarkers(fcg, only.pos = TRUE)
  cluster_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>%
    ungroup() -> top5

  saveRDS(cluster_markers,
          file.path(config$fcg$base_dir, "doubletless_clustermarkers.rds"))
  saveRDS(fcg, file.path(config$fcg$base_dir, "doubletless.rds"))
  cat("Saved: doubletless.rds, doubletless_clustermarkers.rds\n\n")

  #---------------------------------------------------------------------------
  # 1.5 Cell type annotation
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("1.5 Cell type annotation\n")
  cat(strrep("-", 60), "\n\n")

  # FCG cluster-to-celltype mapping (based on marker expression)
  # Clusters: 0, 1, 2, 3, 4, 5, 6, 7, 8
  fcg_cluster_mapping <- c(
    "Hepatocytes",      # Cluster 0
    "Hepatocytes",      # Cluster 1
    "Endothelial cells", # Cluster 2
    "Immune cells",     # Cluster 3
    "Stellate cells",   # Cluster 4
    "Immune cells",     # Cluster 5
    "Immune cells",     # Cluster 6
    "Hepatocytes",      # Cluster 7
    "Cholangiocytes"    # Cluster 8
  )

  fcg <- annotate_cell_types(fcg, fcg_cluster_mapping)

  # Visualize
  if (has_SCpubr) {
    SCpubr::do_DimPlot(fcg, pt.size = 0.03, legend.position = "right", label = TRUE)
    ggsave(file.path(config$fcg$base_dir, "annotated_dimplot.pdf"),
           width = 8, height = 8)

    SCpubr::do_DotPlot(sample = fcg, features = cell_type_markers, flip = TRUE)
    ggsave(file.path(config$fcg$base_dir, "annotated_celltypemarkersdotplot.pdf"),
           width = 10, height = 8)
  }

  saveRDS(fcg, file.path(config$fcg$base_dir, "annotated.rds"))
  cat("Saved: annotated.rds\n\n")

  #---------------------------------------------------------------------------
  # 1.6 Pseudobulking
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("1.6 Pseudobulking\n")
  cat(strrep("-", 60), "\n\n")

  # Exclude wild-type controls for main analysis
  fcg_main <- subset(x = fcg, group != "XYT-WT")

  fcg_pseudobulks <- create_pseudobulk(
    fcg_main,
    group_by = c("sample_name", "cell_type", "allosome_type")
  )

  saveRDS(fcg_pseudobulks,
          file.path(config$fcg$base_dir, "pseudobulks_scc.rds"))
  cat("Saved: pseudobulks_scc.rds\n\n")

  #---------------------------------------------------------------------------
  # 1.7 Hepatocyte zonation annotation
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("1.7 Hepatocyte zonation annotation\n")
  cat(strrep("-", 60), "\n\n")

  # Subset hepatocytes
  fcg_heps <- subset(x = fcg, subset = cell_type == "Hepatocytes")
  cat("Hepatocytes:", ncol(fcg_heps), "\n")
  rm(fcg)

  saveRDS(fcg_heps, file.path(config$fcg$base_dir, "heps.rds"))

  # Annotate zonation
  fcg_heps <- annotate_zonation(fcg_heps)

  saveRDS(fcg_heps, file.path(config$fcg$base_dir, "heps_zonationsubsetted.rds"))
  cat("Saved: heps_zonationsubsetted.rds\n\n")

  # Create zonation pseudobulk
  fcg_heps_pseudobulks <- create_pseudobulk(
    fcg_heps,
    group_by = c("sample_name", "zonation_type", "gonadal_type")
  )

  saveRDS(fcg_heps_pseudobulks,
          file.path(config$fcg$base_dir, "heps_zonationsubsetted_pseudobulks_gonad.rds"))
  cat("Saved: heps_zonationsubsetted_pseudobulks_gonad.rds\n\n")

  rm(fcg_heps)
  gc()
}

cat("\nFCG pre-processing section complete.\n")
cat("(Set if(FALSE) to if(TRUE) to run with actual data)\n\n")


################################################################################
# 2. GDX (GONADECTOMY HORMONE REPLACEMENT) PRE-PROCESSING
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("2. GDX (GONADECTOMY HORMONE REPLACEMENT) PRE-PROCESSING\n")
cat(strrep("=", 80), "\n\n")

cat("This section processes the GDX validation cohort.\n")
cat("Gonadectomized mice with or without hormone replacement.\n\n")

if (FALSE) {  # Set to TRUE to run GDX processing

  #---------------------------------------------------------------------------
  # 2.1 Import raw data
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("2.1 Importing raw data\n")
  cat(strrep("-", 60), "\n\n")

  setwd(config$gdx$base_dir)

  # Load metadata
  metadata <- readxl::read_xlsx(path = config$gdx$metadata_file, sheet = "metadata") %>%
    dplyr::filter(!is.na(sequencing_name)) %>%
    dplyr::mutate(sample_name = sequencing_name) %>%
    dplyr::arrange(sample_name)

  cat("Metadata loaded:", nrow(metadata), "samples\n")
  print(table(metadata$group))

  # Import samples
  gdxhr <- import_samples_seurat_raw(
    folder = config$gdx$alignment_dir,
    samples = unique(metadata$sample_name),
    gtf = config$fcg$gtf_file,
    metadata = metadata %>% dplyr::distinct(sample_name, group, pump_status)
  )

  saveRDS(gdxhr, file.path(config$gdx$base_dir, "imported.rds"))
  cat("Saved: imported.rds\n\n")

  #---------------------------------------------------------------------------
  # 2.2 Quality control: UMI filtering
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("2.2 UMI filtering\n")
  cat(strrep("-", 60), "\n\n")

  gdxhr <- filter_umi_outliers(
    gdxhr,
    lower_threshold = qc_params$lower_umi_threshold,
    upper_sd_mult = qc_params$upper_umi_sd_multiplier
  )

  saveRDS(gdxhr, file.path(config$gdx$base_dir, "filtered.rds"))
  cat("Saved: filtered.rds\n\n")

  #---------------------------------------------------------------------------
  # 2.3 Clustering
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("2.3 Clustering\n")
  cat(strrep("-", 60), "\n\n")

  gdxhr <- run_clustering_workflow(
    gdxhr,
    pca_dims = qc_params$pca_dims,
    resolution = qc_params$cluster_resolution
  )

  # Generate diagnostic plots
  if (has_SCpubr) {
    SCpubr::do_DimPlot(gdxhr, pt.size = 0.03, legend.position = "right")
    ggsave(file.path(config$gdx$base_dir, "clustered_dimplot.pdf"),
           width = 8, height = 8)

    SCpubr::do_DimPlot(gdxhr, pt.size = 0.03, legend.position = "right",
                        group.by = "group")
    ggsave(file.path(config$gdx$base_dir, "clustered_dimplot_group.pdf"),
           width = 9, height = 8)
  }

  saveRDS(gdxhr, file.path(config$gdx$base_dir, "clustered.rds"))
  cat("Saved: clustered.rds\n\n")

  #---------------------------------------------------------------------------
  # 2.4 Doublet detection and removal
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("2.4 Doublet removal\n")
  cat(strrep("-", 60), "\n\n")

  gdxhr <- remove_doublets(gdxhr, doublet_rate = qc_params$doublet_rate)

  # Visualize cell type markers
  if (has_SCpubr) {
    SCpubr::do_DotPlot(sample = gdxhr, features = cell_type_markers, flip = TRUE)
    ggsave(file.path(config$gdx$base_dir, "doubletless_celltypemarkersdotplot.pdf"),
           width = 10, height = 8)
  }

  # Find cluster markers
  cluster_markers <- FindAllMarkers(gdxhr, only.pos = TRUE)
  cluster_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>%
    ungroup() -> top5

  saveRDS(gdxhr, file.path(config$gdx$base_dir, "doubletless.rds"))
  cat("Saved: doubletless.rds\n\n")

  #---------------------------------------------------------------------------
  # 2.5 Cell type annotation
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("2.5 Cell type annotation\n")
  cat(strrep("-", 60), "\n\n")

  # GDX cluster-to-celltype mapping (based on marker expression)
  # Note: Different cluster arrangement than FCG
  gdx_cluster_mapping <- c(
    "Hepatocytes",       # Cluster 0
    "Endothelial cells", # Cluster 1
    "Immune cells",      # Cluster 2
    "Stellate cells",    # Cluster 3
    "Hepatocytes",       # Cluster 4
    "Cholangiocytes"     # Cluster 5
  )

  gdxhr <- annotate_cell_types(gdxhr, gdx_cluster_mapping)

  # Visualize
  if (has_SCpubr) {
    SCpubr::do_DimPlot(gdxhr, pt.size = 0.03, legend.position = "right", label = TRUE)
    ggsave(file.path(config$gdx$base_dir, "annotated_dimplot.pdf"),
           width = 8, height = 8)

    SCpubr::do_DotPlot(sample = gdxhr, features = cell_type_markers, flip = TRUE)
    ggsave(file.path(config$gdx$base_dir, "annotated_celltypemarkersdotplot.pdf"),
           width = 10, height = 8)
  }

  saveRDS(gdxhr, file.path(config$gdx$base_dir, "annotated.rds"))
  cat("Saved: annotated.rds\n\n")

  #---------------------------------------------------------------------------
  # 2.6 Pseudobulking
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("2.6 Pseudobulking\n")
  cat(strrep("-", 60), "\n\n")

  # Without pump status subsetting
  gdxhr_pseudobulks <- create_pseudobulk(
    gdxhr,
    group_by = c("sample_name", "cell_type", "group")
  )

  saveRDS(gdxhr_pseudobulks,
          file.path(config$gdx$base_dir, "pseudobulks_withoutpumpstatussubsetting.rds"))
  cat("Saved: pseudobulks_withoutpumpstatussubsetting.rds\n")

  # With pump status subsetting (only samples with pump = "yes")
  gdxhr_pump <- subset(x = gdxhr, pump_status == "yes")
  gdxhr_pump_pseudobulks <- create_pseudobulk(
    gdxhr_pump,
    group_by = c("sample_name", "cell_type", "group")
  )

  saveRDS(gdxhr_pump_pseudobulks,
          file.path(config$gdx$base_dir, "pseudobulks_withpumpstatussubsetting.rds"))
  cat("Saved: pseudobulks_withpumpstatussubsetting.rds\n\n")

  rm(gdxhr_pump)

  #---------------------------------------------------------------------------
  # 2.7 Hepatocyte zonation annotation
  #---------------------------------------------------------------------------
  cat(strrep("-", 60), "\n")
  cat("2.7 Hepatocyte zonation annotation\n")
  cat(strrep("-", 60), "\n\n")

  # Subset hepatocytes
  gdxhr_heps <- subset(x = gdxhr, subset = cell_type == "Hepatocytes")
  cat("Hepatocytes:", ncol(gdxhr_heps), "\n")
  rm(gdxhr)

  saveRDS(gdxhr_heps, file.path(config$gdx$base_dir, "heps.rds"))

  # Annotate zonation
  gdxhr_heps <- annotate_zonation(gdxhr_heps)

  saveRDS(gdxhr_heps, file.path(config$gdx$base_dir, "heps_zonationsubsetted.rds"))
  cat("Saved: heps_zonationsubsetted.rds\n\n")

  # Create zonation pseudobulk (without pump subsetting)
  gdxhr_heps_pseudobulks <- create_pseudobulk(
    gdxhr_heps,
    group_by = c("sample_name", "zonation_type", "group")
  )

  saveRDS(gdxhr_heps_pseudobulks,
          file.path(config$gdx$base_dir, "heps_zonationsubsetted_pseudobulks_withoutpumpstatussubsetting.rds"))
  cat("Saved: heps_zonationsubsetted_pseudobulks_withoutpumpstatussubsetting.rds\n")

  # With pump subsetting
  gdxhr_heps_pump <- subset(x = gdxhr_heps, pump_status == "yes")
  gdxhr_heps_pump_pseudobulks <- create_pseudobulk(
    gdxhr_heps_pump,
    group_by = c("sample_name", "zonation_type", "group")
  )

  saveRDS(gdxhr_heps_pump_pseudobulks,
          file.path(config$gdx$base_dir, "heps_zonationsubsetted_pseudobulks_withpumpstatussubsetting.rds"))
  cat("Saved: heps_zonationsubsetted_pseudobulks_withpumpstatussubsetting.rds\n\n")

  rm(gdxhr_heps)
  rm(gdxhr_heps_pump)
  gc()
}

cat("\nGDX pre-processing section complete.\n")
cat("(Set if(FALSE) to if(TRUE) to run with actual data)\n\n")


################################################################################
# WORKFLOW COMPLETE
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("WORKFLOW 0: snRNA-seq PRE-PROCESSING COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("This workflow documents the pre-processing pipeline that produces:\n")
cat("  - FCG pseudobulk data (main study cohort)\n")
cat("  - GDX pseudobulk data (gonadectomy validation cohort)\n")
cat("  - Zonation-annotated hepatocyte pseudobulk data\n\n")

cat("Output files from FCG processing:\n")
cat("  - imported.rds              Raw imported Seurat object\n")
cat("  - filtered.rds              UMI-filtered object\n")
cat("  - clustered.rds             Clustered object with UMAP\n")
cat("  - doubletless.rds           Doublet-removed object\n")
cat("  - annotated.rds             Cell type annotated object\n")
cat("  - pseudobulks_scc.rds       Main pseudobulk data\n")
cat("  - heps.rds                  Hepatocytes subset\n")
cat("  - heps_zonationsubsetted.rds  Zonation-annotated hepatocytes\n")
cat("  - heps_zonationsubsetted_pseudobulks_gonad.rds  Zonation pseudobulk\n\n")

cat("Output files from GDX processing:\n")
cat("  - imported.rds through annotated.rds  (same as FCG)\n")
cat("  - pseudobulks_withoutpumpstatussubsetting.rds\n")
cat("  - pseudobulks_withpumpstatussubsetting.rds\n")
cat("  - heps_zonationsubsetted_pseudobulks_*.rds\n\n")

cat("Next steps:\n")
cat("  1. Run workflow 01: Data preparation and QC\n")
cat("     source('scripts/workflows/01_data_preparation_qc.R')\n\n")

cat("NOTE: This workflow requires access to raw alignment data.\n")
cat("      For users starting with pre-processed data, begin with workflow 01.\n\n")

cat(strrep("=", 80), "\n\n")
