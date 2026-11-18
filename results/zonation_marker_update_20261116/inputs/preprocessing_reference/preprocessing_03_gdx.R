# ### 3.0 GDX pre-processing ###

# # preprocessing_01_setup.R provides import_samples_seurat_raw() and the shared
# # packages. A bare relative source() resolves against the working directory,
# # which under LSF is not the directory holding this script, so give the path
# # explicitly and fail loudly if it ever moves.
# setup_file <- "/home/j083d/scripts/fcgl_paper/preprocessing_01_setup.R"
# stopifnot(file.exists(setup_file))
# source(setup_file)

# ## Import data ##

# metadata <- readxl::read_xlsx(
#   path = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/sci_metadata.xlsx",
#   sheet = "metadata"
# ) %>%
#   dplyr::filter(!is.na(sequencing_name)) %>%
#   dplyr::mutate(sample_name = sequencing_name) %>%
#   dplyr::filter(!(group %in% c("f_g_e","f_g_t","m_g_e","m_g_t"))) %>%   # <- remove groups here
#   dplyr::arrange(sample_name)

# gdx <- import_samples_seurat_raw(
#   folder = "/omics/odcf/analysis/OE0538_projects/DO-0012/processed_sci-seq_data/41064_gdxhr/41064/alignment",
#   samples = unique(metadata$sample_name),
#   gtf = "/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf",
#   metadata = metadata %>% dplyr::distinct(sample_name, group)
# )

# saveRDS(gdx, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/imported.rds")

# ## Remove UMI outliers ##

# lower_umi_threshold = 500
# gdx <- subset(gdx, subset = nCount_RNA > lower_umi_threshold)
# upper_umi_threshold <- base::mean(gdx$nCount_RNA) + (3 * stats::sd(gdx$nCount_RNA))
# gdx <- subset(gdx, subset = nCount_RNA < upper_umi_threshold)

# mean(gdx$nCount_RNA)
# mean(gdx$nFeature_RNA)

# ## Clustering ##

# gdx <- NormalizeData(gdx)
# gdx <- FindVariableFeatures(gdx)
# gdx <- ScaleData(gdx)
# gdx <- RunPCA(gdx, seed.use = 123)
# gdx <- FindNeighbors(gdx, dims = 1:20)
# gdx <- FindClusters(gdx, resolution = 0.03, random.seed = 123)
# gdx <- RunUMAP(gdx, dims = 1:20, seed.use = 123)

# SCpubr::do_DimPlot(gdx, pt.size = 0.03, legend.position = "right")
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/clustered_dimplot.pdf", width = 8, height = 8)
# SCpubr::do_DimPlot(gdx, pt.size = 0.03, legend.position = "right", group.by = "group")
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/clustered_dimplot_group.pdf", width = 9, height = 8)

# cell_type_markers <- list("Hepatocytes" = c("Alb","Apoa2","Hnf4a"),
#                           "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
#                           "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
#                           "Stellate cells" = c("Lrat","Des","Reln"),
#                           "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
# SCpubr::do_DotPlot(sample = gdx, features = cell_type_markers, flip = TRUE)
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/clustered_celltypemarkersdotplot.pdf", width = 10, height = 8)

# ## Doublet removal ##

# set.seed(123)
# sce <- scDblFinder(GetAssayData(gdx, layer="counts"),
#                    dbr = .02,
#                    clusters=Idents(gdx))
# gdx$scDblFinder.score <- sce$scDblFinder.score
# gdx$scDblFinder.class <- sce$scDblFinder.class
# rm(sce)

# knitr::kable(table(gdx$scDblFinder.class), align = "c")
# doublets <- rownames(gdx[,gdx@meta.data$scDblFinder.class == "doublet"]@meta.data)

# cells_prior <- base::ncol(gdx)
# gdx <- subset(gdx, subset = scDblFinder.class == "singlet")
# cells_post <- base::ncol(gdx)
# cells_lost <- cells_prior - cells_post

# knitr::kable(
#   base::data.frame(
#     cells_prior = cells_prior,
#     cells_post = cells_post,
#     cells_lost = cells_lost,
#     perc_lost = round(cells_lost / cells_prior * 100, 1)
#   )
# )

# SCpubr::do_DimPlot(gdx, pt.size = 0.03, legend.position = "right", label = TRUE)
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/doubletless_dimplot.pdf", width = 8, height = 8)

# cell_type_markers <- list("Hepatocytes" = c("Alb","Apoa2","Hnf4a"),
#                           "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
#                           "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
#                           "Stellate cells" = c("Lrat","Des","Reln"),
#                           "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
# SCpubr::do_DotPlot(sample = gdx, features = cell_type_markers, flip = TRUE)
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/doubletless_celltypemarkersdotplot.pdf", width = 10, height = 8)

# mean(gdx$nCount_RNA)
# mean(gdx$nFeature_RNA)

# ## Broad cell type annotation and pseudobulking ##

# new.cluster.ids <- c("Hepatocytes", "Hepatocytes", "Endothelial cells","Immune cells","Stellate cells","Hepatocytes", "Hepatocytes", "Cholangiocytes")
# names(new.cluster.ids) <- levels(gdx)
# gdx <- RenameIdents(gdx, new.cluster.ids)

# SCpubr::do_DimPlot(gdx, pt.size = 0.03, legend.position = "right", label = TRUE)
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/annotated_dimplot.pdf", width = 8, height = 8)

# # biased marker dotplot
# cell_type_markers <- list("Hepatocytes" = c("Alb","Apoa2","Hnf4a"),
#                           "Immune cells" = c("Ptprc","Vsig4","Adgre1"),
#                           "Endothelial cells" = c("Btnl9","Clec4g","Pecam1"),
#                           "Stellate cells" = c("Lrat","Des","Reln"),
#                           "Cholangiocytes" = c("Krt19","Ankrd1","Epcam"))
# SCpubr::do_DotPlot(sample = gdx, features = cell_type_markers, flip = TRUE)
# ggsave("/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/annotated_celltypemarkersdotplot.pdf", width = 10, height = 8)

# saveRDS(gdx, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/annotated.rds")

# # pseudobulk

# gdx$cell_type <- Seurat::Idents(gdx)
# gdx_pseudobulks <- AggregateExpression(gdx, assays = "RNA", return.seurat = T, group.by = c("sample_name", "cell_type", "group"))
# saveRDS(gdx_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/gdx_pseudobulks.rds")

# ## Hepatocyte zonation annotation and pseudobulking ##

# gdx$cell_type <- Seurat::Idents(gdx)
# gdx_heps <- subset(x = gdx, subset = cell_type == "Hepatocytes")
# gdx_heps

# gdx_heps_pc <- subset(gdx_heps, subset = Cyp2e1 > 1 & Cyp2f2 < 1)
# gdx_heps_pp <- subset(gdx_heps, subset = Cyp2e1 < 1 & Cyp2f2 > 1)
# gdx_heps_pc$zonation_type <- "Pericentral"
# gdx_heps_pp$zonation_type <- "Periportal"
# gdx_heps <- merge(gdx_heps_pc, y = gdx_heps_pp)
# gdx_heps <- JoinLayers(gdx_heps, assay = "RNA")

# saveRDS(gdx_heps, "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/pcppheps.rds")

# # pseudobulk
# gdx_zonation_pseudobulks <- AggregateExpression(gdx_heps, assays = "RNA", return.seurat = T, group.by = c("sample_name", "zonation_type", "group"))
# saveRDS(gdx_zonation_pseudobulks, file = "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci/gdx_zonation_pseudobulks.rds")


## ALTERNATE hepatocyte zonation annotation and pseudobulking ##
#
# This is an ALTERNATE assignment, retained ALONGSIDE the Cyp2e1/Cyp2f2 rule
# above - it does not replace it. The original assignment remains the basis of
# the published analyses; this one exists so the reviewer can see what changes
# when zonation is called from the whole marker panel instead of two genes.
#
# Rationale (Reviewer 3 point 3). In the FCG data Cyp2f2 is male-biased within
# the periportal zone itself (t = +9.57) while Cyp2e1 is close to neutral
# within the pericentral zone (t = -0.47), so the binary rule carries a net sex
# bias of -10.0 where a score over the full marker set carries +0.9. This block
# is identical in preprocessing_02/03/04/05 apart from the cohort constants and
# the single step that loads the hepatocytes, so all four datasets are assigned
# on a common basis with no per-cohort tuning. Whether the binary rule is
# biased in any one cohort is an empirical question, answered by the comparison
# section below rather than assumed here.
#
# Requires UCell (BiocManager::install("UCell")); not loaded in setup.
#
# RUNS STANDALONE. This block loads the hepatocytes itself and reproduces any
# cohort filtering performed above, so it can be executed on its own without
# running any of the preceding sections. It works on its own object,
# `heps_score`, and does not read or modify any object created above.

library(Seurat)
library(dplyr)
library(Matrix)
library(UCell)
set.seed(123)

## ---- cohort constants ------------------------------------------------------
# GDX: gonadectomy/sham x sex, vehicle only; the hormone-treated groups
# (f_g_e, f_g_t, m_g_e, m_g_t) were excluded at import, so no filter is needed.
dir           <- "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci"
id_col        <- "sample_name"
group_col     <- "group"
new_heps_file <- "pcppheps_score.rds"
new_pb_file   <- "gdx_zonation_pseudobulks_score.rds"

## ---- cohort-specific: obtain the hepatocytes as `heps_score` ---------------
annot <- readRDS(file.path(dir, "annotated.rds"))
if (!"cell_type" %in% colnames(annot[[]])) annot$cell_type <- Seurat::Idents(annot)
heps_score <- subset(x = annot, subset = cell_type == "Hepatocytes")
rm(annot); gc()
## ----------------------------------------------------------------------------

stopifnot(all(c(id_col, group_col) %in% colnames(heps_score[[]])))

# AggregateExpression substitutes "_" inside values when it builds profile
# names, so the animal id is normalised here - that way the id held in the
# metadata and the id appearing in the profile name are the same string.
# Idempotent: harmless where the id carries no underscore, or where the
# substitution has already been applied upstream.
new_ids <- gsub("_", "-", as.character(heps_score[[]][[id_col]]))
names(new_ids) <- colnames(heps_score)
heps_score[[id_col]] <- new_ids
rm(new_ids)
stopifnot(!anyNA(heps_score[[]][[id_col]]))

# Standardisation reference for the two signature scores.
#
# FALSE reproduces the assignment behind the numbers in the current reviewer
# response: PC and PP scores are standardised over ALL nuclei, including those
# that will end up unassigned. Nuclei lacking markers at one pole score zero
# there, and periportal markers are detected more often than pericentral ones,
# so the pericentral spike at zero is the larger of the two and the global cut
# sits slightly toward pericentral. This moves where the cut falls, not the
# ranking of samples, so it affects the reported proportions and not the
# between-group test - confirmed empirically in FCG, where the trough cut and
# the zero cut gave identical p (0.0030).
#
# TRUE standardises over scorable nuclei only, which is cleaner but changes the
# reported proportions and would require the response letter to be updated.
scale_on_scorable <- FALSE

if (sum(grepl("^counts", Layers(heps_score, assay = "RNA"))) > 1) {
  heps_score <- JoinLayers(heps_score, assay = "RNA")
}
if (!"data" %in% Layers(heps_score, assay = "RNA")) {
  heps_score <- NormalizeData(heps_score)
}

## Marker sets ##

zonation_gene_file <- "/omics/odcf/analysis/OE0538_projects/DO-0012/james_misc/gene_lists/Zonation_Genes.csv"
stopifnot(file.exists(zonation_gene_file))
zonation_genes <- read.csv(zonation_gene_file, stringsAsFactors = FALSE)
names(zonation_genes)[1:2] <- c("gene", "zone")

pc_markers <- base::intersect(zonation_genes$gene[zonation_genes$zone == "Pericentral"],
                              rownames(heps_score))
pp_markers <- base::intersect(zonation_genes$gene[zonation_genes$zone == "Periportal"],
                              rownames(heps_score))
# A gene listed under both zones would be scored at both poles and mislabelled
# in the detection table, so refuse rather than silently double-count it.
stopifnot(length(base::intersect(pc_markers, pp_markers)) == 0)
# Midlobular markers are excluded, as in the FCG analysis: they carry almost no
# zonation signal while being the most sex-biased group in the marker set.

# Drop markers detected in >95% or <5% of nuclei - no discriminative power.
# Which genes this removes is dataset-specific and is recorded in the exported
# detection table. It is NOT assumed to match the set dropped in the FCG data.
counts_mat <- GetAssayData(heps_score, assay = "RNA", layer = "counts")
marker_detection <- Matrix::rowMeans(counts_mat[c(pc_markers, pp_markers), , drop = FALSE] > 0)

detection_table <- data.frame(
  gene = names(marker_detection),
  detection = as.numeric(marker_detection),
  zone = ifelse(names(marker_detection) %in% pc_markers, "Pericentral", "Periportal")
) %>%
  dplyr::mutate(kept = detection >= 0.05 & detection <= 0.95) %>%
  dplyr::arrange(zone, dplyr::desc(detection))
write.csv(detection_table, file.path(dir, "zonation_score_marker_detection.csv"),
          row.names = FALSE)

pc_markers <- detection_table$gene[detection_table$zone == "Pericentral" & detection_table$kept]
pp_markers <- detection_table$gene[detection_table$zone == "Periportal"  & detection_table$kept]
message("markers used - pericentral: ", length(pc_markers), ", periportal: ", length(pp_markers))
stopifnot(length(pc_markers) >= 3, length(pp_markers) >= 3)

## Signature scores ##

# UCell ranks genes within each nucleus, so the score is robust to dropout and
# is cell-intrinsic: the same threshold means the same thing in every cohort
# without rescoring against each dataset's own background. This is why the
# identical procedure can be applied across the four datasets.
# maxRank is taken from the data because the 1500 default assumes whole cells.
genes_detected <- Matrix::colSums(counts_mat > 0)
ucell_max_rank <- max(1000, round(stats::median(genes_detected)))
message("median genes detected per nucleus: ", stats::median(genes_detected),
        " | UCell maxRank: ", ucell_max_rank)

heps_score <- AddModuleScore_UCell(
  heps_score,
  features = list(PC = pc_markers, PP = pp_markers),
  assay    = "RNA",
  maxRank  = ucell_max_rank,
  name     = "_UCell",
  ncores   = 4
)

## Assignment ##

# Nuclei carrying no marker signal at one pole cannot be placed and are left
# unassigned rather than pushed to an extreme. Scorability is determined before
# standardisation so that it can serve as the standardisation reference.
n_pc_detected <- Matrix::colSums(counts_mat[pc_markers, , drop = FALSE] > 0)
n_pp_detected <- Matrix::colSums(counts_mat[pp_markers, , drop = FALSE] > 0)
heps_score$n_pc_detected <- n_pc_detected
heps_score$n_pp_detected <- n_pp_detected
heps_score$scorable <- n_pc_detected > 0 & n_pp_detected > 0

# Standardise before subtracting, so neither pole dominates through detection
# rate: periportal markers are detected several-fold more often than pericentral
# ones. With scale_on_scorable = FALSE this is identical to scale().
score_ref <- if (scale_on_scorable) heps_score$scorable else rep(TRUE, ncol(heps_score))
stopifnot(sum(score_ref) > 1)
standardise <- function(x) (x - base::mean(x[score_ref])) / stats::sd(x[score_ref])
heps_score$PC_z <- standardise(heps_score$PC_UCell)
heps_score$PP_z <- standardise(heps_score$PP_UCell)
heps_score$zonation_score <- heps_score$PC_z - heps_score$PP_z
heps_score$zonation_score[!heps_score$scorable] <- NA_real_

# Every scorable nucleus is classified, at a single global cut applied
# identically to all samples. A per-sample cut would force equal zone
# proportions in every animal and make the compositional comparison
# untestable. Zero is the natural cut on a standardised difference and gives
# an approximately balanced global split; because the tests below depend only
# on the ranking of samples, shifting this cut changes the proportions reported
# but not the inference.
heps_score$zonation_type_score <- ifelse(
  is.na(heps_score$zonation_score), NA_character_,
  ifelse(heps_score$zonation_score > 0, "Pericentral", "Periportal")
)

message("scorable nuclei: ", sum(heps_score$scorable), " / ", ncol(heps_score),
        " (", round(100 * mean(heps_score$scorable), 1), "%)")
print(table(heps_score$zonation_type_score, useNA = "ifany"))

# Comparison with the original assignment, for the record
zonation_markers <- FetchData(heps_score, vars = c("Cyp2e1", "Cyp2f2"), layer = "data")
heps_score$zonation_type_binary <- ifelse(
  zonation_markers$Cyp2e1 > 1 & zonation_markers$Cyp2f2 < 1, "Pericentral",
  ifelse(zonation_markers$Cyp2e1 < 1 & zonation_markers$Cyp2f2 > 1, "Periportal", NA)
)
print(table(binary = heps_score$zonation_type_binary,
            score  = heps_score$zonation_type_score, useNA = "ifany"))

score_metadata <- heps_score[[]][, base::intersect(
  c(id_col, group_col, "sample_name", "group", "gonadal_type", "allosome_type",
    "diet", "sex", "nCount_RNA", "nFeature_RNA",
    "PC_UCell", "PP_UCell", "PC_z", "PP_z", "zonation_score",
    "n_pc_detected", "n_pp_detected", "scorable",
    "zonation_type_score", "zonation_type_binary"),
  colnames(heps_score[[]])), drop = FALSE]

# write.csv does not compress a .gz filename - it writes plain text under a
# misleading extension - so open the connection explicitly.
score_metadata_con <- gzfile(file.path(dir, "zonation_score_per_nucleus.csv.gz"), "w")
write.csv(score_metadata, score_metadata_con)
close(score_metadata_con)

saveRDS(heps_score, file.path(dir, new_heps_file))

# pseudobulk
heps_score_zoned <- subset(
  heps_score,
  cells = colnames(heps_score)[!is.na(heps_score$zonation_type_score)]
)
heps_score_zoned$zonation_type <- as.character(heps_score_zoned$zonation_type_score)

zonation_pseudobulks_score <- AggregateExpression(
  heps_score_zoned, assays = "RNA", return.seurat = TRUE,
  group.by = c(id_col, "zonation_type", group_col)
)
saveRDS(zonation_pseudobulks_score, file.path(dir, new_pb_file))

rm(heps_score, heps_score_zoned, counts_mat); gc()




## Comparison of the original and alternate zonation assignments ##
#
# Reviewer 3 point 3. Quantifies what changes between the two assignments at
# the two levels that matter:
#
#   A. the nucleus subsets      - basis of the compartment proportions
#   B. the zonation pseudobulks - basis of the zonation-resolved gene analyses
#
# This block is identical in preprocessing_02/03/04/05; only the constants
# below differ. In some cohorts the original assignment merges two subsets, so
# the old object holds only the nuclei the binary rule could call; in others it
# annotates in place and holds all hepatocytes with NA. The code handles both.
#
# RUNS STANDALONE. Reads only the .rds files written above. References no
# object from any earlier section.

library(Seurat)
library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)
set.seed(123)

## ---- cohort constants ------------------------------------------------------
dir           <- "/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdx_sci"
id_col        <- "sample_name"
group_col     <- "group"
old_heps_file <- "pcppheps.rds"
new_heps_file <- "pcppheps_score.rds"
old_pb_file   <- "gdx_zonation_pseudobulks.rds"
new_pb_file   <- "gdx_zonation_pseudobulks_score.rds"

comparison_dir <- file.path(dir, "zonation_assignment_comparison")

dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)

FDR   <- 0.05
ZONES <- c("Pericentral", "Periportal")

# Group-polarized zonation genes (section B2b). A gene counts as polarized
# toward one group if zonation is strong there and absent in the other.
# POLAR_HI / POLAR_LO are the manuscript's FCG definition; POLAR_GATE
# additionally requires the interaction contrast to be significant.
POLAR_HI   <- 5.4
POLAR_LO   <- 1
POLAR_GATE <- FALSE

## ============================ A. Nucleus subsets ============================

heps_old <- readRDS(file.path(dir, old_heps_file))
heps_new <- readRDS(file.path(dir, new_heps_file))

md_old <- heps_old[[]]; md_new <- heps_new[[]]

# Depending on the dataset the original object either holds only the assigned
# nuclei (merged PC + PP) or all hepatocytes with NA for the unassigned ones.
# Either way the new object holds all hepatocytes, so md_new supplies the
# animal and group covariates for every nucleus either assignment can place.
old_call <- setNames(as.character(md_old$zonation_type),       rownames(md_old))
new_call <- setNames(as.character(md_new$zonation_type_score), rownames(md_new))

all_cells <- union(names(old_call), names(new_call))
stopifnot(all(all_cells %in% rownames(md_new)))

cmp <- tibble(
  cell = all_cells,
  old  = unname(old_call[all_cells]),
  new  = unname(new_call[all_cells]),
  id   = as.character(md_new[all_cells, id_col]),
  grp  = as.character(md_new[all_cells, group_col])
)
stopifnot(!anyNA(cmp$id), !anyNA(cmp$grp))
GROUPS <- sort(unique(cmp$grp))
message("groups found: ", paste(GROUPS, collapse = ", "))

cat("\n== A1. assignment cross-tabulation ==\n")
ct <- table(old = cmp$old, new = cmp$new, useNA = "ifany"); print(ct)
write.csv(as.data.frame(ct), file.path(comparison_dir, "A1_crosstab.csv"), row.names = FALSE)

both <- !is.na(cmp$old) & !is.na(cmp$new)
cat(sprintf("\nnuclei assigned: old %d, new %d, both %d\n",
            sum(!is.na(cmp$old)), sum(!is.na(cmp$new)), sum(both)))
cat(sprintf("agreement where both assign: %.2f%%\n", 100 * mean(cmp$old[both] == cmp$new[both])))

cat("\n== A2. where the two differ, split by group ==\n")
trans <- cmp %>%
  mutate(transition = paste0(ifelse(is.na(old), "none", old), " -> ",
                             ifelse(is.na(new), "none", new))) %>%
  count(transition, grp) %>%
  pivot_wider(names_from = grp, values_from = n, values_fill = 0) %>%
  mutate(across(all_of(GROUPS), ~ round(100 * .x / sum(.x), 2), .names = "pct_{.col}")) %>%
  mutate(total = rowSums(across(all_of(GROUPS)))) %>%
  arrange(desc(total))
# a single enrichment ratio is only meaningful with two groups
if (length(GROUPS) == 2) {
  trans$ratio <- round(trans[[paste0("pct_", GROUPS[1])]] /
                       pmax(trans[[paste0("pct_", GROUPS[2])]], 1e-9), 2)
}
print(as.data.frame(trans))
write.csv(trans, file.path(comparison_dir, "A2_transitions_by_group.csv"), row.names = FALSE)

cat("\n== A3. compartment proportions per animal, both assignments ==\n")
props <- bind_rows(
  cmp %>% filter(!is.na(old)) %>% count(id, grp, zone = old) %>% mutate(method = "original"),
  cmp %>% filter(!is.na(new)) %>% count(id, grp, zone = new) %>% mutate(method = "alternate")
)

# count() omits empty combinations, so an animal with no nuclei in one zone
# would have no row for it and would be dropped from the comparison entirely
# instead of contributing a zero. Fill those in - but only when it actually
# happens, so that the common path stays exactly the one that produced the
# numbers already reported, and the fill announces itself when it is used.
missing_zone <- props %>% count(method, id, grp, name = "n_zones") %>%
  filter(n_zones < length(ZONES))
if (nrow(missing_zone) > 0) {
  message("  filling ", nrow(missing_zone), " animal-zone combinations with zero: ",
          paste(unique(missing_zone$id), collapse = ", "))
  props <- props %>%
    group_by(method) %>%
    tidyr::complete(tidyr::nesting(id, grp), zone = ZONES, fill = list(n = 0L)) %>%
    ungroup()
}
props <- props %>% group_by(method, id, grp) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(props, file.path(comparison_dir, "A3_per_sample_proportions.csv"), row.names = FALSE)

# an animal with no assigned nuclei at all under one method is still absent
lopsided <- props %>% distinct(method, id) %>% count(id) %>% filter(n < 2)
if (nrow(lopsided) > 0) {
  message("  animals assigned under only one method: ", paste(lopsided$id, collapse = ", "))
}

pc_frac <- props %>% filter(zone == "Pericentral")

group_summary <- pc_frac %>%
  group_by(method, grp) %>%
  summarise(n_animals = dplyr::n(), mean_frac = mean(frac), sd_frac = sd(frac), .groups = "drop")
print(as.data.frame(group_summary))
write.csv(group_summary, file.path(comparison_dir, "A3_group_summary.csv"), row.names = FALSE)

# all pairwise comparisons; with two groups this is the single test reported
# for the FCG data, with more it keeps every contrast on the same footing
safe_wilcox <- function(a, b) {
  tryCatch(suppressWarnings(wilcox.test(a, b)$p.value), error = function(e) NA_real_)
}
group_pairs <- utils::combn(GROUPS, 2, simplify = FALSE)
pairwise <- lapply(c("original", "alternate"), function(m) {
  d <- pc_frac %>% filter(method == m)
  lapply(group_pairs, function(p) {
    a <- d$frac[d$grp == p[1]]; b <- d$frac[d$grp == p[2]]
    tibble(method = m, group_a = p[1], group_b = p[2],
           n_a = length(a), n_b = length(b),
           mean_a = mean(a), mean_b = mean(b),
           difference = mean(a) - mean(b),
           relative_pct = 100 * (mean(a) / mean(b) - 1),
           p_wilcox = safe_wilcox(a, b))
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  # more than one pair means more than one test
  group_by(method) %>% mutate(p_BH = p.adjust(p_wilcox, "BH")) %>% ungroup()
print(as.data.frame(pairwise))
write.csv(pairwise, file.path(comparison_dir, "A3_pairwise_tests.csv"), row.names = FALSE)

ggsave(file.path(comparison_dir, "A3_proportions.pdf"),
       ggplot(pc_frac, aes(grp, frac, fill = grp)) +
         geom_boxplot(outlier.shape = NA, alpha = .5) +
         geom_jitter(width = .15, size = 1.2) +
         facet_wrap(~method) + ylim(0, 1) +
         labs(y = "pericentral fraction", x = NULL) + theme_bw() +
         theme(legend.position = "none"),
       width = 7, height = 4)

rm(heps_old, heps_new); gc()

## ========================== B. Zonation pseudobulks =========================

pb_old <- readRDS(file.path(dir, old_pb_file))

# The alternate zonation pseudobulk must be split by zone. Verify, and rebuild
# from the per-nucleus object if what is on disk is not zone-split - an earlier
# file aggregated by animal only would give profile names with no zone field
# and make the zonation contrasts impossible to fit.
pb_new_path <- file.path(dir, new_pb_file)

zone_split <- function(obj) {
  md <- obj[[]]
  if ("zonation_type" %in% colnames(md) &&
      all(ZONES %in% as.character(md$zonation_type))) return(TRUE)
  nm <- colnames(GetAssayData(obj, assay = "RNA", layer = "counts"))
  all(ZONES %in% unlist(strsplit(nm, "_", fixed = TRUE)))
}

pb_new <- if (file.exists(pb_new_path)) readRDS(pb_new_path) else NULL

if (is.null(pb_new) || !zone_split(pb_new)) {
  message("alternate pseudobulk is not split by zone - rebuilding from ", new_heps_file)
  heps_rebuild <- readRDS(file.path(dir, new_heps_file))
  if (sum(grepl("^counts", Layers(heps_rebuild, assay = "RNA"))) > 1) {
    heps_rebuild <- JoinLayers(heps_rebuild, assay = "RNA")
  }
  heps_rebuild <- subset(
    heps_rebuild,
    cells = colnames(heps_rebuild)[!is.na(heps_rebuild$zonation_type_score)]
  )
  heps_rebuild$zonation_type <- as.character(heps_rebuild$zonation_type_score)
  stopifnot(all(ZONES %in% heps_rebuild$zonation_type))

  pb_new <- AggregateExpression(
    heps_rebuild, assays = "RNA", return.seurat = TRUE,
    group.by = c(id_col, "zonation_type", group_col)
  )
  saveRDS(pb_new, pb_new_path)
  rm(heps_rebuild); gc()
  stopifnot(zone_split(pb_new))
}

# Pull counts and rebuild the design from the profile metadata. AggregateExpression
# replaces "_" with "-" inside values, so both spellings of a group label are
# accepted when falling back to the profile names.
GROUP_TOKENS <- unique(c(GROUPS, gsub("_", "-", GROUPS)))

prep <- function(obj) {
  m  <- GetAssayData(obj, assay = "RNA", layer = "counts")
  md <- obj[[]]
  need <- c(id_col, "zonation_type", group_col)

  if (all(need %in% colnames(md))) {
    smp  <- as.character(md[[id_col]])
    zone <- as.character(md$zonation_type)
    grp  <- as.character(md[[group_col]])
  } else {
    # AggregateExpression does not reliably rebuild the group.by variables as
    # metadata columns, and the profile-name format varies by Seurat version
    # and by whether underscores inside values were substituted. Rather than
    # assume a field count or ordering, find the zone and group tokens by
    # value and treat whatever remains as the animal id.
    message("  metadata columns absent; rebuilding from profile names")
    message("  example profile names: ",
            paste(utils::head(colnames(m), 3), collapse = " | "))
    parts <- strsplit(colnames(m), "_", fixed = TRUE)

    pick <- function(x, allowed) {
      hit <- x[x %in% allowed]
      if (length(hit) == 1L) hit else NA_character_
    }
    zone <- vapply(parts, pick, "", allowed = ZONES)
    grp  <- vapply(parts, pick, "", allowed = GROUP_TOKENS)
    smp  <- mapply(function(x, z, g) paste(setdiff(x, c(z, g)), collapse = "_"),
                   parts, zone, grp, USE.NAMES = FALSE)

    if (anyNA(zone) || anyNA(grp)) {
      bad <- colnames(m)[is.na(zone) | is.na(grp)]
      no_zone <- !any(unlist(parts) %in% ZONES)
      stop("could not parse zone/group from these profile names: ",
           paste(utils::head(bad, 5), collapse = " | "),
           if (no_zone)
             "\n  NO zone field is present in any profile name: this pseudobulk was aggregated without zonation_type in group.by and must be rebuilt"
           else
             paste0("\n  expected one of {", paste(ZONES, collapse = ", "),
                    "} and one of {", paste(GROUP_TOKENS, collapse = ", "),
                    "} among the underscore-separated fields"))
    }
  }

  # the two pseudobulk sets must label groups the same way to be comparable
  grp <- gsub("_", "-", grp)

  stopifnot(all(zone %in% ZONES), length(smp) == ncol(m))
  message("  ", ncol(m), " profiles | ", length(unique(smp)), " animals | zones: ",
          paste(names(table(zone)), table(zone), sep = "=", collapse = ", "),
          " | groups: ", paste(sort(unique(grp)), collapse = ", "))

  list(counts = as.matrix(m),
       id     = smp,
       zone   = factor(zone, levels = c("Periportal", "Pericentral")),
       grp    = grp)
}

# limma with the design used in the Methods: zonation, group and their
# interaction, blocking on animal identity with duplicateCorrelation because
# each animal contributes paired zonation pseudobulks. The contrast set is
# built programmatically so it adapts to however many groups a dataset has;
# with two groups it reduces exactly to the set used for the FCG data.
fit_pb <- function(d) {
  # keep only animals contributing both compartments, so the block structure is
  # balanced and duplicateCorrelation is estimable
  keep <- d$id %in% names(which(table(d$id) == 2))
  if (!all(keep)) message("  dropping ", sum(!keep), " unpaired profiles")
  d$counts <- d$counts[, keep, drop = FALSE]
  d$id     <- d$id[keep]
  d$zone   <- droplevels(d$zone[keep])
  d$grp    <- d$grp[keep]

  # group labels go into design column names, so they must be syntactically
  # valid - a "-" would otherwise be read as subtraction by makeContrasts
  glab <- make.names(d$grp)
  stopifnot(!anyDuplicated(unique(data.frame(d$grp, glab))$glab))   # no label collision

  grp_cell <- factor(paste(as.character(d$zone), glab, sep = "_"))
  design <- model.matrix(~ 0 + grp_cell)
  colnames(design) <- levels(grp_cell)

  y <- DGEList(d$counts); y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  v  <- voom(y, design)
  cr <- duplicateCorrelation(v, design, block = d$id)
  v  <- voom(y, design, block = d$id, correlation = cr$consensus)
  cr <- duplicateCorrelation(v, design, block = d$id)
  fit <- lmFit(v, design, block = d$id, correlation = cr$consensus)

  G    <- sort(unique(glab))
  cell <- function(z, g) paste(z, g, sep = "_")
  stopifnot(all(c(cell("Pericentral", G), cell("Periportal", G)) %in% colnames(design)))

  # makeContrasts uses non-standard evaluation, so the contrasts are built as
  # strings and passed through `contrasts =`
  ctr <- c(zonation = sprintf("(%s)/%d - (%s)/%d",
                              paste(cell("Periportal",  G), collapse = " + "), length(G),
                              paste(cell("Pericentral", G), collapse = " + "), length(G)))
  for (g in G) {
    ctr[paste0("zonation_", g)] <- sprintf("%s - %s", cell("Periportal", g), cell("Pericentral", g))
  }
  for (p in utils::combn(G, 2, simplify = FALSE)) {
    a <- p[1]; b <- p[2]; tag <- paste0(b, "_vs_", a)
    ctr[tag] <- sprintf("(%s + %s)/2 - (%s + %s)/2",
                        cell("Pericentral", b), cell("Periportal", b),
                        cell("Pericentral", a), cell("Periportal", a))
    ctr[paste0(tag, "_in_PC")] <- sprintf("%s - %s", cell("Pericentral", b), cell("Pericentral", a))
    ctr[paste0(tag, "_in_PP")] <- sprintf("%s - %s", cell("Periportal",  b), cell("Periportal",  a))
    ctr[paste0(tag, "_interaction")] <- sprintf("(%s - %s) - (%s - %s)",
                                                cell("Periportal", b), cell("Pericentral", b),
                                                cell("Periportal", a), cell("Pericentral", a))
  }

  cm <- makeContrasts(contrasts = ctr, levels = design)
  colnames(cm) <- names(ctr)

  f2 <- eBayes(contrasts.fit(fit, cm))
  setNames(lapply(colnames(cm),
                  function(k) topTable(f2, coef = k, number = Inf, sort.by = "none")),
           colnames(cm))
}

res_old <- fit_pb(prep(pb_old))
res_new <- fit_pb(prep(pb_new))

contrasts_shared <- intersect(names(res_old), names(res_new))
if (!setequal(names(res_old), names(res_new))) {
  warning("contrast sets differ between assignments; comparing the ",
          length(contrasts_shared), " shared contrasts")
}

cat("\n== B1. per-contrast comparison ==\n")
b1 <- lapply(contrasts_shared, function(k) {
  a <- res_old[[k]]; b <- res_new[[k]]
  g <- intersect(rownames(a), rownames(b))
  sa <- rownames(a)[a$adj.P.Val < FDR]; sb <- rownames(b)[b$adj.P.Val < FDR]
  tibble(contrast = k,
         n_tested_old = nrow(a), n_tested_new = nrow(b), n_shared_tested = length(g),
         n_sig_old = length(sa), n_sig_new = length(sb),
         n_sig_both = length(intersect(sa, sb)),
         jaccard = round(length(intersect(sa, sb)) / max(length(union(sa, sb)), 1), 3),
         r_t = round(cor(a[g, "t"], b[g, "t"]), 3),
         r_logFC = round(cor(a[g, "logFC"], b[g, "logFC"]), 3),
         sign_agreement = round(mean(sign(a[g, "t"]) == sign(b[g, "t"])), 3))
}) %>% bind_rows()
print(as.data.frame(b1))
write.csv(b1, file.path(comparison_dir, "B1_contrast_comparison.csv"), row.names = FALSE)

# per-gene t-statistics side by side, for the supplementary table
tt <- lapply(contrasts_shared, function(k) {
  a <- res_old[[k]]; b <- res_new[[k]]; g <- intersect(rownames(a), rownames(b))
  tibble(gene = g, contrast = k,
         t_original = a[g, "t"], t_alternate = b[g, "t"],
         fdr_original = a[g, "adj.P.Val"], fdr_alternate = b[g, "adj.P.Val"])
}) %>% bind_rows()
write.csv(tt, gzfile(file.path(comparison_dir, "B1_per_gene_tstats.csv.gz")), row.names = FALSE)

ggsave(file.path(comparison_dir, "B1_tstat_scatter.pdf"),
       ggplot(tt, aes(t_original, t_alternate)) +
         geom_point(size = .2, alpha = .15) +
         geom_abline(slope = 1, intercept = 0, colour = "red", linetype = 2) +
         facet_wrap(~contrast, scales = "free") + theme_bw() +
         labs(x = "t, original assignment", y = "t, alternate assignment"),
       width = 11, height = 8)

## ---- B2. zonation agreement between groups ---------------------------------
# How similar the zonation axis is from one group to the next, under each
# assignment. If the alternate assignment were adding noise rather than
# removing it, these correlations would fall.

zon_keys <- grep("^zonation_", names(res_old), value = TRUE)

cat("\n== B2. between-group correlation of the zonation axis ==\n")
b2 <- lapply(c(original = "original", alternate = "alternate"), function(nm) {
  r <- if (nm == "original") res_old else res_new
  lapply(utils::combn(zon_keys, 2, simplify = FALSE), function(p) {
    g <- intersect(rownames(r[[p[1]]]), rownames(r[[p[2]]]))
    tibble(method = nm, group_a = sub("^zonation_", "", p[1]),
           group_b = sub("^zonation_", "", p[2]),
           n_genes = length(g),
           r = round(cor(r[[p[1]]][g, "t"], r[[p[2]]][g, "t"]), 3))
  }) %>% bind_rows()
}) %>% bind_rows()
print(as.data.frame(b2))
write.csv(b2, file.path(comparison_dir, "B2_zonation_agreement.csv"), row.names = FALSE)

## ---- B2b. group-polarized zonation genes -----------------------------------
# Only has a clear reading with two groups, so it is computed only in that case.
if (length(zon_keys) == 2) {
  int_key <- grep("_interaction$", names(res_old), value = TRUE)
  stopifnot(length(int_key) == 1)
  ga <- sub("^zonation_", "", zon_keys[1]); gb <- sub("^zonation_", "", zon_keys[2])

  polarized <- function(res) {
    a <- res[[zon_keys[1]]]; b <- res[[zon_keys[2]]]; i <- res[[int_key]]
    g <- Reduce(intersect, list(rownames(a), rownames(b), rownames(i)))
    out <- tibble(gene = g,
                  t_a = a[g, "t"], fdr_a = a[g, "adj.P.Val"],
                  t_b = b[g, "t"], fdr_b = b[g, "adj.P.Val"],
                  t_int = i[g, "t"], fdr_int = i[g, "adj.P.Val"])
    keep_a <- abs(out$t_a) > POLAR_HI & abs(out$t_b) < POLAR_LO
    keep_b <- abs(out$t_b) > POLAR_HI & abs(out$t_a) < POLAR_LO
    if (POLAR_GATE) {
      keep_a <- keep_a & out$fdr_int < FDR
      keep_b <- keep_b & out$fdr_int < FDR
    }
    out$class <- NA_character_
    out$class[keep_a] <- paste(ga, ifelse(out$t_a[keep_a] > 0, "periportal", "pericentral"))
    out$class[keep_b] <- paste(gb, ifelse(out$t_b[keep_b] > 0, "periportal", "pericentral"))
    out[!is.na(out$class), ]
  }
  pol_old <- polarized(res_old); pol_new <- polarized(res_new)

  cat("\n== B2b. group-polarized zonation genes ==\n")
  shared_counts <- inner_join(pol_old[, c("gene", "class")], pol_new[, c("gene", "class")],
                              by = c("gene", "class")) %>% count(class, name = "shared")
  b2b <- full_join(count(pol_old, class, name = "original"),
                   count(pol_new, class, name = "alternate"), by = "class") %>%
    full_join(shared_counts, by = "class") %>%
    mutate(across(c(original, alternate, shared), ~tidyr::replace_na(.x, 0L)))
  print(as.data.frame(b2b))
  cat(sprintf("total: original %d, alternate %d, shared %d\n",
              nrow(pol_old), nrow(pol_new), length(intersect(pol_old$gene, pol_new$gene))))
  write.csv(b2b, file.path(comparison_dir, "B2b_polarized_counts.csv"), row.names = FALSE)

  rename_pol <- function(d, suffix) {
    setNames(d, c("gene", paste0(c(paste0("t_", ga), paste0("fdr_", ga),
                                   paste0("t_", gb), paste0("fdr_", gb),
                                   "t_int", "fdr_int", "class"), "_", suffix)))
  }
  pol_genes <- full_join(rename_pol(pol_old, "original"),
                         rename_pol(pol_new, "alternate"), by = "gene") %>%
    mutate(status = dplyr::case_when(
      !is.na(class_original) & !is.na(class_alternate) ~ "both",
      !is.na(class_original)                           ~ "original only",
      TRUE                                             ~ "alternate only"))
  write.csv(pol_genes, file.path(comparison_dir, "B2b_polarized_genes.csv"), row.names = FALSE)
} else {
  message("more than two groups - the polarized-gene classification is skipped")
}

## ---- B3. overlap of group-biased and zonation-biased genes -----------------

cat("\n== B3. group-biased genes that are also zonation-biased ==\n")
grp_keys <- grep("_vs_", names(res_old), value = TRUE)
grp_keys <- grp_keys[!grepl("_(in_PC|in_PP|interaction)$", grp_keys)]

b3 <- lapply(c(original = "original", alternate = "alternate"), function(nm) {
  r <- if (nm == "original") res_old else res_new
  g <- rownames(r$zonation)
  zonb <- g[r$zonation[g, "adj.P.Val"] < FDR]
  lapply(grp_keys, function(k) {
    grpb <- g[r[[k]][g, "adj.P.Val"] < FDR]
    tibble(method = nm, contrast = k,
           n_group_biased = length(grpb), n_zonation_biased = length(zonb),
           n_also_zonation_biased = length(intersect(grpb, zonb)),
           pct = round(100 * length(intersect(grpb, zonb)) / max(length(grpb), 1), 1))
  }) %>% bind_rows()
}) %>% bind_rows()
print(as.data.frame(b3))
write.csv(b3, file.path(comparison_dir, "B3_group_zonation_overlap.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(comparison_dir, "sessionInfo.txt"))
message("\nComparison outputs written to: ", comparison_dir)
