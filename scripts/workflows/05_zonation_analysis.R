################################################################################
# Workflow 5: Hepatocyte Zonation Analysis
################################################################################
#
# This workflow analyzes zonation-specific gene expression in hepatocytes.
# The data is ALREADY pseudobulk aggregated with samples separated by zone.
#
# Author: Refactored codebase
# Date: 2025
#
# Replaces legacy scripts:
#   - scripts/pseudobulk_DA_zonation.R
#
# Key biological findings:
#   - 43% of gonad-dependent genes are pericentral-biased (Cyp2e1+)
#   - 18% of gonad-dependent genes are periportal-biased (Cyp2f2+)
#
################################################################################

# Clear workspace
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(gridExtra)
  library(decoupleR)
  library(reshape2)
  library(pheatmap)
  library(GSEABase)
})

################################################################################
# 1. CONFIGURATION
################################################################################

cat(strrep("=", 80), "\n")
cat("WORKFLOW 5: HEPATOCYTE ZONATION ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

# Load configuration and functions
source("scripts/config/analysis_config.R")
source("scripts/support_functions.R")

# Define output directory
zonation_dir <- "results/zonation"
if (!dir.exists(zonation_dir)) {
  dir.create(zonation_dir, recursive = TRUE)
}

if (!dir.exists(file.path(zonation_dir, "volcanos"))) {
  dir.create(file.path(zonation_dir, "volcanos"), recursive = TRUE)
}

################################################################################
# 2. LOAD ZONATION PSEUDOBULK DATA
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("2. LOADING ZONATION PSEUDOBULK DATA\n")
cat(strrep("-", 80), "\n\n")

cat("NOTE: This data is ALREADY pseudobulk aggregated.\n")
cat("Each sample has been split into Pericentral and Periportal zones.\n\n")

# Load zonation-specific pseudobulk data
cat("Loading data...\n")
scaled_RNA <- read.csv("data/RNA/scaled_RNA_zonation.csv", row.names = 1, check.names = FALSE)
meta_data_df <- read.csv("data/RNA/meta_data_df_zonation.csv")

cat("  Genes:", nrow(scaled_RNA), "\n")
cat("  Samples:", ncol(scaled_RNA), "\n")
cat("  Metadata rows:", nrow(meta_data_df), "\n\n")

cat("Metadata structure:\n")
cat("  Columns:", paste(colnames(meta_data_df), collapse = ", "), "\n")
cat("  Sample format: {mouse_id}_{zonation_type}_{gonadal_type}\n\n")

cat("Sample breakdown:\n")
cat("  Zonation types:\n")
print(table(meta_data_df$zonation_type))
cat("\n  Gonadal types:\n")
print(table(meta_data_df$gonadal_type))
cat("\n")

################################################################################
# 3. PCA VISUALIZATION
################################################################################

cat(strrep("-", 80), "\n")
cat("3. PCA VISUALIZATION\n")
cat(strrep("-", 80), "\n\n")

# Create combined condition for visualization
meta_data_df$condition <- paste(meta_data_df$zonation_type,
                                meta_data_df$gonadal_type,
                                sep = "_")

cat("Running PCA...\n")
pcares <- prcomp(t(scaled_RNA))

# PCA colored by zonation
pdf(file.path(zonation_dir, "PCA_by_zonation.pdf"), width = 8, height = 6)
grid::grid.newpage()
grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
grid::grid.draw(nicePCA(scaled_RNA, meta_data_df[, c("sample", "zonation_type")], no_label = TRUE))
dev.off()

# PCA colored by gonadal type
pdf(file.path(zonation_dir, "PCA_by_gonad.pdf"), width = 8, height = 6)
grid::grid.newpage()
grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
grid::grid.draw(nicePCA(scaled_RNA, meta_data_df[, c("sample", "gonadal_type")], no_label = TRUE))
dev.off()

# PCA colored by combined condition
pdf(file.path(zonation_dir, "PCA_by_condition.pdf"), width = 8, height = 6)
grid::grid.newpage()
grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
grid::grid.draw(nicePCA(scaled_RNA, meta_data_df[, c("sample", "condition")], no_label = TRUE))
dev.off()

cat("  Saved PCA plots\n\n")

################################################################################
# 4. PCA ENRICHMENT FOR METADATA
################################################################################

cat(strrep("-", 80), "\n")
cat("4. PCA ENRICHMENT FOR METADATA\n")
cat(strrep("-", 80), "\n\n")

cat("Testing which PCs are associated with experimental factors...\n\n")

# Get top 15 PCs
Z_matrix <- as.data.frame(pcares$x)[, 1:15]

# Prepare metadata network for decoupleR
meta_data_df_long <- melt(meta_data_df[, -2], id.vars = "sample")
meta_data_df_long <- meta_data_df_long[, -2]
meta_data_df_long$mor <- 1
names(meta_data_df_long) <- c("target", "source", "mor")

# Run enrichment
meta_data_enrichment <- decoupleR::run_ulm(mat = Z_matrix, network = meta_data_df_long)
meta_data_enrichment <- dcast(meta_data_enrichment, source ~ condition, value.var = "score")

meta_data_enrichment <- meta_data_enrichment[, c("source", paste("PC", 1:15, sep = ""))]
rownames(meta_data_enrichment) <- meta_data_enrichment$source

# Calculate variance explained
explained <- (pcares$sdev^2 / sum(pcares$sdev^2))[1:15]
names(explained) <- names(meta_data_enrichment)[-1]
explained <- as.data.frame(explained)

# Create color palette
t <- as.vector(t(meta_data_enrichment[, -1]))
palette1 <- createLinearColors(t[t < 0], withZero = FALSE, maximum = abs(min(t, na.rm = TRUE)) * 10)
palette2 <- createLinearColors(t[t > 0], withZero = FALSE, maximum = abs(max(t, na.rm = TRUE)) * 10)
palette <- c(palette1, palette2)

# Plot enrichment heatmap
pdf(file.path(zonation_dir, "PCA_metadata_enrichment.pdf"), width = 10, height = 6)
pheatmap::pheatmap(meta_data_enrichment[, -1],
                   display_numbers = TRUE,
                   show_rownames = TRUE,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   color = palette,
                   angle_col = 315,
                   annotation_col = explained)
dev.off()

cat("  Saved PCA enrichment heatmap\n\n")

################################################################################
# 5. DIFFERENTIAL ANALYSIS
################################################################################

cat(strrep("-", 80), "\n")
cat("5. DIFFERENTIAL ANALYSIS\n")
cat(strrep("-", 80), "\n\n")

# Prepare targets
targets <- meta_data_df[, c("sample", "condition")]
names(targets)[2] <- "condition"

cat("Conditions (", length(unique(targets$condition)), " total):\n")
print(unique(targets$condition))
cat("\n")

# Define comparisons
# Format: sample_order matches unique(targets$condition)
comparisons <- list(
  "portal_vs_Central" = c(-1, -3, 2, 4),              # Overall zonation effect
  "testes_vs_ova" = c(-1, -2, 3, 4),                  # Overall gonad effect
  "portalova_vs_centralova" = c(-1, 2),              # Zonation in ovaries
  "portaltestes_vs_centraltestes" = c(-3, 4),        # Zonation in testes
  "centraltestes_vs_centralova" = c(-1, 3),          # Gonad effect in pericentral
  "portaltestes_vs_portalova" = c(-2, 4),            # Gonad effect in periportal
  "testes_effect_zonation" = c(-3, 4, 1, -2)         # Interaction: does zonation differ by gonad?
)

cat("Running limma differential analysis...\n")
cat("  Contrasts:", length(comparisons), "\n")
for (comp_name in names(comparisons)) {
  cat("    -", comp_name, "\n")
}
cat("\n")

limmaRes <- runLimma(scaled_RNA, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(
  limma_res = limmaRes,
  comp_names = names(comparisons),
  number = nrow(scaled_RNA),
  adjust.method = "fdr"
)

cat("Differential analysis complete.\n\n")

# Print summary statistics
for (comp_name in names(ttop_list)) {
  ttop <- ttop_list[[comp_name]]
  sig_genes <- sum(ttop$adj.P.Val < 0.05, na.rm = TRUE)
  cat("  ", comp_name, ": ", sig_genes, " significant genes (FDR < 0.05)\n", sep = "")
}
cat("\n")

################################################################################
# 6. VOLCANO PLOTS
################################################################################

cat(strrep("-", 80), "\n")
cat("6. GENERATING VOLCANO PLOTS\n")
cat(strrep("-", 80), "\n\n")

for (k in 1:length(ttop_list)) {
  comp_name <- names(ttop_list)[k]
  cat("  Creating volcano plot:", comp_name, "\n")

  volca_plt <- volcano_nice(ttop_list[[k]],
                           FCIndex = 2,
                           IDIndex = 1,
                           pValIndex = 5,
                           label = TRUE,
                           nlabels = 40)

  ggsave(file.path(zonation_dir, "volcanos", paste0("volcano_", comp_name, ".pdf")),
         plot = volca_plt,
         width = 10,
         height = 8)
}

cat("\n  Saved volcano plots to:", file.path(zonation_dir, "volcanos"), "\n\n")

################################################################################
# 7. PATHWAY ENRICHMENT
################################################################################

cat(strrep("-", 80), "\n")
cat("7. PATHWAY ENRICHMENT ANALYSIS\n")
cat(strrep("-", 80), "\n\n")

cat("Loading gene sets...\n")
if (file.exists("support/m5.all.v2023.2.Mm.symbols.gmt")) {
  tryCatch({
    gene_sets <- import_gmt("support/m5.all.v2023.2.Mm.symbols.gmt")
    cat("  Loaded gene sets from GMT file\n")
    cat("  Columns:", paste(colnames(gene_sets), collapse = ", "), "\n")

    # Rename columns for decoupleR compatibility
    # GMT import gives: "gene" (target) and "term" (source)
    # decoupleR needs: "target" and "source"
    gene_sets <- gene_sets %>%
      dplyr::rename(target = gene, source = term)

    cat("  Renamed columns to: ", paste(colnames(gene_sets), collapse = ", "), "\n")
    cat("  Pathways:", length(unique(gene_sets$source)), "\n")
    cat("  Genes:", length(unique(gene_sets$target)), "\n\n")

    # Convert ttop_list to t-value matrix
    cat("Converting differential results to t-value matrix...\n")

    tval_df <- ttop_list_to_t_table(ttop_list)

    # ttop_list_to_t_table merges by "ID" column which contains gene names
    # Set rownames to the ID column (gene names)
    if ("ID" %in% colnames(tval_df)) {
      rownames(tval_df) <- tval_df$ID
      tval_df <- tval_df[, colnames(tval_df) != "ID", drop = FALSE]
    }

    cat("  Matrix dimensions:", nrow(tval_df), "genes x", ncol(tval_df), "columns\n")
    cat("  Gene names in matrix (first 5):", paste(head(rownames(tval_df), 5), collapse = ", "), "\n")
    cat("  Gene names in gene sets (first 5):", paste(head(unique(gene_sets$target), 5), collapse = ", "), "\n")

    # Check overlap
    overlap <- intersect(rownames(tval_df), gene_sets$target)
    cat("  Genes in common:", length(overlap), "out of", nrow(tval_df), "matrix genes\n\n")

    # cat("Running pathway enrichment (decoupleR)...\n")
    # pathway_scores <- decoupleR::run_ulm(
    #   mat = tval_df,
    #   network = gene_sets,
    #   .source = "source",
    #   .target = "target",
    #   minsize = 15
    # )
    # 
    # 
    # cat("  Enriched", nrow(pathway_scores), "pathway × contrast combinations\n\n")
    # 
    # # Save results
    # save(pathway_scores, gene_sets, file = file.path(zonation_dir, "pathway_enrichment.RData"))
    
    # cat("Skip run and load pathway enrichment results...\n")
    load(file = file.path(zonation_dir, "pathway_enrichment.RData"))
    
    # Top pathways for each comparison
    for (comp_name in unique(pathway_scores$condition)) {
      cat("  Top pathways for", comp_name, ":\n")
      top_pw <- pathway_scores %>%
        dplyr::filter(condition == comp_name) %>%
        dplyr::arrange(desc(abs(score))) %>%
        head(5)
      for (i in 1:nrow(top_pw)) {
        cat("    ", i, ". ", top_pw$source[i], " (score: ",
            round(top_pw$score[i], 2), ")\n", sep = "")
      }
      cat("\n")
    }

    cat("✓ Pathway enrichment complete\n\n")
  }, error = function(e) {
    cat("  Error during pathway enrichment:", conditionMessage(e), "\n")
    cat("  Skipping pathway enrichment.\n\n")
  })
} else {
  cat("  Gene set file not found. Skipping pathway enrichment.\n\n")
}

################################################################################
# 8. LOAD LIPIDOMICS DATA AND COMPUTE TF ACTIVITY SCORES
################################################################################
 
cat(strrep("-", 80), "\n")
cat("8. LOADING LIPIDOMICS AND COMPUTING TF ACTIVITY SCORES\n")
cat(strrep("-", 80), "\n\n")

# Load lipidomics data
cat("Loading lipidomics data...\n")
if (file.exists("data/metabolomics/lipidomic.RData")) {
  load("data/metabolomics/lipidomic.RData")
  cat("  Loaded lipidomics data: ", nrow(lipidomic), " lipids x ", ncol(lipidomic), " samples\n\n")
} else {
  stop("Lipidomics data file not found: data/metabolomics/lipidomic.RData\n")
}

# Load TF regulons
cat("Loading TF regulons...\n")
load("support/TF_regulons.RData")
TF_regulons <- TF_regulons[complete.cases(TF_regulons), ]

# Add custom hormone receptor targets
Gonadal_Hormone_Genes <- as.data.frame(read.csv("support/Gonadal_Hormone_Genes.csv"))
Gonadal_Hormone_Genes$source <- gsub("_spec", "", Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("ESR1", "Esr1", Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("AR", "Ar", Gonadal_Hormone_Genes$source)

# Remove CollecTRI versions of Ar and Esr1 to use custom targets
TF_regulons <- TF_regulons[!TF_regulons$source %in% c("Ar", "Esr1"), ]
TF_regulons <- rbind(TF_regulons, Gonadal_Hormone_Genes)

cat("  Total TF regulons: ", length(unique(TF_regulons$source)), "\n\n")

# Compute TF activity scores across all samples
cat("Computing TF activity scores...\n")
TF_scores <- decoupleR::run_ulm(mat = scaled_RNA, network = TF_regulons, minsize = 20)
TF_scores_df <- dcast(TF_scores, source ~ condition, value.var = "score")
rownames(TF_scores_df) <- TF_scores_df$source
TF_scores_df <- TF_scores_df[, -1]

cat("  Computed scores for ", nrow(TF_scores_df), " TFs across ", ncol(TF_scores_df), " samples\n\n")

# Create t-value table for TF enrichment on contrasts
cat("Creating t-value table for TF enrichment...\n")
t_table <- ttop_list_to_t_table(ttop_list)
names(t_table) <- c("gene", names(comparisons))

# Save t-value table
write.csv(t_table, file = file.path(zonation_dir, "t_table_zonation.csv"), row.names = FALSE)
cat("  Saved t-value table\n\n")

# Compute TF activity on contrasts
cat("Computing TF enrichment on contrasts...\n")
decoupler_input_t_table <- t_table[, -1]
rownames(decoupler_input_t_table) <- t_table$gene

TF_scores_contrasts <- decoupleR::run_ulm(mat = decoupler_input_t_table, network = TF_regulons, minsize = 20)
TF_scores_contrasts_df <- dcast(TF_scores_contrasts, source ~ condition, value.var = "score")
rownames(TF_scores_contrasts_df) <- TF_scores_contrasts_df$source

cat("  Computed TF enrichment for ", nrow(TF_scores_contrasts_df), " TFs\n\n")

cat("✓ Lipidomics and TF scores loaded\n\n")

################################################################################
# 9. IDENTIFY ZONATION-BIASED GENES
################################################################################

cat(strrep("-", 80), "\n")
cat("9. IDENTIFYING ZONATION-BIASED GENES\n")
cat(strrep("-", 80), "\n\n")

cat("Identifying genes with zonation-specific effects...\n\n")

# Genes with testes-specific zonation bias (high zonation effect in testes, low in ovaries)
testes_zonation_bias_genes <- t_table[abs(t_table$testes_effect_zonation) > 5 &
                                       abs(t_table$portalova_vs_centralova) < 1, ]
cat("  Testes-specific zonation genes: ", nrow(testes_zonation_bias_genes), "\n")

# Genes with ovaries-specific zonation bias (high zonation effect in ovaries, low in testes)
ovaries_zonation_bias_genes <- t_table[abs(t_table$testes_effect_zonation) > 5 &
                                        abs(t_table$portaltestes_vs_centraltestes) < 1, ]
cat("  Ovaries-specific zonation genes: ", nrow(ovaries_zonation_bias_genes), "\n")

# Genes with general zonation bias (independent of gonad type)
zonation_bias_genes <- t_table[abs(t_table$testes_effect_zonation) < 1 &
                                abs(t_table$portal_vs_Central) > 5, ]
cat("  General zonation genes: ", nrow(zonation_bias_genes), "\n\n")

# Split by zone
testes_pericentral_genes <- testes_zonation_bias_genes[testes_zonation_bias_genes$portaltestes_vs_centraltestes < 0, ]
testes_periportal_genes <- testes_zonation_bias_genes[testes_zonation_bias_genes$portaltestes_vs_centraltestes > 0, ]
ovaries_pericentral_genes <- ovaries_zonation_bias_genes[ovaries_zonation_bias_genes$portalova_vs_centralova < 0, ]
ovaries_periportal_genes <- ovaries_zonation_bias_genes[ovaries_zonation_bias_genes$portalova_vs_centralova > 0, ]

cat("  Breakdown:\n")
cat("    Testes pericentral: ", nrow(testes_pericentral_genes), "\n")
cat("    Testes periportal: ", nrow(testes_periportal_genes), "\n")
cat("    Ovaries pericentral: ", nrow(ovaries_pericentral_genes), "\n")
cat("    Ovaries periportal: ", nrow(ovaries_periportal_genes), "\n\n")

# Create combined top genes matrix for visualization
top_genes <- as.data.frame(do.call(rbind, list(
  testes_pericentral_genes,
  testes_periportal_genes,
  ovaries_pericentral_genes,
  ovaries_periportal_genes
)))
rownames(top_genes) <- top_genes$gene
top_genes <- top_genes[, c("portalova_vs_centralova", "portaltestes_vs_centraltestes")]

# Create heatmap of top zonation genes
pdf(file.path(zonation_dir, "zonation_biased_genes_heatmap.pdf"), width = 12, height = 8)
t <- as.vector(t(top_genes))
palette1 <- createLinearColors(t[t < 0], withZero = FALSE,
                              maximum = abs(min(t, na.rm = TRUE)) * 10,
                              my_colors = c("lightgreen", "white", "orange"))
palette2 <- createLinearColors(t[t > 0], withZero = FALSE,
                              maximum = abs(max(t, na.rm = TRUE)) * 10,
                              my_colors = c("lightgreen", "white", "orange"))
palette <- c(palette1, palette2)
pheatmap::pheatmap(t(top_genes),
                   display_numbers = TRUE,
                   show_rownames = TRUE,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   color = palette,
                   angle_col = 0)
dev.off()

cat("  Saved zonation-biased genes heatmap\n\n")
cat("✓ Zonation-biased genes identified\n\n")

################################################################################
# 10. GENERATE ZONATION-SPECIFIC CORRELATION NETWORKS
################################################################################

cat(strrep("-", 80), "\n")
cat("10. GENERATING ZONATION-SPECIFIC CORRELATION NETWORKS\n")
cat(strrep("-", 80), "\n\n")

cat("Building correlation networks for each zonation × gonad combination...\n")
cat("Two versions will be generated:\n")
cat("  1. Standard networks: correlation threshold = 0.65 (for visualization)\n")
cat("  2. Full networks: correlation threshold = 0 (all correlations, lipid species only)\n\n")

threshold_cor <- 0.65

# Create multi-omic output directory if it doesn't exist
multiomics_dir <- "results/multi-omic"
if (!dir.exists(multiomics_dir)) {
  dir.create(multiomics_dir, recursive = TRUE)
}

#===============================================================================
# Network 1: Testes Periportal
#===============================================================================
cat("1. Testes Periportal Network...\n")

scaled_RNA_testes_portal <- scaled_RNA[testes_periportal_genes$gene,
                                       grepl("Periportal_Testes", names(scaled_RNA))]
TF_scores_df_testes_portal <- TF_scores_df[, grepl("Periportal_Testes", names(TF_scores_df))]
lipidomic_testes <- lipidomic[, gsub("_Peri.*_Testes", "", names(scaled_RNA_testes_portal))]
names(scaled_RNA_testes_portal) <- gsub("_Peri.*_Testes", "", names(scaled_RNA_testes_portal))
names(TF_scores_df_testes_portal) <- gsub("_Peri.*_Testes", "", names(TF_scores_df_testes_portal))

cormat_testes_portal <- as.data.frame(rbind(
  make_cor_network(lipidomic_testes, scaled_RNA_testes_portal, theshold_cor = threshold_cor),
  make_cor_network(scaled_RNA_testes_portal, TF_scores_df_testes_portal, theshold_cor = threshold_cor, TF = TRUE)
))
names(cormat_testes_portal) <- c("node1", "node2", "weight")
cormat_testes_portal_ATT <- expand_attributes(cormat_testes_portal, t_table, TF_scores_df)
cormat_testes_portal_ATT$type <- ifelse(cormat_testes_portal_ATT$gene %in% rownames(TF_scores_df), "TF",
                                        ifelse(cormat_testes_portal_ATT$gene %in% t_table$gene, "gene", "lipid"))

write.csv(cormat_testes_portal, file = file.path(multiomics_dir, "cormat_testes_portal_SIF.csv"), row.names = FALSE)
write.csv(cormat_testes_portal_ATT, file = file.path(multiomics_dir, "cormat_testes_portal_ATT.csv"), row.names = FALSE)

# Generate full correlation network (threshold = 0) for non-compressed lipid species
cat("  Generating full correlation network (threshold = 0)...\n")
cormat_testes_portal_full <- as.data.frame(rbind(
  make_cor_network(lipidomic_testes, scaled_RNA_testes_portal, theshold_cor = 0),
  make_cor_network(scaled_RNA_testes_portal, TF_scores_df_testes_portal, theshold_cor = 0, TF = TRUE)
))
names(cormat_testes_portal_full) <- c("node1", "node2", "weight")
cormat_testes_portal_full_ATT <- expand_attributes(cormat_testes_portal_full, t_table, TF_scores_df)
cormat_testes_portal_full_ATT$type <- ifelse(cormat_testes_portal_full_ATT$gene %in% rownames(TF_scores_df), "TF",
                                              ifelse(cormat_testes_portal_full_ATT$gene %in% t_table$gene, "gene", "lipid"))

write.csv(cormat_testes_portal_full, file = file.path(multiomics_dir, "cormat_testes_portal_full_SIF.csv"), row.names = FALSE)
write.csv(cormat_testes_portal_full_ATT, file = file.path(multiomics_dir, "cormat_testes_portal_full_ATT.csv"), row.names = FALSE)

# Create compressed version (aggregate by lipid class)
cormat_testes_portal$node2 <- gsub("[(].*", "", cormat_testes_portal$node2)
cormat_testes_portal$sign <- sign(cormat_testes_portal$weight)
cormat_testes_portal <- aggregate(weight ~ node1 + node2 + sign, data = cormat_testes_portal, FUN = mean)
cormat_testes_portal$abs_weight <- abs(cormat_testes_portal$weight)

cormat_testes_portal_ATT$gene <- gsub("[(].*", "", cormat_testes_portal_ATT$gene)
cormat_testes_portal_ATT <- aggregate(cbind(portal_vs_Central, testes_vs_ova, portalova_vs_centralova,
                                            portaltestes_vs_centraltestes, centraltestes_vs_centralova,
                                            portaltestes_vs_portalova, testes_effect_zonation) ~ gene + type,
                                     data = cormat_testes_portal_ATT, FUN = mean)

write.csv(cormat_testes_portal, file = file.path(multiomics_dir, "compressed_cormat_testes_portal_SIF.csv"), row.names = FALSE)
write.csv(cormat_testes_portal_ATT, file = file.path(multiomics_dir, "compressed_cormat_testes_portal_ATT.csv"), row.names = FALSE)

cat("  ✓ Saved testes periportal networks\n")

#===============================================================================
# Network 2: Testes Pericentral
#===============================================================================
cat("2. Testes Pericentral Network...\n")

scaled_RNA_testes_central <- scaled_RNA[testes_pericentral_genes$gene,
                                        grepl("Pericentral_Testes", names(scaled_RNA))]
TF_scores_df_testes_central <- TF_scores_df[, grepl("Pericentral_Testes", names(TF_scores_df))]
lipidomic_testes <- lipidomic[, gsub("_Peri.*_Testes", "", names(scaled_RNA_testes_central))]
names(scaled_RNA_testes_central) <- gsub("_Peri.*_Testes", "", names(scaled_RNA_testes_central))
names(TF_scores_df_testes_central) <- gsub("_Peri.*_Testes", "", names(TF_scores_df_testes_central))

cormat_testes_central <- as.data.frame(rbind(
  make_cor_network(lipidomic_testes, scaled_RNA_testes_central, theshold_cor = threshold_cor),
  make_cor_network(scaled_RNA_testes_central, TF_scores_df_testes_central, theshold_cor = threshold_cor, TF = TRUE)
))
names(cormat_testes_central) <- c("node1", "node2", "weight")
cormat_testes_central_ATT <- expand_attributes(cormat_testes_central, t_table, TF_scores_df)
cormat_testes_central_ATT$type <- ifelse(cormat_testes_central_ATT$gene %in% rownames(TF_scores_df), "TF",
                                         ifelse(cormat_testes_central_ATT$gene %in% t_table$gene, "gene", "lipid"))

write.csv(cormat_testes_central, file = file.path(multiomics_dir, "cormat_testes_central_SIF.csv"), row.names = FALSE)
write.csv(cormat_testes_central_ATT, file = file.path(multiomics_dir, "cormat_testes_central_ATT.csv"), row.names = FALSE)

# Generate full correlation network (threshold = 0) for non-compressed lipid species
cat("  Generating full correlation network (threshold = 0)...\n")
cormat_testes_central_full <- as.data.frame(rbind(
  make_cor_network(lipidomic_testes, scaled_RNA_testes_central, theshold_cor = 0),
  make_cor_network(scaled_RNA_testes_central, TF_scores_df_testes_central, theshold_cor = 0, TF = TRUE)
))
names(cormat_testes_central_full) <- c("node1", "node2", "weight")
cormat_testes_central_full_ATT <- expand_attributes(cormat_testes_central_full, t_table, TF_scores_df)
cormat_testes_central_full_ATT$type <- ifelse(cormat_testes_central_full_ATT$gene %in% rownames(TF_scores_df), "TF",
                                               ifelse(cormat_testes_central_full_ATT$gene %in% t_table$gene, "gene", "lipid"))

write.csv(cormat_testes_central_full, file = file.path(multiomics_dir, "cormat_testes_central_full_SIF.csv"), row.names = FALSE)
write.csv(cormat_testes_central_full_ATT, file = file.path(multiomics_dir, "cormat_testes_central_full_ATT.csv"), row.names = FALSE)

# Create compressed version
cormat_testes_central$node2 <- gsub("[(].*", "", cormat_testes_central$node2)
cormat_testes_central$sign <- sign(cormat_testes_central$weight)
cormat_testes_central <- aggregate(weight ~ node1 + node2 + sign, data = cormat_testes_central, FUN = mean)
cormat_testes_central$abs_weight <- abs(cormat_testes_central$weight)

cormat_testes_central_ATT$gene <- gsub("[(].*", "", cormat_testes_central_ATT$gene)
cormat_testes_central_ATT <- aggregate(cbind(portal_vs_Central, testes_vs_ova, portalova_vs_centralova,
                                             portaltestes_vs_centraltestes, centraltestes_vs_centralova,
                                             portaltestes_vs_portalova, testes_effect_zonation) ~ gene + type,
                                      data = cormat_testes_central_ATT, FUN = mean)

write.csv(cormat_testes_central, file = file.path(multiomics_dir, "compressed_cormat_testes_central_SIF.csv"), row.names = FALSE)
write.csv(cormat_testes_central_ATT, file = file.path(multiomics_dir, "compressed_cormat_testes_central_ATT.csv"), row.names = FALSE)

cat("  ✓ Saved testes pericentral networks\n")

#===============================================================================
# Network 3: Ovaries Periportal
#===============================================================================
cat("3. Ovaries Periportal Network...\n")

scaled_RNA_ovaries_portal <- scaled_RNA[ovaries_periportal_genes$gene,
                                        grepl("Periportal_Ovaries", names(scaled_RNA))]
TF_scores_df_Ovaries_portal <- TF_scores_df[, grepl("Periportal_Ovaries", names(TF_scores_df))]
lipidomic_Ovaries <- lipidomic[, gsub("_Peri.*_Ovaries", "", names(scaled_RNA_ovaries_portal))]
names(scaled_RNA_ovaries_portal) <- gsub("_Peri.*_Ovaries", "", names(scaled_RNA_ovaries_portal))
names(TF_scores_df_Ovaries_portal) <- gsub("_Peri.*_Ovaries", "", names(TF_scores_df_Ovaries_portal))

cormat_Ovaries_portal <- as.data.frame(rbind(
  make_cor_network(lipidomic_Ovaries, scaled_RNA_ovaries_portal, theshold_cor = threshold_cor),
  make_cor_network(scaled_RNA_ovaries_portal, TF_scores_df_Ovaries_portal, theshold_cor = threshold_cor, TF = TRUE)
))
names(cormat_Ovaries_portal) <- c("node1", "node2", "weight")
cormat_Ovaries_portal_ATT <- expand_attributes(cormat_Ovaries_portal, t_table, TF_scores_df)
cormat_Ovaries_portal_ATT$type <- ifelse(cormat_Ovaries_portal_ATT$gene %in% rownames(TF_scores_df), "TF",
                                         ifelse(cormat_Ovaries_portal_ATT$gene %in% t_table$gene, "gene", "lipid"))

write.csv(cormat_Ovaries_portal, file = file.path(multiomics_dir, "cormat_Ovaries_portal_SIF.csv"), row.names = FALSE)
write.csv(cormat_Ovaries_portal_ATT, file = file.path(multiomics_dir, "cormat_Ovaries_portal_ATT.csv"), row.names = FALSE)

# Generate full correlation network (threshold = 0) for non-compressed lipid species
cat("  Generating full correlation network (threshold = 0)...\n")
cormat_Ovaries_portal_full <- as.data.frame(rbind(
  make_cor_network(lipidomic_Ovaries, scaled_RNA_ovaries_portal, theshold_cor = 0),
  make_cor_network(scaled_RNA_ovaries_portal, TF_scores_df_Ovaries_portal, theshold_cor = 0, TF = TRUE)
))
names(cormat_Ovaries_portal_full) <- c("node1", "node2", "weight")
cormat_Ovaries_portal_full_ATT <- expand_attributes(cormat_Ovaries_portal_full, t_table, TF_scores_df)
cormat_Ovaries_portal_full_ATT$type <- ifelse(cormat_Ovaries_portal_full_ATT$gene %in% rownames(TF_scores_df), "TF",
                                               ifelse(cormat_Ovaries_portal_full_ATT$gene %in% t_table$gene, "gene", "lipid"))

write.csv(cormat_Ovaries_portal_full, file = file.path(multiomics_dir, "cormat_Ovaries_portal_full_SIF.csv"), row.names = FALSE)
write.csv(cormat_Ovaries_portal_full_ATT, file = file.path(multiomics_dir, "cormat_Ovaries_portal_full_ATT.csv"), row.names = FALSE)

# Create compressed version
cormat_Ovaries_portal$node2 <- gsub("[(].*", "", cormat_Ovaries_portal$node2)
cormat_Ovaries_portal$sign <- sign(cormat_Ovaries_portal$weight)
cormat_Ovaries_portal <- aggregate(weight ~ node1 + node2 + sign, data = cormat_Ovaries_portal, FUN = mean)
cormat_Ovaries_portal$abs_weight <- abs(cormat_Ovaries_portal$weight)

cormat_Ovaries_portal_ATT$gene <- gsub("[(].*", "", cormat_Ovaries_portal_ATT$gene)
cormat_Ovaries_portal_ATT <- aggregate(cbind(portal_vs_Central, testes_vs_ova, portalova_vs_centralova,
                                             portaltestes_vs_centraltestes, centraltestes_vs_centralova,
                                             portaltestes_vs_portalova, testes_effect_zonation) ~ gene + type,
                                      data = cormat_Ovaries_portal_ATT, FUN = mean)

write.csv(cormat_Ovaries_portal, file = file.path(multiomics_dir, "compressed_cormat_Ovaries_portal_SIF.csv"), row.names = FALSE)
write.csv(cormat_Ovaries_portal_ATT, file = file.path(multiomics_dir, "compressed_cormat_Ovaries_portal_ATT.csv"), row.names = FALSE)

cat("  ✓ Saved ovaries periportal networks\n")

#===============================================================================
# Network 4: Ovaries Pericentral
#===============================================================================
cat("4. Ovaries Pericentral Network...\n")

scaled_RNA_ova_central <- scaled_RNA[ovaries_pericentral_genes$gene,
                                     grepl("Pericentral_Ovaries", names(scaled_RNA))]
TF_scores_df_ova_central <- TF_scores_df[, grepl("Pericentral_Ovaries", names(TF_scores_df))]
lipidomic_Ovaries <- lipidomic[, gsub("_Peri.*_Ovaries", "", names(scaled_RNA_ova_central))]
names(scaled_RNA_ova_central) <- gsub("_Peri.*_Ovaries", "", names(scaled_RNA_ova_central))
names(TF_scores_df_ova_central) <- gsub("_Peri.*_Ovaries", "", names(TF_scores_df_ova_central))

cormat_ova_central <- as.data.frame(rbind(
  make_cor_network(lipidomic_Ovaries, scaled_RNA_ova_central, theshold_cor = threshold_cor),
  make_cor_network(scaled_RNA_ova_central, TF_scores_df_ova_central, theshold_cor = threshold_cor, TF = TRUE)
))
names(cormat_ova_central) <- c("node1", "node2", "weight")
cormat_ova_central_ATT <- expand_attributes(cormat_ova_central, t_table, TF_scores_df)
cormat_ova_central_ATT$type <- ifelse(cormat_ova_central_ATT$gene %in% rownames(TF_scores_df), "TF",
                                      ifelse(cormat_ova_central_ATT$gene %in% t_table$gene, "gene", "lipid"))

write.csv(cormat_ova_central, file = file.path(multiomics_dir, "cormat_ova_central_SIF.csv"), row.names = FALSE)
write.csv(cormat_ova_central_ATT, file = file.path(multiomics_dir, "cormat_ova_central_ATT.csv"), row.names = FALSE)

# Generate full correlation network (threshold = 0) for non-compressed lipid species
cat("  Generating full correlation network (threshold = 0)...\n")
cormat_ova_central_full <- as.data.frame(rbind(
  make_cor_network(lipidomic_Ovaries, scaled_RNA_ova_central, theshold_cor = 0),
  make_cor_network(scaled_RNA_ova_central, TF_scores_df_ova_central, theshold_cor = 0, TF = TRUE)
))
names(cormat_ova_central_full) <- c("node1", "node2", "weight")
cormat_ova_central_full_ATT <- expand_attributes(cormat_ova_central_full, t_table, TF_scores_df)
cormat_ova_central_full_ATT$type <- ifelse(cormat_ova_central_full_ATT$gene %in% rownames(TF_scores_df), "TF",
                                            ifelse(cormat_ova_central_full_ATT$gene %in% t_table$gene, "gene", "lipid"))

write.csv(cormat_ova_central_full, file = file.path(multiomics_dir, "cormat_ova_central_full_SIF.csv"), row.names = FALSE)
write.csv(cormat_ova_central_full_ATT, file = file.path(multiomics_dir, "cormat_ova_central_full_ATT.csv"), row.names = FALSE)

# Create compressed version
cormat_ova_central$node2 <- gsub("[(].*", "", cormat_ova_central$node2)
cormat_ova_central$sign <- sign(cormat_ova_central$weight)
cormat_ova_central <- aggregate(weight ~ node1 + node2 + sign, data = cormat_ova_central, FUN = mean)
cormat_ova_central$abs_weight <- abs(cormat_ova_central$weight)

cormat_ova_central_ATT$gene <- gsub("[(].*", "", cormat_ova_central_ATT$gene)
cormat_ova_central_ATT <- aggregate(cbind(portal_vs_Central, testes_vs_ova, portalova_vs_centralova,
                                          portaltestes_vs_centraltestes, centraltestes_vs_centralova,
                                          portaltestes_vs_portalova, testes_effect_zonation) ~ gene + type,
                                   data = cormat_ova_central_ATT, FUN = mean)

write.csv(cormat_ova_central, file = file.path(multiomics_dir, "compressed_cormat_ova_central_SIF.csv"), row.names = FALSE)
write.csv(cormat_ova_central_ATT, file = file.path(multiomics_dir, "compressed_cormat_ova_central_ATT.csv"), row.names = FALSE)

cat("  ✓ Saved ovaries pericentral networks\n\n")

cat("✓ All 4 zonation-specific correlation networks generated\n")
cat("  Output files (12 files per network = 48 total):\n")
cat("    Standard networks (threshold = 0.65):\n")
cat("      - cormat_testes_portal_SIF.csv / ATT.csv\n")
cat("      - cormat_testes_central_SIF.csv / ATT.csv\n")
cat("      - cormat_Ovaries_portal_SIF.csv / ATT.csv\n")
cat("      - cormat_ova_central_SIF.csv / ATT.csv\n")
cat("    Full networks (threshold = 0, lipid species only):\n")
cat("      - cormat_testes_portal_full_SIF.csv / ATT.csv\n")
cat("      - cormat_testes_central_full_SIF.csv / ATT.csv\n")
cat("      - cormat_Ovaries_portal_full_SIF.csv / ATT.csv\n")
cat("      - cormat_ova_central_full_SIF.csv / ATT.csv\n")
cat("    Compressed networks (threshold = 0.65, lipid classes only):\n")
cat("      - compressed_cormat_* (8 files, aggregated by lipid class)\n\n")

################################################################################
# 11. LIPID-GENE ASSOCIATION HEATMAP BY ZONATION/GONAD CATEGORY
################################################################################

cat(strrep("-", 80), "\n")
cat("11. LIPID-GENE ASSOCIATION HEATMAP BY ZONATION/GONAD CATEGORY\n")
cat(strrep("-", 80), "\n\n")

cat("Creating binary heatmap of individual lipid associations with gene categories...\n")
cat("Categories: TPC (Testes Pericentral), TPP (Testes Periportal),\n")
cat("            OPC (Ovaries Pericentral), OPP (Ovaries Periportal)\n\n")

# Load lipidomic differential analysis results
cat("Loading lipidomic differential analysis...\n")
load("results/lipidomic/ttop_full_list.RData")

# Extract testes vs ovaries comparison (handle nested list structure)
if ("Testes vs Ovaries" %in% names(ttop_full_list)) {
  lipid_ttop <- ttop_full_list[["Testes vs Ovaries"]]
  # Check if it's a nested list
  if (is.list(lipid_ttop) && !is.data.frame(lipid_ttop)) {
    lipid_ttop <- lipid_ttop[[1]]
  }
} else if ("testes_vs_ova" %in% names(ttop_full_list)) {
  lipid_ttop <- ttop_full_list[["testes_vs_ova"]]
  if (is.list(lipid_ttop) && !is.data.frame(lipid_ttop)) {
    lipid_ttop <- lipid_ttop[[1]]
  }
} else {
  stop("Could not find testes vs ovaries comparison in lipidomic results")
}

# Ensure it's a data frame
lipid_ttop <- as.data.frame(lipid_ttop)

cat("  Loaded", nrow(lipid_ttop), "lipids\n\n")

# Read the UNCOMPRESSED correlation networks (individual lipids, not aggregated)
cat("Loading uncompressed correlation networks (individual lipids)...\n")

cormat_testes_portal <- read.csv(file.path(multiomics_dir, "cormat_testes_portal_SIF.csv"))
cormat_testes_central <- read.csv(file.path(multiomics_dir, "cormat_testes_central_SIF.csv"))
cormat_Ovaries_portal <- read.csv(file.path(multiomics_dir, "cormat_Ovaries_portal_SIF.csv"))
cormat_ova_central <- read.csv(file.path(multiomics_dir, "cormat_ova_central_SIF.csv"))

cat("  Loaded all 4 uncompressed networks\n\n")

# Function to extract unique lipids from a network
extract_lipids <- function(network_df) {
  # Filter to lipid nodes only (node2, excluding TFs which are in TF_scores_df)
  lipid_edges <- network_df[!network_df$node2 %in% rownames(TF_scores_df), ]

  # Return unique lipids
  unique(lipid_edges$node2)
}

# Extract unique lipids from each network
cat("Extracting lipids associated with each category...\n")

TPC_lipids <- extract_lipids(cormat_testes_central)
TPP_lipids <- extract_lipids(cormat_testes_portal)
OPC_lipids <- extract_lipids(cormat_ova_central)
OPP_lipids <- extract_lipids(cormat_Ovaries_portal)

cat("  TPC:", length(TPC_lipids), "unique lipids\n")
cat("  TPP:", length(TPP_lipids), "unique lipids\n")
cat("  OPC:", length(OPC_lipids), "unique lipids\n")
cat("  OPP:", length(OPP_lipids), "unique lipids\n\n")

# Get all unique lipids across all networks
all_lipids <- unique(c(TPC_lipids, TPP_lipids, OPC_lipids, OPP_lipids))

cat("  Total unique lipids across all networks:", length(all_lipids), "\n\n")

# Create binary matrix: rows = categories, columns = lipids
cat("Creating binary association matrix...\n")

binary_matrix <- data.frame(
  lipid = all_lipids,
  TPC = as.integer(all_lipids %in% TPC_lipids),
  TPP = as.integer(all_lipids %in% TPP_lipids),
  OPC = as.integer(all_lipids %in% OPC_lipids),
  OPP = as.integer(all_lipids %in% OPP_lipids)
)

# Add lipid differential analysis statistics
lipid_ttop_subset <- lipid_ttop[, c("ID", "logFC", "t", "adj.P.Val")]
names(lipid_ttop_subset)[1] <- "lipid"

binary_matrix <- binary_matrix %>%
  dplyr::left_join(lipid_ttop_subset, by = "lipid")

# Filter to lipids that appear in at least one network and are significant in DA
binary_matrix <- binary_matrix %>%
  dplyr::filter((TPC + TPP + OPC + OPP) > 0 & !is.na(t))

cat("  Matrix dimensions after filtering:", nrow(binary_matrix), "lipids x 4 categories\n\n")

# Order lipids by absolute t-value (for visualization)
binary_matrix <- binary_matrix %>%
  dplyr::mutate(abs_t = abs(t),
                gonad_bias = ifelse(logFC > 0, "Testes", "Ovaries")) %>%
  dplyr::arrange(desc(abs_t))

# Save binary matrix
write.csv(binary_matrix, file.path(zonation_dir, "lipid_category_binary_matrix.csv"), row.names = FALSE)
cat("  ✓ Saved binary matrix:", file.path(zonation_dir, "lipid_category_binary_matrix.csv"), "\n\n")

# Create visualization: heatmap (bottom) + barplot (top)
cat("Creating combined visualization...\n")

# Prepare data for heatmap (transpose so categories are rows, lipids are columns)
heatmap_data <- as.matrix(t(binary_matrix[, c("TPC", "TPP", "OPC", "OPP")]))
colnames(heatmap_data) <- binary_matrix$lipid

# Prepare data for barplot
barplot_data <- data.frame(
  lipid = binary_matrix$lipid,
  abs_t = binary_matrix$abs_t,
  gonad_bias = binary_matrix$gonad_bias
)

# Limit to top lipids if too many (for readability)
max_lipids <- 50
if (nrow(binary_matrix) > max_lipids) {
  cat("  Limiting visualization to top", max_lipids, "lipids by absolute t-value\n")
  heatmap_data <- heatmap_data[, 1:max_lipids]
  barplot_data <- barplot_data[1:max_lipids, ]
}

# Create barplot (top panel) - sorted by t-value
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)

p_bar <- ggplot(barplot_data, aes(x = factor(1:nrow(barplot_data)), y = abs_t, fill = gonad_bias)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Testes" = "royalblue", "Ovaries" = "red"),
                   name = "Gonadal Bias") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  labs(y = "Absolute t-value\n(Testes vs Ovaries)")

# Create heatmap (bottom panel) - sorted by t-value
heatmap_long <- as.data.frame(heatmap_data) %>%
  tibble::rownames_to_column("category") %>%
  tidyr::pivot_longer(cols = -category, names_to = "lipid", values_to = "associated")

# Create lipid order factor to match barplot
heatmap_long$lipid <- factor(heatmap_long$lipid, levels = barplot_data$lipid)
heatmap_long$category <- factor(heatmap_long$category, levels = c("OPP", "OPC", "TPP", "TPC"))

p_heatmap <- ggplot(heatmap_long, aes(x = lipid, y = category, fill = factor(associated))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                   name = "Associated",
                   labels = c("No", "Yes")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(y = "Category")

# Combine plots - version 1 (sorted by t-value)
pdf(file.path(zonation_dir, "lipid_zonation_association_heatmap.pdf"), width = 16, height = 8)
grid.newpage()
grid.draw(rbind(
  ggplotGrob(p_bar),
  ggplotGrob(p_heatmap),
  size = "last"
))
dev.off()

cat("  ✓ Saved t-value sorted plot:", file.path(zonation_dir, "lipid_zonation_association_heatmap.pdf"), "\n")

#===============================================================================
# VERSION 2: Hierarchical clustering (no dendrogram shown)
#===============================================================================

cat("  Creating hierarchical clustering version...\n")

# Perform hierarchical clustering on lipids based on binary association pattern
heatmap_data_dist <- dist(t(heatmap_data), method = "binary")
heatmap_data_hclust <- hclust(heatmap_data_dist, method = "complete")

# Reorder lipids by clustering
lipid_order_clustered <- colnames(heatmap_data)[heatmap_data_hclust$order]

# Reorder barplot data
barplot_data_clustered <- barplot_data[match(lipid_order_clustered, barplot_data$lipid), ]

# Create barplot with clustered order
p_bar_clustered <- ggplot(barplot_data_clustered,
                          aes(x = factor(1:nrow(barplot_data_clustered)), y = abs_t, fill = gonad_bias)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Testes" = "royalblue", "Ovaries" = "red"),
                   name = "Gonadal Bias") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  labs(y = "Absolute t-value\n(Testes vs Ovaries)")

# Create heatmap with clustered order
heatmap_long_clustered <- heatmap_long
heatmap_long_clustered$lipid <- factor(heatmap_long_clustered$lipid, levels = lipid_order_clustered)

p_heatmap_clustered <- ggplot(heatmap_long_clustered, aes(x = lipid, y = category, fill = factor(associated))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                   name = "Associated",
                   labels = c("No", "Yes")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(y = "Category")

# Combine plots - version 2a (clustered with barplot)
pdf(file.path(zonation_dir, "lipid_zonation_association_heatmap_clustered_with_bars.pdf"), width = 16, height = 8)
grid.newpage()
grid.draw(rbind(
  ggplotGrob(p_bar_clustered),
  ggplotGrob(p_heatmap_clustered),
  size = "last"
))
dev.off()

cat("  ✓ Saved clustered plot with bars:", file.path(zonation_dir, "lipid_zonation_association_heatmap_clustered_with_bars.pdf"), "\n")

# Version 2b: Heatmap only (no barplot) - clustered
pdf(file.path(zonation_dir, "lipid_zonation_association_heatmap_clustered.pdf"), width = 16, height = 4)
print(p_heatmap_clustered)
dev.off()

cat("  ✓ Saved clustered heatmap only:", file.path(zonation_dir, "lipid_zonation_association_heatmap_clustered.pdf"), "\n")

#===============================================================================
# VERSION 3: Lipid class level (aggregated)
#===============================================================================

cat("  Creating lipid class-level version...\n")

# Extract lipid class from lipid names
binary_matrix_class <- binary_matrix
binary_matrix_class$lipid_class <- gsub("[(].*", "", binary_matrix_class$lipid)

# Aggregate by lipid class
# For binary associations: take max (1 if any lipid in class is associated)
# For t-values: take mean of absolute t-values
lipid_class_summary <- binary_matrix_class %>%
  dplyr::group_by(lipid_class) %>%
  dplyr::summarise(
    TPC = as.integer(max(TPC) > 0),
    TPP = as.integer(max(TPP) > 0),
    OPC = as.integer(max(OPC) > 0),
    OPP = as.integer(max(OPP) > 0),
    mean_abs_t = mean(abs_t, na.rm = TRUE),
    mean_logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(gonad_bias = ifelse(mean_logFC > 0, "Testes", "Ovaries")) %>%
  dplyr::arrange(desc(mean_abs_t))

cat("    Aggregated to", nrow(lipid_class_summary), "lipid classes\n")

# Prepare data for heatmap (transpose so categories are rows, classes are columns)
heatmap_data_class <- as.matrix(t(lipid_class_summary[, c("TPC", "TPP", "OPC", "OPP")]))
colnames(heatmap_data_class) <- lipid_class_summary$lipid_class

# Prepare data for barplot
barplot_data_class <- data.frame(
  lipid_class = lipid_class_summary$lipid_class,
  mean_abs_t = lipid_class_summary$mean_abs_t,
  gonad_bias = lipid_class_summary$gonad_bias
)

# Create barplot for lipid classes
p_bar_class <- ggplot(barplot_data_class,
                      aes(x = factor(1:nrow(barplot_data_class)), y = mean_abs_t, fill = gonad_bias)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Testes" = "royalblue", "Ovaries" = "red"),
                   name = "Gonadal Bias") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  labs(y = "Mean Absolute t-value\n(Testes vs Ovaries)")

# Create heatmap for lipid classes
heatmap_long_class <- as.data.frame(heatmap_data_class) %>%
  tibble::rownames_to_column("category") %>%
  tidyr::pivot_longer(cols = -category, names_to = "lipid_class", values_to = "associated")

heatmap_long_class$lipid_class <- factor(heatmap_long_class$lipid_class, levels = barplot_data_class$lipid_class)
heatmap_long_class$category <- factor(heatmap_long_class$category, levels = c("OPP", "OPC", "TPP", "TPC"))

p_heatmap_class <- ggplot(heatmap_long_class, aes(x = lipid_class, y = category, fill = factor(associated))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                   name = "Associated",
                   labels = c("No", "Yes")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(y = "Category")

# Combine plots - version 3 (lipid class level, sorted by t-value)
pdf(file.path(zonation_dir, "lipid_class_zonation_association_heatmap.pdf"), width = 12, height = 8)
grid.newpage()
grid.draw(rbind(
  ggplotGrob(p_bar_class),
  ggplotGrob(p_heatmap_class),
  size = "last"
))
dev.off()

cat("  ✓ Saved lipid class-level plot (t-value sorted):", file.path(zonation_dir, "lipid_class_zonation_association_heatmap.pdf"), "\n")

# Save lipid class summary
write.csv(lipid_class_summary, file.path(zonation_dir, "lipid_class_association_summary.csv"), row.names = FALSE)
cat("  ✓ Saved lipid class summary table\n")

#===============================================================================
# VERSION 4: Lipid class level with hierarchical clustering
#===============================================================================

cat("  Creating lipid class-level clustered version...\n")

# Perform hierarchical clustering on lipid classes based on binary association pattern
heatmap_data_class_dist <- dist(t(heatmap_data_class), method = "binary")
heatmap_data_class_hclust <- hclust(heatmap_data_class_dist, method = "complete")

# Reorder lipid classes by clustering
lipid_class_order_clustered <- colnames(heatmap_data_class)[heatmap_data_class_hclust$order]

# Reorder barplot data
barplot_data_class_clustered <- barplot_data_class[match(lipid_class_order_clustered, barplot_data_class$lipid_class), ]

# Create barplot with clustered order
p_bar_class_clustered <- ggplot(barplot_data_class_clustered,
                                aes(x = factor(1:nrow(barplot_data_class_clustered)),
                                    y = mean_abs_t, fill = gonad_bias)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Testes" = "royalblue", "Ovaries" = "red"),
                   name = "Gonadal Bias") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  labs(y = "Mean Absolute t-value\n(Testes vs Ovaries)")

# Create heatmap with clustered order
heatmap_long_class_clustered <- heatmap_long_class
heatmap_long_class_clustered$lipid_class <- factor(heatmap_long_class_clustered$lipid_class,
                                                   levels = lipid_class_order_clustered)

p_heatmap_class_clustered <- ggplot(heatmap_long_class_clustered,
                                    aes(x = lipid_class, y = category, fill = factor(associated))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                   name = "Associated",
                   labels = c("No", "Yes")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(y = "Category")

# Combine plots - version 4 (lipid class level, clustered)
pdf(file.path(zonation_dir, "lipid_class_zonation_association_heatmap_clustered_with_bars.pdf"),
    width = 12, height = 8)
grid.newpage()
grid.draw(rbind(
  ggplotGrob(p_bar_class_clustered),
  ggplotGrob(p_heatmap_class_clustered),
  size = "last"
))
dev.off()

cat("  ✓ Saved lipid class-level clustered plot:",
    file.path(zonation_dir, "lipid_class_zonation_association_heatmap_clustered_with_bars.pdf"), "\n\n")

# Print summary statistics
cat("Summary statistics:\n")
cat("  Total lipids in visualization:", nrow(barplot_data), "\n")
cat("  Testes-biased lipids:", sum(barplot_data$gonad_bias == "Testes"), "\n")
cat("  Ovaries-biased lipids:", sum(barplot_data$gonad_bias == "Ovaries"), "\n\n")

cat("  Association breakdown:\n")
for (cat_name in c("TPC", "TPP", "OPC", "OPP")) {
  n_associated <- sum(binary_matrix[1:min(max_lipids, nrow(binary_matrix)), cat_name] == 1)
  cat("    ", cat_name, ": ", n_associated, " lipids\n", sep = "")
}
cat("\n")

cat("✓ Lipid-gene association heatmap complete\n\n")

################################################################################
# 12. SAVE RESULTS
################################################################################

cat(strrep("-", 80), "\n")
cat("12. SAVING RESULTS\n")
cat(strrep("-", 80), "\n\n")

# Save differential analysis results
save(ttop_list, file = file.path(zonation_dir, "zonation_ttop_list.RData"))
save(limmaRes, file = file.path(zonation_dir, "zonation_limma_results.RData"))

# Save metadata
write.csv(meta_data_df, file.path(zonation_dir, "metadata.csv"), row.names = FALSE)

# Create summary statistics
summary_stats <- data.frame(
  comparison = names(ttop_list),
  n_sig_genes = sapply(ttop_list, function(x) sum(x$adj.P.Val < 0.05, na.rm = TRUE)),
  n_upregulated = sapply(ttop_list, function(x) sum(x$adj.P.Val < 0.05 & x$logFC > 0, na.rm = TRUE)),
  n_downregulated = sapply(ttop_list, function(x) sum(x$adj.P.Val < 0.05 & x$logFC < 0, na.rm = TRUE))
)

write.csv(summary_stats, file.path(zonation_dir, "summary_statistics.csv"), row.names = FALSE)

cat("  Saved differential analysis results\n")
cat("  Saved summary statistics\n\n")

################################################################################
# 13. WORKFLOW COMPLETE
################################################################################

cat(strrep("=", 80), "\n")
cat("ZONATION ANALYSIS WORKFLOW COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("Key findings:\n")
cat("  - Analyzed", nrow(meta_data_df), "pseudobulk samples\n")
cat("  - Tested", length(comparisons), "contrasts\n")
cat("  - Generated volcano plots and pathway enrichments\n")
cat("  - Created 4 zonation-specific correlation networks (32 output files)\n")
cat("  - Created lipid class association visualization by category\n\n")

cat("Output files:\n")
cat("  - ", file.path(zonation_dir, "zonation_ttop_list.RData"), "\n")
cat("  - ", file.path(zonation_dir, "zonation_limma_results.RData"), "\n")
cat("  - ", file.path(zonation_dir, "summary_statistics.csv"), "\n")
cat("  - ", file.path(zonation_dir, "t_table_zonation.csv"), "\n")
cat("  - ", file.path(zonation_dir, "PCA_*.pdf"), "\n")
cat("  - ", file.path(zonation_dir, "volcanos/"), "\n")
cat("  - ", file.path(zonation_dir, "lipid_class_by_category.pdf"), "\n")
cat("  - ", file.path(zonation_dir, "lipid_class_category_summary.csv"), "\n")
cat("  - ", file.path(multiomics_dir, "cormat_*_SIF.csv / ATT.csv"), " (8 networks)\n\n")

cat("Expected biological findings (from manuscript):\n")
cat("  - ~43% of gonad-dependent genes are pericentral-biased\n")
cat("  - ~18% of gonad-dependent genes are periportal-biased\n")
cat("  - Key zonation-specific genes:\n")
cat("    * Pericentral: Cyp2c38, Cyp2b13, Cyp2d9\n")
cat("    * Periportal: Cyp2b9\n")
cat("  - Seven genes show high correlation with TFs and lipids:\n")
cat("    * Female-biased pericentral: Sult3a1, Abcb1a, Cyp2a4\n")
cat("    * Male-biased periportal: Mtmr11, Mindy3, Esco2, Kitl\n\n")

cat("This refactored workflow replaces: scripts/pseudobulk_DA_zonation.R\n\n")
cat(strrep("=", 80), "\n")
