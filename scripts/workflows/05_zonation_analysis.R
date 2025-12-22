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

# Load required libraries for plotting
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)

#===============================================================================
# Prepare lipid class level data for visualization
#===============================================================================

cat("  Preparing lipid class-level data...\n")

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

# Save lipid class summary
write.csv(lipid_class_summary, file.path(zonation_dir, "lipid_class_association_summary.csv"), row.names = FALSE)
cat("  ✓ Saved lipid class summary table\n")

#===============================================================================
# Lipid class level visualization with jitter plot
# Uses the same lipid classes and ordering as 03_fcg_lipidomic_analysis.R
#===============================================================================

cat("  Creating lipid class-level jitter plot version...\n")

# Load lipid class enrichment scores (same as used in 03_fcg_lipidomic_analysis.R)
cat("    Loading lipid class enrichment scores...\n")
load("results/lipidomic/enrichment/lipid_class_enrichment_scores.RData")

# Extract enrichment scores for gonadal comparison (testes vs ovaries)
# Condition names in enrichment results are formatted as "Testes vs Ovaries_tval_Lipids"
lipid_enrichment_gonade <- lipid_class_results %>%
  dplyr::filter(grepl("^Testes vs Ovaries", condition)) %>%
  dplyr::select(source, enrichment = score) %>%
  dplyr::mutate(direction = ifelse(enrichment < 0, "negative", "positive"))

# Apply enrichment threshold (same as 03_fcg_lipidomic_analysis.R)
lipid_enrichment_gonade$direction <- ifelse(
  abs(lipid_enrichment_gonade$enrichment) > 1.8,
  lipid_enrichment_gonade$direction,
  "non-significant"
)

cat("    Found enrichment scores for", nrow(lipid_enrichment_gonade), "lipid classes\n")

# Create jitter data from FULL lipidomic analysis (same approach as 03_fcg_lipidomic_analysis.R)
# This uses lipid_ttop which was loaded earlier (line ~753)
lipid_classes_map <- data.frame(
  ID = lipid_ttop$ID,
  source = gsub("[(].*", "", lipid_ttop$ID),  # Extract lipid class from name
  testes_vs_ova = lipid_ttop$t,
  stringsAsFactors = FALSE
)

# Join individual lipid t-values with enrichment scores (same as 03_fcg_lipidomic_analysis.R)
# This uses ALL lipid classes with enrichment scores - NO filtering to heatmap classes
jitter_data <- lipid_classes_map %>%
  dplyr::inner_join(lipid_enrichment_gonade, by = "source")

cat("    Jitter data: ", nrow(jitter_data), " individual lipids across ",
    length(unique(jitter_data$source)), " lipid classes\n")

# Order lipid classes by enrichment score (descending, like 03_fcg_lipidomic_analysis.R)
# Use ALL classes with enrichment - same as workflow 03
classes_with_enrichment <- lipid_enrichment_gonade %>%
  dplyr::arrange(desc(enrichment)) %>%
  dplyr::pull(source)

cat("    Classes ordered by enrichment:", length(classes_with_enrichment), "\n")

# Apply enrichment-based order to jitter data
jitter_data$source <- factor(jitter_data$source, levels = classes_with_enrichment)

#===============================================================================
# Helper function and sample subsets (used by threshold analysis)
#===============================================================================

# Helper function to get lipids correlated with genes
get_correlated_lipids <- function(gene_list, rna_data, lipid_data, cor_threshold = 0.65) {
  if (length(gene_list) == 0) return(character(0))

  # Filter to available genes
  available_genes <- intersect(gene_list, rownames(rna_data))
  if (length(available_genes) == 0) return(character(0))

  # Compute correlations
  cor_mat <- cor(t(rna_data[available_genes, , drop = FALSE]), t(lipid_data), use = "pairwise.complete.obs")

  # Get lipids with high correlation to any gene
  max_cors <- apply(abs(cor_mat), 2, max, na.rm = TRUE)
  lipids_passing <- names(max_cors[max_cors >= cor_threshold])

  return(lipids_passing)
}

# Get sample subsets for each condition
testes_portal_samples <- grep("Periportal_Testes", colnames(scaled_RNA), value = TRUE)
testes_central_samples <- grep("Pericentral_Testes", colnames(scaled_RNA), value = TRUE)
ovaries_portal_samples <- grep("Periportal_Ovaries", colnames(scaled_RNA), value = TRUE)
ovaries_central_samples <- grep("Pericentral_Ovaries", colnames(scaled_RNA), value = TRUE)

# Get matching lipidomic samples
lip_testes_samples <- gsub("_Peri.*_Testes", "", testes_portal_samples)
lip_ovaries_samples <- gsub("_Peri.*_Ovaries", "", ovaries_portal_samples)

#===============================================================================
# Threshold sensitivity analysis
# Shows how lipid class proportions change across zonation categories
# as the interaction threshold varies from 0 to 10
#===============================================================================

cat("  Creating threshold sensitivity analysis...\n")

# Define threshold range
threshold_range <- seq(0, 10, by = 0.1)
cat("    Testing", length(threshold_range), "thresholds from 0 to 10\n")

# Function to compute zonation assignments for a given threshold
compute_zonation_at_threshold <- function(thresh, t_table, scaled_RNA, lipidomic,
                                          testes_portal_samples, testes_central_samples,
                                          ovaries_portal_samples, ovaries_central_samples,
                                          lip_testes_samples, lip_ovaries_samples,
                                          cor_threshold = 0.65) {

  # Get genes at this threshold
  testes_genes <- t_table[abs(t_table$testes_effect_zonation) > thresh &
                           abs(t_table$portalova_vs_centralova) < 1, ]
  ovaries_genes <- t_table[abs(t_table$testes_effect_zonation) > thresh &
                            abs(t_table$portaltestes_vs_centraltestes) < 1, ]

  # Split by zone
  tpc_genes <- testes_genes[testes_genes$portaltestes_vs_centraltestes < 0, ]$gene
  tpp_genes <- testes_genes[testes_genes$portaltestes_vs_centraltestes > 0, ]$gene
  opc_genes <- ovaries_genes[ovaries_genes$portalova_vs_centralova < 0, ]$gene
  opp_genes <- ovaries_genes[ovaries_genes$portalova_vs_centralova > 0, ]$gene

  # Get correlated lipids for each category
  tpc_lipids <- get_correlated_lipids(tpc_genes, scaled_RNA[, testes_central_samples],
                                       lipidomic[, lip_testes_samples], cor_threshold)
  tpp_lipids <- get_correlated_lipids(tpp_genes, scaled_RNA[, testes_portal_samples],
                                       lipidomic[, lip_testes_samples], cor_threshold)
  opc_lipids <- get_correlated_lipids(opc_genes, scaled_RNA[, ovaries_central_samples],
                                       lipidomic[, lip_ovaries_samples], cor_threshold)
  opp_lipids <- get_correlated_lipids(opp_genes, scaled_RNA[, ovaries_portal_samples],
                                       lipidomic[, lip_ovaries_samples], cor_threshold)

  # Determine zonation for each lipid
  all_lipids <- rownames(lipidomic)
  periportal_lipids <- unique(c(tpp_lipids, opp_lipids))
  pericentral_lipids <- unique(c(tpc_lipids, opc_lipids))

  zonation <- dplyr::case_when(
    all_lipids %in% periportal_lipids & all_lipids %in% pericentral_lipids ~ "both",
    all_lipids %in% periportal_lipids ~ "periportal",
    all_lipids %in% pericentral_lipids ~ "pericentral",
    TRUE ~ "no_association"
  )

  data.frame(
    lipid = all_lipids,
    lipid_class = gsub("[(].*", "", all_lipids),
    zonation = zonation,
    stringsAsFactors = FALSE
  )
}

# Compute zonation at each threshold
cat("    Computing zonation assignments across thresholds...\n")
sensitivity_results <- lapply(threshold_range, function(thresh) {
  result <- compute_zonation_at_threshold(
    thresh, t_table, scaled_RNA, lipidomic,
    testes_portal_samples, testes_central_samples,
    ovaries_portal_samples, ovaries_central_samples,
    lip_testes_samples, lip_ovaries_samples
  )
  result$threshold <- thresh
  result
})

sensitivity_df <- do.call(rbind, sensitivity_results)
cat("    Generated", nrow(sensitivity_df), "lipid × threshold combinations\n")

# Calculate proportions by lipid class and zonation category for each threshold
cat("    Calculating proportions by lipid class...\n")

proportion_by_class <- sensitivity_df %>%
  dplyr::group_by(threshold, lipid_class, zonation) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(threshold, lipid_class) %>%
  dplyr::mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  dplyr::ungroup()

# Order zonation factor (bottom to top: green, purple, orange, grey)
proportion_by_class$zonation <- factor(
  proportion_by_class$zonation,
  levels = c("periportal", "both", "pericentral", "no_association")
)

# Get lipid classes ordered by enrichment (same as jitter plot)
proportion_by_class$lipid_class <- factor(
  proportion_by_class$lipid_class,
  levels = classes_with_enrichment
)

# Filter to classes that exist in our data
proportion_by_class <- proportion_by_class %>%
  dplyr::filter(!is.na(lipid_class))

cat("    Lipid classes:", length(unique(proportion_by_class$lipid_class)), "\n")

# Note: Per-class stacked bar plot will be created after optimal_threshold is determined

# Also create a summary plot showing total proportions across all lipids
proportion_total <- sensitivity_df %>%
  dplyr::group_by(threshold, zonation) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(threshold) %>%
  dplyr::mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  dplyr::ungroup()

proportion_total$zonation <- factor(
  proportion_total$zonation,
  levels = c("periportal", "both", "pericentral", "no_association")
)

# Note: p_sensitivity_total will be created and saved after optimal_threshold is determined
# so we can add the threshold line to the plot

# Save sensitivity data
write.csv(proportion_by_class,
          file.path(zonation_dir, "threshold_sensitivity_by_class.csv"),
          row.names = FALSE)
write.csv(proportion_total,
          file.path(zonation_dir, "threshold_sensitivity_total.csv"),
          row.names = FALSE)

cat("  ✓ Saved sensitivity data to CSV files\n\n")

#===============================================================================
# Version 5: Optimal threshold detection (Entropy + Stability)
# Combines information content (entropy) with stability analysis
#===============================================================================

cat("  Finding optimal threshold using entropy + stability...\n")

# Function to calculate Shannon entropy (normalized to 0-1)
calc_entropy <- function(proportions) {
  # Remove zeros to avoid log(0)
  p <- proportions[proportions > 0]
  if (length(p) <= 1) return(0)
  # Shannon entropy, normalized by max possible entropy (log of number of categories)
  entropy <- -sum(p * log(p))
  max_entropy <- log(length(proportions))
  return(entropy / max_entropy)
}

# Calculate entropy at each threshold (excluding no_association from entropy calc
# since we want to maximize info in the ASSOCIATED categories)
entropy_by_threshold <- proportion_total %>%
  dplyr::group_by(threshold) %>%
  dplyr::summarise(
    # Entropy of all categories
    entropy_all = calc_entropy(proportion),
    # Entropy of associated categories only (periportal, pericentral, both)
    entropy_associated = {
      assoc_props <- proportion[zonation != "no_association"]
      if (sum(assoc_props) > 0) {
        calc_entropy(assoc_props / sum(assoc_props))
      } else {
        0
      }
    },
    # Proportion of lipids with any association
    prop_associated = 1 - proportion[zonation == "no_association"],
    # Proportion of "both" lipids (among all lipids) - we want to minimize this
    prop_both = ifelse(any(zonation == "both"), proportion[zonation == "both"], 0),
    .groups = "drop"
  )

# Calculate stability: how much do classifications change between adjacent thresholds?
# Lower change = more stable
stability_data <- sensitivity_df %>%
  dplyr::arrange(lipid, threshold) %>%
  dplyr::group_by(lipid) %>%
  dplyr::mutate(
    prev_zonation = dplyr::lag(zonation),
    changed = zonation != prev_zonation
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(prev_zonation)) %>%
  dplyr::group_by(threshold) %>%
  dplyr::summarise(
    n_changed = sum(changed),
    prop_changed = mean(changed),
    .groups = "drop"
  )

# Merge entropy and stability data
optimal_metrics <- entropy_by_threshold %>%
  dplyr::left_join(stability_data, by = "threshold") %>%
  dplyr::mutate(
    # Stability score: inverse of change rate (1 = perfectly stable, 0 = everything changes)
    stability = 1 - prop_changed,
    # Specificity score: penalize having many "both" lipids (we want clear zonation)
    specificity = 1 - prop_both,
    # Combined score: balance entropy, stability, coverage, and specificity
    # We want HIGH entropy (diverse associations), HIGH stability, and LOW "both" proportion
    # Using cubed specificity to strongly penalize shared lipids
    combined_score = entropy_associated * stability * prop_associated * specificity^3
  )

# Find optimal threshold (maximum combined score)
optimal_idx <- which.max(optimal_metrics$combined_score)
optimal_threshold <- optimal_metrics$threshold[optimal_idx]

cat("    Optimal threshold:", optimal_threshold, "\n")
cat("    At this threshold:\n")
cat("      - Entropy (associated):", round(optimal_metrics$entropy_associated[optimal_idx], 3), "\n")
cat("      - Stability:", round(optimal_metrics$stability[optimal_idx], 3), "\n")
cat("      - Proportion associated:", round(optimal_metrics$prop_associated[optimal_idx], 3), "\n")
cat("      - Proportion 'both':", round(optimal_metrics$prop_both[optimal_idx], 3), "\n")
cat("      - Specificity (1-both):", round(optimal_metrics$specificity[optimal_idx], 3), "\n")
cat("      - Combined score:", round(optimal_metrics$combined_score[optimal_idx], 3), "\n")

# Create visualization of the optimization
p_optimal <- ggplot(optimal_metrics, aes(x = threshold)) +
  # Entropy of associated categories
  geom_line(aes(y = entropy_associated, color = "Entropy (associated)"), linewidth = 1) +
  # Stability
  geom_line(aes(y = stability, color = "Stability"), linewidth = 1) +
  # Proportion associated
  geom_line(aes(y = prop_associated, color = "Prop. associated"), linewidth = 1, linetype = "dashed") +
  # Specificity (1 - prop_both)
  geom_line(aes(y = specificity, color = "Specificity (1-both)"), linewidth = 1, linetype = "dashed") +
  # Combined score
  geom_line(aes(y = combined_score, color = "Combined score"), linewidth = 1.5) +
  # Mark optimal threshold
  geom_vline(xintercept = optimal_threshold, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = optimal_threshold + 0.3, y = 0.95,
           label = paste0("Optimal: ", optimal_threshold),
           hjust = 0, color = "red", size = 4) +
  scale_color_manual(
    name = "Metric",
    values = c("Entropy (associated)" = "blue", "Stability" = "darkgreen",
               "Prop. associated" = "grey50", "Specificity (1-both)" = "orange",
               "Combined score" = "red")
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Interaction threshold (|t| for zonation × gonad)",
    y = "Score (normalized 0-1)",
    title = "Optimal threshold detection: Entropy × Stability × Coverage × Specificity",
    subtitle = paste0("Optimal threshold = ", optimal_threshold,
                      " (penalizes shared 'both' lipids)")
  )

ggsave(file.path(zonation_dir, "optimal_threshold_detection.pdf"),
       plot = p_optimal,
       width = 10, height = 6)

cat("  ✓ Saved optimal threshold plot:",
    file.path(zonation_dir, "optimal_threshold_detection.pdf"), "\n")

# Save the metrics data
write.csv(optimal_metrics,
          file.path(zonation_dir, "optimal_threshold_metrics.csv"),
          row.names = FALSE)

cat("  ✓ Saved optimal threshold metrics to CSV\n\n")

# Now create the total sensitivity plot with the optimal threshold line
p_sensitivity_total <- ggplot(proportion_total,
                              aes(x = threshold, y = proportion, fill = zonation)) +
  geom_area(alpha = 0.8) +
  geom_vline(xintercept = optimal_threshold, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = optimal_threshold + 0.2, y = 0.5,
           label = paste0("Optimal: ", optimal_threshold),
           hjust = 0, color = "red", size = 4, fontface = "bold") +
  scale_fill_manual(
    name = "Zonation Association",
    values = c("periportal" = "lightgreen", "both" = "mediumpurple",
               "pericentral" = "orange", "no_association" = "grey70"),
    labels = c("Periportal", "Both", "Pericentral", "No association")
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Interaction threshold (|t| for zonation × gonad)",
    y = "Proportion of all lipid species",
    title = "Overall lipid zonation associations across interaction thresholds"
  )

ggsave(file.path(zonation_dir, "lipid_zonation_threshold_sensitivity_total.pdf"),
       plot = p_sensitivity_total,
       width = 10, height = 6)

cat("  ✓ Saved total sensitivity plot with optimal threshold:",
    file.path(zonation_dir, "lipid_zonation_threshold_sensitivity_total.pdf"), "\n\n")

# Print summary of zonation at optimal threshold
cat("  Zonation distribution at optimal threshold (", optimal_threshold, "):\n", sep = "")
optimal_zonation <- sensitivity_df %>%
  dplyr::filter(threshold == optimal_threshold) %>%
  dplyr::count(zonation) %>%
  dplyr::mutate(percent = round(100 * n / sum(n), 1))

for (i in 1:nrow(optimal_zonation)) {
  cat("    ", optimal_zonation$zonation[i], ": ", optimal_zonation$n[i],
      " (", optimal_zonation$percent[i], "%)\n", sep = "")
}
cat("\n")

#===============================================================================
# Lipid class proportions at optimal threshold (stacked bar chart)
# Split by gonad bias (Testes vs Ovaries)
#===============================================================================

cat("  Creating lipid class proportions at optimal threshold...\n")

# Get gonad bias for each lipid from lipid_ttop
lipid_gonad_bias <- data.frame(
  lipid = lipid_ttop$ID,
  gonad_bias = ifelse(lipid_ttop$logFC > 0, "Testes", "Ovaries"),
  stringsAsFactors = FALSE
)

# Get zonation for each lipid at optimal threshold and add gonad bias
lipid_zonation_optimal <- sensitivity_df %>%
  dplyr::filter(threshold == optimal_threshold) %>%
  dplyr::left_join(lipid_gonad_bias, by = "lipid")

# Calculate proportions by lipid class, gonad bias, and zonation
class_proportions_by_gonad <- lipid_zonation_optimal %>%
  dplyr::group_by(lipid_class, gonad_bias, zonation) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(lipid_class, gonad_bias) %>%
  dplyr::mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  dplyr::ungroup()

# Order lipid classes by enrichment (same as jitter plot)
# Reverse order for facet_grid with coord_flip to match jitter plot orientation
class_proportions_by_gonad$lipid_class <- factor(
  class_proportions_by_gonad$lipid_class,
  levels = rev(classes_with_enrichment)
)

# Filter to classes that exist
class_proportions_by_gonad <- class_proportions_by_gonad %>%
  dplyr::filter(!is.na(lipid_class))

# Order zonation factor
class_proportions_by_gonad$zonation <- factor(
  class_proportions_by_gonad$zonation,
  levels = c("periportal", "both", "pericentral", "no_association")
)

# Order gonad_bias factor
class_proportions_by_gonad$gonad_bias <- factor(
  class_proportions_by_gonad$gonad_bias,
  levels = c("Testes", "Ovaries")
)

cat("    Lipid classes with data:", length(unique(class_proportions_by_gonad$lipid_class)), "\n")
cat("    Testes-biased lipids:", sum(lipid_gonad_bias$gonad_bias == "Testes"), "\n")
cat("    Ovaries-biased lipids:", sum(lipid_gonad_bias$gonad_bias == "Ovaries"), "\n")

# Create stacked bar chart with gonad bias stacked vertically within each lipid class
# Use faceting by lipid class so T/O bars are touching within each class
p_class_proportions <- ggplot(class_proportions_by_gonad,
                               aes(x = gonad_bias, y = proportion, fill = zonation)) +
  geom_bar(stat = "identity", width = 1, position = "stack") +
  facet_grid(lipid_class ~ ., switch = "y") +
  scale_fill_manual(
    name = "Zonation\nAssociation",
    values = c("periportal" = "lightgreen", "both" = "mediumpurple",
               "pericentral" = "orange", "no_association" = "grey70"),
    labels = c("Periportal", "Both", "Pericentral", "No association")
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Testes" = "T", "Ovaries" = "O")) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0.15, "lines"),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
    strip.placement = "outside"
  ) +
  labs(
    x = NULL,
    y = "Proportion of lipid species",
    title = paste0("Lipid class zonation proportions at optimal threshold (", optimal_threshold, ")"),
    subtitle = "T = Testes-biased, O = Ovaries-biased"
  )

ggsave(file.path(zonation_dir, "lipid_class_zonation_proportions.pdf"),
       plot = p_class_proportions,
       width = 4, height = 8)

cat("  ✓ Saved lipid class proportions plot:",
    file.path(zonation_dir, "lipid_class_zonation_proportions.pdf"), "\n\n")

#===============================================================================
# Standalone vertical jitter plot (like workflow 03)
# Uses the OPTIMAL THRESHOLD determined above
# Duplicates points that are in "both" categories into periportal and pericentral
#===============================================================================

cat("  Creating standalone vertical jitter plot using optimal threshold...\n")
cat("    Using optimal interaction threshold:", optimal_threshold, "\n")

# Get gene sets at optimal threshold
v3_testes_zonation_genes <- t_table[abs(t_table$testes_effect_zonation) > optimal_threshold &
                                     abs(t_table$portalova_vs_centralova) < 1, ]
v3_ovaries_zonation_genes <- t_table[abs(t_table$testes_effect_zonation) > optimal_threshold &
                                      abs(t_table$portaltestes_vs_centraltestes) < 1, ]

# Split by zone
v3_testes_pericentral <- v3_testes_zonation_genes[v3_testes_zonation_genes$portaltestes_vs_centraltestes < 0, ]$gene
v3_testes_periportal <- v3_testes_zonation_genes[v3_testes_zonation_genes$portaltestes_vs_centraltestes > 0, ]$gene
v3_ovaries_pericentral <- v3_ovaries_zonation_genes[v3_ovaries_zonation_genes$portalova_vs_centralova < 0, ]$gene
v3_ovaries_periportal <- v3_ovaries_zonation_genes[v3_ovaries_zonation_genes$portalova_vs_centralova > 0, ]$gene

cat("    Gene counts at optimal threshold (", optimal_threshold, "):\n", sep = "")
cat("      Testes pericentral:", length(v3_testes_pericentral), "\n")
cat("      Testes periportal:", length(v3_testes_periportal), "\n")
cat("      Ovaries pericentral:", length(v3_ovaries_pericentral), "\n")
cat("      Ovaries periportal:", length(v3_ovaries_periportal), "\n")

# Compute lipid correlations at optimal threshold
v3_TPC_lipids <- get_correlated_lipids(
  v3_testes_pericentral,
  scaled_RNA[, testes_central_samples],
  lipidomic[, lip_testes_samples]
)
v3_TPP_lipids <- get_correlated_lipids(
  v3_testes_periportal,
  scaled_RNA[, testes_portal_samples],
  lipidomic[, lip_testes_samples]
)
v3_OPC_lipids <- get_correlated_lipids(
  v3_ovaries_pericentral,
  scaled_RNA[, ovaries_central_samples],
  lipidomic[, lip_ovaries_samples]
)
v3_OPP_lipids <- get_correlated_lipids(
  v3_ovaries_periportal,
  scaled_RNA[, ovaries_portal_samples],
  lipidomic[, lip_ovaries_samples]
)

cat("    Lipid associations at optimal threshold:\n")
cat("      TPC lipids:", length(v3_TPC_lipids), "\n")
cat("      TPP lipids:", length(v3_TPP_lipids), "\n")
cat("      OPC lipids:", length(v3_OPC_lipids), "\n")
cat("      OPP lipids:", length(v3_OPP_lipids), "\n")

# Create zonation mapping
v3_lipid_zonation <- data.frame(
  lipid = unique(c(v3_TPC_lipids, v3_TPP_lipids, v3_OPC_lipids, v3_OPP_lipids)),
  stringsAsFactors = FALSE
)
v3_lipid_zonation$is_periportal <- v3_lipid_zonation$lipid %in% c(v3_TPP_lipids, v3_OPP_lipids)
v3_lipid_zonation$is_pericentral <- v3_lipid_zonation$lipid %in% c(v3_TPC_lipids, v3_OPC_lipids)
v3_lipid_zonation$zonation_direction <- dplyr::case_when(
  v3_lipid_zonation$is_periportal & v3_lipid_zonation$is_pericentral ~ "both",
  v3_lipid_zonation$is_periportal ~ "periportal",
  v3_lipid_zonation$is_pericentral ~ "pericentral",
  TRUE ~ "no_association"
)

cat("    Zonation breakdown:\n")
cat("      Periportal only:", sum(v3_lipid_zonation$zonation_direction == "periportal"), "\n")
cat("      Pericentral only:", sum(v3_lipid_zonation$zonation_direction == "pericentral"), "\n")
cat("      Both:", sum(v3_lipid_zonation$zonation_direction == "both"), "\n")

# Create jitter data with zonation coloring
jitter_data_v3 <- jitter_data %>%
  dplyr::left_join(
    v3_lipid_zonation[, c("lipid", "zonation_direction")],
    by = c("ID" = "lipid")
  ) %>%
  dplyr::mutate(
    zonation_direction = ifelse(is.na(zonation_direction), "no_network_data", zonation_direction)
  )

cat("    Total lipids colored:", sum(jitter_data_v3$zonation_direction != "no_network_data"), "\n")

# Split "both" into two rows: one periportal, one pericentral
both_rows <- jitter_data_v3 %>% dplyr::filter(zonation_direction == "both")

if (nrow(both_rows) > 0) {
  # Duplicated points go into their respective zonation groups but get "both" color
  both_as_periportal <- both_rows %>% dplyr::mutate(zonation_direction = "periportal", color_group = "both")
  both_as_pericentral <- both_rows %>% dplyr::mutate(zonation_direction = "pericentral", color_group = "both")

  # Non-duplicated points: color matches their zonation
  jitter_data_v3 <- jitter_data_v3 %>%
    dplyr::filter(zonation_direction != "both") %>%
    dplyr::mutate(color_group = as.character(zonation_direction)) %>%
    dplyr::bind_rows(both_as_periportal, both_as_pericentral)

  cat("    Duplicated", nrow(both_rows), "lipids from 'both' into periportal and pericentral (colored purple)\n")
} else {
  jitter_data_v3$color_group <- as.character(jitter_data_v3$zonation_direction)
}

# Reorder zonation groups: green (periportal), grey (no_network_data), orange (pericentral)
jitter_data_v3$zonation_direction <- factor(
  jitter_data_v3$zonation_direction,
  levels = c("periportal", "no_network_data", "pericentral")
)

# Reorder color groups for legend
jitter_data_v3$color_group <- factor(
  jitter_data_v3$color_group,
  levels = c("periportal", "no_network_data", "pericentral", "both")
)

# Set point sizes: orange/green = 2, purple (both) = 1, grey = 1
jitter_data_v3$point_size <- dplyr::case_when(
  jitter_data_v3$color_group %in% c("periportal", "pericentral") ~ 2,
  jitter_data_v3$color_group == "both" ~ 1,
  TRUE ~ 1  # no_network_data
)

# Create alternating background shading for lipid classes
n_classes <- length(levels(jitter_data_v3$source))
shade_data <- data.frame(
  xmin = seq(0.5, n_classes - 0.5, by = 2),
  xmax = seq(1.5, n_classes + 0.5, by = 2)
)
# Ensure we don't exceed the number of classes
shade_data <- shade_data[shade_data$xmin < n_classes + 0.5, ]
shade_data$xmax <- pmin(shade_data$xmax, n_classes + 0.5)

# Create vertical jitter plot with coord_flip (like workflow 03)
# Use dodged positions for side-by-side groups (3 groups), 4 colors
# Duplicated points (originally "both") are purple
p_jitter_vertical <- ggplot(jitter_data_v3,
                            aes(x = source, y = testes_vs_ova,
                                fill = color_group, color = color_group,
                                group = zonation_direction)) +
  # Add alternating background shading (before points so it's behind)
  geom_rect(data = shade_data,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "grey95", alpha = 0.8) +
  geom_point(
    aes(size = point_size),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.6
  ) +
  scale_size_identity() +  # Use actual size values, no legend
  scale_fill_manual(
    name   = "Zonation\nAssociation",
    values = c("periportal" = "lightgreen", "no_network_data" = "grey70",
               "pericentral" = "orange", "both" = "mediumpurple"),
    labels = c("Periportal only", "No network data", "Pericentral only", "Both zones")
  ) +
  scale_color_manual(
    name   = "Zonation\nAssociation",
    values = c("periportal" = "darkgreen", "no_network_data" = "grey50",
               "pericentral" = "darkorange", "both" = "purple4"),
    labels = c("Periportal only", "No network data", "Pericentral only", "Both zones")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(hjust = 1),
    legend.position = "right",
    panel.grid.major.y = element_blank(),  # Remove horizontal grid lines (cleaner with shading)
    plot.margin = margin(10, 10, 10, 10, unit = "pt")
  ) +
  labs(x = "Lipid Class", y = "t-value (Testes vs Ovaries)",
       title = paste0("Lipid zonation associations (threshold = ", optimal_threshold, ")"))

ggsave(file.path(zonation_dir, "lipid_class_zonation_jitter_vertical.pdf"),
       plot = p_jitter_vertical,
       width = 5, height = 8)

cat("  ✓ Saved standalone vertical jitter plot:",
    file.path(zonation_dir, "lipid_class_zonation_jitter_vertical.pdf"), "\n\n")

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

cat("✓ Lipid-gene association visualization complete\n\n")

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
cat("  - ", file.path(zonation_dir, "lipid_class_zonation_association_heatmap*.pdf"), " (5 versions)\n")
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
