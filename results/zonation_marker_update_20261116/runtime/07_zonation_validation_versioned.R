################################################################################
# Workflow 7: Validation of Gene-Lipid Correlations with Knockdown Experiments
################################################################################
#
# This workflow validates the gene-lipid correlations identified in the MOFA
# integration (Workflow 4) and zonation analysis (Workflow 5) using in vitro
# knockdown experiments in primary hepatocytes.
#
# Author: Refactored codebase
# Date: 2025
#
# Replaces legacy scripts:
#   - scripts/legacy/lipidomic_gonad_validation_2.R
#   - scripts/legacy/lipidomic_gonad_sszb_validation_2.R
#
# Key biological insights:
#   - Validates that genes identified as correlated with lipids in FCGL mice
#     actually regulate those lipids when knocked down
#   - Tests both global hepatocyte correlations and zonation-specific effects
#   - Quantifies consistency between prediction (correlation) and validation (knockdown)
#
# Genes validated:
#   - Global hepatocyte network: Cyp2d9, Cyp2b9, Cyp2c38
#   - Zonation-specific network: Sult3a1, Mtmr11
#
################################################################################

# Clear workspace
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(limma)
  library(reshape2)
  library(decoupleR)
  library(grid)
  library(gridExtra)
  library(scales)
  library(cowplot)
  library(patchwork)
  library(broom)
})

################################################################################
# 1. CONFIGURATION
################################################################################

cat(strrep("=", 80), "\n")
cat("WORKFLOW 7: VALIDATION OF GENE-LIPID CORRELATIONS\n")
cat(strrep("=", 80), "\n\n")

# Load configuration
source("scripts/config/analysis_config.R")
source("scripts/config/datasets.R")

# Load core functions
source("scripts/support_functions.R")

# Define output directory
validation_dir <- "/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation_marker_update_20261116/fcg07"
if (!dir.exists(validation_dir)) {
  dir.create(validation_dir, recursive = TRUE)
}

################################################################################
t_table <- read.csv("results/validation/knockdown_t_table.csv",check.names=FALSE)
rownames(t_table) <- t_table$ID
comparisons <- setNames(as.list(rep(NA,ncol(t_table)-1L)),names(t_table)[-1L])
available_global <- character(0)
mean_cl_normal <- zmu_mean_cl_normal
global_xlim <- c(-1.1,1.1); global_ylim <- c(-25,25)
# 8. VALIDATION: ZONATION-SPECIFIC NETWORK
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("8. VALIDATION: ZONATION-SPECIFIC NETWORK\n")
cat(strrep("-", 80), "\n\n")

cat("Validating Sult3a1, Mtmr11 against zonation-specific gene-lipid correlations...\n")
cat("  These genes show spatially biased expression (periportal vs pericentral)\n\n")

# Load zonation-specific networks
cat("Loading zonation-specific networks:\n")

network_files <- list(
  testes_portal = "/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation_marker_update_20261116/fcg05/networks/cormat_testes_portal_full_SIF.csv",
  testes_central = "/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation_marker_update_20261116/fcg05/networks/cormat_testes_central_full_SIF.csv",
  ova_portal = "/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation_marker_update_20261116/fcg05/networks/cormat_Ovaries_portal_full_SIF.csv",
  ova_central = "/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation_marker_update_20261116/fcg05/networks/cormat_ova_central_full_SIF.csv"
)

# Check which files exist
available_networks <- list()
for (net_name in names(network_files)) {
  if (file.exists(network_files[[net_name]])) {
    cat("  Loading", net_name, "...\n")
    net <- read_csv(network_files[[net_name]], show_col_types = FALSE)

    # Filter to lipid targets
    net <- net[grepl("[(]", net$node2), ]

    # Add network ID to gene names
    net$node1 <- paste(net$node1, net_name, sep = "_")

    available_networks[[net_name]] <- net
  } else {
    cat("  WARNING: File not found:", network_files[[net_name]], "\n")
  }
}

if (length(available_networks) == 0) {
  cat("\nWARNING: No zonation-specific networks found. Skipping this section.\n")
} else {
  # Combine all networks
  sszb_network <- do.call(rbind, available_networks)
  names(sszb_network) <- c("source", "target", "mor")

  cat("\n  Combined network:", nrow(sszb_network), "edges\n\n")

  # Run decoupleR validation
  cat("Running decoupleR ULM for zonation-specific validation...\n")

  sszb_validation_scores <- decoupleR::run_ulm(
    mat = as.matrix(t_table[, -1]),
    network = sszb_network,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = 2
  )

  # Parse gonad/zonation info from source names
  sszb_validation_scores$gonade_zonation <- sub("^[^_]*_", "", sszb_validation_scores$source)
  sszb_validation_scores$source <- tolower(gsub("_.*", "", sszb_validation_scores$source))

  # Filter to on-target validation
  validation_set_only_sszb <- sszb_validation_scores %>%
    dplyr::filter(source == tolower(condition))

  cat("\nZonation-specific validation scores (on-target):\n")
  print(validation_set_only_sszb %>%
          dplyr::select(source, condition, gonade_zonation, score, p_value) %>%
          dplyr::arrange(source, dplyr::desc(abs(score))))

  # Save validation scores
  write.csv(
    sszb_validation_scores,
    file = file.path(validation_dir, "zonation_network_validation_scores.csv"),
    row.names = FALSE
  )

  ############################################################################
  # 9. SCATTER PLOTS: ZONATION-SPECIFIC VALIDATION
  ############################################################################

  cat("\n")
  cat(strrep("-", 80), "\n")
  cat("9. SCATTER PLOTS: ZONATION-SPECIFIC VALIDATION\n")
  cat(strrep("-", 80), "\n\n")

  cat("Creating scatter plots for zonation-biased genes...\n")

  # Identify which knockdowns to validate against zonation network
  zonation_genes <- c("sult3a1", "mtmr11")
  available_zonation <- intersect(tolower(names(comparisons)), zonation_genes)

  # Curated labels for zonation genes (from manuscript figures)
  zonation_curated_labels <- list(
    "mtmr11" = c("BMP(22:6_22:6)", "PC(18:0_18:1)", "PI(36:1)",
                  "PS(18:0_18:1)", "PI(18:0_18:1)")
  )

  for (gene_name in available_zonation) {
    # Assign point color based on zonation x gonad bias
    # mtmr11 = testes periportal, sult3a1 = ovaries pericentral
    point_color <- switch(gene_name,
      "mtmr11"  = "#7ba7ff",   # Testes Periportal
      "sult3a1" = "#bf2e2e",   # Ovaries Pericentral
      "black"
    )
    cat("  Plotting", gene_name, "...\n")

    # Find column index in t_table
    KO_index <- which(tolower(names(t_table)) == tolower(gene_name))

    if (length(KO_index) == 0) next

    # Identify best zonation/gonad context for this gene
    best_context <- validation_set_only_sszb %>%
      dplyr::filter(tolower(source) == tolower(gene_name)) %>%
      dplyr::arrange(dplyr::desc(abs(score))) %>%
      dplyr::slice(1) %>%
      dplyr::pull(gonade_zonation)

    if (length(best_context) == 0) {
      cat("    WARNING: No zonation context found for", gene_name, "\n")
      next
    }

    cat("    Best context:", best_context, "\n")

    # Filter network to this context
    sszb_network_sub <- sszb_network %>%
      dplyr::filter(grepl(best_context, source))

    # Merge with knockdown data
    lipid_scoring_Scatter <- merge(
      t_table[, c(1, KO_index)],
      sszb_network_sub,
      by.x = "ID",
      by.y = "target"
    )

    # Clean up source names
    lipid_scoring_Scatter$source <- gsub("_.*", "", lipid_scoring_Scatter$source)

    # Filter to on-target predictions
    lipid_scoring_Scatter$mor <- ifelse(
      tolower(lipid_scoring_Scatter$source) == names(t_table)[KO_index],
      lipid_scoring_Scatter$mor,
      0
    )

    names(lipid_scoring_Scatter) <- c("ID", "KO_tval", "KO_target", "predicted_association_weight")

    # Calculate consistency
    lipid_scoring_Scatter$consistency <-
      lipid_scoring_Scatter$KO_tval * lipid_scoring_Scatter$predicted_association_weight * -1

    # Rescale consistency to 0-1 range (0.5 = no consistency, symmetric around 0.5)
    max_abs_consistency <- max(abs(lipid_scoring_Scatter$consistency), na.rm = TRUE)
    lipid_scoring_Scatter$consistency_scaled <- 0.5 + lipid_scoring_Scatter$consistency / (2 * max_abs_consistency)

    lipid_scoring_Scatter <- lipid_scoring_Scatter %>%
      dplyr::arrange(dplyr::desc(consistency))

    # Use curated labels if available, otherwise label top 20
    if (gene_name %in% names(zonation_curated_labels)) {
      lipid_scoring_Scatter$label <- ifelse(
        lipid_scoring_Scatter$ID %in% zonation_curated_labels[[gene_name]],
        lipid_scoring_Scatter$ID,
        NA
      )
    } else {
      lipid_scoring_Scatter$label <- c(
        lipid_scoring_Scatter$ID[1:min(20, nrow(lipid_scoring_Scatter))],
        rep(NA, max(0, nrow(lipid_scoring_Scatter) - 20))
      )
      lipid_scoring_Scatter$label <- ifelse(
        lipid_scoring_Scatter$predicted_association_weight != 0,
        lipid_scoring_Scatter$label,
        NA
      )
    }

    # Format context for plot
    context_label <- paste0("(", gsub("_", " ", best_context), " specific)")

    # Create scatter plot with flat color based on zonation x gonad bias
    p <- ggplot(
      lipid_scoring_Scatter,
      aes(x = predicted_association_weight, y = KO_tval, label = label)
    ) +
      geom_smooth(
        method      = "lm",
        formula     = y ~ x
      ) +
      geom_point(colour = point_color, size = 2) +
      stat_summary(fun.data = mean_cl_normal, colour = point_color) +
      coord_cartesian(xlim = global_xlim, ylim = global_ylim) +
      ggtitle(paste0(toupper(gene_name), " Knockdown Validation\n", context_label)) +
      xlab("Predicted Association\n(Zonation-specific correlation in FCGL mice)") +
      ylab("Observed Lipid Change\n(t-statistic in knockdown experiment)") +
      geom_text_repel(max.overlaps = 20, size = 2.5,
                      segment.color = "grey30", segment.size = 0.4,
                      min.segment.length = 0, box.padding = 0.5, point.padding = 0.3,
                      arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold")
      )

    ggsave(
      filename = file.path(validation_dir, paste0(gene_name, "_zonation_validation_scatter.pdf")),
      plot = p,
      width = 6,
      height = 6
    )
  }
}
  cat("\nZonation-specific scatter plots saved to: results/validation/\n")

################################################################################
# 11. STRUCTURAL ANALYSIS OF KO RESPONSE CONSISTENCY
################################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("11. STRUCTURAL ANALYSIS OF KO RESPONSE CONSISTENCY\n")
cat(strrep("-", 80), "\n\n")

cat("Analyzing links between KO response consistency and lipid structural features.\n")
cat("Consistency = predicted_association (correlation) x KO t-value x (-1)\n")
cat("Positive consistency = KO result matches FCGL prediction.\n\n")

# Load structural analysis functions for lipid name parsing
source("scripts/core/structural_analysis.R")

# Focus lipid classes (matching structural analysis workflow 06)
focus_classes <- c("SM", "BMP", "PG", "TG")

# Output directory
structural_consistency_dir <- file.path(validation_dir, "structural_consistency")
if (!dir.exists(structural_consistency_dir)) {
  dir.create(structural_consistency_dir, recursive = TRUE)
}

#-------------------------------------------------------------------------------
# Helper: create 3-panel structural consistency plot for one gene x class subset
#-------------------------------------------------------------------------------
create_consistency_structural_plot <- function(plot_data, cls_label, gene_label,
                                               chain_limits, db_limits,
                                               consistency_limits) {

  n_pts <- nrow(plot_data)

  # Run regression: consistency ~ chain_length + db_position + chain_length:db_position
  chain_p <- NA; db_p <- NA; interaction_p <- NA; r_sq <- NA
  if (n_pts >= 10) {
    model <- tryCatch({
      lm(consistency ~ chain_length + db_position + chain_length:db_position,
         data = plot_data)
    }, error = function(e) NULL)

    if (!is.null(model)) {
      model_summary <- summary(model)
      coefs <- broom::tidy(model)
      chain_row <- coefs[coefs$term == "chain_length", ]
      db_row <- coefs[coefs$term == "db_position", ]
      int_row <- coefs[coefs$term == "chain_length:db_position", ]

      if (nrow(chain_row) > 0) chain_p <- chain_row$p.value
      if (nrow(db_row) > 0) db_p <- db_row$p.value
      if (nrow(int_row) > 0) interaction_p <- int_row$p.value
      r_sq <- round(model_summary$r.squared, 3)
    }
  }

  # Format p-values
  format_pval <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return(sprintf("%.2e", p))
    if (p < 0.01) return(sprintf("%.3f", p))
    return(sprintf("%.2f", p))
  }

  chain_sig <- !is.na(chain_p) && chain_p < 0.05
  db_sig <- !is.na(db_p) && db_p < 0.05
  interaction_sig <- !is.na(interaction_p) && interaction_p < 0.05

  chain_color <- if(chain_sig) "#00BA38" else "grey50"
  db_color <- if(db_sig) "#00BA38" else "grey50"
  interaction_color <- if(interaction_sig) "#00BA38" else "grey50"

  # Panel 1: Chain Length vs Consistency
  p1 <- ggplot(plot_data, aes(x = chain_length, y = consistency)) +
    geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB", shape = 15) +
    geom_smooth(method = "lm", se = TRUE, color = chain_color, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    coord_cartesian(xlim = chain_limits, ylim = consistency_limits) +
    labs(
      title = paste0(cls_label, " (n=", n_pts, ", R2=", r_sq, ")"),
      subtitle = paste0("Chain Length: p = ", format_pval(chain_p)),
      x = "Chain Length (C)",
      y = paste0("KO Response Consistency\n(", gene_label, ")")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = chain_color, face = "bold", size = 11),
      axis.title.x = element_text(color = chain_color,
                                   face = if(chain_sig) "bold" else "plain")
    )

  # Panel 2: Double Bonds vs Consistency
  p2 <- ggplot(plot_data, aes(x = db_position, y = consistency)) +
    geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB", shape = 15) +
    geom_smooth(method = "lm", se = TRUE, color = db_color, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    coord_cartesian(xlim = db_limits, ylim = consistency_limits) +
    labs(
      title = "",
      subtitle = paste0("Double Bonds: p = ", format_pval(db_p)),
      x = "Double Bonds (DB)",
      y = paste0("KO Response Consistency\n(", gene_label, ")")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.subtitle = element_text(color = db_color, face = "bold", size = 11),
      axis.title.x = element_text(color = db_color,
                                   face = if(db_sig) "bold" else "plain")
    )

  # Panel 3: C vs DB colored by Consistency (2D interaction plot)
  p3 <- ggplot(plot_data, aes(x = chain_length, y = db_position, color = consistency)) +
    geom_point(size = 3, alpha = 1, shape = 15) +
    scale_color_gradient2(low = "#e54141", mid = "white", high = "#4173d8",
                          midpoint = 0, name = "Consistency",
                          limits = c(-1, 1),
                          oob = scales::squish) +
    coord_cartesian(xlim = chain_limits, ylim = db_limits) +
    labs(
      title = "",
      subtitle = paste0("C x DB Interaction: p = ", format_pval(interaction_p)),
      x = "Chain Length (C)",
      y = "Double Bonds (DB)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.subtitle = element_text(color = interaction_color, face = "bold", size = 11),
      legend.position = "right"
    )

  # Combine
  combined <- p1 + p2 + p3 + plot_layout(ncol = 3)
  return(list(
    plot = combined,
    stats = data.frame(
      gene = gene_label,
      lipid_class = cls_label,
      n_lipids = n_pts,
      r_squared = r_sq,
      chain_length_p = chain_p,
      db_position_p = db_p,
      interaction_p = interaction_p,
      chain_length_sig = chain_sig,
      db_position_sig = db_sig,
      interaction_sig = interaction_sig,
      stringsAsFactors = FALSE
    )
  ))
}

#-------------------------------------------------------------------------------
# Prepare consistency data for each validated gene
#-------------------------------------------------------------------------------
cat("Preparing consistency data for each validated gene...\n\n")

all_gene_data <- list()
all_regression_stats <- list()

# Global network genes
for (gene_name in available_global) {
  KO_index <- which(tolower(names(t_table)) == tolower(gene_name))
  if (length(KO_index) == 0) next

  # Filter network to this gene's edges only
  gene_edges <- lipid_network[tolower(lipid_network$source) == gene_name, ]
  if (nrow(gene_edges) == 0) {
    cat("  ", gene_name, ": no network edges found, skipping\n")
    next
  }

  # Merge with KO t-values
  merged <- merge(
    t_table[, c(1, KO_index)],
    gene_edges[, c("target", "mor")],
    by.x = "ID", by.y = "target"
  )
  colnames(merged) <- c("ID", "KO_tval", "predicted_association_weight")

  # Compute consistency = correlation x t-value x (-1)
  merged$consistency <- merged$KO_tval * merged$predicted_association_weight * -1

  # Parse structural features from lipid names
  parsed <- parse_lipid_names_to_df(merged$ID, verbose = FALSE)
  merged$lipid_class <- parsed$lipid_class
  merged$chain_length <- parsed$chain_length
  merged$db_position <- parsed$db_position

  # Remove lipids with missing structural data
  merged <- merged[!is.na(merged$chain_length) & !is.na(merged$db_position), ]

  all_gene_data[[gene_name]] <- merged
  cat("  ", gene_name, "(global):", nrow(merged), "lipids with structural data\n")
}

# Zonation network genes
if (exists("sszb_network") && exists("available_zonation") && exists("validation_set_only_sszb")) {
  for (gene_name in available_zonation) {
    KO_index <- which(tolower(names(t_table)) == tolower(gene_name))
    if (length(KO_index) == 0) next

    # Find best zonation context
    best_context <- validation_set_only_sszb %>%
      dplyr::filter(tolower(source) == tolower(gene_name)) %>%
      dplyr::arrange(dplyr::desc(abs(score))) %>%
      dplyr::slice(1) %>%
      dplyr::pull(gonade_zonation)

    if (length(best_context) == 0) {
      cat("  ", gene_name, ": no zonation context found, skipping\n")
      next
    }

    # Filter network to this gene in its best context
    sszb_sub <- sszb_network[grepl(best_context, sszb_network$source), ]
    sszb_sub$source_clean <- tolower(gsub("_.*", "", sszb_sub$source))
    gene_edges <- sszb_sub[sszb_sub$source_clean == gene_name, ]

    if (nrow(gene_edges) == 0) {
      cat("  ", gene_name, ": no network edges in context", best_context, "\n")
      next
    }

    # Merge with KO t-values
    merged <- merge(
      t_table[, c(1, KO_index)],
      gene_edges[, c("target", "mor")],
      by.x = "ID", by.y = "target"
    )
    colnames(merged) <- c("ID", "KO_tval", "predicted_association_weight")

    # Compute consistency
    merged$consistency <- merged$KO_tval * merged$predicted_association_weight * -1

    # Parse structural features
    parsed <- parse_lipid_names_to_df(merged$ID, verbose = FALSE)
    merged$lipid_class <- parsed$lipid_class
    merged$chain_length <- parsed$chain_length
    merged$db_position <- parsed$db_position

    # Remove lipids with missing structural data
    merged <- merged[!is.na(merged$chain_length) & !is.na(merged$db_position), ]

    all_gene_data[[gene_name]] <- merged
    cat("  ", gene_name, "(zonation,", best_context, "):", nrow(merged),
        "lipids with structural data\n")
  }
}

if (length(all_gene_data) == 0) {
  cat("\nWARNING: No gene data available for structural consistency analysis.\n")
} else {

  #-----------------------------------------------------------------------------
  # Use the SAME axis limits as workflow 06 structural analysis plots
  # Computed from lipid_structure filtered to focus classes (SM, BMP, PG, TG)
  #-----------------------------------------------------------------------------
  cat("Loading axis limits from structural analysis (workflow 06)...\n")
  lipid_structure_ref <- read.csv("results/lipidomic/structural_model/lipids_with_structure.csv")
  lipid_structure_ref <- lipid_structure_ref[
    lipid_structure_ref$lipid_class %in% focus_classes &
    !is.na(lipid_structure_ref$chain_length) &
    !is.na(lipid_structure_ref$db_position) &
    !is.na(lipid_structure_ref$logfc), ]

  ref_chain_range <- range(lipid_structure_ref$chain_length)
  ref_chain_pad <- diff(ref_chain_range) * 0.05
  global_chain_lim <- c(ref_chain_range[1] - ref_chain_pad,
                         ref_chain_range[2] + ref_chain_pad)

  ref_db_range <- range(lipid_structure_ref$db_position)
  ref_db_pad <- diff(ref_db_range) * 0.05
  global_db_lim <- c(ref_db_range[1] - ref_db_pad,
                      ref_db_range[2] + ref_db_pad)

  rm(lipid_structure_ref)

  cat("\n  Global axis limits (matching workflow 06):\n")
  cat("    Chain Length:", round(global_chain_lim[1], 1), "to",
      round(global_chain_lim[2], 1), "\n")
  cat("    Double Bonds:", round(global_db_lim[1], 1), "to",
      round(global_db_lim[2], 1), "\n\n")

  #-----------------------------------------------------------------------------
  # Generate plots for each gene
  #-----------------------------------------------------------------------------
  cat("Generating structural consistency plots...\n\n")

  for (gene_name in names(all_gene_data)) {
    gene_data <- all_gene_data[[gene_name]]
    gene_label <- toupper(gene_name)

    # Per-gene consistency range (for y-axis and color scale)
    gene_consistency_range <- range(gene_data$consistency, na.rm = TRUE)
    cons_padding <- diff(gene_consistency_range) * 0.05
    gene_consistency_lim <- c(gene_consistency_range[1] - cons_padding,
                               gene_consistency_range[2] + cons_padding)

    cat("  ", gene_label, ":\n")

    plot_list <- list()
    stats_list <- list()

    for (cls in c(focus_classes, "All Lipids")) {
      if (cls == "All Lipids") {
        class_data <- gene_data
      } else {
        class_data <- gene_data[gene_data$lipid_class == cls, ]
      }

      if (nrow(class_data) < 3) {
        cat("    ", cls, ": n =", nrow(class_data), "(skipped, < 3 lipids)\n")
        next
      }

      result <- create_consistency_structural_plot(
        plot_data = class_data,
        cls_label = cls,
        gene_label = gene_label,
        chain_limits = global_chain_lim,
        db_limits = global_db_lim,
        consistency_limits = gene_consistency_lim
      )

      plot_list[[cls]] <- result$plot
      stats_list[[cls]] <- result$stats

      sig_terms <- c()
      if (result$stats$chain_length_sig) sig_terms <- c(sig_terms, "C")
      if (result$stats$db_position_sig) sig_terms <- c(sig_terms, "DB")
      if (result$stats$interaction_sig) sig_terms <- c(sig_terms, "C:DB")

      sig_text <- if (length(sig_terms) > 0) paste0(" [sig: ", paste(sig_terms, collapse = ", "), "]") else ""
      cat("    ", cls, ": n =", result$stats$n_lipids,
          ", R2 =", result$stats$r_squared, sig_text, "\n")
    }

    # Collect regression stats
    if (length(stats_list) > 0) {
      all_regression_stats[[gene_name]] <- do.call(rbind, stats_list)
    }

    # Save combined PDF (one page per gene, rows = lipid classes)
    if (length(plot_list) > 0) {
      n_rows <- length(plot_list)

      pdf(file.path(structural_consistency_dir,
                    paste0(gene_name, "_structural_consistency.pdf")),
          width = 14, height = 5 * n_rows)

      combined_plot <- wrap_plots(plot_list, ncol = 1) +
        plot_annotation(
          title = paste0(gene_label, " KO Response Consistency vs Lipid Structure"),
          subtitle = paste("Consistency = correlation (FCGL network) x KO t-value x (-1).",
                           "Green = significant at p < 0.05",
                           "(model: consistency ~ C + DB + C:DB)"),
          theme = theme(
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 10, color = "grey30")
          )
        )

      print(combined_plot)
      dev.off()

      cat("    Saved:", gene_name, "_structural_consistency.pdf\n\n")
    }
  }

  #-----------------------------------------------------------------------------
  # Save regression summary across all genes and classes
  #-----------------------------------------------------------------------------
  if (length(all_regression_stats) > 0) {
    regression_summary <- do.call(rbind, all_regression_stats)
    rownames(regression_summary) <- NULL

    write.csv(
      regression_summary,
      file = file.path(structural_consistency_dir,
                       "consistency_structural_regression_summary.csv"),
      row.names = FALSE
    )

    cat("  Regression summary saved:\n")
    cat("    results/validation/structural_consistency/consistency_structural_regression_summary.csv\n")

    # Print significant findings
    sig_findings <- regression_summary[
      regression_summary$chain_length_sig |
      regression_summary$db_position_sig |
      regression_summary$interaction_sig, ]

    if (nrow(sig_findings) > 0) {
      cat("\n  Significant structural associations found:\n")
      for (i in seq_len(nrow(sig_findings))) {
        terms <- c()
        if (sig_findings$chain_length_sig[i]) {
          terms <- c(terms, paste0("C (p=", sprintf("%.3f", sig_findings$chain_length_p[i]), ")"))
        }
        if (sig_findings$db_position_sig[i]) {
          terms <- c(terms, paste0("DB (p=", sprintf("%.3f", sig_findings$db_position_p[i]), ")"))
        }
        if (sig_findings$interaction_sig[i]) {
          terms <- c(terms, paste0("C:DB (p=", sprintf("%.3f", sig_findings$interaction_p[i]), ")"))
        }
        cat("    ", sig_findings$gene[i], "-", sig_findings$lipid_class[i], ":",
            paste(terms, collapse = ", "), "\n")
      }
    } else {
      cat("\n  No significant structural associations found at p < 0.05.\n")
    }
  }
}

cat("\nStructural consistency analysis complete.\n")
cat("Output directory: results/validation/structural_consistency/\n")
