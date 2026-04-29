## Quick analysis: Correlation between lipid insaturation (double bonds) and
## gonadal t-values, plus stress pathway extraction
## Author: Claude Code analysis
## Date: 2026-03-20

library(ggplot2)
library(dplyr)
library(ggrepel)

cat("=== PART 1: Insaturation vs gonadal t-values ===\n\n")

# --- Load gonadal ttop (has t-values per lipid) ---
ttop <- read.csv("results/lipidomic/gonade_type_ttop.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(ttop), "lipids from gonadal ttop\n")

# --- Load all bubble files to get DB (double bond) info per lipid ---
bubble_dir <- "data/metabolomics/for_bubbles_with_names/"
bubble_files <- list.files(bubble_dir, pattern = "_with_name\\.csv$", full.names = TRUE)

all_struct <- do.call(rbind, lapply(bubble_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  class_name <- gsub("_with_name\\.csv$", "", basename(f))
  df$lipid_class <- class_name
  df
}))
cat("Loaded structural data for", nrow(all_struct), "lipid entries\n")
cat("Columns:", paste(colnames(all_struct), collapse = ", "), "\n")

# --- Merge ttop with structural data ---
# ttop has column "ID" with lipid names, bubble data has "lipid"
merged <- merge(ttop, all_struct, by.x = "ID", by.y = "lipid", all.x = FALSE, all.y = FALSE)
cat("Merged:", nrow(merged), "lipids with both t-values and structural data\n\n")

# --- Global correlation: DB vs t-value ---
cor_test <- cor.test(merged$DB, merged$t, method = "pearson")
cat("--- Global correlation: Double bonds (DB) vs t-value ---\n")
cat("  Pearson r =", round(cor_test$estimate, 4), "\n")
cat("  p-value =", format(cor_test$p.value, digits = 4), "\n")
cat("  95% CI: [", round(cor_test$conf.int[1], 4), ",", round(cor_test$conf.int[2], 4), "]\n\n")

# Also Spearman (more robust)
cor_spearman <- cor.test(merged$DB, merged$t, method = "spearman")
cat("  Spearman rho =", round(cor_spearman$estimate, 4), "\n")
cat("  p-value =", format(cor_spearman$p.value, digits = 4), "\n\n")

# --- Linear model: t ~ DB ---
lm_global <- lm(t ~ DB, data = merged)
cat("Linear model (t ~ DB):\n")
print(summary(lm_global))
cat("\n")

# --- Linear model controlling for lipid class: t ~ DB + lipid_class ---
lm_class <- lm(t ~ DB + lipid_class, data = merged)
cat("Linear model (t ~ DB + lipid_class):\n")
cat("  DB coefficient:", round(coef(lm_class)["DB"], 4), "\n")
cat("  DB p-value:", format(summary(lm_class)$coefficients["DB", "Pr(>|t|)"], digits = 4), "\n\n")

# --- Also control for chain length ---
lm_full <- lm(t ~ DB + C + lipid_class, data = merged)
cat("Linear model (t ~ DB + C + lipid_class):\n")
cat("  DB coefficient:", round(coef(lm_full)["DB"], 4), "\n")
cat("  DB p-value:", format(summary(lm_full)$coefficients["DB", "Pr(>|t|)"], digits = 4), "\n")
cat("  C coefficient:", round(coef(lm_full)["C"], 4), "\n")
cat("  C p-value:", format(summary(lm_full)$coefficients["C", "Pr(>|t|)"], digits = 4), "\n\n")

# --- Per-class correlation ---
cat("--- Per-class correlations: DB vs t-value ---\n")
class_cors <- merged %>%
  group_by(lipid_class) %>%
  dplyr::filter(n() >= 10) %>%
  dplyr::summarise(
    n = n(),
    pearson_r = cor(DB, t, method = "pearson"),
    p_value = cor.test(DB, t, method = "pearson")$p.value,
    mean_t = mean(t),
    .groups = "drop"
  ) %>%
  dplyr::arrange(p_value)

print(as.data.frame(class_cors), digits = 3)
cat("\n")

# --- Plot 1: Global scatter DB vs t-value ---
p1 <- ggplot(merged, aes(x = DB, y = t)) +
  geom_point(aes(color = lipid_class), alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 1) +
  labs(
    x = "Number of double bonds (insaturations)",
    y = "t-value (Testes vs Ovaries)",
    title = paste0("Lipid insaturation vs gonadal effect (r = ",
                   round(cor_test$estimate, 3), ", p = ",
                   format(cor_test$p.value, digits = 3), ")"),
    subtitle = "Positive t = male-biased, Negative t = female-biased"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7))

ggsave("results/lipidomic/insaturation_vs_tvalue_global.pdf", p1, width = 10, height = 7)
cat("Saved: results/lipidomic/insaturation_vs_tvalue_global.pdf\n")

# --- Plot 2: Per-class faceted ---
# Focus on classes with >= 20 lipids
big_classes <- merged %>%
  group_by(lipid_class) %>%
  dplyr::filter(n() >= 20) %>%
  ungroup()

# Add correlation annotation per class
class_labels <- big_classes %>%
  group_by(lipid_class) %>%
  dplyr::summarise(
    r = round(cor(DB, t), 3),
    p = cor.test(DB, t)$p.value,
    label = paste0("r=", round(cor(DB, t), 2),
                   ifelse(cor.test(DB, t)$p.value < 0.05, "*", "")),
    .groups = "drop"
  )

p2 <- ggplot(big_classes, aes(x = DB, y = t)) +
  geom_point(alpha = 0.5, size = 1.2, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = class_labels, aes(label = label),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3, color = "red") +
  facet_wrap(~ lipid_class, scales = "free_x", ncol = 5) +
  labs(
    x = "Number of double bonds",
    y = "t-value (Testes vs Ovaries)",
    title = "Insaturation vs gonadal effect by lipid class"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

ggsave("results/lipidomic/insaturation_vs_tvalue_by_class.pdf", p2, width = 14, height = 10)
cat("Saved: results/lipidomic/insaturation_vs_tvalue_by_class.pdf\n\n")


cat("=== PART 2: Stress-related pathways in hepatocyte transcriptomics ===\n\n")

# --- Load pathway enrichment results ---
load("results/pseudobulk/pathways/pathway_scores_df_all.RData")
cat("Loaded pathway scores. Object name(s):", ls(pattern = "pathway|score|df"), "\n")

# The object should be pathway_scores_df_all or similar
# Let's check what we have
for (obj_name in ls(pattern = "pathway|score|df")) {
  obj <- get(obj_name)
  cat("Object '", obj_name, "':", class(obj), "\n")
  if (is.data.frame(obj)) {
    cat("  Dims:", nrow(obj), "x", ncol(obj), "\n")
    cat("  Columns:", paste(head(colnames(obj), 10), collapse = ", "), "\n")
    cat("  First few rows:\n")
    print(head(obj, 3))
  }
}

# --- Search for stress-related pathways ---
# Define stress-related keywords
stress_keywords <- c(
  "OXIDAT", "REACTIVE_OXYGEN", "PEROXID", "RADICAL",
  "ANTIOXIDANT", "GLUTATHIONE", "THIOREDOXIN",
  "STRESS", "UNFOLDED_PROTEIN", "HEAT_SHOCK",
  "APOPTOT", "CELL_DEATH", "NECROPT",
  "INFLAMMAT", "NF.?KB", "CYTOKINE",
  "DNA_DAMAGE", "DNA_REPAIR",
  "XENOBIOTIC", "DETOXIF",
  "LIPID_PEROXI", "FERROPTOSIS",
  "AUTOPHAGY",
  "TUMOR", "CARCINO", "NEOPLAS"
)

# Try to find the right data structure
# Pathway scores might be in wide format (pathways as rows, contrasts as columns)
# or long format. Let's handle both.

# Get the pathway data
pathway_df <- NULL
for (obj_name in ls(pattern = "pathway|score")) {
  obj <- get(obj_name)
  if (is.data.frame(obj) && nrow(obj) > 100) {
    pathway_df <- obj
    cat("\nUsing object:", obj_name, "\n")
    break
  }
}

if (!is.null(pathway_df)) {
  cat("Pathway data dimensions:", nrow(pathway_df), "x", ncol(pathway_df), "\n")
  cat("Column names:", paste(colnames(pathway_df), collapse = ", "), "\n\n")

  # Check if it's in wide or long format
  # If wide: rows = pathways, columns = contrasts/cell types
  # If long: has source/condition/score columns

  if ("source" %in% colnames(pathway_df)) {
    # Long format from decoupleR
    cat("Format: long (decoupleR output)\n")

    # Filter for hepatocytes and gonadal contrast
    hepa_stress <- pathway_df %>%
      dplyr::filter(grepl("Hepato|hepa", condition, ignore.case = TRUE) |
                      grepl("gonad|testes.*ova", condition, ignore.case = TRUE))

    cat("Hepatocyte entries:", nrow(hepa_stress), "\n")

    # Search stress pathways
    pattern <- paste(stress_keywords, collapse = "|")
    stress_pathways <- pathway_df %>%
      dplyr::filter(grepl(pattern, source, ignore.case = TRUE))

    cat("Stress-related pathway entries:", nrow(stress_pathways), "\n")
    print(head(stress_pathways, 20))

  } else {
    # Likely wide format - pathways as rownames or first column
    cat("Format: wide\n")
    cat("First 3 column names:", paste(head(colnames(pathway_df), 3), collapse = ", "), "\n")
    cat("Sample row names:", paste(head(rownames(pathway_df), 3), collapse = ", "), "\n")

    # Get pathway names (either rownames or first column)
    if (is.character(pathway_df[[1]])) {
      pathway_names <- pathway_df[[1]]
    } else {
      pathway_names <- rownames(pathway_df)
    }

    # Search stress pathways
    pattern <- paste(stress_keywords, collapse = "|")
    stress_idx <- grepl(pattern, pathway_names, ignore.case = TRUE)

    stress_results <- pathway_df[stress_idx, ]
    cat("\nFound", nrow(stress_results), "stress-related pathways\n\n")

    # Print hepatocyte-relevant columns if they exist
    hepa_cols <- grep("Hepato|hepa|Testes.*Ova|gonad", colnames(pathway_df),
                      ignore.case = TRUE, value = TRUE)
    cat("Hepatocyte/gonadal columns:", paste(hepa_cols, collapse = ", "), "\n\n")

    if (length(hepa_cols) > 0) {
      cat("--- Stress pathways in hepatocytes (gonadal contrast) ---\n")
      # Get the gonadal contrast for hepatocytes
      tvo_col <- grep("Testes.*Ova.*Hepato|Hepato.*Testes|gonad.*Hepato",
                      hepa_cols, ignore.case = TRUE, value = TRUE)

      if (length(tvo_col) == 0) {
        tvo_col <- hepa_cols[1]
      }

      cat("Using column(s):", paste(tvo_col, collapse = ", "), "\n\n")

      # Sort by score
      if (is.character(pathway_df[[1]])) {
        stress_out <- data.frame(
          pathway = stress_results[[1]],
          score = as.numeric(stress_results[[tvo_col[1]]]),
          stringsAsFactors = FALSE
        )
      } else {
        stress_out <- data.frame(
          pathway = rownames(stress_results),
          score = as.numeric(stress_results[[tvo_col[1]]]),
          stringsAsFactors = FALSE
        )
      }
      stress_out <- stress_out[order(-abs(stress_out$score)), ]

      cat("Top 40 stress-related pathways (sorted by |score|):\n")
      print(head(stress_out, 40), row.names = FALSE)
    }
  }
} else {
  cat("Could not find pathway data frame. Checking RData contents...\n")
  cat("All objects in environment:", paste(ls(), collapse = ", "), "\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
