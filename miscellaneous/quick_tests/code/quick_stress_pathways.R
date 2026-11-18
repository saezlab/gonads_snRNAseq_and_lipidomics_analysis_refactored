## Focused extraction of stress-related pathways from hepatocyte data
library(dplyr)

load("results/pseudobulk/pathways/pathway_scores_df_all.RData")

# Check condition names
cat("Unique conditions:\n")
print(unique(pathway_results$condition))

# Filter for hepatocytes + Testes vs Ovaries (gonadal contrast)
hepa_gonadal <- pathway_results %>%
  dplyr::filter(grepl("Testes.*Ova.*Hepato|Hepato.*Testes", condition, ignore.case = TRUE))

cat("\nHepatocyte gonadal entries:", nrow(hepa_gonadal), "\n")

if (nrow(hepa_gonadal) == 0) {
  # Try broader search
  cat("Trying broader search...\n")
  hepa_conds <- unique(pathway_results$condition)[grepl("Hepato", unique(pathway_results$condition))]
  cat("Hepatocyte conditions:", paste(hepa_conds, collapse = "\n  "), "\n")

  gonadal_conds <- unique(pathway_results$condition)[grepl("Testes|Ova", unique(pathway_results$condition))]
  cat("Gonadal conditions:", paste(gonadal_conds, collapse = "\n  "), "\n")

  # Use the gonadal hepatocyte condition
  target_cond <- gonadal_conds[grepl("Hepato", gonadal_conds)]
  if (length(target_cond) == 0) {
    target_cond <- gonadal_conds[1]
  }
  cat("Using condition:", target_cond, "\n")
  hepa_gonadal <- pathway_results %>%
    dplyr::filter(condition == target_cond[1])
}

cat("Working with", nrow(hepa_gonadal), "pathways\n\n")

# Define stress-related keywords (comprehensive)
stress_keywords <- c(
  "OXIDAT", "REACTIVE_OXYGEN", "ROS", "PEROXID", "SUPEROXID", "RADICAL",
  "ANTIOXIDANT", "GLUTATHIONE", "THIOREDOXIN", "CATALASE",
  "STRESS_RESPONSE", "RESPONSE_TO_STRESS", "CELLULAR_RESPONSE_TO_STRESS",
  "UNFOLDED_PROTEIN", "HEAT_SHOCK", "ENDOPLASMIC_RETICULUM_STRESS",
  "APOPTOT", "PROGRAMMED_CELL_DEATH", "NECROPT", "FERROPTOSIS",
  "INFLAMMAT", "NF_KAPPA_B", "NFKB", "CYTOKINE",
  "DNA_DAMAGE", "DNA_REPAIR",
  "XENOBIOTIC", "DETOXIF",
  "LIPID_PEROXI", "LIPID_OXIDAT",
  "AUTOPHAGY",
  "TUMOR", "CARCINO", "NEOPLAS",
  "HYDROGEN_PEROXIDE", "NITRIC_OXIDE",
  "IRON_ION", "METAL_ION_STRESS",
  "HYPOXIA", "OXYGEN_LEVEL",
  "P53", "TP53"
)

pattern <- paste(stress_keywords, collapse = "|")

stress_hepa <- hepa_gonadal %>%
  dplyr::filter(grepl(pattern, source, ignore.case = TRUE)) %>%
  dplyr::arrange(desc(abs(score)))

cat("=== STRESS-RELATED PATHWAYS IN HEPATOCYTES (Testes vs Ovaries) ===\n")
cat("Positive score = male-biased (testes), Negative score = female-biased (ovaries)\n\n")

# Print top 60
cat("Top 60 stress pathways by |score|:\n")
stress_print <- stress_hepa %>%
  dplyr::select(source, score, p_value) %>%
  head(60)

for (i in seq_len(nrow(stress_print))) {
  row <- stress_print[i, ]
  sig <- ifelse(row$p_value < 0.001, "***",
           ifelse(row$p_value < 0.01, "**",
             ifelse(row$p_value < 0.05, "*", "")))
  direction <- ifelse(row$score > 0, "MALE", "FEMALE")
  cat(sprintf("  %6.2f %s %-6s %s\n",
              row$score, sig, direction, row$source))
}

# Now categorize by stress type
cat("\n\n=== ORGANIZED BY STRESS CATEGORY ===\n\n")

categories <- list(
  "OXIDATIVE STRESS" = "OXIDAT|REACTIVE_OXYGEN|PEROXID|SUPEROXID|RADICAL|ANTIOXIDANT|HYDROGEN_PEROXIDE",
  "GLUTATHIONE/REDOX" = "GLUTATHIONE|THIOREDOXIN",
  "ER STRESS" = "ENDOPLASMIC_RETICULUM_STRESS|UNFOLDED_PROTEIN",
  "APOPTOSIS/CELL DEATH" = "APOPTOT|PROGRAMMED_CELL_DEATH|NECROPT|FERROPTOSIS",
  "INFLAMMATION" = "INFLAMMAT|NF_KAPPA_B|NFKB|CYTOKINE",
  "DNA DAMAGE" = "DNA_DAMAGE|DNA_REPAIR|P53|TP53",
  "XENOBIOTIC/DETOX" = "XENOBIOTIC|DETOXIF",
  "HYPOXIA" = "HYPOXIA|OXYGEN_LEVEL",
  "AUTOPHAGY" = "AUTOPHAGY",
  "LIPID OXIDATION" = "LIPID_PEROXI|LIPID_OXIDAT",
  "IRON/METAL" = "IRON_ION|METAL_ION"
)

for (cat_name in names(categories)) {
  cat_pattern <- categories[[cat_name]]
  cat_results <- hepa_gonadal %>%
    dplyr::filter(grepl(cat_pattern, source, ignore.case = TRUE)) %>%
    dplyr::arrange(desc(abs(score)))

  if (nrow(cat_results) > 0) {
    cat(sprintf("--- %s (%d pathways) ---\n", cat_name, nrow(cat_results)))

    for (i in seq_len(min(nrow(cat_results), 10))) {
      row <- cat_results[i, ]
      sig <- ifelse(row$p_value < 0.001, "***",
               ifelse(row$p_value < 0.01, "**",
                 ifelse(row$p_value < 0.05, "*", "")))
      direction <- ifelse(row$score > 0, "MALE", "FEMALE")
      cat(sprintf("  %7.2f %s %-6s %s\n",
                  row$score, sig, direction, row$source))
    }
    cat("\n")
  }
}

cat("\n=== DONE ===\n")
