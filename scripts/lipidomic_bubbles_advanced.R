# new/extra packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(scatterpie)

# --- options controlling pie rendering by regulation ---
only_regulated  <- T   # <- set FALSE to show pies regardless of LogFC
positive_or_negative <- 1      # <- e.g. 1 or -1
min_abs_logfc <- 0

# --- your existing loading/merging code (with a tiny fix: by.y) ---
zonation_lipids_files <- list.files(".", recursive = TRUE) 
zonation_lipids_files <- zonation_lipids_files[grepl("results/multi-omic/cormat_.*_SIF", zonation_lipids_files)]

correlation_tables <- list()
for (filename in zonation_lipids_files) {
  correlation_tables[[filename]] <- read_csv(filename)
  correlation_tables[[filename]]$condition <- str_extract(filename, "(?<=cormat_).*(?=_SIF)")
}
correlation_table_merged <- as.data.frame(do.call(rbind, correlation_tables))

df <- read.csv("data/metabolomics/TG.csv")
df_with_names <- read.csv("data/metabolomics/TG_with_names.csv")

df_with_names <- merge(df_with_names, correlation_table_merged, by.x = "lipid", by.y = "node2")
df_with_names$CDB_ID <- paste(df_with_names$C, df_with_names$DB, sep = "_")
df$CDB_ID <- paste(df$C, df$DB, sep = "_")

# --- build presence (0/1) per condition for each CDB_ID ---
cond_levels <- c("ova_central", "Ovaries_portal", "testes_central", "testes_portal")

presence_wide <- df_with_names %>%
  distinct(CDB_ID, condition) %>%
  filter(condition %in% cond_levels) %>%
  mutate(value = 1L) %>%
  tidyr::pivot_wider(
    names_from = condition, values_from = value,
    values_fill = 0L
  )

for (cc in cond_levels) {
  if (!cc %in% names(presence_wide)) presence_wide[[cc]] <- 0L
}

# --- join onto df; compute hits and regulation status ---
df_anno <- df %>%
  left_join(presence_wide, by = "CDB_ID") %>%
  mutate(across(all_of(cond_levels), ~replace_na(., 0L))) %>%
  mutate(
    any_hit   = rowSums(across(all_of(cond_levels))) > 0,
    regulated = !is.na(LogFC) & LogFC * positive_or_negative >= 0 & LogFC != 0
  )

# choose which rows get pies
pie_data <- df_anno %>%
  filter(any_hit, if (only_regulated) regulated else TRUE)

# pie radius in data units
pie_r <- 0.6

# --- plot: background grey for all, pies for selected rows ---
p <- ggplot() +
  # baseline: light grey (covers both "no condition" and "not regulated")
  geom_point(
    data = df_anno,
    aes(x = C, y = DB),
    shape = 21, fill = "grey85", color = "grey55",
    size = 3, stroke = 0.4
  ) +
  # pies only for (hits) ∧ (regulated if requested)
  geom_scatterpie(
    data = pie_data,
    aes(x = C, y = DB, r = pie_r),
    cols = cond_levels,
    color = "grey15", linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "ova_central"    = "firebrick",
      "Ovaries_portal" = "darkorange",
      "testes_central" = "steelblue",
      "testes_portal"  = "forestgreen"
    ),
    breaks = cond_levels,
    name = "Condition"
  ) +
  labs(
    title = if (only_regulated)
      paste0("TG — pies only for |LogFC| > ", min_abs_logfc, " (testes biais) // (grey = no zonation/testes biais)")
    else
      "TG — pies for any condition (grey = none)",
    x = "C", y = "DB"
  ) +
  ylim(c(0, 15)) +
  xlim(c(30, 65)) +
  coord_fixed() +
  theme_minimal() +
  theme(legend.position = "right")

p


positive_or_negative <- -1      # <- e.g. 1 or -1

df <- read.csv("data/metabolomics/SM.csv")
df_with_names <- read.csv("data/metabolomics/SM_with_names.csv")

df_with_names <- merge(df_with_names, correlation_table_merged, by.x = "lipid", by.y = "node2")
df_with_names$CDB_ID <- paste(df_with_names$C, df_with_names$DB, sep = "_")
df$CDB_ID <- paste(df$C, df$DB, sep = "_")

# --- build presence (0/1) per condition for each CDB_ID ---
cond_levels <- c("ova_central", "Ovaries_portal", "testes_central", "testes_portal")

presence_wide <- df_with_names %>%
  distinct(CDB_ID, condition) %>%
  filter(condition %in% cond_levels) %>%
  mutate(value = 1L) %>%
  tidyr::pivot_wider(
    names_from = condition, values_from = value,
    values_fill = 0L
  )

for (cc in cond_levels) {
  if (!cc %in% names(presence_wide)) presence_wide[[cc]] <- 0L
}

# --- join onto df; compute hits and regulation status ---
df_anno <- df %>%
  left_join(presence_wide, by = "CDB_ID") %>%
  mutate(across(all_of(cond_levels), ~replace_na(., 0L))) %>%
  mutate(
    any_hit   = rowSums(across(all_of(cond_levels))) > 0,
    regulated = !is.na(LogFC) & LogFC * positive_or_negative >= 0 & LogFC != 0
  )

# choose which rows get pies
pie_data <- df_anno %>%
  filter(any_hit, if (only_regulated) regulated else TRUE)

# pie radius in data units
pie_r <- 0.6

# --- plot: background grey for all, pies for selected rows ---
p <- ggplot() +
  # baseline: light grey (covers both "no condition" and "not regulated")
  geom_point(
    data = df_anno,
    aes(x = C, y = DB),
    shape = 21, fill = "grey85", color = "grey55",
    size = 3, stroke = 0.4
  ) +
  # pies only for (hits) ∧ (regulated if requested)
  geom_scatterpie(
    data = pie_data,
    aes(x = C, y = DB, r = pie_r),
    cols = cond_levels,
    color = "grey15", linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "ova_central"    = "firebrick",
      "Ovaries_portal" = "darkorange",
      "testes_central" = "steelblue",
      "testes_portal"  = "forestgreen"
    ),
    breaks = cond_levels,
    name = "Condition"
  ) +
  labs(
    title = if (only_regulated)
      paste0("SM — pies only for |LogFC| < ", min_abs_logfc, " (ovary biais) // (grey = no zonation/ovary biais)")
    else
      "SM — pies for any condition (grey = none)",
    x = "C", y = "DB"
  ) +
  ylim(c(0, 15)) +
  xlim(c(30, 65)) +
  coord_fixed() +
  theme_minimal() +
  theme(legend.position = "right")

p
