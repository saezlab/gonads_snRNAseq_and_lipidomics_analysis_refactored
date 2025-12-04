# ── libraries ────────────────────────────────────────────────────────────────
library(tidyverse)    # readr, dplyr, ggplot2, etc.

# ── import ───────────────────────────────────────────────────────────────────
df1 <- read_csv("results/in_vitro_validation/validation_scores_df_1.csv")
df2 <- read_csv("results/in_vitro_validation/validation_scores_df_1_sszb.csv")

# ── harmonise & merge ────────────────────────────────────────────────────────
df_comb <- df1 %>%
  mutate(gonade_zonation = "global") %>%            # add the missing column
  bind_rows(df2) %>%                                # stack the two dataframes
  mutate(
    significant = p_value < 0.05,                   # flag p‑value < 0.05
    condition   = fct_reorder(condition, score)     # order by score for plotting
  )

# ── plot ─────────────────────────────────────────────────────────────────────
plot_scores <- ggplot(df_comb,
                      aes(x = condition, y = score, fill = significant)) +
  geom_col(width = 0.75, show.legend = FALSE) +      # bar‑style
  geom_hline(yintercept = 0, linetype = "dashed") +  # 0 reference
  geom_text(aes(label = ifelse(significant, "*", "")),
            vjust = ifelse(df_comb$score > 0, -0.6, 1.4), size = 6) +
  scale_fill_manual(values = c(`TRUE`  = "steelblue",
                               `FALSE` = "grey80")) +
  facet_wrap(~ gonade_zonation, scales = "free_x") + # separate panels if desired
  labs(title = "Validation scores ( *p* < 0.05 flagged )",
       x = NULL, y = "Score") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ── display or save ─────────────────────────────────────────────────────────
print(plot_scores)
# ggsave("validation_scores_plot.png", plot_scores, width = 8, height = 5, dpi = 300)

# ── plot with harmonised bar widths ──────────────────────────────────────────
plot_scores <- ggplot(df_comb,
                      aes(x = condition, y = score, fill = significant)) +
  geom_col(width = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = ifelse(significant, "*", "")),
            vjust = ifelse(df_comb$score > 0, -0.6, 1.4), size = 6) +
  scale_fill_manual(values = c(`TRUE`  = "steelblue",
                               `FALSE` = "grey80")) +
  facet_grid(~ gonade_zonation,
             scales = "free_x",
             space  = "free_x") +
  labs(title = "Validation scores ( *p* < 0.05 flagged )",
       x = NULL, y = "Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1),  # rotate condition names a bit more
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)  # rotate group labels
  )

print(plot_scores)

first_res <- df_comb[c(1:4),]
plot_scores <- ggplot(first_res,
                      aes(x = condition, y = score, fill = significant)) +
  geom_col(width = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = ifelse(significant, "*", "")),
            vjust = ifelse(first_res$score > 0, -0.6, 1.4), size = 6) +
  scale_fill_manual(values = c(`TRUE`  = "steelblue",
                               `FALSE` = "grey80")) +
  facet_grid(~ gonade_zonation,
             scales = "free_x",
             space  = "free_x") +
  labs(title = "Validation scores ( *p* < 0.05 flagged )",
       x = NULL, y = "Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1),  # rotate condition names a bit more
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)  # rotate group labels
  )
print(plot_scores)


second_res <- df_comb[c(5:11),]
plot_scores <- ggplot(second_res,
                      aes(x = condition, y = score, fill = significant)) +
  geom_col(width = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = ifelse(significant, "*", "")),
            vjust = ifelse(second_res$score > 0, -0.6, 1.4), size = 6) +
  scale_fill_manual(values = c(`TRUE`  = "steelblue",
                               `FALSE` = "grey80")) +
  facet_grid(~ gonade_zonation,
             scales = "free_x",
             space  = "free_x") +
  labs(title = "Validation scores ( *p* < 0.05 flagged )",
       x = NULL, y = "Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1),  # rotate condition names a bit more
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)  # rotate group labels
  )
print(plot_scores)
