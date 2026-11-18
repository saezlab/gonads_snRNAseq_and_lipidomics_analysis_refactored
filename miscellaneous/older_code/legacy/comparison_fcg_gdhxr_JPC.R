library(readr)
library(ggplot2)
library(ggrepel)

gdhxr_groups_tval_df <- as.data.frame(read_csv("results/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv"))

tval_full_df <- as.data.frame(read_csv("results/pseudobulk/tval_full_df.csv"))

t_val_compararison <- merge(gdhxr_groups_tval_df[,c(1,2)],tval_full_df[,c(1,2)])

plot(t_val_compararison$testes_vs_ova_tval_Hepatocytes, t_val_compararison$testes_gonadosectomy)

t_val_compararison$labels <- ifelse(
  abs(t_val_compararison$testes_gonadosectomy) > 5 & 
    abs(t_val_compararison$testes_vs_ova_tval_Hepatocytes) > 17,
  t_val_compararison$ID,
  NA)

names(t_val_compararison)[c(2,3)] <- c("gdhxr_tval","testes_vs_ova_tval")

# Compute Pearson correlation
cor_value <- cor(t_val_compararison$testes_vs_ova_tval, t_val_compararison$gdhxr_tval, method = "pearson")
cor_text <- paste0("Pearson r = ", round(cor_value, 2))

t_val_compararison$color <- ifelse(t_val_compararison$testes_vs_ova_tval > 0, "Higher in testes group","Higer in ovary group")

fontsize <- 12
# Plot with annotation in top right
ggplot(t_val_compararison, aes(x = testes_vs_ova_tval, y = gdhxr_tval, label = labels, color = color)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  xlab("FCG gonad effect size (t-value)") +
  ylab("GDX gonad effect size (t-value)") +
  theme_minimal() +
  theme(axis.title = element_text(size = fontsize),
        axis.text = element_text(size = fontsize),
        plot.title = element_text(size = fontsize, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  geom_label_repel(size = 5, max.overlaps = Inf) +
  scale_color_manual(values = c("red", "royalblue")) +
  annotate("text", x = Inf, y = Inf, label = cor_text, hjust = 1.1, vjust = 1.5, size = 6)

# ggsave("Fig_panels/Fig2_fcg_vs_gdx.pdf", width = 10, height = 9, device = cairo_pdf)
