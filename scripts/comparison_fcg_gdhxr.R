library(readr)
library(ggplot2)
library(ggrepel)

gdhxr_groups_tval_df <- as.data.frame(read_csv("results/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv"))

tval_full_df <- as.data.frame(read_csv("results/pseudobulk/tval_full_df.csv"))

t_val_compararison <- merge(gdhxr_groups_tval_df[,c(1,2)],tval_full_df[,c(1,2)])

plot(t_val_compararison$testes_vs_ova_tval_Hepatocytes, t_val_compararison$testes_gonadosectomy)

t_val_compararison$labels <- ifelse(
  abs(t_val_compararison$testes_gonadosectomy) > 3.5 & 
    abs(t_val_compararison$testes_vs_ova_tval_Hepatocytes) > 15,
  t_val_compararison$ID,
  NA)

names(t_val_compararison)[c(2,3)] <- c("gdhxr_tval","testes_vs_ova_tval")

ggplot(t_val_compararison, aes(x = testes_vs_ova_tval, y = gdhxr_tval, label = labels)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = T, color = "royalblue") +
  geom_label_repel() +
  theme_minimal(base_family = "Arial")
  
