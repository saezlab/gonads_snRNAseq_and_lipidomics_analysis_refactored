library(readr)
library(readxl)
library(pheatmap)
library(limma)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(cowplot)
library(gridExtra)
library(decoupleR)
library(grid)
source("scripts/support_functions.R")


#from https://github.com/saezlab/liver_sex_differences
lipidomic <- read_csv("data/metabolomics/lipidomics_gdxhr.csv")

lipidomic <- lipidomic[,c(4,35:72)]

lipidomic <- as.data.frame(clean_lipid_data(lipidomic))

row.names(lipidomic) <- lipidomic$Lipid
lipidomic <- as.data.frame(lipidomic[,-1])

targets <- as.data.frame(cbind(names(lipidomic),names(lipidomic)))
targets$V2 <- sub("^[^_]*_", "", targets$V1)
names(targets) <- c("sample","condition")

targets_gdx <- targets[grepl("_v$",targets$condition),]


grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(lipidomic,targets_gdx, no_label = T))


unique(targets$condition)

comparisons <- list("testes_gonadosectomy" = c(7,-8),
                    "testes_oestro" = c(5,-7),
                    "testes_testo" = c(6,-7),
                    "ova_gonadosectomy" = c(3,-4),
                    "ova_oestro" = c(1,-3),
                    "ova_testo" = c(2,-3))

limmaRes <- runLimma(lipidomic, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(lipidomic[,1]),
                                    adjust.method = "fdr")

t_table <- merge(ttop_list[[1]][,c(1,4)],ttop_list[[2]][,c(1,4)], by = "ID")
for(i in 3:length(ttop_list))
{
  t_table <- merge(t_table,ttop_list[[i]][,c(1,4)], by = "ID")
}
names(t_table)[-1] <- names(ttop_list)

write_csv(t_table, file = "results/lipidomic/gdxhr_t_table.csv")

#CONNECT LIPID CLASSE WITH ENZYMES

lipid_classes <- as.data.frame(cbind(t_table$ID,gsub("[(].*","",t_table$ID)))
names(lipid_classes) <- c("target","source")
barplot(table(lipid_classes))
lipid_classes$mor <- 1

row.names(t_table) <- t_table$ID

lipid_classe_enrichment <- decoupleR::run_ulm(mat = t_table[,-1], network = lipid_classes, minsize = 10)
lipid_classe_enrichment <- dcast(lipid_classe_enrichment, source~condition, value.var = "score")
row.names(lipid_classe_enrichment) <- lipid_classe_enrichment$source

lipid_counts_in_classes <- table(lipid_classes$source)
lipid_counts_in_classes <- lipid_counts_in_classes[lipid_counts_in_classes >= 10]

row.names(lipid_classe_enrichment) <- paste(row.names(lipid_classe_enrichment), as.numeric(lipid_counts_in_classes), sep = " n=")

t <- as.vector(t(lipid_classe_enrichment[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10, my_colors = c("darkgrey","white","darkgrey"))
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10, my_colors = c("darkgrey","white","darkgrey"))
palette <- c(palette1, palette2)
pheatmap::pheatmap(lipid_classe_enrichment[,-1], display_numbers = T,
                   show_rownames = T, cluster_cols = T, cluster_rows = T,color = palette, angle_col = 315)

lipid_classe_enrichment$source

lipid_classification <- as.data.frame(read_csv("support/lipid_classification.csv"))

lipid_enrichment_gonade <- lipid_classe_enrichment[,c(1,5)]
lipid_enrichment_gonade$n <- row.names(lipid_enrichment_gonade)

lipid_enrichment_gonade <- merge(lipid_enrichment_gonade, lipid_classification, by.x = "source", by.y = "Lipid")
lipid_enrichment_gonade <- lipid_enrichment_gonade[order(abs(lipid_enrichment_gonade$testes_gonadosectomy), decreasing = F),]
lipid_enrichment_gonade$n <- factor(lipid_enrichment_gonade$n, levels = lipid_enrichment_gonade$n)

ggplot(lipid_enrichment_gonade, aes(x = n, y = testes_gonadosectomy, color = testes_gonadosectomy)) + 
  geom_point(size = 4) + 
  geom_segment(aes(x=n, xend=n, y=0, yend=testes_gonadosectomy), size = 1) + 
  coord_flip() + 
  theme_minimal() +
  scale_color_gradient2(
    midpoint = 0,    # zero is the midpoint
    low = "red",     # negative values
    mid = "grey",   # zero
    high = "royalblue"    # positive values
  )

gonade_type_ttop <- as.data.frame(
  read_csv("results/lipidomic/gonade_type_ttop.csv"))

t_val_compararison <- merge(gonade_type_ttop[,c(1,4)],t_table[,c(1,2)])
names(t_val_compararison) <- c("ID","fcgl","testes_gonadosectomy")

plot(t_val_compararison$fcgl, t_val_compararison$testes_gonadosectomy)

t_val_compararison$labels <- ifelse(
  abs(t_val_compararison$fcgl) > 13 & 
    abs(t_val_compararison$testes_gonadosectomy) > 3,
  t_val_compararison$ID,
  NA)

t_val_compararison$color <- ifelse(t_val_compararison$fcgl > 0, "testes","ovary")
ggplot(t_val_compararison, aes(x = fcgl, y = testes_gonadosectomy, label = labels, color = color)) + 
  geom_point(alpha = 0.5) + 
  scale_color_manual(values = c("red", "royalblue")) +
  geom_smooth(method = "lm", se = T, color = "black") +
  geom_label_repel() +
  theme_minimal()



