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

lipidomic <- as.data.frame(
  read_csv("data/metabolomics/lipidomics_20250319transfections.csv"))

lipidomic <- as.data.frame(clean_lipid_data(lipidomic))

targets <- as.data.frame(cbind(names(lipidomic),names(lipidomic)))[-1,]
targets[,2] <- gsub("_.*","",targets[,2])
names(targets) <- c("sample","condition")

row.names(lipidomic) <- lipidomic$Lipid
lipidomic <- lipidomic[,-1]

grid.newpage()  # Start a new page
grid.rect(gp = gpar(fill = "white", col = NA))  # Fill the entire background with white
grid.draw(nicePCA(lipidomic, targets, no_label = T))

unique(targets$condition)

comparisons <- list("cyp2d9" = c(5,-1),
                    "cyp2b9" = c(2,-1),
                    "cyp2b13" = c(3,-1),
                    "cyp2c38" = c(4,-1))

limmaRes <- runLimma(lipidomic, targets, comparisons = comparisons)

ttop_list <- limma_res_to_ttop_list(limma_res = limmaRes,
                                    comp_names = names(comparisons),
                                    number = length(lipidomic[,1]),
                                    adjust.method = "fdr")

t_table <- ttop_list_to_t_table(ttop_list)
row.names(t_table) <- t_table$ID

t <- as.vector(t(t_table[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(t_table[,-1], show_rownames = F, cluster_cols = T, cluster_rows = T,color = palette, angle_col = 315, display_numbers = F)
dev.off()

lipid_network <- read_csv(file = "results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full_largenet.csv")
names(lipid_network) <- c("target","source","mor")
lipid_network <- lipid_network[which(grepl("[(]",lipid_network$target)),]
# lipid_network <- lipid_network[which(abs(lipid_network$mor) > 0.7),]
# View(lipid_network[lipid_network$source == "Cyp2b13" & lipid_network$target %in% t_table$ID,])

validation_scores <- decoupleR::run_ulm(t_table[,-1], network = lipid_network, minsize = 2)
validation_set_only <- validation_scores[tolower(validation_scores$source) == validation_scores$condition,]
validation_scores_df <- dcast(validation_scores, source~condition, value.var = "score")

write_csv(validation_set_only, file = "results/in_vitro_validation/validation_scores_df_1.csv")

KO_index <- 3
lipid_scoring_Scatter <- merge(t_table[,c(1,KO_index)], lipid_network, by.x = "ID", by.y = "target")
lipid_scoring_Scatter$mor <- ifelse(tolower(lipid_scoring_Scatter$source) == names(lipid_scoring_Scatter)[2],lipid_scoring_Scatter$mor,0)

names(lipid_scoring_Scatter) <- c("ID","KO_tval","KO_target","predicted_association_weight")
lipid_scoring_Scatter$ID <- ifelse(lipid_scoring_Scatter$predicted_association_weight != 0, lipid_scoring_Scatter$ID, NA)

ggplot(lipid_scoring_Scatter, aes(x = predicted_association_weight, y = KO_tval, label = ID)) + 
  geom_point() + 
  geom_label_repel() +
  ggtitle(names(t_table)[KO_index]) +
  theme_minimal() + 
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x)

