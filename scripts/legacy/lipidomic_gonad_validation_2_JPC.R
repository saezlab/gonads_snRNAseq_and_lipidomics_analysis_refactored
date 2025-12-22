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
library(scales)
library(Cairo)

source("scripts/support_functions.R")

lipidomic <- as.data.frame(
  read_csv("data/metabolomics/lipidomics_20250423transfections.csv"))

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

comparisons <- list("cyp2d9" = c(4,-1),
                    "cyp2b9" = c(2,-1),
                    "cyp2c38" = c(3,-1),
                    "mtmr11" = c(6,-1),
                    "sult3a1" = c(5,-1))

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

KO_index <- 4
lipid_scoring_Scatter <- merge(t_table[,c(1,KO_index)], lipid_network, by.x = "ID", by.y = "target")
lipid_scoring_Scatter$mor <- ifelse(tolower(lipid_scoring_Scatter$source) == names(lipid_scoring_Scatter)[2],lipid_scoring_Scatter$mor,0)

names(lipid_scoring_Scatter) <- c("ID","KO_tval","KO_target","predicted_association_weight")
lipid_scoring_Scatter$ID <- ifelse(lipid_scoring_Scatter$predicted_association_weight != 0, lipid_scoring_Scatter$ID, NA)

lipid_scoring_Scatter$consistency <- lipid_scoring_Scatter$KO_tval * lipid_scoring_Scatter$predicted_association_weight * -1
# lipid_scoring_Scatter$consistency <- (lipid_scoring_Scatter$consistency + abs(min(lipid_scoring_Scatter$consistency)))
lipid_scoring_Scatter <- lipid_scoring_Scatter[order(lipid_scoring_Scatter$consistency, decreasing = T),]

lipid_scoring_Scatter$label <- c(lipid_scoring_Scatter$ID[1:10],rep(NA,length(lipid_scoring_Scatter$ID)-10))

ggplot(lipid_scoring_Scatter,
       aes(predicted_association_weight, KO_tval, label = label, colour = consistency)
) +
  geom_point(size=3) +
  geom_smooth(
    inherit.aes = FALSE,
    mapping     = aes(x = predicted_association_weight, y = KO_tval),
    method      = "lm",
    formula     = y ~ x,
    color       = "grey70",   # <- make the line grey
    se          = TRUE       # optional: hide confidence band
  ) +
  stat_summary(fun.data = mean_cl_normal) +
  scale_colour_gradient(
    name   = "Consistency",
    low    = "grey90",
    high   = "#fbb42c",
    limits = c(-5, 5),
    oob    = squish,
    breaks = c(-5, 5)) +
  xlab("In vivo correlation between Mtmr11 exp. and lipid abund. (corr. coeff.)") +
  ylab("In vitro effect of Mtmr11 KD on lipid abundance (t-value)") +
  theme_classic(base_family = "Arial", base_size = 20) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20,face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = c(1, 1),        # top right corner
        legend.justification = c(1, 1)    # anchor the legend's top right corner
  ) +
  geom_label_repel(
    size = 7,
    max.overlaps = 100,     # limit number plotted, allows ggrepel to push them apart
    box.padding = 1.0,      # add space around the text boxes
    point.padding = 0.5,    # add space between points and labels
    force = 2,              # make the repulsion stronger
    min.segment.length = 0  # always draw connecting lines
  )

ggsave("Fig_panels/Fig4_cyp2d9.pdf", width = 10, height = 9, device = cairo_pdf)

