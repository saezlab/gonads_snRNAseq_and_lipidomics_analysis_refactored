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

comparisons <- list("sult3a1" = c(5,-1),
                    "mtmr11" = c(6,-1))

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

testes_portal <- as.data.frame(read_csv("results/multi-omic/cormat_testes_portal_SIF.csv")) 
testes_portal <- testes_portal[grepl("[()]",testes_portal$node2),]
testes_portal$node1 <- paste(testes_portal$node1,"testes_portal",sep = "_")

testes_central <- as.data.frame(read_csv("results/multi-omic/cormat_testes_central_SIF.csv")) 
testes_central <- testes_central[grepl("[()]",testes_central$node2),]
testes_central$node1 <- paste(testes_central$node1,"testes_central",sep = "_")

ova_portal <- as.data.frame(read_csv("results/multi-omic/cormat_Ovaries_portal_SIF.csv")) 
ova_portal <- ova_portal[grepl("[()]",ova_portal$node2),]
ova_portal$node1 <- paste(ova_portal$node1,"ova_portal",sep = "_")

ova_central <- as.data.frame(read_csv("results/multi-omic/cormat_ova_central_SIF.csv")) 
ova_central <- ova_central[grepl("[()]",ova_central$node2),]
ova_central$node1 <- paste(ova_central$node1,"ova_central",sep = "_")

sszb_network <- as.data.frame(do.call(rbind,list(testes_portal,testes_central,ova_portal,ova_central)))
names(sszb_network) <- c("source","target","mor")

sszb_validation_scores <- decoupleR::run_ulm(t_table[,-1], network = sszb_network, minsize = 2)
sszb_validation_scores$gonade_zonation <- sub("^[^_]*_", "", sszb_validation_scores$source)
sszb_validation_scores$source <- tolower(gsub("_.*","",sszb_validation_scores$source))

validation_set_only <- sszb_validation_scores[sszb_validation_scores$source == sszb_validation_scores$condition,]

KO_index <- 3
zzsb <- as.character(validation_set_only[validation_set_only$source == names(t_table)[KO_index],"gonade_zonation"])
sszb_network_sub <- sszb_network[which(grepl(zzsb, sszb_network$source)),]
lipid_scoring_Scatter <- merge(t_table[,c(1,KO_index)], sszb_network_sub, by.x = "ID", by.y = "target")
lipid_scoring_Scatter$source <- gsub("_.*","",lipid_scoring_Scatter$source)
lipid_scoring_Scatter$mor <- ifelse(tolower(lipid_scoring_Scatter$source) == names(lipid_scoring_Scatter)[2],lipid_scoring_Scatter$mor,0)

names(lipid_scoring_Scatter) <- c("ID","KO_tval","KO_target","predicted_association_weight")
lipid_scoring_Scatter$ID <- ifelse(lipid_scoring_Scatter$predicted_association_weight != 0, lipid_scoring_Scatter$ID, NA)

lipid_scoring_Scatter$label <- ifelse(lipid_scoring_Scatter$KO_tval < -5 & lipid_scoring_Scatter$predicted_association_weight > 0.2, lipid_scoring_Scatter$ID, NA)
lipid_scoring_Scatter$color <- ifelse(!is.na(lipid_scoring_Scatter$label),"darkgreen","black")

zzsb <- paste(paste("(specific expression bias in ",gsub("_"," ",zzsb),sep = ""),")",sep = "")

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

ggsave("Fig_panels/Fig5_mtmr11.pdf", width = 10, height = 9, device = cairo_pdf)

