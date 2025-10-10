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

lipid_scoring_Scatter$label <- c(lipid_scoring_Scatter$ID[1:20],rep(NA,length(lipid_scoring_Scatter$ID)-20))

ggplot(lipid_scoring_Scatter,
       aes(predicted_association_weight, KO_tval, label = label, colour = consistency)
) +
  geom_point() +
  geom_smooth(
    inherit.aes = FALSE,
    mapping     = aes(x = predicted_association_weight, y = KO_tval),
    method      = "lm",
    formula     = y ~ x
  ) +
  stat_summary(fun.data = mean_cl_normal) +
  scale_colour_gradient(
    name   = "Consistency",
    low    = "black",    # darkest end
    high   = "chartreuse3",# brightest end
    limits = c(-5, 5),
    oob = squish
  ) +
  ggtitle(paste(names(t_table)[KO_index], zzsb, sep = " ")) +
  xlab("Correlation between hepatocyte gene expression\nand lipid abundance in FCGL mice liver") +
  ylab("Deregulation of lipid species (t-value)\nin knock-down experiment") +
  geom_label_repel() +
  theme_minimal()

#for slideck
ggplot(lipid_scoring_Scatter,
       aes(predicted_association_weight, KO_tval, label = label, colour = consistency)
) +
  
  ggtitle(paste(names(t_table)[KO_index], zzsb, sep = " ")) +
  xlab("Correlation between hepatocyte gene expression\nand lipid abundance in FCGL mice liver") +
  ylab("Deregulation of lipid species (t-value)\nin knock-down experiment") +
  theme_minimal()

#for slideck
ggplot(lipid_scoring_Scatter,
       aes(predicted_association_weight, KO_tval, label = label, colour = consistency)
) +
  geom_point() +
  geom_smooth(
    inherit.aes = FALSE,
    mapping     = aes(x = predicted_association_weight, y = KO_tval),
    method      = "lm",
    formula     = y ~ x
  ) +
  stat_summary(fun.data = mean_cl_normal) +
  scale_colour_gradient(
    name   = "Consistency",
    low    = "black",    # darkest end
    high   = "chartreuse3"      # brightest end
  ) +
  ggtitle(paste(names(t_table)[KO_index], zzsb, sep = " ")) +
  xlab("Correlation between hepatocyte gene expression\nand lipid abundance in FCGL mice liver") +
  ylab("Deregulation of lipid species (t-value)\nin knock-down experiment") +
  theme_minimal()


# # install.packages("plotly")
# library(plotly)
# 
# # 1) Basic 3D scatter, coloured by consistency
# p <- plot_ly(
#   data = lipid_scoring_Scatter,
#   x    = ~predicted_association_weight,
#   y    = ~KO_tval,
#   z    = ~consistency,
#   text = lipid_scoring_Scatter$label,  # pulls it directly by vector
#   type = "scatter3d", mode = "markers",
#   marker = list(
#     size      = 4,
#     color     = ~consistency,
#     colorscale= list(c(0,1), c("black","chartreuse3")),
#     showscale = TRUE,
#     colorbar  = list(title="Consistency")
#   )
# )
# 
# # 2) (Optional) Add a regression plane of consistency ~ x + y
# mod <- lm(consistency ~ predicted_association_weight + KO_tval,
#           data = lipid_scoring_Scatter)
# 
# # make a grid over x,y
# grid <- expand.grid(
#   predicted_association_weight = seq(min(lipid_scoring_Scatter$predicted_association_weight),
#                                      max(lipid_scoring_Scatter$predicted_association_weight),
#                                      length = 30),
#   KO_tval = seq(min(lipid_scoring_Scatter$KO_tval),
#                 max(lipid_scoring_Scatter$KO_tval),
#                 length = 30)
# )
# grid$consistency <- predict(mod, newdata = grid)
# 
# # 3) Tidy up axis labels and title
# p <- p %>% layout(
#   title = paste(names(t_table)[KO_index], zzsb),
#   scene = list(
#     xaxis = list(title="Correlation between hepatocyte gene expression\nand lipid abundance in FCGL mice liver"),
#     yaxis = list(title="Deregulation of lipid species (t-value)\nin knock-down experiment"),
#     zaxis = list(title="Consistency")
#   )
# )
# 
# p  # print the widget
# 



