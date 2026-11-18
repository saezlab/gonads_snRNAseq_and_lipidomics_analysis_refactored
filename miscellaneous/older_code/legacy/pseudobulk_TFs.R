library(decoupleR)
library(readr)
library(GSEABase)
library(reshape2)
library(pheatmap)

source("scripts/support_functions.R")

load(file = "results/pseudobulk/tval_full_df.RData")

# TF_regulons <- decoupleR::get_collectri(organism = "mouse")
# save(TF_regulons, file = "support/TF_regulons.RData")
load("support/TF_regulons.RData")
TF_regulons <- TF_regulons[complete.cases(TF_regulons),]


tval_full_matrix <- tval_full_df
row.names(tval_full_matrix) <- tval_full_matrix$ID
tval_full_matrix <- tval_full_matrix[,-1]

TF_scores <- decoupleR::run_ulm(mat = tval_full_matrix, network = TF_regulons[,], minsize = 20)
TF_scores_df <- dcast(TF_scores, source~condition, value.var = "score")

TF_scores_df_gonade <- TF_scores_df[,c(1,which(grepl("testes",names(TF_scores_df))))]
TF_scores_df_allosome <- TF_scores_df[,c(1,which(grepl("XX_vs_XY",names(TF_scores_df))))]
TF_scores_df_group <- TF_scores_df[,c(1,which(grepl("XXO",names(TF_scores_df))))]


#Just to make some plresenation
# library(ggplot2)
# 
# TF_int <- "Dbp"
# 
# reg_sub <- TF_regulons[TF_regulons$source == TF_int,]
# 
# reg_sub <- merge(reg_sub,tval_full_df, by.x = "target", by.y = "ID", all.y = T)
# 
# reg_sub[is.na(reg_sub$mor),"mor"] <- 0
# 
# ggplot(reg_sub, aes(x = mor, y = testes_vs_ova_tval_Hepatocytes)) + geom_point() +
#   stat_summary(fun.data= mean_cl_normal) + 
#   geom_smooth(method='lm') + theme_minimal()

##

TF_scores_df_gonade_top <- TF_scores_df_gonade[apply(TF_scores_df_gonade[,-1],1,function(x){max(abs(x)) > 4}),]
row.names(TF_scores_df_gonade_top) <- TF_scores_df_gonade_top$source

t <- as.vector(t(TF_scores_df_gonade_top[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(TF_scores_df_gonade_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T)
pheatmap::pheatmap(TF_scores_df_gonade_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, filename = "results/pseudobulk/TFs/TF_scores_df_gonade_top.pdf", width = 11, height = 7.5)
dev.off()
##

TF_scores_df_allosome_top <- TF_scores_df_allosome[apply(TF_scores_df_allosome[,-1],1,function(x){max(abs(x)) > 2}),]
row.names(TF_scores_df_allosome_top) <- TF_scores_df_allosome_top$source

t <- as.vector(t(TF_scores_df_allosome_top[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(TF_scores_df_allosome_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T)
pheatmap::pheatmap(TF_scores_df_allosome_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, filename = "results/pseudobulk/TFs/TF_scores_df_allosome_top.pdf", width = 10, height = 5.5)
dev.off()
##

TF_scores_df_group_top <- TF_scores_df_group[apply(TF_scores_df_group[,-1],1,function(x){max(abs(x)) > 4}),]
row.names(TF_scores_df_group_top) <- TF_scores_df_group_top$source

t <- as.vector(t(TF_scores_df_group_top[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(TF_scores_df_group_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T)
pheatmap::pheatmap(TF_scores_df_group_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, filename = "results/pseudobulk/TFs/TF_scores_df_group_top.pdf", width = 8, height = 7.5)
dev.off()

TF_scores_df_all <- merge(TF_scores_df_allosome,TF_scores_df_gonade, by = "source")
TF_scores_df_all <- merge(TF_scores_df_all,TF_scores_df_group, by = "source")

save(TF_scores_df_all, file = "results/pseudobulk/TFs/TF_scores_df_all.RData")
