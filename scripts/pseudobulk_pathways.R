library(decoupleR)
library(readr)
library(GSEABase)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)            # colourRamp2 lives here

source("scripts/support_functions.R")

load("results/pseudobulk/ttop_full_list.RData")

tval_df_list <- list()
for(grouping_type in names(ttop_full_list))
{
  ttop_list <- lapply(ttop_full_list[[grouping_type]], function(x){return(x[,c(1,4)])})
  tval_df <- merge(ttop_list[[1]],ttop_list[[2]])
  for(i in 3:length(ttop_list))
  {
    tval_df <- merge(tval_df,ttop_list[[i]])
  }
  names(tval_df)[-1] <- paste(paste(grouping_type,"tval",sep = "_"),c(names(ttop_full_list[[grouping_type]])),sep = "_")
  tval_df_list[[grouping_type]] <- tval_df
}

tval_full_df <- merge(tval_df_list[[1]], tval_df_list[[2]])
tval_full_df <- merge(tval_full_df, tval_df_list[[3]])

save(tval_full_df, file = "results/pseudobulk/tval_full_df.RData")
write_csv(tval_full_df, file = "results/pseudobulk/tval_full_df.csv")

pathways <- import_gmt("support/m5.all.v2023.2.Mm.symbols.gmt")
names(pathways) <- c("target","source")
pathways$mor <- 1

tval_full_matrix <- tval_full_df
row.names(tval_full_matrix) <- tval_full_matrix$ID
tval_full_matrix <- tval_full_matrix[,-1]

pathway_scores <- decoupleR::run_ulm(tval_full_matrix, network = pathways)
pathway_scores_df <- dcast(pathway_scores, source~condition, value.var = "score")

pathway_scores_df_gonade <- pathway_scores_df[,c(1,which(grepl("testes",names(pathway_scores_df))))]
pathway_scores_df_allosome <- pathway_scores_df[,c(1,which(grepl("XY_vs",names(pathway_scores_df))))]
pathway_scores_df_group <- pathway_scores_df[,c(1,which(grepl("XXT_vs",names(pathway_scores_df))))]

##

pathway_scores_df_gonade_top <- pathway_scores_df_gonade[apply(pathway_scores_df_gonade[,-1],1,function(x){max(abs(x)) > 12}),]
row.names(pathway_scores_df_gonade_top) <- pathway_scores_df_gonade_top$source
row.names(pathway_scores_df_gonade_top) <- gsub("_"," ",gsub("GO[CMB][CFP]","",row.names(pathway_scores_df_gonade_top)))

names(pathway_scores_df_gonade_top) <- gsub(".*_","",names(pathway_scores_df_gonade_top))
row.names(pathway_scores_df_gonade_top) <- tolower(row.names(pathway_scores_df_gonade_top))
pathway_scores_df_gonade_top <- pathway_scores_df_gonade_top[,c(1,4,3,5,6,2)]
t <- as.vector(t(pathway_scores_df_gonade_top[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(pathway_scores_df_gonade_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T)
pheatmap::pheatmap(pathway_scores_df_gonade_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, filename = "results/pseudobulk/pathways/pathway_scores_df_gonade_top.pdf", width = 11, height = 7.5)
dev.off()
##

mat <- pathway_scores_df_gonade_top[c(2,8,9,11,12),-1]
pheatmap::pheatmap(df, show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T)
dev.off()


ht <- Heatmap(
  mat,
  name               = "logFC",
  col                = col_fun,          # symmetric colour map you already built
  cluster_rows       = TRUE,
  cluster_columns    = FALSE,
  column_names_rot   = 45,
  heatmap_legend_param = list(
    at     = ticks,
    labels = ticks
  )
)

## ---- DRAW, with extra padding ---------------------------------------------
## padding = unit(c(bottom, left, top, right), "mm")
draw(ht,
     heatmap_legend_side = "left",
     padding = unit(c(4, 4, 4, 25), "mm"))  # ← tweak the 14 mm until nothing overlaps

pathway_scores_df_allosome_top <- pathway_scores_df_allosome[apply(pathway_scores_df_allosome[,-1],1,function(x){max(abs(x)) > 11}),]
row.names(pathway_scores_df_allosome_top) <- pathway_scores_df_allosome_top$source

t <- as.vector(t(pathway_scores_df_allosome_top[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(pathway_scores_df_allosome_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T)
pheatmap::pheatmap(pathway_scores_df_allosome_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, filename = "results/pseudobulk/pathways/pathway_scores_df_allosome_top.pdf", width = 10, height = 5.5)
dev.off()
##

pathway_scores_df_group_top <- pathway_scores_df_group[apply(pathway_scores_df_group[,-1],1,function(x){max(abs(x)) > 11}),]
row.names(pathway_scores_df_group_top) <- pathway_scores_df_group_top$source

t <- as.vector(t(pathway_scores_df_group_top[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(pathway_scores_df_group_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T)
pheatmap::pheatmap(pathway_scores_df_group_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, filename = "results/pseudobulk/pathways/pathway_scores_df_group_top.pdf", width = 8, height = 7.5)
dev.off()

pathway_scores_df_all <- merge(pathway_scores_df_allosome,pathway_scores_df_gonade, by = "source")
pathway_scores_df_all <- merge(pathway_scores_df_all,pathway_scores_df_group, by = "source")

save(pathway_scores_df_all, file = "results/pseudobulk/pathways/pathway_scores_df_all.RData")
