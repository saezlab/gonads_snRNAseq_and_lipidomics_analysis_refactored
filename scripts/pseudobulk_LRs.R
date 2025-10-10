library(decoupleR)
library(readr)
library(GSEABase)
library(reshape2)
library(pheatmap)
library(liana)
library(dplyr)

source("scripts/support_functions.R")

load(file = "results/pseudobulk/tval_full_df.RData")

# # Here, we will convert LIANA's Consensus resource to murine symbols
# op_resource <- select_resource("Consensus")[[1]]
# 
# # Generate orthologous resource
# ortholog_resource <- generate_homologs(op_resource = op_resource,
#                                        target_organism = 10090) # mouse
# 
# save(ortholog_resource, file = "support/ortholog_resource.RData")
load("support/ortholog_resource.RData")

ligrec_regulons <- liana_to_decoupleR(ortholog_resource)


tval_full_matrix <- tval_full_df
row.names(tval_full_matrix) <- tval_full_matrix$ID
tval_full_matrix <- tval_full_matrix[,-1]

#"Ova_vs_testes_tval_","XX_vs_XY_tval_","XXO_vs_XXT_tval_"
comparison_group <- "testes_vs_ova_tval_"

tval_matrix <- tval_full_matrix[,which(grepl(comparison_group,names(tval_full_matrix)))]

cell_types <- gsub(".*_","",names(tval_matrix))

#let's make an inde carousel so that each sender cell can match all the possible receiver cells
cell_type_carousel <- list()
cell_type_carousel[[1]] <- c(1:length(cell_types))
for(i in 1:(length(cell_types)-1))
{
  cell_type_carousel[[i+1]] <- c(cell_type_carousel[[i]][-1],cell_type_carousel[[i]][1])
}

#this will match each sender matrix with each receiver matrix, in every possible pairwise combinations (including self-self)
LR_decoupler_mat <- list()
for(i in 1:length(cell_types))
{
  sender_mat <- tval_matrix
  receiver_mat <- tval_matrix[,cell_type_carousel[[i]]]
  names(sender_mat) <- paste(names(tval_matrix),names(receiver_mat),sep = "_toward_")
  names(receiver_mat) <- paste(names(tval_matrix),names(receiver_mat),sep = "_toward_")
  
  row.names(sender_mat) <- paste("ligand_",row.names(sender_mat), sep = "")
  row.names(receiver_mat) <- paste("receptor_",row.names(receiver_mat), sep ="")
  
  decoupler_mat <- as.data.frame(rbind(sender_mat,receiver_mat))
  
  names(decoupler_mat) <- gsub(comparison_group,"",names(decoupler_mat))
  
  LR_decoupler_mat[[i]] <- decoupler_mat
}

LR_decoupler_mat <- as.data.frame(do.call(cbind,LR_decoupler_mat))

LR_decoupler_mat <- LR_decoupler_mat[which(row.names(LR_decoupler_mat) %in% ligrec_regulons$target),]

LR_scores <- decoupleR::run_ulm(mat = LR_decoupler_mat, network = ligrec_regulons, minsize = 2)

LR_scores_df <- dcast(LR_scores,formula = source~condition,value.var = "score")
row.names(LR_scores_df) <- LR_scores_df$source

ligrec_scores_df_top <- LR_scores_df[apply(LR_scores_df[,-1],1,function(x){max(abs(x)) > 5}),]

#
t <- as.vector(t(ligrec_scores_df_top[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
palette <- c(palette1, palette2)
pheatmap::pheatmap(ligrec_scores_df_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, gaps_col = c(5,10,15,20))
name_of_file <- paste("results/pseudobulk/ligrecs/ligrec_scores_df_top",paste(comparison_group,".pdf",sep = ""),sep = "_")
pheatmap::pheatmap(ligrec_scores_df_top[,-1], show_rownames = T, cluster_cols = F, cluster_rows = T,color = palette, angle_col = 315, display_numbers = T, gaps_col = c(5,10,15,20), filename = name_of_file, width = 17, height = 8)
dev.off()
##
