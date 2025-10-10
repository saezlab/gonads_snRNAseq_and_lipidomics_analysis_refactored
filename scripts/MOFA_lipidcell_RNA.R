library(reshape2)
library(decoupleR)
library(GSEABase)
library(readr)
library(MOFA2)
library(pheatmap)

library(ComplexHeatmap)
library(circlize)
library(grid)

source("scripts/support_functions.R")

#I sploit the data in three parts so that the file size doesn't exceed github limit
scaled_RNA_1 <- as.data.frame(read_csv("data/RNA/scaled_RNA_1.csv"))
scaled_RNA_2 <- as.data.frame(read_csv("data/RNA/scaled_RNA_2.csv"))
scaled_RNA_3 <- as.data.frame(read_csv("data/RNA/scaled_RNA_3.csv"))

scaled_RNA <- merge(scaled_RNA_1, scaled_RNA_2, by = "gene")
scaled_RNA <- merge(scaled_RNA, scaled_RNA_3, by = "gene")

row.names(scaled_RNA) <- scaled_RNA$gene
scaled_RNA <- scaled_RNA[,-1]

rm(scaled_RNA_1)
rm(scaled_RNA_2)
rm(scaled_RNA_3)

targets <- data.frame(sample = names(scaled_RNA), 
                      condition = names(scaled_RNA))
targets$condition <- gsub(".*_","",targets$condition)
targets$sample <- gsub("_[A-Z].*","",targets$sample)

names(scaled_RNA) <- gsub("_Testes","",names(scaled_RNA))
names(scaled_RNA) <- gsub("_Ovaries","",names(scaled_RNA))


load("data/metabolomics/lipidomic.RData")

scaled_RNA$feature <- row.names(scaled_RNA)

RNA_MOFA <- melt(scaled_RNA)
RNA_MOFA$variable <- as.character(RNA_MOFA$variable)
RNA_MOFA$view <- gsub(".*_","",RNA_MOFA$variable)
RNA_MOFA$variable <- gsub("_[A-Z].*","",RNA_MOFA$variable)
names(RNA_MOFA) <- c("feature","sample","value","view")
RNA_MOFA <- RNA_MOFA[,c(2,1,4,3)]
RNA_MOFA <- RNA_MOFA[RNA_MOFA$sample %in% names(lipidomic),]

lipidomic_MOFA <- lipidomic
lipidomic_MOFA$feature <- row.names(lipidomic_MOFA)
lipidomic_MOFA <- melt(lipidomic_MOFA)
lipidomic_MOFA$variable <- as.character(lipidomic_MOFA$variable)
lipidomic_MOFA$view <- "lipidomic"
names(lipidomic_MOFA) <- c("feature","sample","value","view")

MOFA_data <- as.data.frame(rbind(RNA_MOFA, lipidomic_MOFA))
MOFA_data$feature <- paste(MOFA_data$feature, MOFA_data$view, sep = "_")

write_csv(MOFA_data, file = "data/mofa/MOFA_RNA_data.csv")

# only run to retrain mofa model
# for(i in c(5:15)){
#   wd <- getwd()
#   write_MOFA_options_file(input_file = '/data/mofa/MOFA_RNA_data.csv', output_file = paste0('/results/mofa/RNA/mofa_res_',i,'factor.hdf5'), factors = i, convergence_mode = "slow")
#   system(paste0(paste('python3', wd, sep = " "), paste('/scripts/mofa/mofa2.py', "options_MOFA.csv")))
# }
# 
# Load MOFA output
for(i in c(5:15)){
  filepath <- paste0('results/mofa/RNA/mofa_res_',i,'factor.hdf5')
  model <- load_model(filepath)
  # meta_data <- read_csv("data/metadata/RNA_metadata_cluster.csv")[,c(1,2)]
  # colnames(meta_data) <- c("sample","cluster")
  # mofa_metadata <- merge(samples_metadata(model), meta_data, by = "sample")
  # names(mofa_metadata) <- c("sample","group","cluster")
  # samples_metadata(model) <- mofa_metadata
  assign(paste0("model",i),model)
}

compare_factors(list(model7,model8,model9,model10,model11,model12,model13,model15), cluster_rows = F, cluster_cols = F,)

model <- model12

plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
calculate_variance_explained(model)

r2_per_factor <- as.data.frame(model@cache$variance_explained$r2_per_factor[[1]])
# names(r2_per_factor)[c(1:5)] <- paste(names(r2_per_factor)[c(1:5)], "RNA")

r2_per_factor <- r2_per_factor[,c(3,2,4,5,1,6)]
names(r2_per_factor)[1:5] <- paste(names(r2_per_factor)[1:5],"RNA")
names(r2_per_factor)[6] <- "Liver lipidomic"

dim(r2_per_factor)
names(r2_per_factor)
row.names(r2_per_factor)

pheatmap(r2_per_factor, display_numbers = T, angle_col = "315", 
         legend_labels = c("0","10", "20", "30", "40", "Variance\n\n"), 
         legend = T, main = "", 
         legend_breaks = c(0,10, 20, 30, 40, max(r2_per_factor)), 
         cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white","grey45"))(100), fontsize_number = 10)




Z_matrix <- as.data.frame(model@expectations$Z$single_group)

targets$mor <- 1
names(targets) <- c("target","source","mor")
targets <- unique(targets)

gonade_enrichment <- decoupleR::run_ulm(Z_matrix, targets)
gonade_enrichment <- reshape2::dcast(gonade_enrichment, source~condition, value.var = "score")
row.names(gonade_enrichment) <- gonade_enrichment$source
gonade_enrichment <- as.data.frame(t(gonade_enrichment[,-1]))
gonade_enrichment <- gonade_enrichment[c(1,3,4,5,6,7,8,9,10,2),]

dim(gonade_enrichment)
names(gonade_enrichment)
row.names(gonade_enrichment)

limit <- max(abs(as.numeric(as.matrix(gonade_enrichment))))
# Define a color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Set breaks from -1 to 1
breaks <- seq(-limit, limit, length.out = 101)
pheatmap::pheatmap(gonade_enrichment,  angle_col = "0",color = my_palette, breaks = breaks, cluster_rows = F, cluster_cols = F, display_numbers = T)


mat_r2    <- as.matrix(r2_per_factor)
mat_gonad <- abs(as.matrix(gonade_enrichment))

# Row-annotation: horizontal, side-by-side bars with opacity
gonad_ann <- rowAnnotation(
  "Gonad\nEnrichment" = anno_barplot(
    mat_gonad,
    gp           = gpar(
      fill  = c("red", "royalblue"),  # ovaries, testes
      alpha = 0.7                      # 70% opacity
    ),
    border       = FALSE,
    axis         = TRUE,
    beside       = TRUE,
    orientation  = "horizontal",
    bar_width    = unit(0.6, "cm")
  ),
  annotation_width = unit(4, "cm")
)

# Grey ramp for heatmap
col_fun <- colorRamp2(
  c(min(mat_r2), median(mat_r2), max(mat_r2)),
  c("white", "grey70", "black")
)

# Build heatmap, printing values in each cell
ht <- Heatmap(
  mat_r2,
  name             = "R²",
  col              = col_fun,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  right_annotation = gonad_ann,
  row_names_side   = "right",
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.1f", mat_r2[i, j]), x, y,
      gp = gpar(
        fontsize = 9,
        col      = ifelse(mat_r2[i, j] > median(mat_r2), "white", "black")
      )
    )
  }
)

# Draw with padding so nothing is clipped
draw(ht, padding = unit(c(4, 4, 4, 4), "cm"))


##Weight analysis
weights <- get_weights(model, views = "all", factors = "all")

load("support/TF_regulons.RData")
TF_regulons <- TF_regulons[complete.cases(TF_regulons),]

Gonadal_Hormone_Genes <- as.data.frame(read_csv("support/Gonadal_Hormone_Genes.csv"))
Gonadal_Hormone_Genes$source <- gsub("_spec","",Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("ESR1","Esr1",Gonadal_Hormone_Genes$source)
Gonadal_Hormone_Genes$source <- gsub("AR","Ar",Gonadal_Hormone_Genes$source)

TF_regulons <- TF_regulons[!TF_regulons$source %in% c("Ar", "Esr1"),]

TF_regulons <- rbind(TF_regulons, Gonadal_Hormone_Genes)

RNA_weight_input <- weights$Hepatocytes
row.names(RNA_weight_input) <- gsub("_.*","",row.names(RNA_weight_input))

TF_scores <- decoupleR::run_ulm(mat = RNA_weight_input, network = TF_regulons, minsize = 20)
TF_scores_factors <- dcast(TF_scores, source~condition, value.var = "score")
row.names(TF_scores_factors) <- TF_scores_factors$source

row.names(TF_scores_factors) <- paste(row.names(TF_scores_factors),"_Hepatocytes_TFs",sep = "")

weights$Hepatocytes_TFs <- as.matrix(TF_scores_factors[,-1])

save(weights, file = "results/mofa/main_model_RNA_lipid_weight.RData")
