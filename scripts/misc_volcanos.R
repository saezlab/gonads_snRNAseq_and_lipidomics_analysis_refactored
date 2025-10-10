library(readxl)
library(ggplot2)
library(ggrepel)

source("scripts/support_functions.R")

Doublets_and_Sex_DEG_results_all_hepatocytes <- as.data.frame(
  read_excel("results/Doublets_and_Sex_DEG_results_all_hepatocytes.xlsx", 
                                                           sheet = "All_Hepatocytes_Sex"))

Doublets_and_Sex_DEG_results_all_hepatocytes$pvals <- ifelse(Doublets_and_Sex_DEG_results_all_hepatocytes$pvals == 0, 10e-301,Doublets_and_Sex_DEG_results_all_hepatocytes$pvals)
volcano_nice(Doublets_and_Sex_DEG_results_all_hepatocytes, IDIndex = 1, FCIndex = 3, pValIndex = 4, nlabels = 100, label = T)


DGE_hepatocytes_results <- as.data.frame(
  read_excel("results/DGE_hepatocytes_results.xlsx", 
                                      sheet = "Portal_vs_Central"))
DGE_hepatocytes_results$pvals <- ifelse(DGE_hepatocytes_results$pvals == 0, 10e-299,DGE_hepatocytes_results$pvals)

volcano_nice(DGE_hepatocytes_results, IDIndex = 1, FCIndex = 3, pValIndex = 4, nlabels = 100, label = T)

