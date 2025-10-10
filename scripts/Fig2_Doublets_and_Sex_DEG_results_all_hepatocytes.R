library(readxl)
library(readr)
library(ggplot2)
library(ggrepel)

source("scripts/support_functions.R")

ttop <- as.data.frame(
  read_excel("results/Doublets_and_Sex_DEG_results_all_hepatocytes.xlsx", , 
             sheet = "All_Hepatocytes_Sex"))

X9WZDEU_mMulti_100g_gene_list <- as.data.frame(read_csv("support/FCGL_Xenium_Panel.csv"))
X9WZDEU_mMulti_100g_gene_list <- X9WZDEU_mMulti_100g_gene_list[X9WZDEU_mMulti_100g_gene_list$Reason != "cell type marker",]

ttop <- ttop[ttop$names %in% X9WZDEU_mMulti_100g_gene_list$Gene,]
names(ttop)[1] <- "ID"

c("Slco1a1", "C8a", "Cyp2b9","Nipal1")

ttop$pvals <- ifelse(ttop$pvals < 10^-295, 10^-300, ttop$pvals)

labels <- ifelse((-log10(ttop$pvals) < 300 & -log10(ttop$pvals) > 100) | ttop$ID %in% c("Slco1a1", "C8a", "Cyp2b9","Nipal1"), ttop$ID, NA)

volcano_nice <- function (df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                          label = FALSE, straight = FALSE, nlabels, manual_labels = NA, col_high_neut_low = c("red","grey30","royalblue"))
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log10(hAss)
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log10(df[, pValIndex]))))
  }
  else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log10(df[, pValIndex])))
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log10(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log10(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  if (straight) {
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$adj.P.Val <=
                           hAssOri, "1", "0")
  }
  else {
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                      as.numeric(x[pValIndex]), vAss))
  }
  
  df <- df[order(df$adj.P.Val, decreasing = F), ]
  df$condLabel <- df[, IDIndex]
  # df[df$couleur == "0", "condLabel"] <- NA
  labels_to_keep <- c(manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  
  print(table(df$couleur))
  
  if (label) {
    a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = couleur)) +
      geom_point(alpha = 0.5) + geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) + ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
      scale_colour_manual(values = c(col_high_neut_low[2], col_high_neut_low[3],
                                     col_high_neut_low[1])) + theme_minimal() + theme(legend.position = "none")
  }
  else {
    if (straight) {
      a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        geom_vline(xintercept = -vAss, color = col_high_neut_low[3]) +
        geom_vline(xintercept = vAss, color = col_high_neut_low[3]) +
        ylim(c(0, ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = col_high_neut_low[1]) +
        scale_colour_manual(values = c(col_high_neut_low[2], col_high_neut_low[1],
                                       col_high_neut_low[3])) + theme_minimal() + theme(legend.position = "none")
    }
    else {
      a <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        stat_function(fun = xneg, xlim = c(-xlimAbs,
                                           -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                          ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
        scale_colour_manual(values = c(col_high_neut_low[2], col_high_neut_low[1],
                                       col_high_neut_low[3])) + theme_minimal() + theme(legend.position = "none")
    }
  }
  return(a)
}

volcano_nice(ttop, IDIndex = 1, FCIndex = 3, pValIndex = 4, nlabels = 100, label = T, manual_labels = na.omit(labels)) + ylab("-log10(p-value)") + coord_flip()


