library(ggplot2)
df <- read.csv("data/metabolomics/TG.csv")

df$abs_logFC <- df$X.LogFC.

ggplot(df, aes(x = C, y = DB)) +
  geom_point(size = 2,aes(color = LogFC)) +
  scale_color_gradient2(low = "royalblue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 6)) +
  labs(
    title = "Lipid class: TG, Comparison: O vs. T",
    x = "C",
    y = "DB",
    color = "LogFC"
  ) +
  ylim(c(0,15)) +
  xlim(c(30,65)) +
  theme_minimal()


df <- read.csv("data/metabolomics/SM.csv")

df$abs_logFC <- df$X.LogFC.

ggplot(df, aes(x = C, y = DB)) +
  geom_point(size = 2,aes(color = LogFC)) +
  scale_color_gradient2(low = "royalblue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 6)) +
  labs(
    title = "Lipid class: SM, Comparison: O vs. T",
    x = "C",
    y = "DB",
    color = "LogFC"
  ) +
  ylim(c(0,15)) +
  xlim(c(30,65)) +
  theme_minimal()

all_lipid_class_files <- list.files("data/metabolomics/for_bubbles/")

all_lipid_class_files <- paste("data/metabolomics/for_bubbles/",all_lipid_class_files,sep = "")

for(lipid_class_file in all_lipid_class_files)
{
  df <- read.csv(lipid_class_file)
  
  lipid_class <- gsub("data/.*/","",gsub(".csv$","",lipid_class_file))
  
  df$abs_logFC <- df$X.LogFC.
  
  bubbleplot <- ggplot(df, aes(x = C, y = DB)) +
    geom_point(size = 2,aes(color = LogFC)) +
    scale_color_gradient2(low = "royalblue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(range = c(1, 6)) +
    labs(
      title = "Lipid class: SM, Comparison: O vs. T",
      x = "C",
      y = "DB",
      color = "LogFC"
    ) +
    ylim(c(0,25)) +
    xlim(c(0,65)) +
    theme_minimal()
  
  ggsave(filename = paste("results/lipidomic/bubbles_DBC/",paste(lipid_class,".pdf", sep = ""), sep = ""),bubbleplot, width = 8.42, height = 5.53)

  # Save plot data for testing
  if (Sys.getenv("TESTING_MODE") == "TRUE") {
    write.csv(df, paste0("tests/outputs/lipidomic_bubbles_", lipid_class, "_data.csv"), row.names = FALSE)
  }
}
