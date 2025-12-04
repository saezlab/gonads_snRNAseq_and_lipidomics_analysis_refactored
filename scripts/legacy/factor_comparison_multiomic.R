library(ggplot2)
library(readr)

PCA_sample_coord_lipid <- as.data.frame(
  read_csv("results/lipidomic/PCA_sample_coord_lipid.csv"))

PCA_sample_coord_RNA_HEPA <- as.data.frame(
  read_csv("results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv"))
targets <- as.data.frame(cbind(PCA_sample_coord_RNA_HEPA$samples, gsub(".*_Hepatocytes_","",PCA_sample_coord_RNA_HEPA$samples)))
targets$V1 <-  gsub("_Hepatocytes_.*","",targets$V1)
PCA_sample_coord_RNA_HEPA$samples <- gsub("_Hepatocytes_.*","",PCA_sample_coord_RNA_HEPA$samples)


n_factor <- 1
merge_components <- merge(PCA_sample_coord_lipid[,c(n_factor,48)], PCA_sample_coord_RNA_HEPA[,c(n_factor,44)], by = "samples")
names(merge_components) <- c("samples","lipid_factor","RNA_factor")
merge_components <- merge(merge_components, targets, by.x = "samples", by.y = "V1")

ggplot(merge_components, aes(x = lipid_factor, y = RNA_factor)) + geom_point()

data <- merge_components
# Create the scatter plot
ggplot(data, aes(x = lipid_factor, y = RNA_factor, color = V2)) +
  geom_point(size = 4) +
  
  # Linear regression for points lower than 0 (both axes)
  geom_smooth(data = subset(data, lipid_factor < 0 & RNA_factor < 0),
              method = "lm", se = FALSE) +
  
  # Linear regression for points higher than 0 (both axes)
  geom_smooth(data = subset(data, lipid_factor > 0 & RNA_factor > 0),
              method = "lm", se = FALSE) +
  
  # Labels and theme
  labs(title = "First component of PCA from RNA and lipid data",
       x = "Lipid Factor 1",
       y = "Hepatocytes RNA Factor 1") +
  theme_minimal()

# Subset the data for the two groups
group_below_0 <- subset(data, lipid_factor < 0 & RNA_factor < 0)
group_above_0 <- subset(data, lipid_factor > 0 & RNA_factor > 0)

# Calculate the correlation coefficients for both groups
cor_below_0 <- cor(group_below_0$lipid_factor, group_below_0$RNA_factor)
cor_above_0 <- cor(group_above_0$lipid_factor, group_above_0$RNA_factor)

cor.test(group_below_0$lipid_factor, group_below_0$RNA_factor)
cor.test(group_above_0$lipid_factor, group_above_0$RNA_factor)

# Print the correlation coefficients
cat("Correlation coefficient for group below 0:", cor_below_0, "\n")
cat("Correlation coefficient for group above 0:", cor_above_0, "\n")
