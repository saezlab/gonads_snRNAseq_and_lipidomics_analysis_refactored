library(readr)
library(readxl)
library(uwot)
library(cluster)



source("scripts/support_functions.R")


Suppl_Table_1_OverviewData <- as.data.frame(read_excel("data/Suppl_Table_1_OverviewData.xlsx"))
Suppl_Table_1_OverviewData <- Suppl_Table_1_OverviewData[,-4]

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

scaled_RNA <- scaled_RNA[,grepl("Hepa",names(scaled_RNA))]

# Function for UMAP grid search
umap_grid_search <- function(data_matrix, eval_function, 
                             n_neighbors_range = 5:50, 
                             min_dist_range = seq(0.01, 0.5, by = 0.05),
                             n_components_range = 2:15, 
                             max_configs = 200, 
                             seed = 42) {
  
  set.seed(seed) # For reproducibility
  
  # Generate all possible combinations of parameters
  param_grid <- expand.grid(
    n_neighbors = n_neighbors_range,
    min_dist = min_dist_range,
    n_components = n_components_range
  )
  
  # Randomly sample up to max_configs configurations
  if (nrow(param_grid) > max_configs) {
    param_grid <- param_grid[sample(1:nrow(param_grid), max_configs), ]
  }
  
  # Initialize results storage
  results <- data.frame(
    n_neighbors = integer(),
    min_dist = numeric(),
    n_components = integer(),
    eval_score = numeric()
  )
  
  # Perform UMAP with each parameter set
  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]
    message(sprintf("Running UMAP with config %d/%d: n_neighbors=%d, min_dist=%.2f, n_components=%d",
                    i, nrow(param_grid), params$n_neighbors, params$min_dist, params$n_components))
    
    # Run UMAP
    umap_result <- umap(
      data_matrix,
      n_neighbors = params$n_neighbors,
      min_dist = params$min_dist,
      n_components = params$n_components
    )
    
    # Evaluate result
    score <- eval_function(umap_result, data_matrix)
    
    # Save results
    results <- rbind(results, cbind(params, eval_score = score))
  }
  
  # Find the best configuration
  best_config <- results[which.max(results$eval_score), ]
  
  list(
    best_config = best_config,
    all_results = results
  )
}


# Example evaluation function: Mean silhouette width
eval_function <- function(umap_result, original_data) {
  cluster_labels <- kmeans(umap_result, centers = 3)$cluster
  silhouette_score <- mean(silhouette(cluster_labels, dist(umap_result))[, 3])
  return(silhouette_score)
}

# Example data
set.seed(42)

# Perform grid search
grid_results <- umap_grid_search(
  data_matrix = t(scaled_RNA),
  eval_function = eval_function,
  n_neighbors_range = 5:20, 
  min_dist_range = seq(0.05, 0.3, by = 0.05),
  n_components_range = 2:10, 
  max_configs = 200
)

# View the best configuration
print(grid_results$best_config)

# View all results
head(grid_results$all_results)
# Plot UMAP
plot(umap_coords, main = "UMAP Projection")



# Compute correlations
correlations <- cor(t(scaled_RNA), umap_coords)

