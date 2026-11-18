# Synthetic plotting smoke test. These fixture values are not study results.
if (dir.exists(".local/R-library")) .libPaths(c(".local/R-library", .libPaths()))
source("scripts/core/cf453_figure6_plots.R")
set.seed(1406)
metadata <- expand.grid(replicate = seq_len(7), sex = c("female", "male"), diet = c("SD", "WD"), stringsAsFactors = FALSE)
metadata$mouse_id <- paste0("fixture_", seq_len(nrow(metadata)))
metadata$group <- paste(metadata$sex, metadata$diet, sep = "_")
genes <- paste0("Fixture", LETTERS[1:12])
zonation <- expand.grid(gene = genes, diet = c("SD", "WD"), stringsAsFactors = FALSE)
zonation$t_female <- stats::rnorm(nrow(zonation), 0, 5)
zonation$t_male <- stats::rnorm(nrow(zonation), 0, 5)
zonation$effect_female <- zonation$t_female / 4
zonation$effect_male <- zonation$t_male / 4
zonation$q_interaction <- seq(.0001, .2, length.out = nrow(zonation))
zonation$category <- rep(c("female_PP", "female_PC", "male_PP", "male_PC", "shared", "unselected"), 4)
classes <- c("PC", "PE", "TG")
class_summary <- expand.grid(diet = c("SD", "WD"), lipid_class = classes,
                             category = c("female_PP", "female_PC", "male_PP", "male_PC", "shared", "unassigned"), stringsAsFactors = FALSE)
class_summary$n <- 2
class_summary$total <- 12
class_summary$proportion <- class_summary$n / 12
representatives <- expand.grid(sex = c("female", "male"), zone = c("periportal", "pericentral"), stringsAsFactors = FALSE)
representatives$gene <- genes[seq_len(4)]
representatives$rank <- 1
representatives$stable_degree <- 6
representatives$q_interaction <- .001
representative_edges <- merge(representatives[, c("sex", "zone", "gene")],
                              expand.grid(partner_rank = seq_len(6), diet = c("SD", "WD"), stringsAsFactors = FALSE))
representative_edges$lipid <- paste0("PE(O-", c("16:0", "18:0", "18:1", "20:1", "22:1", "24:2")[representative_edges$partner_rank], "_22:6)")
representative_edges$r <- rep(c(.78, -.72, .68, -.61, .88, -.77), length.out = nrow(representative_edges))
representative_edges$n <- 7
representative_edges$lipid_class <- "PE O-"
representative_edges$stable_display <- TRUE
examples <- merge(representative_edges[, c("sex", "zone", "gene", "lipid", "diet")],
                  data.frame(replicate = seq_len(7)))
examples$mouse_id <- paste0("fixture_", examples$replicate)
examples$lipid_label <- examples$lipid
examples$rna <- stats::rnorm(nrow(examples))
examples$lipid_log2 <- 18 + examples$rna + stats::rnorm(nrow(examples), 0, .5)
zonation_diet_effects <- expand.grid(gene = genes, sex = c("female", "male"), stringsAsFactors = FALSE)
zonation_diet_effects$effect <- stats::rnorm(nrow(zonation_diet_effects))
zonation_diet_effects$se <- .3
zonation_diet_effects$low <- zonation_diet_effects$effect - .6
zonation_diet_effects$high <- zonation_diet_effects$effect + .6
zonation_diet_effects$p <- stats::runif(nrow(zonation_diet_effects))
zonation_diet_effects$q <- rep(c(.01, .3), length.out = nrow(zonation_diet_effects))
lipid_results <- data.frame(lipid = paste0("PC(", seq_len(20), ":0_18:1)"), lipid_class = "PC",
                            contrast = "male_minus_female_SD", effect = stats::rnorm(20),
                            q = stats::runif(20), stringsAsFactors = FALSE)
qc_pca <- metadata
qc_pca$kind <- "lipid"
qc_pca$profile_id <- qc_pca$mouse_id
qc_pca$zone <- NA_character_
qc_pca$PC1 <- stats::rnorm(28)
qc_pca$PC2 <- stats::rnorm(28)
qc_pca$variance_PC1 <- .22
qc_pca$variance_PC2 <- .13
thresholds <- expand.grid(diet = c("SD", "WD"), threshold = c(.5, .65, .8), lipid_class = classes, stringsAsFactors = FALSE)
thresholds$assigned <- 2
thresholds$total <- 10
class_pair_density <- expand.grid(diet = c("SD", "WD"), sex = c("female", "male"),
                                  zone = c("PP", "PC"), lipid_class = classes, stringsAsFactors = FALSE)
class_pair_density$total_pairs <- 120L
class_pair_density$display_pairs <- seq_len(nrow(class_pair_density))
class_pair_density$stable_pairs <- class_pair_density$display_pairs - 1L
class_pair_density$fdr05_pairs <- 0L
class_pair_density$display_fraction <- class_pair_density$display_pairs / class_pair_density$total_pairs
class_pair_density$stable_fraction <- class_pair_density$stable_pairs / class_pair_density$total_pairs
class_pair_density$fdr05_fraction <- 0
quality_sensitivity <- expand.grid(diet = c("SD", "WD"), quality = c("keep_broad", "keep_primary", "keep_fill"),
                                   lipid_class = classes, stringsAsFactors = FALSE)
quality_sensitivity$assigned <- 7L
quality_sensitivity$total <- 10L
primary_interactions <- data.frame(gene = genes, contrast = "sex_diet_zone", q = seq(.0001, .2, length.out = 12))
fixture <- list(metadata = metadata, zonation = zonation, class_summary = class_summary,
                representatives = representatives, representative_edges = representative_edges,
                examples = examples, zonation_diet_effects = zonation_diet_effects,
                lipid_results = lipid_results, qc_pca = qc_pca, thresholds = thresholds,
                class_pair_density = class_pair_density, quality_sensitivity = quality_sensitivity,
                primary_interactions = primary_interactions)
output <- Sys.getenv("CF6_PLOT_SMOKE_OUTPUT", file.path(tempdir(), "cf453_plot_smoke"))
manifest <- cf6_make_panels(fixture, output)
stopifnot(all(file.exists(manifest$file)), all(file.info(manifest$file)$size > 100),
          all(table(manifest$panel) == 3L), all(file.exists(manifest$source_table)))
network <- manifest[manifest$family == "D", ]
stopifnot(length(unique(network$width_in)) == 1L, length(unique(network$height_in)) == 1L)
empty <- cf6_plot_network(representative_edges[FALSE, ], "female", "PC", "SD")
stopifnot(nrow(empty$source) == 0L, inherits(empty$plot, "ggplot"))
cat("CF453 plotting smoke passed:", length(unique(manifest$panel)), "panels in", output, "\n")
