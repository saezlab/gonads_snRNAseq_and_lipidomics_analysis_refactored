#!/usr/bin/env Rscript
# Factor 1 hepatocyte TF -> measured lipid MOON pilot.
# Run from the repository root after prepare_factor1_inputs.R.
options(stringsAsFactors = FALSE)
source("scripts/support_functions.R")
required <- c("cosmosR", "decoupleR", "igraph", "jsonlite")
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
input_dir <- "results/cosmos_factor1_pilot/inputs"
use_compression <- !"--no-compression" %in% commandArgs(TRUE)
out <- if (use_compression) "results/cosmos_factor1_pilot" else
  "results/cosmos_factor1_pilot/uncompressed_sensitivity"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if ("--validate-setup" %in% commandArgs(TRUE)) {
  cat("Runner parsed; installed packages:", paste(required, collapse = ", "), "\n")
  quit(save = "no")
}
source("scripts/moon/moon_safe_helpers.R")
input <- readRDS(file.path(input_dir, "factor1_inputs.rds"))
needed <- c("tf_activity", "lipid_weight", "hepatocyte_gene_weight", "tf_regulon")
stopifnot(all(needed %in% names(input)))
validate_vector <- function(x) {
  stopifnot(is.numeric(x), length(x) > 0, !is.null(names(x)),
            !anyDuplicated(names(x)), all(nzchar(names(x))), all(is.finite(x)))
}
invisible(lapply(input[needed[1:3]], validate_vector))
pkn_path <- "results/cosmos_mouse_pkn/pkn_lipid_rules.rds"
pkn <- cosmosR::meta_network_cleanup(readRDS(pkn_path))
pkn <- pkn[, c("source", "interaction", "target")]
stopifnot(!anyDuplicated(pkn[, c("source", "target")]),
          all(pkn$interaction %in% c(-1, 1)), all(pkn$source != pkn$target))
node_ids <- function(e) unique(c(e$source, e$target))
edge_key <- function(e) paste(e$source, e$interaction, e$target, sep = "\t")
write_csv <- function(x, f) utils::write.csv(x, file.path(out, f), row.names = FALSE, na = "")
reach <- function(e, seeds, direction = "out", depth = Inf) {
  if (direction == "in") e <- data.frame(source = e$target, target = e$source)
  seen <- intersect(seeds, node_ids(e)); frontier <- seen; step <- 0L
  while (length(frontier) && step < depth) {
    new <- setdiff(unique(e$target[e$source %in% frontier]), seen)
    seen <- c(seen, new); frontier <- new; step <- step + 1L
  }
  seen
}
n_layers <- 6L
tf_all <- input$tf_activity
tf <- tf_all[abs(tf_all) > 2 & names(tf_all) %in% node_ids(pkn)]
lipid <- input$lipid_weight[input$lipid_weight != 0 & names(input$lipid_weight) %in% node_ids(pkn)]
rna <- input$hepatocyte_gene_weight
reg <- input$tf_regulon[, c("source", "target", "mor")]
stopifnot(length(tf) > 0, length(lipid) > 2,
          all(grepl("^Gene__", names(tf))), all(grepl("^Metab__", names(lipid))),
          all(is.finite(reg$mor)), all(reg$mor %in% c(-1, 1)))
cat("PKN:", nrow(pkn), "edges; TF seeds:", length(tf),
    "; mapped lipid node observations:", length(lipid), "\n")
# Match moon(): observations are terminal. Do this BEFORE reachability pruning.
net <- pkn[!pkn$source %in% names(lipid), ]
pruning <- list()
repeat {
  old <- nrow(net)
  forward <- reach(net, names(tf), depth = n_layers)
  net <- net[net$source %in% forward & net$target %in% forward, ]
  backward <- reach(net, names(lipid), direction = "in", depth = n_layers)
  net <- net[net$source %in% backward & net$target %in% backward, ]
  tf <- tf[names(tf) %in% node_ids(net)]
  lipid <- lipid[names(lipid) %in% node_ids(net)]
  pruning[[length(pruning) + 1L]] <- data.frame(iteration = length(pruning) + 1L,
      edges = nrow(net), nodes = length(node_ids(net)), upstream = length(tf), downstream = length(lipid))
  if (old == nrow(net)) break
}
stopifnot(nrow(net) > 0, length(tf) > 0, length(lipid) > 2)
write_csv(do.call(rbind, pruning), "pruning_ledger.csv")
write_csv(data.frame(source = names(tf), score = unname(tf)), "upstream_input.csv")
write_csv(data.frame(source = names(lipid), score = unname(lipid)), "downstream_input.csv")
saveRDS(net, file.path(out, "pkn_pruned_uncompressed.rds"))
# Preserve all identifiers touched by the RNA coherence test, plus inputs.
protected <- unique(c(names(tf), names(lipid), reg$source, reg$target))
compressed <- moon_safe_compress(net, if (use_compression) protected else node_ids(net))
score_net <- compressed$pkn_scoring
map <- compressed$node_map
write_csv(map, "compression_mapping.csv")
saveRDS(compressed, file.path(out, "compression.rds"))
cat("Pruned:", nrow(net), "edges /", length(node_ids(net)), "nodes; compressed:",
    nrow(score_net), "edges /", length(node_ids(score_net)), "nodes\n")
score_once <- function(e) {
  ans <- cosmosR::moon(upstream_input = tf, downstream_input = lipid,
       meta_network = e, n_layers = n_layers, statistic = "ulm")
  ans <- as.data.frame(ans[, c("source", "score", "level")])
  stopifnot(!anyDuplicated(ans$source), all(is.finite(ans$score)),
            all(names(lipid) %in% ans$source[ans$level == 0]))
  ans
}
dir.create(file.path(out, "iterations"), showWarnings = FALSE)
ledger <- list(); converged <- FALSE; removed <- list()
for (iteration in seq_len(10L)) {
  used <- score_net
  cat("MOON coherence iteration", iteration, "with", nrow(used), "edges\n")
  scored <- score_once(used)
  filtered <- cosmosR::filter_incohrent_TF_target(scored, reg, used, rna)
  delta <- used[!edge_key(used) %in% edge_key(filtered), ]
  if (nrow(delta)) { delta$iteration <- iteration; removed[[length(removed) + 1L]] <- delta }
  ledger[[iteration]] <- data.frame(iteration = iteration, edges_scored = nrow(used),
     edges_removed = nrow(delta), scored_nodes = nrow(scored),
     upstream_scored = sum(names(tf) %in% scored$source), max_level = max(scored$level))
  saveRDS(list(network = used, scores = scored),
          file.path(out, "iterations", sprintf("iteration_%02d.rds", iteration)))
  score_net <- filtered
  if (identical(edge_key(used), edge_key(filtered))) { converged <- TRUE; break }
}
if (!converged) scored <- score_once(score_net)
write_csv(do.call(rbind, ledger), "coherence_ledger.csv")
write_csv(if (length(removed)) do.call(rbind, removed) else
  data.frame(source = character(), interaction = numeric(), target = character(), iteration = integer()),
  "coherence_removed_edges.csv")
saveRDS(score_net, file.path(out, "pkn_final_compressed.rds"))
write_csv(scored, "moon_scores_compressed.csv")
# Use an explicit complete mapping (also handles the zero-compression case).
scores <- merge(map, scored, by.x = "compressed", by.y = "source")
scores <- scores[, c("original", "score", "level")]
names(scores)[1] <- "source"
stopifnot(!anyDuplicated(scores$source), all(is.finite(scores$score)))
lookup <- setNames(map$compressed, map$original)
mapped <- net
mapped$source <- unname(lookup[net$source]); mapped$target <- unname(lookup[net$target])
final <- net[edge_key(mapped) %in% edge_key(score_net), ]
stopifnot(!any(final$source %in% names(lipid)))
saveRDS(final, file.path(out, "pkn_final_uncompressed.rds"))
write_csv(scores, "moon_scores.csv")
catalog <- utils::read.csv("results/cosmos_mouse_pkn/nodes.csv", check.names = FALSE)
score_annotate <- function(x) {
  x$type <- ifelse(x$source %in% names(tf), "upstream_input",
                   ifelse(x$level == 0, "level0", "other"))
  x$RNA_input <- unname(rna[x$source])
  x$upstream_activity <- unname(tf_all[x$source])
  merge(x, catalog, by.x = "source", by.y = "node_id", all.x = TRUE, sort = FALSE)
}
write_csv(score_annotate(scores), "moon_scores_annotated.csv")
settings <- data.frame(name = c("relaxed", "primary", "strict"),
                       primary = c(1, 1.7, 2), secondary = c(0, 1, 1.5))
sensitivity <- list()
for (k in seq_len(nrow(settings))) {
  setting <- settings[k, ]
  sol <- moon_extract_paths(scores, final, tf, setting$primary, setting$secondary)
  sif <- sol$SIF; att <- score_annotate(sol$ATT)
  # Independent invariant checks: all retained edges and nodes support a signed path.
  sm <- setNames(scores$score, scores$source)
  stopifnot(all(edge_key(sif) %in% edge_key(final)),
            all(sign(sm[sif$source] * sm[sif$target]) == sif$interaction),
            all(node_ids(sif) %in% reach(sif, names(tf))),
            all(node_ids(sif) %in% reach(sif, names(lipid), direction = "in")),
            !any(sif$source %in% names(lipid)),
            all(sign(att$score[att$type == "upstream_input"]) == sign(att$upstream_activity[att$type == "upstream_input"])))
  write_csv(sif, paste0("solution_", setting$name, "_SIF.csv"))
  write_csv(att, paste0("solution_", setting$name, "_ATT.csv"))
  sensitivity[[k]] <- data.frame(setting, nodes = nrow(att), edges = nrow(sif),
       upstream_tfs = sum(att$type == "upstream_input"), downstream_nodes = sum(att$level == 0))
}
write_csv(do.call(rbind, sensitivity), "threshold_sensitivity.csv")
upstream_audit <- data.frame(source = names(tf_all), input_activity = unname(tf_all),
  passes_abs2 = abs(tf_all) > 2, in_pkn = names(tf_all) %in% node_ids(pkn),
  reachable = names(tf_all) %in% names(tf), scored = names(tf_all) %in% scores$source,
  moon_score = scores$score[match(names(tf_all), scores$source)])
write_csv(upstream_audit, "upstream_coverage.csv")
write_csv(as.data.frame(table(scores$level)), "scoring_levels.csv")
manifest_paths <- c(pkn_path, "scripts/moon/run_factor1_moon.R", "scripts/moon/moon_safe_helpers.R",
                     file.path(input_dir, "factor1_inputs.rds"))
manifest <- list(created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  command = paste("Rscript scripts/moon/run_factor1_moon.R", if (use_compression) "" else "--no-compression"),
  compression_enabled = use_compression, cosmosR = as.character(utils::packageVersion("cosmosR")),
  decoupleR = as.character(utils::packageVersion("decoupleR")), statistic = "ulm", n_layers = n_layers,
  tf_cutoff = "abs(activity) > 2", expression_filter = "none: RNA weights are not expression presence measurements",
  downstream_scaling = "none; level0 weights exempt from inferred-score secondary cutoff, not primary seeds",
  convergence = converged, iterations = length(ledger), upstream_inputs = length(tf), downstream_inputs = length(lipid),
  scored_upstream = sum(names(tf) %in% scores$source), scored_nodes = nrow(scores),
  orientation = input$orientation, input_and_script_md5 = as.list(tools::md5sum(manifest_paths)),
  template_url = "https://github.com/saezlab/COSMOS_basic/blob/7108c6b849c712b1282440dc122c9d783f563da8/scripts/net_compr_MOON.R",
  threshold_sensitivity = do.call(rbind, sensitivity))
jsonlite::write_json(manifest, file.path(out, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Completed MOON pilot. Converged:", converged, "\n")
print(do.call(rbind, sensitivity), row.names = FALSE)
