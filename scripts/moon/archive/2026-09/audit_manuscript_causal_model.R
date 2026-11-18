#!/usr/bin/env Rscript
# Audit frozen manuscript results after relocation; no MOON/model refitting.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
base <- "results/causal_modelling"
review <- file.path(base, "cosmos_factor1_compressed_review")
run <- file.path(base, "cosmos_factor1_parameter_search/runs/radius6_depth4_compressed")
candidate <- file.path(review, "threshold_search/C03")
out <- file.path(base, "manuscript_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
edge_key <- function(x) paste(x$source, x$interaction, x$target, sep = "\t")
node_ids <- function(x) unique(c(x$source, x$target))

pkn_path <- file.path(base, "cosmos_mouse_pkn/pkn_lipid_rules.rds")
input_path <- file.path(base, "cosmos_factor1_pilot/inputs/factor1_inputs.rds")
pkn <- readRDS(pkn_path)
inputs <- readRDS(input_path)
state <- readRDS(file.path(run, "run_state.rds"))
compression_path <- file.path(base, "cosmos_factor1_parameter_search/topology_radius_6/compressed/compression.rds")
comp <- readRDS(compression_path)
scoring <- readRDS(file.path(run, "pkn_final_scoring.rds"))
final <- readRDS(file.path(run, "pkn_final_uncompressed.rds"))
manifest <- jsonlite::read_json(file.path(run, "run_manifest.json"), simplifyVector = TRUE)
sif <- read_csv(file.path(candidate, "network_SIF.csv"))
att <- read_csv(file.path(candidate, "network_ATT.csv"))
cat("PKN class:", paste(class(pkn), collapse = ","), "; dimensions:", nrow(pkn), "x", ncol(pkn), "\n")
cat("Frozen run:", state$run_id, "; candidate nodes/edges:", nrow(att), "/", nrow(sif), "\n")

stopifnot(nrow(pkn) == 254589L, length(node_ids(pkn)) == 52755L,
          nrow(sif) == 388L, nrow(att) == 153L,
          identical(state$run_id, "radius6_depth4_compressed"),
          state$n_layers == 4L, state$topology_radius == 6L,
          state$converged, state$iterations == 3L,
          length(state$tf) == 49L, length(state$lipid) == 74L,
          !anyDuplicated(edge_key(sif)), !anyDuplicated(att$source),
          all(pkn$interaction %in% c(-1, 1)), !anyNA(pkn),
          !any(pkn$source == pkn$target))

# Historical file names identify hashes; relocated files are the actual inputs.
expected <- manifest$input_and_script_md5
stopifnot(unname(tools::md5sum(pkn_path)) == expected[["results/cosmos_mouse_pkn/pkn_lipid_rules.rds"]],
          unname(tools::md5sum(input_path)) == expected[["results/cosmos_factor1_pilot/inputs/factor1_inputs.rds"]],
          unname(tools::md5sum("scripts/moon/moon_safe_helpers.R")) == expected[["scripts/moon/moon_safe_helpers.R"]])

# Original-ID extraction graph must project exactly to the compressed scored graph.
mv <- stats::setNames(comp$node_map$compressed, comp$node_map$original)
projected <- final
projected$source <- unname(mv[final$source])
projected$target <- unname(mv[final$target])
stopifnot(!anyNA(projected),
          identical(moon_check_pkn(projected), moon_check_pkn(scoring)),
          identical(moon_check_pkn(final), moon_check_pkn(state$pkn_final)),
          all(mv[names(state$tf)] == names(state$tf)),
          all(mv[names(state$lipid)] == names(state$lipid)),
          comp$n_merged_nodes > 0L)

# Re-extract from frozen scores without changing the scientific scoring state.
replay <- moon_extract_paths(state$scores, final, state$tf,
                            primary_thresh = 1.85, secondary_thresh = 0.55,
                            RNA_input = inputs$hepatocyte_gene_weight)
stopifnot(setequal(edge_key(replay$SIF), edge_key(sif)),
          nrow(replay$SIF) == nrow(sif), setequal(replay$ATT$source, att$source),
          replay$audit$all_edges_consistent,
          replay$audit$all_nodes_on_seed_to_observation_walk,
          replay$audit$retained_upstream == 4L, replay$audit$retained_observed == 3L)
ii <- match(att$source, state$scores$source)
stopifnot(!anyNA(ii), max(abs(att$score - state$scores$score[ii])) < 1e-10,
          identical(as.integer(att$level), as.integer(state$scores$level[ii])))

all_scores <- read_csv(file.path(review, "score_retention/radius6_depth4_compressed_scores_all.csv"))
focus <- all_scores[match(paste0("Gene__", c("Esr1", "Nfic", "Dbp")), all_scores$source), ]
stopifnot(nrow(all_scores) == 2238L, sum(all_scores$returned_by_moon) == 2207L,
          sum(all_scores$retention_reason == "upstream_sign_conflict") == 31L,
          all(focus$score > 0), all(focus$upstream_activity < 0),
          all(!focus$returned_by_moon), all(focus$source %in% node_ids(scoring)))
grid_baseline <- read_csv(file.path(review, "threshold_search/grid_complete.csv"))
grid_sensitivity <- read_csv(file.path(review, "sensitivity_threshold_search/grid_complete.csv"))
stopifnot(nrow(grid_baseline) == 2376L, nrow(grid_sensitivity) == 5187L)
matched <- read_csv(file.path(review, "sensitivity_threshold_search/matched_threshold_sensitivity.csv"))
same <- matched[matched$run_id == "radius6_depth4_compressed_inputs_only" &
                matched$primary == 1.85 & matched$secondary == 0.55, ]
stopifnot(nrow(same) == 1L, same$edges == 388L,
          same$exact_edge_set_equals_baseline_C03, same$edge_jaccard_to_baseline_C03 == 1)

artifact_paths <- c(pkn_path, input_path, compression_path, file.path(run, c("run_state.rds",
  "pkn_final_scoring.rds", "pkn_final_uncompressed.rds", "run_manifest.json")),
  file.path(candidate, c("network_SIF.csv", "network_ATT.csv", "network.graphml")),
  "scripts/moon/moon_safe_helpers.R", "scripts/moon/audit_manuscript_causal_model.R")
utils::write.csv(data.frame(path = artifact_paths, md5 = unname(tools::md5sum(artifact_paths))),
                 file.path(out, "artifact_manifest.csv"), row.names = FALSE)
summary <- list(audited_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  operation = "frozen-state audit and exact solution re-extraction; no model or MOON refitting",
  run_id = state$run_id, candidate = "cosmos_factor1_compressed_review/threshold_search/C03",
  prior = list(edges = nrow(pkn), nodes = length(node_ids(pkn)), archived_hash_matches = TRUE),
  inputs = list(tf_inputs = length(state$tf), lipid_inputs = length(state$lipid), archived_hash_matches = TRUE),
  compression = list(initial_edges = nrow(comp$pkn_scoring), merged_original_nodes = comp$n_merged_nodes,
                     groups = comp$n_groups, final_projection_matches = TRUE),
  scoring = list(pruning_distance = 6L, n_layers = 4L, iterations = state$iterations,
                 final_scoring_edges = nrow(scoring), final_original_edges = nrow(final)),
  extraction = replay$audit, exact_reextraction_matches_saved = TRUE,
  all_scores = list(before_sign_filter = nrow(all_scores), returned = sum(all_scores$returned_by_moon),
                    upstream_sign_conflicts = sum(all_scores$retention_reason == "upstream_sign_conflict")),
  focus_tf_scores = focus[, c("source", "score", "upstream_activity", "retention_reason")],
  threshold_evaluations = nrow(grid_baseline) + nrow(grid_sensitivity),
  stronger_compression_reproduces_selected_graph = TRUE,
  current_versions = list(R = as.character(getRversion()), cosmosR = as.character(utils::packageVersion("cosmosR")),
                          decoupleR = as.character(utils::packageVersion("decoupleR"))),
  relocation_note = "Historical scripts/manifests/RDS path strings may reference results/cosmos_*; this audit resolves relocated artifacts explicitly. This is not an end-to-end replay of historical runners.")
jsonlite::write_json(summary, file.path(out, "audit_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines(c("# Manuscript causal-model audit", "",
  "Run from the repository root: `Rscript scripts/moon/audit_manuscript_causal_model.R`.", "",
  "All checks passed. The relocated PKN and input RDS hashes match the archived run manifest; the original-ID final graph projects exactly to the compressed scored graph. Re-extracting the selected network from its frozen scores reproduces all 388 signed edges and 153 nodes, with four TFs and three measured DGs. Input signs, score-retention counts, threshold-grid counts and the saved stronger-compression comparison also agree.", "",
  "This audit does not refit MOFA, recompute MOON scores, modify the archived inputs, or demonstrate end-to-end replay of historical runners. Some runners and RDS metadata retain pre-move paths and require migration before such a replay. `audit_summary.json` and `artifact_manifest.csv` record the audited values and files."),
  file.path(out, "README.md"))
cat("All manuscript frozen-state checks passed; exact 388-edge network reproduced.\n")
