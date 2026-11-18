#!/usr/bin/env Rscript
# Audit actual shared pruning/compression; never run MOON scoring or mutate fits.
options(stringsAsFactors = FALSE, digits = 17)
script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script_arg)) stop("Run with Rscript")
script_path <- normalizePath(script_arg)
root <- normalizePath(file.path(dirname(script_path), "../../.."))
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/audit_normalized_rules_pruning.R --output=NEW_EMPTY_DIRECTORY [--analysis=ANALYSIS_DIRECTORY] [--topology=TOPOLOGY_AUDIT_DIRECTORY]\n",
      "Run the Python normalized-rules topology audit first. This script runs pruning/compression only, never MOON scoring.\n")
  quit(save = "no")
}
stopifnot(all(grepl("^--(output|analysis|topology)=.+$", args)))
option <- function(key, default = NULL) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}
output_arg <- option("--output")
if (is.null(output_arg)) stop("An explicit --output=NEW_EMPTY_DIRECTORY is required")
analysis_dir <- normalizePath(option("--analysis", "results/causal_modelling/factor1_lipid_to_tf"))
topology_dir <- normalizePath(option("--topology", file.path(analysis_dir, "normalized_rules_topology")))
source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")
stopifnot(requireNamespace("digest", quietly = TRUE))
analysis_manifest_path <- file.path(analysis_dir, "run_manifest.json")
analysis_manifest <- jsonlite::read_json(analysis_manifest_path, simplifyVector = TRUE)
stopifnot(identical(analysis_manifest$config$analysis, "factor1_lipid_to_tf"))
contract_path <- analysis_manifest$config$dependency_contract
contract <- moon_read_config(contract_path)
out <- moon_new_output(output_arg)
input <- readRDS(contract$inputs)
raw <- moon_check_pkn(readRDS(contract$pkn))
tf <- input$tf_activity[input$tf_activity != 0 & names(input$tf_activity) %in% moon_nodes(raw)]
panels <- list(all_tfs = tf, strong_tfs = tf[abs(tf) > 2])
stopifnot(length(input$lipid_weight) == 155L, length(tf) == 353L, length(panels$strong_tfs) == 49L)
graph_paths <- c(
  prebuilt_curated_ph73 = file.path(analysis_dir, "coverage", "pkn_ph73_lipid_precursor.csv"),
  projected_rules_ph73 = file.path(topology_dir, "pkn_lipid_rules_ph73_precursor.csv"))
graphs <- lapply(graph_paths, function(path) moon_check_pkn(utils::read.csv(path)))
canonical_edges <- function(net) {
  net <- net[order(net$source, net$target, net$interaction), , drop = FALSE]
  paste(paste(net$source, net$interaction, net$target, sep = "\t"), collapse = "\n")
}
edge_hash <- function(net) digest::digest(canonical_edges(net), algo = "sha256", serialize = FALSE)
vector_hash <- function(vec) {
  vec <- vec[order(names(vec))]
  digest::digest(paste(paste(names(vec), sprintf("%.17g", vec), sep = "\t"), collapse = "\n"),
                 algo = "sha256", serialize = FALSE)
}
comparison <- ledgers <- list()
existing_dependency_paths <- character()
for (panel in names(panels)) for (steps in c(6L, 8L, 12L)) {
  message("Pruning/compression equivalence: ", panel, ", steps=", steps)
  prepared <- lapply(graphs, function(pkn) {
    upstream <- input$lipid_weight[names(input$lipid_weight) %in% moon_nodes(pkn)]
    downstream <- panels[[panel]][names(panels[[panel]]) %in% moon_nodes(pkn)]
    moon_prepare(pkn, upstream, downstream, steps)
  })
  base <- prepared$prebuilt_curated_ph73
  rules <- prepared$projected_rules_ph73
  stopifnot(identical(base$net, rules$net), identical(base$upstream, rules$upstream),
            identical(base$downstream, rules$downstream))
  protected <- unique(c(names(base$upstream), names(base$downstream),
                       input$tf_regulon$source, input$tf_regulon$target))
  compressed <- lapply(prepared, function(x) moon_safe_compress(x$net, protected))
  stopifnot(identical(compressed$prebuilt_curated_ph73$pkn_scoring,
                      compressed$projected_rules_ph73$pkn_scoring),
            identical(compressed$prebuilt_curated_ph73$node_map,
                      compressed$projected_rules_ph73$node_map))
  id <- paste0(panel, "_pruning", steps)
  for (name in names(prepared)) {
    x <- prepared[[name]]
    x$ledger$graph <- name
    x$ledger$tf_panel <- panel
    x$ledger$pruning_steps <- steps
    ledgers[[length(ledgers) + 1L]] <- x$ledger
    utils::write.csv(x$net, file.path(out, paste0("pruned_", name, "_", id, ".csv")), row.names = FALSE)
  }
  existing_dirs <- file.path(analysis_dir, "runs",
    paste0("ph73_lipid_precursor_", panel, "_pruning", steps, "_layers", c(4, 6, 8, 12)))
  existing_paths <- file.path(existing_dirs, "prepared_topology.rds")
  existing_paths <- existing_paths[file.exists(existing_paths)]
  existing_dependency_paths <- c(existing_dependency_paths, existing_paths)
  existing_equal <- vapply(existing_paths, function(path) {
    saved <- readRDS(path)
    identical(saved$net, base$net) && identical(saved$upstream, base$upstream) &&
      identical(saved$downstream, base$downstream)
  }, logical(1))
  stopifnot(all(existing_equal))
  comparison[[length(comparison) + 1L]] <- data.frame(
    tf_panel = panel, pruning_steps = steps,
    base_pruned_edges = nrow(base$net), rules_pruned_edges = nrow(rules$net),
    base_pruned_nodes = length(moon_nodes(base$net)), rules_pruned_nodes = length(moon_nodes(rules$net)),
    upstream_lipids = length(base$upstream), downstream_tfs = length(base$downstream),
    exact_signed_edges_identical = identical(base$net, rules$net),
    exact_upstream_vectors_identical = identical(base$upstream, rules$upstream),
    exact_downstream_vectors_identical = identical(base$downstream, rules$downstream),
    base_edge_sha256 = edge_hash(base$net), rules_edge_sha256 = edge_hash(rules$net),
    upstream_sha256 = vector_hash(base$upstream), downstream_sha256 = vector_hash(base$downstream),
    compressed_edges = nrow(compressed$prebuilt_curated_ph73$pkn_scoring),
    compressed_signed_edges_identical = identical(compressed$prebuilt_curated_ph73$pkn_scoring,
                                                  compressed$projected_rules_ph73$pkn_scoring),
    compression_map_identical = identical(compressed$prebuilt_curated_ph73$node_map,
                                         compressed$projected_rules_ph73$node_map),
    compressed_edge_sha256 = edge_hash(compressed$prebuilt_curated_ph73$pkn_scoring),
    matched_existing_prepared_states = length(existing_paths),
    all_existing_prepared_states_identical = all(existing_equal),
    scoring_performed = FALSE)
}
result <- do.call(rbind, comparison)
utils::write.csv(result, file.path(out, "pruning_equivalence.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, ledgers), file.path(out, "pruning_equivalence_ledgers.csv"), row.names = FALSE)
dependencies <- unique(c(script_path, unname(graph_paths), contract$inputs, contract$pkn,
  contract_path, analysis_manifest_path, existing_dependency_paths, "scripts/moon/moon_safe_helpers.R",
  "scripts/moon/lib/workflow.R"))
manifest <- data.frame(path = dependencies, bytes = file.info(dependencies)$size,
  sha256 = vapply(dependencies, function(path) digest::digest(file = path, algo = "sha256"), character(1)))
utils::write.csv(manifest, file.path(out, "pruning_equivalence_dependency_hashes.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "pruning_equivalence_sessionInfo.txt"))
writeLines(c(
  "# Actual MOON preparation and compression equivalence",
  "",
  "Run `Rscript scripts/moon/ancillary/audit_normalized_rules_pruning.R --analysis=ANALYSIS_DIRECTORY --topology=TOPOLOGY_AUDIT_DIRECTORY --output=NEW_EMPTY_DIRECTORY` after `audit_normalized_rules_topology.py`. Occupied output directories are refused.",
  "",
  "This audit calls the same shared moon_prepare and moon_safe_compress functions as the lipid-to-TF runner, under its frozen software/input contract. It performs no MOON scoring and does not overwrite any existing fit.",
  "",
  "For both TF observation panels (all 353 mapped TFs and the 49 with absolute activity above two), and pruning steps 6, 8 and 12, the normalized curated and rule-expanded graphs produce identical canonical signed pruned edges and identical upstream/downstream numeric vectors. The compressed scoring edge sets and complete compression maps are also identical. The pruned states match every available saved normalized-curated prepared state.",
  "",
  "The 149 restored measured-species inputs have no route to a protein or observed TF under the lipid-only precursor policy. The shared alternating forward/backward pruning removes those disconnected additions before scoring. Restoring them therefore cannot change these fits: the identical scoring inputs, protected nodes, compression and shared RNA/regulon evidence imply the same deterministic MOON/coherence computation for each fixed scoring-depth/threshold configuration. This is equivalence of the computation supplied to scoring, not a claim that a second fit was performed.",
  "",
  "pruning_equivalence.csv records exact equality checks, SHA256 hashes, edge/node/input counts and comparisons against saved states. Edge hashes use canonical sorted source/sign/target tab-delimited lines without a trailing newline. Numeric-vector hashes sort by node and use 17 significant digits. pruning_equivalence_ledgers.csv preserves the actual alternating-pruning iterations for both starting graphs. Dependency hashes and sessionInfo record the implementation and environment.",
  "",
  "The full rule-expanded normalized graph remains a derived diagnostic copy, and the 155 accepted measured identifiers are unchanged. No original PKN or fit was modified."
), file.path(out, "pruning_equivalence.md"))
print(result[, c("tf_panel", "pruning_steps", "base_pruned_edges", "upstream_lipids",
                 "downstream_tfs", "compressed_edges", "exact_signed_edges_identical",
                 "compression_map_identical", "matched_existing_prepared_states")])
