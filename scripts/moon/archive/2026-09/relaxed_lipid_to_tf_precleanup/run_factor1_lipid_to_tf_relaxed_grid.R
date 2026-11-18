#!/usr/bin/env Rscript
# Score-first relaxed Factor 1 lipid -> TF hypothesis screen.
# Candidate chemical identities are topology seeds; measurement weights are
# assignment-level annotations checked after scoring, never ULM inputs.
options(stringsAsFactors = FALSE, digits = 17)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
setwd(normalizePath(file.path(dirname(script), "../..")))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R ",
    "[--output=NEW_DIRECTORY] [--config=PATH] [--validate-only]\n",
    "Scores candidate chemicals from TF activities on the direct chemical-to-protein ",
    "plus protein-to-protein overlay. Factor 1 lipid weights are retained only as ",
    "post-score mapping/sign-consistency annotations.\n"
  ))
  quit(save = "no")
}
stopifnot(all(grepl("^--(output|config)=.+$|^--validate-only$", args)))
option <- function(key, default) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}

source("scripts/moon/moon_safe_helpers.R")
source("scripts/moon/lib/workflow.R")

config_path <- option("--config", "scripts/moon/config/factor1_lipid_to_tf_relaxed.json")
cfg <- jsonlite::read_json(config_path, simplifyVector = TRUE)
stopifnot(identical(cfg$analysis, "factor1_lipid_to_tf_relaxed_regulatory"),
          cfg$max_iterations >= 1L, cfg$strong_tf_cutoff > 0,
          cfg$maximum_solution_edges_exclusive > 0L,
          cfg$minimum_upstream_candidates >= 1L,
          cfg$minimum_downstream_tfs >= 2L)

write_json <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE,
                       na = "null", null = "null", digits = NA)
}
protein_node <- function(x) grepl("^(Gene|UniProt)__", x)
truth <- function(x) tolower(as.character(x)) %in% c("true", "1", "yes")
required_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

# Preserve the frozen factor/software contract while using a separate PKN and mapping ledger.
contract <- moon_read_config(cfg$dependency_contract)
stopifnot(identical(cfg$inputs, contract$inputs))
input <- readRDS(cfg$inputs)
overlay <- moon_check_pkn(readRDS(cfg$pkn))
assignments <- utils::read.csv(cfg$mapping, check.names = FALSE)
required_columns(assignments,
  c("assignment_id", "measurement_feature", "candidate_measurement_value", "pkn_source_node",
    "pkn_source_chebi", "default_relaxed_candidate"),
  "Relaxed mapping ledger")
assignments <- assignments[truth(assignments$default_relaxed_candidate), , drop = FALSE]
if (!nrow(assignments)) stop("No default relaxed assignments in mapping ledger")
assignments$candidate_measurement_value <- as.numeric(assignments$candidate_measurement_value)
if (any(!is.finite(assignments$candidate_measurement_value)) ||
    any(assignments$candidate_measurement_value == 0)) {
  stop("Relaxed mapping ledger has non-finite or zero Factor 1 measurements")
}
if (anyDuplicated(assignments$assignment_id)) stop("Assignment IDs must be unique")
if (!all(grepl("^Metab__", assignments$pkn_source_node))) {
  stop("Relaxed candidates must be Metab__ PKN source nodes")
}

# A source may receive multiple feature assignments. Its unit score is calculated once;
# topology seed values are placeholders because moon() does not use them for ULM scoring.
candidate_nodes <- sort(unique(assignments$pkn_source_node))
if (!all(candidate_nodes %in% moon_nodes(overlay))) {
  stop("Mapping ledger contains candidate node(s) absent from directional overlay")
}
topology_seed <- stats::setNames(rep(1, length(candidate_nodes)), candidate_nodes)
moon_vector_check(topology_seed)

# Explicitly omit every reaction/metabolite-to-metabolite transition: this is a direct
# chemical-regulation screen, not a lipid-remodelling or precursor-path analysis.
direct_pkn <- overlay[
  (grepl("^Metab__", overlay$source) & protein_node(overlay$target)) |
    (protein_node(overlay$source) & protein_node(overlay$target)),
  c("source", "target", "interaction"), drop = FALSE
]
direct_pkn <- moon_check_pkn(direct_pkn)
if (!all(candidate_nodes %in% moon_nodes(direct_pkn))) {
  stop("Candidate node(s) have no direct regulatory representation in direct signaling graph")
}

moon_vector_check(input$tf_activity)
tf <- input$tf_activity[input$tf_activity != 0 &
  names(input$tf_activity) %in% moon_nodes(direct_pkn)]
if (!length(tf)) stop("No Factor 1 TF activities are present in the direct signaling PKN")
panels <- list(
  all_tfs = tf,
  strong_tfs = tf[abs(tf) > cfg$strong_tf_cutoff]
)
if (any(vapply(panels, length, integer(1)) == 0L)) stop("A configured TF panel is empty")

# moon() filters upstream values only after calculating scores. This local fit preserves
# scores for all candidate nodes, while moon_prepare still uses the candidate IDs to
# retain the directed candidate-to-TF corridor.
moon_fit_unconstrained_upstream <- function(data, pruning_steps, n_layers,
                                            max_iterations = 10L) {
  prepared <- moon_prepare(data$pkn, data$topology_seed, data$downstream, pruning_steps)
  protected <- unique(c(names(prepared$upstream), names(prepared$downstream),
                         data$regulon$source, data$regulon$target))
  compression <- moon_safe_compress(prepared$net, protected)
  scoring <- compression$pkn_scoring
  ledger <- removed <- list()
  converged <- FALSE
  for (iteration in seq_len(max_iterations)) {
    message("MOON iteration ", iteration, ": ", nrow(scoring), " compressed edges")
    scores <- moon_score(scoring, NULL, prepared$downstream, n_layers)
    filtered <- cosmosR::filter_incohrent_TF_target(scores, data$regulon, scoring, data$rna)
    stopifnot(all(moon_edge_keys(filtered) %in% moon_edge_keys(scoring)))
    delta <- scoring[!moon_edge_keys(scoring) %in% moon_edge_keys(filtered), , drop = FALSE]
    ledger[[iteration]] <- data.frame(
      iteration = iteration, edges_scored = nrow(scoring), edges_removed = nrow(delta),
      scored_nodes = nrow(scores), topology_seed_nodes = length(prepared$upstream),
      max_level = max(scores$level), stringsAsFactors = FALSE
    )
    if (nrow(delta)) {
      delta$iteration <- iteration
      removed[[length(removed) + 1L]] <- delta
    }
    if (identical(moon_edge_keys(scoring), moon_edge_keys(filtered))) {
      converged <- TRUE
      break
    }
    scoring <- filtered
  }
  if (!converged) stop("RNA coherence did not converge in ", max_iterations, " iterations")
  node_map <- stats::setNames(compression$node_map$compressed, compression$node_map$original)
  projected <- prepared$net
  projected$source <- unname(node_map[projected$source])
  projected$target <- unname(node_map[projected$target])
  final <- prepared$net[moon_edge_keys(projected) %in% moon_edge_keys(scoring), , drop = FALSE]
  stopifnot(setequal(
    moon_edge_keys(projected[moon_edge_keys(projected) %in% moon_edge_keys(scoring), , drop = FALSE]),
    moon_edge_keys(scoring)
  ))
  list(
    prepared = prepared, compression = compression, pkn_scoring = scoring,
    pkn_final = final, scores_compressed = scores,
    scores = moon_safe_decompress_scores(scores, compression$node_map),
    coherence_ledger = do.call(rbind, ledger),
    removed_edges = if (length(removed)) do.call(rbind, removed) else data.frame(
      source = character(), target = character(), interaction = numeric(), iteration = integer()),
    pruning_steps = pruning_steps, n_layers = n_layers, converged = converged,
    scoring_mode = "unconstrained_candidate_scores"
  )
}

candidate_status <- function(assignments, fit) {
  score_map <- stats::setNames(fit$scores$score, fit$scores$source)
  level_map <- stats::setNames(fit$scores$level, fit$scores$source)
  retained_nodes <- moon_nodes(fit$prepared$net)
  out <- assignments
  out$candidate_in_prepared_topology <- out$pkn_source_node %in% retained_nodes
  out$candidate_scored <- out$pkn_source_node %in% names(score_map)
  out$moon_score <- unname(score_map[out$pkn_source_node])
  out$moon_level <- unname(level_map[out$pkn_source_node])
  out$assignment_sign_agrees <- out$candidate_scored &
    sign(out$moon_score) == sign(out$candidate_measurement_value)
  out$assignment_status <- ifelse(!out$candidate_in_prepared_topology, "pruned_no_candidate_to_TF_corridor",
    ifelse(!out$candidate_scored, "not_scored_within_realized_layers",
      ifelse(out$assignment_sign_agrees, "scored_sign_concordant", "scored_sign_discordant")))
  out
}

candidate_node_summary <- function(status) {
  grouped <- split(seq_len(nrow(status)), status$pkn_source_node)
  rows <- lapply(names(grouped), function(node) {
    x <- status[grouped[[node]], , drop = FALSE]
    present <- x$candidate_scored
    values <- sort(unique(sign(x$candidate_measurement_value)))
    data.frame(
      pkn_source_node = node,
      pkn_source_chebi = x$pkn_source_chebi[1],
      n_measurement_assignments = nrow(x),
      measurement_features = paste(sort(unique(x$measurement_feature)), collapse = ";"),
      assigned_loading_signs = paste(values, collapse = ";"),
      assignment_signs_are_consistent = length(values) == 1L,
      candidate_in_prepared_topology = all(x$candidate_in_prepared_topology),
      candidate_scored = all(present),
      moon_score = if (all(present)) x$moon_score[1] else NA_real_,
      moon_level = if (all(present)) x$moon_level[1] else NA_real_,
      assignments_sign_concordant = sum(x$assignment_sign_agrees),
      assignments_sign_discordant = sum(present & !x$assignment_sign_agrees),
      all_assignments_sign_concordant = all(present) && all(x$assignment_sign_agrees),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

seed_from_scores <- function(summary, mode) {
  if (mode == "topology") {
    keep <- summary$candidate_scored & summary$moon_level > 0 & summary$moon_score != 0
    return(stats::setNames(sign(summary$moon_score[keep]), summary$pkn_source_node[keep]))
  }
  if (mode == "sign_consensus") {
    keep <- summary$all_assignments_sign_concordant & summary$candidate_scored &
      summary$moon_level > 0 & summary$moon_score != 0
    return(stats::setNames(sign(summary$moon_score[keep]), summary$pkn_source_node[keep]))
  }
  stop("Unknown seed mode: ", mode)
}

selection_is_eligible <- function(solution, cfg) {
  nrow(solution$SIF) > 0L &&
    nrow(solution$SIF) < cfg$maximum_solution_edges_exclusive &&
    solution$audit$retained_upstream >= cfg$minimum_upstream_candidates &&
    solution$audit$retained_observed >= cfg$minimum_downstream_tfs
}

if ("--validate-only" %in% args) {
  cat("Validated ", nrow(assignments), " relaxed feature-to-chemical assignments from ",
      length(unique(assignments$measurement_feature)), " measured lipids; ",
      length(candidate_nodes), " unique candidate chemical nodes; ", nrow(direct_pkn),
      " direct signaling edges; TF panels: ",
      paste(names(panels), vapply(panels, length, integer(1)), sep = "=", collapse = ", "), "\n", sep = "")
  quit(save = "no")
}

out <- moon_new_output(option("--output", cfg$output))
utils::write.csv(assignments, file.path(out, "relaxed_candidate_mapping.csv"), row.names = FALSE)
utils::write.csv(data.frame(pkn_source_node = names(topology_seed), topology_seed_value = unname(topology_seed)),
                 file.path(out, "candidate_topology_seeds.csv"), row.names = FALSE)
utils::write.csv(direct_pkn, file.path(out, "pkn_direct_signaling.csv"), row.names = FALSE)
write_json(list(
  user_direction = "lipids -> TFs",
  factor = "Factor1",
  orientation = input$orientation,
  mapping = "Relaxed feature-to-chemical assignments; all source-node collisions retained as separate rows.",
  upstream = "Candidate chemical node IDs are topology seeds only. Their numeric seed values are not MOON measurements.",
  downstream = "Hepatocyte Factor 1 TF footprint activities are terminal observations.",
  score_construction = "ULM scores are calculated with upstream_input=NULL. Factor 1 lipid loadings are joined afterwards for assignment-level sign concordance.",
  pkn_policy = "Direct chemical-to-protein and protein-to-protein edges from the directional-STITCH exploratory overlay; no reaction traversal.",
  timing = "Cross-sectional signatures; no temporal ordering or direct hepatocyte ligand exposure inferred.",
  candidate_mapping_policy = "Exact identity, explicit acid/base, same-formula/same-charge structural descendant, or compatible-more-specific SwissLipids candidate only; broad/class mappings excluded.",
  limitations = cfg$limitations
), file.path(out, "input_decision.json"))

run_rows <- grid_rows <- status_rows <- node_rows <- list()
grid <- expand.grid(primary = cfg$primary_thresholds, secondary = cfg$secondary_thresholds)
grid <- grid[grid$primary >= grid$secondary, , drop = FALSE]
for (panel_name in names(panels)) {
  for (steps in cfg$pruning_steps) {
    for (layers in cfg$n_layers) {
      run_id <- paste("direct", panel_name, paste0("pruning", steps), paste0("layers", layers), sep = "_")
      data <- list(pkn = direct_pkn, topology_seed = topology_seed, downstream = panels[[panel_name]],
                   regulon = input$tf_regulon[, c("source", "target", "mor")],
                   rna = input$hepatocyte_gene_weight)
      prepared <- tryCatch(moon_prepare(data$pkn, data$topology_seed, data$downstream, steps), error = identity)
      status <- if (inherits(prepared, "error")) "no_directed_candidate_to_TF_topology" else
        if (length(prepared$downstream) < cfg$minimum_downstream_tfs) "insufficient_TF_observations" else "ready"
      row <- data.frame(run_id = run_id, tf_panel = panel_name, pruning_steps = steps,
        n_layers = layers, status = status, topology_seed_candidates = length(topology_seed),
        retained_candidates = if (inherits(prepared, "error")) 0L else length(prepared$upstream),
        retained_tfs = if (inherits(prepared, "error")) 0L else length(prepared$downstream),
        prepared_edges = if (inherits(prepared, "error")) 0L else nrow(prepared$net),
        converged = FALSE, iterations = NA_integer_, scoring_edges = NA_integer_, error = "",
        stringsAsFactors = FALSE)
      fit <- NULL
      if (status == "ready") {
        message("Fitting ", run_id)
        fit <- tryCatch(moon_fit_unconstrained_upstream(data, steps, layers, cfg$max_iterations), error = identity)
        if (inherits(fit, "error")) {
          row$status <- "scoring_failed"
          row$error <- conditionMessage(fit)
        } else {
          row$status <- "scored"
          row$converged <- fit$converged
          row$iterations <- nrow(fit$coherence_ledger)
          row$scoring_edges <- nrow(fit$pkn_scoring)
          run_dir <- file.path(out, "runs", run_id)
          dir.create(run_dir, recursive = TRUE)
          moon_save_fit(fit, run_dir)
          moon_write_csv(fit$scores, file.path(run_dir, "moon_scores_unconstrained.csv"))
          assignment_status <- candidate_status(assignments, fit)
          node_status <- candidate_node_summary(assignment_status)
          moon_write_csv(assignment_status, file.path(run_dir, "candidate_assignment_status.csv"))
          moon_write_csv(node_status, file.path(run_dir, "candidate_node_status.csv"))
          status_rows[[length(status_rows) + 1L]] <- cbind(run_id = run_id, assignment_status)
          node_rows[[length(node_rows) + 1L]] <- cbind(run_id = run_id, node_status)
          for (seed_mode in c("topology", "sign_consensus")) {
            seeds <- seed_from_scores(node_status, seed_mode)
            for (j in seq_len(nrow(grid))) {
              solution <- if (length(seeds)) moon_extract_paths(
                fit$scores, fit$pkn_final, seeds, grid$primary[j], grid$secondary[j], data$rna
              ) else list(SIF = data.frame(), ATT = data.frame(), audit = list(
                retained_upstream = 0L, retained_observed = 0L, final_edges = 0L, final_nodes = 0L,
                all_nodes_on_seed_to_observation_walk = TRUE))
              if (length(seeds)) moon_validate_solution(solution, fit, grid$primary[j], grid$secondary[j])
              grid_rows[[length(grid_rows) + 1L]] <- data.frame(
                run_id = run_id, tf_panel = panel_name, pruning_steps = steps, n_layers = layers,
                seed_mode = seed_mode, primary = grid$primary[j], secondary = grid$secondary[j],
                edges = nrow(solution$SIF), nodes = nrow(solution$ATT),
                upstream_candidates = solution$audit$retained_upstream,
                downstream_tfs = solution$audit$retained_observed,
                eligible = selection_is_eligible(solution, cfg), stringsAsFactors = FALSE)
            }
          }
        }
      }
      run_rows[[length(run_rows) + 1L]] <- row
      moon_write_csv(do.call(rbind, run_rows), file.path(out, "run_summary.csv"))
    }
  }
}

run_summary <- do.call(rbind, run_rows)
thresholds <- if (length(grid_rows)) do.call(rbind, grid_rows) else data.frame(
  run_id = character(), tf_panel = character(), pruning_steps = integer(), n_layers = integer(),
  seed_mode = character(), primary = numeric(), secondary = numeric(), edges = integer(), nodes = integer(),
  upstream_candidates = integer(), downstream_tfs = integer(), eligible = logical())
moon_write_csv(thresholds, file.path(out, "threshold_grid.csv"))
moon_write_csv(thresholds[thresholds$eligible, , drop = FALSE], file.path(out, "eligible_candidates.csv"))
if (length(status_rows)) moon_write_csv(do.call(rbind, status_rows), file.path(out, "candidate_assignment_status_all_runs.csv"))
if (length(node_rows)) moon_write_csv(do.call(rbind, node_rows), file.path(out, "candidate_node_status_all_runs.csv"))

# Materialize one transparent compact sign-consistent representative, if a grid row meets
# the predeclared topology/size criteria. This is a selection convenience, not model selection.
selected <- thresholds[thresholds$eligible & thresholds$seed_mode == "sign_consensus", , drop = FALSE]
if (nrow(selected)) {
  selected <- selected[order(-selected$primary, -selected$secondary, selected$edges, selected$run_id), , drop = FALSE]
  chosen <- selected[1, , drop = FALSE]
  fit <- readRDS(file.path(out, "runs", chosen$run_id, "run_state.rds"))
  assignment_status <- utils::read.csv(file.path(out, "runs", chosen$run_id, "candidate_assignment_status.csv"), check.names = FALSE)
  node_status <- candidate_node_summary(assignment_status)
  seeds <- seed_from_scores(node_status, "sign_consensus")
  solution <- moon_extract_paths(fit$scores, fit$pkn_final, seeds, chosen$primary, chosen$secondary,
                                 input$hepatocyte_gene_weight)
  moon_validate_solution(solution, fit, chosen$primary, chosen$secondary)
  selected_dir <- file.path(out, "selected_sign_consensus")
  dir.create(selected_dir, recursive = TRUE)
  moon_write_csv(solution$SIF, file.path(selected_dir, "network_SIF.csv"))
  att <- solution$ATT
  candidate_annotations <- assignment_status[, intersect(c("assignment_id", "measurement_feature",
    "pkn_source_node", "candidate_measurement_value", "primary_mapping_policy", "assignment_sign_agrees"),
    names(assignment_status)), drop = FALSE]
  names(candidate_annotations)[names(candidate_annotations) == "pkn_source_node"] <- "source"
  att <- merge(att, candidate_annotations, by = "source", all.x = TRUE, sort = FALSE)
  moon_write_csv(att, file.path(selected_dir, "network_ATT.csv"))
  write_json(as.list(chosen), file.path(selected_dir, "selection.json"))
}

sources <- c(config_path, cfg$dependency_contract, cfg$inputs, cfg$pkn, cfg$mapping,
  "scripts/moon/run_factor1_lipid_to_tf_relaxed.R", "scripts/moon/lib/workflow.R",
  "scripts/moon/moon_safe_helpers.R")
moon_write_csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
               file.path(out, "dependency_manifest.csv"))
write_json(list(
  analysis = cfg$analysis, runs = nrow(run_summary), runs_scored = sum(run_summary$status == "scored"),
  status_counts = as.list(table(run_summary$status)), threshold_evaluations = nrow(thresholds),
  eligible_topology_candidates = sum(thresholds$eligible & thresholds$seed_mode == "topology"),
  eligible_sign_consensus_candidates = sum(thresholds$eligible & thresholds$seed_mode == "sign_consensus"),
  interpretation = "Candidate chemicals are scored once from TF observations. Mapping rows are retained as post-score sign-consistency hypotheses and are not independent causal tests."
), file.path(out, "analysis_summary.json"))
write_json(list(created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), config = cfg,
  command_args = args, R_version = R.version.string, packages = contract$packages,
  cosmos_commit = contract$cosmos_commit, completed = TRUE), file.path(out, "run_manifest.json"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Completed relaxed lipid-to-TF score-first workflow: ", sum(run_summary$status == "scored"),
    " scored configurations; ", sum(thresholds$eligible & thresholds$seed_mode == "sign_consensus"),
    " compact sign-consensus extractions. Results: ", out, "\n", sep = "")
