#!/usr/bin/env Rscript
# Compressed topology sensitivities; accepted Factor 1 mapping is unchanged.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
args <- commandArgs(TRUE)
get_arg <- function(key, default) {
  x <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(x)) sub(paste0("^", key, "="), "", tail(x, 1)) else default
}
radii <- as.integer(strsplit(get_arg("--radii", "4,5,7,12"), ",", fixed = TRUE)[[1]])
depth <- as.integer(get_arg("--depth", "8"))
prescore <- "--prescore-rna" %in% args
stopifnot(all(is.finite(radii)), all(radii >= 1), is.finite(depth), depth >= 1)
base <- "results/cosmos_factor1_compressed_review"
out <- file.path(base, if (prescore) "prescore_sensitivity" else "topology_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- readRDS("results/cosmos_factor1_pilot/inputs/factor1_inputs.rds")
pkn <- moon_check_pkn(readRDS("results/cosmos_mouse_pkn/pkn_lipid_rules.rds"))
nodes <- function(e) unique(c(e$source, e$target))
keys <- function(e) paste(e$source, e$interaction, e$target, sep = "\t")
tf_all <- x$tf_activity[abs(x$tf_activity) > 2 & names(x$tf_activity) %in% nodes(pkn)]
lipid_all <- x$lipid_weight[x$lipid_weight != 0 & names(x$lipid_weight) %in% nodes(pkn)]
rna <- x$hepatocyte_gene_weight
reg <- x$tf_regulon[, c("source", "target", "mor")]
wc <- function(x, name) utils::write.csv(x, name, row.names = FALSE, na = "")
reach <- function(e, seed, radius, reverse = FALSE) {
  from <- if (reverse) e$target else e$source
  to <- if (reverse) e$source else e$target
  seen <- intersect(seed, nodes(e)); frontier <- seen
  for (i in seq_len(radius)) {
    frontier <- setdiff(unique(to[from %in% frontier]), seen)
    seen <- union(seen, frontier)
    if (!length(frontier)) break
  }
  seen
}
prune <- function(net, radius) {
  tf <- tf_all; lipid <- lipid_all; ledger <- list()
  repeat {
    old <- nrow(net)
    keep <- reach(net, names(tf), radius)
    net <- net[net$source %in% keep & net$target %in% keep, , drop = FALSE]
    keep <- reach(net, names(lipid), radius, TRUE)
    net <- net[net$source %in% keep & net$target %in% keep, , drop = FALSE]
    tf <- tf[names(tf) %in% nodes(net)]; lipid <- lipid[names(lipid) %in% nodes(net)]
    ledger[[length(ledger) + 1L]] <- data.frame(iteration = length(ledger)+1L,
      edges=nrow(net), nodes=length(nodes(net)), tfs=length(tf), lipids=length(lipid))
    if (nrow(net)==old) break
  }
  list(net=moon_check_pkn(net), tf=tf, lipid=lipid, ledger=do.call(rbind,ledger))
}
terminal <- pkn[!pkn$source %in% names(lipid_all), , drop = FALSE]
if (prescore) {
  # Generated from exact PKN edge evidence, accepting only grn provenance.
  grn_path <- file.path(base, "tf_audit/exact_transcription_evidence.csv")
  stopifnot(file.exists(grn_path))
  grn <- utils::read.csv(grn_path, check.names = FALSE)
  r <- unique(reg); names(r)[3] <- "interaction"
  allowed <- intersect(keys(grn), keys(r))
  t <- unname(tf_all[terminal$source])
  w <- unname(rna[terminal$target])
  bad <- keys(terminal) %in% allowed & is.finite(t) & is.finite(w) &
    sign(t * w * terminal$interaction) < 0
  removed <- terminal[bad, , drop = FALSE]
  wc(removed, file.path(out, "prescore_removed_edges.csv"))
  cat("Prescore: removed", nrow(removed), "exact signed GRN/regulon edges using ALL measured TF signs.\n")
  terminal <- terminal[!bad, , drop = FALSE]
}
focus <- paste0("Gene__",c("Esr1","Nfic","Dbp"))
rows <- list()
for (radius in radii) {
  run_id <- paste0("radius",radius,"_depth",depth,"_compressed",
                   if (prescore) "_prescore_rna" else "")
  run_dir <- file.path(out,"runs",run_id)
  dir.create(run_dir, recursive=TRUE, showWarnings=FALSE)
  state_file <- file.path(run_dir,"run_state.rds")
  if (file.exists(state_file) && !"--force" %in% args) {
    state <- readRDS(state_file)
    rows[[length(rows)+1L]] <- utils::read.csv(file.path(run_dir,"focus_tf_audit.csv"))
    cat("Cached:",run_id,"\n")
    next
  }
  pr <- prune(terminal,radius); original <- pr$net
  wc(pr$ledger,file.path(run_dir,"pruning_ledger.csv"))
  comp <- moon_safe_compress(original,unique(c(names(pr$tf),names(pr$lipid),reg$source,reg$target)))
  stopifnot(comp$n_merged_nodes > 0)
  score_net <- comp$pkn_scoring; removed <- list(); ledger <- list(); converged <- FALSE
  for (iteration in 1:10) {
    cat(run_id,"iteration",iteration,"edges",nrow(score_net),"\n")
    used <- score_net
    scores <- as.data.frame(cosmosR::moon(pr$tf,pr$lipid,used,n_layers=depth,statistic="ulm"))
    filtered <- cosmosR::filter_incohrent_TF_target(scores,reg,used,rna)
    delta <- used[!keys(used) %in% keys(filtered),,drop=FALSE]
    if (nrow(delta)) {
      delta$iteration <- iteration; removed[[length(removed)+1L]] <- delta
    }
    ledger[[iteration]] <- data.frame(iteration=iteration,edges_scored=nrow(used),
      edges_removed=nrow(delta),scored_nodes=nrow(scores),max_level=max(scores$level))
    score_net <- filtered
    if (identical(keys(used),keys(filtered))) {converged <- TRUE;break}
  }
  stopifnot(converged)
  mapping <- comp$node_map
  full <- moon_safe_decompress_scores(scores,mapping)
  unfiltered <- as.data.frame(cosmosR::moon(NULL,pr$lipid,score_net,n_layers=depth,statistic="ulm"))
  raw <- moon_safe_decompress_scores(unfiltered,mapping)
  mv <- stats::setNames(mapping$compressed,mapping$original)
  projected <- original
  projected$source <- unname(mv[projected$source]); projected$target <- unname(mv[projected$target])
  final <- original[keys(projected) %in% keys(score_net),,drop=FALSE]
  stopifnot(!anyDuplicated(full$source),all(names(pr$lipid) %in% full$source[full$level==0]),
            !any(final$source %in% names(pr$lipid)))
  state <- list(run_id=run_id,run_dir=run_dir,compression="compressed",n_layers=depth,
    topology_radius=radius,converged=converged,iterations=length(ledger),tf=pr$tf,lipid=pr$lipid,
    scores=full,pkn_final=final,realized_max_level=max(full$level),scored_nodes=nrow(full),
    final_scoring_edges=nrow(score_net),final_original_edges=nrow(final),
    prescore_rna=prescore)
  saveRDS(state,state_file); saveRDS(comp,file.path(run_dir,"compression.rds"))
  saveRDS(original,file.path(run_dir,"pkn_pruned_uncompressed.rds"))
  saveRDS(score_net,file.path(run_dir,"pkn_final_scoring.rds"))
  saveRDS(final,file.path(run_dir,"pkn_final_uncompressed.rds"))
  wc(full,file.path(run_dir,"moon_scores.csv"))
  wc(raw,file.path(run_dir,"moon_scores_before_upstream_sign_filter.csv"))
  wc(mapping,file.path(run_dir,"compression_mapping.csv"))
  wc(do.call(rbind,ledger),file.path(run_dir,"coherence_ledger.csv"))
  if(length(removed)) wc(do.call(rbind,removed),file.path(run_dir,"coherence_removed_edges.csv"))
  tf_audit <- data.frame(run_id=run_id,tf=focus,topology_radius=radius,n_layers=depth,
    input_activity=unname(tf_all[focus]),in_prepared_pkn=focus %in% nodes(original),
    in_final_pkn=focus %in% nodes(final),inferred_before_sign_filter=raw$score[match(focus,raw$source)],
    score_level=raw$level[match(focus,raw$source)],retained_score=focus %in% full$source)
  tf_audit$sign_agrees <- with(tf_audit,sign(input_activity)==sign(inferred_before_sign_filter))
  wc(tf_audit,file.path(run_dir,"focus_tf_audit.csv"));rows[[length(rows)+1L]] <- tf_audit
  jsonlite::write_json(list(run_id=run_id,converged=converged,compression="same_signed_children",
    radius=radius,depth=depth,prescore_rna=prescore,initial_edges=nrow(original),
    compressed_edges=nrow(comp$pkn_scoring),merged_nodes=comp$n_merged_nodes,groups=comp$n_groups,
    final_scoring_edges=nrow(score_net),final_original_edges=nrow(final),
    n_tfs=length(pr$tf),n_lipids=length(pr$lipid),input=normalizePath("results/cosmos_factor1_pilot/inputs/factor1_inputs.rds"),
    input_and_script_md5=as.list(tools::md5sum(c("scripts/moon/run_factor1_compressed_sensitivities.R",
      "scripts/moon/moon_safe_helpers.R","results/cosmos_factor1_pilot/inputs/factor1_inputs.rds",
      "results/cosmos_mouse_pkn/pkn_lipid_rules.rds")))),
    file.path(run_dir,"run_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
  capture.output(sessionInfo(),file=file.path(run_dir,"sessionInfo.txt"))
  print(tf_audit)
}
wc(do.call(rbind,rows),file.path(out,"focus_tf_audit.csv"))
