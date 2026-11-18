#!/usr/bin/env Rscript
# More compression without losing original-ID TF/RNA coherence checks.
source("scripts/support_functions.R")
source("scripts/moon/moon_safe_helpers.R")
args <- commandArgs(TRUE)
depth <- as.integer(if (any(grepl("^--depth=",args))) sub("^--depth=","",grep("^--depth=",args,value=TRUE)[1]) else 8)
stopifnot(depth>=1)
root <- "results/cosmos_factor1_compressed_review/inputs_only_sensitivity"
run_id <- paste0("radius6_depth",depth,"_compressed_inputs_only")
out <- file.path(root,"runs",run_id)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
x <- readRDS("results/cosmos_factor1_pilot/inputs/factor1_inputs.rds")
pr <- readRDS("results/cosmos_factor1_parameter_search/topology_radius_6/prepared_topology.rds")
reg <- x$tf_regulon[,c("source","target","mor")]
keys <- function(e) paste(e$source,e$interaction,e$target,sep="\t")
nodes <- function(e) unique(c(e$source,e$target))
wc <- function(x,file) utils::write.csv(x,file.path(out,file),row.names=FALSE,na="")
original <- moon_check_pkn(pr$net); tf <- pr$tf; lipid <- pr$lipid
ledger <- list(); removed <- list(); converged <- FALSE
for(i in 1:15) {
  used_original <- original
  tf <- tf[names(tf) %in% nodes(original)]
  lipid <- lipid[names(lipid) %in% nodes(original)]
  comp <- moon_safe_compress(original,unique(c(names(tf),names(lipid))))
  score_net <- comp$pkn_scoring
  cat(run_id,"iteration",i,":",nrow(original),"original ->",nrow(score_net),"scoring edges\n")
  scores <- as.data.frame(cosmosR::moon(tf,lipid,score_net,n_layers=depth,statistic="ulm"))
  full <- moon_safe_decompress_scores(scores,comp$node_map)
  # Always filter in the original domain; merged regulon IDs are never missed.
  filtered <- cosmosR::filter_incohrent_TF_target(full,reg,original,x$hepatocyte_gene_weight)
  delta <- original[!keys(original) %in% keys(filtered),,drop=FALSE]
  ledger[[i]] <- data.frame(iteration=i,edges_scored=nrow(score_net),
    original_edges=nrow(original),edges_removed=nrow(delta),scored_nodes=nrow(full),
    merged_nodes=comp$n_merged_nodes,groups=comp$n_groups,inputs_tfs=length(tf),inputs_lipids=length(lipid))
  if(nrow(delta)) {delta$iteration <- i;removed[[length(removed)+1L]] <- delta}
  saveRDS(list(pkn_original=original,compression=comp,scores_compressed=scores),
    file.path(out,paste0("iteration_",i,".rds")))
  original <- filtered
  if(identical(keys(used_original),keys(filtered))) {converged <- TRUE;break}
}
stopifnot(converged)
# Projection excludes original edges canceled by cleanup after compression.
mv <- stats::setNames(comp$node_map$compressed,comp$node_map$original)
projected <- original
projected$source <- unname(mv[projected$source]);projected$target <- unname(mv[projected$target])
final <- original[keys(projected) %in% keys(score_net),,drop=FALSE]
stopifnot(!any(final$source %in% names(lipid)),!anyDuplicated(full$source))
raw <- moon_safe_decompress_scores(as.data.frame(cosmosR::moon(NULL,lipid,score_net,
  n_layers=depth,statistic="ulm")),comp$node_map)
state <- list(run_id=run_id,run_dir=out,compression="compressed",n_layers=depth,
  topology_radius=6L,converged=converged,iterations=i,tf=tf,lipid=lipid,scores=full,
  pkn_final=final,realized_max_level=max(full$level),scored_nodes=nrow(full),
  final_scoring_edges=nrow(score_net),final_original_edges=nrow(final),
  compression_policy="inputs_only; recompress current original graph after each original-ID coherence iteration")
saveRDS(state,file.path(out,"run_state.rds"))
saveRDS(comp,file.path(out,"compression.rds"))
saveRDS(score_net,file.path(out,"pkn_final_scoring.rds"))
saveRDS(final,file.path(out,"pkn_final_uncompressed.rds"))
wc(full,"moon_scores.csv");wc(raw,"moon_scores_before_upstream_sign_filter.csv")
wc(comp$node_map,"compression_mapping.csv");wc(do.call(rbind,ledger),"coherence_ledger.csv")
if(length(removed)) wc(do.call(rbind,removed),"coherence_removed_edges.csv")
focus <- paste0("Gene__",c("Esr1","Nfic","Dbp"))
a <- data.frame(run_id=run_id,tf=focus,input_activity=unname(pr$tf[focus]),
  in_initial_pkn=focus %in% nodes(pr$net),in_final_pkn=focus %in% nodes(final),
  inferred_before_sign_filter=raw$score[match(focus,raw$source)],
  score_level=raw$level[match(focus,raw$source)],retained_score=focus %in% full$source)
a$sign_agrees <- sign(a$input_activity)==sign(a$inferred_before_sign_filter)
wc(a,"focus_tf_audit.csv")
jsonlite::write_json(list(run_id=run_id,converged=converged,depth=depth,radius=6,
  policy=state$compression_policy,iteration_count=i,initial_tfs=length(pr$tf),
  final_tfs=length(tf),initial_lipids=length(pr$lipid),final_lipids=length(lipid),
  initial_edges=nrow(pr$net),final_edges=nrow(final),final_scoring_edges=nrow(score_net),
  final_merged_nodes=comp$n_merged_nodes,final_groups=comp$n_groups,
  md5=as.list(tools::md5sum(c("scripts/moon/run_factor1_inputs_only_compression.R",
    "scripts/moon/moon_safe_helpers.R","results/cosmos_factor1_pilot/inputs/factor1_inputs.rds",
    "results/cosmos_factor1_parameter_search/topology_radius_6/prepared_topology.rds")))),
  file.path(out,"run_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(a)

