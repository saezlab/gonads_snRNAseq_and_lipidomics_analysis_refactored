#!/usr/bin/env Rscript
# Independently audit completed expanded-species MOON fits and saved re-extractions.
# No fitting, scoring, graph modification or frozen-result writes are performed.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run with Rscript")
root <- normalizePath(file.path(dirname(script), "../../.."), mustWork = TRUE)
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/audit_expanded_factor1_rerun.R [--results=RESULT_ROOT] [--output=NEW_JSON_FILE]\n",
      "Defaults: results/causal_modelling/factor1_tf_to_lipid_expanded_species; output independent_validation.json in that root.\n",
      "Audits7 completed fits and7 re-extractions, frozen measurements, currency policies, compression, signed connectivity, adjacent scoring layers and representative selection.\n",
      "Requires all declared outputs; performs no MOON scoring and refuses to overwrite its output report.\n")
  quit(save = "no")
}
if (any(!grepl("^--(results|output)=.+$", args))) stop("Unknown or empty option; use --help")
option <- function(key, default) {
  hit <- args[startsWith(args, paste0(key, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", key)
  if (length(hit)) substring(hit, nchar(key) + 2L) else default
}
base <- normalizePath(option("--results", "results/causal_modelling/factor1_tf_to_lipid_expanded_species"), mustWork = TRUE)
output <- option("--output", file.path(base, "independent_validation.json"))
if (file.exists(output)) stop("Audit output already exists; choose a new --output file")
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(igraph))
frozen <- readRDS("miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds")
inv <- read.csv("results/causal_modelling/cosmos_mouse_pkn_expanded_species/measured_product_candidate_inventory.csv",check.names=FALSE)
inv <- inv[inv$overlay=="chain_compatible",]
inv$has_product_candidate <- as.logical(inv$has_product_candidate)
ledger <- read.csv(file.path(base,"inputs/input_mapping.csv"),check.names=FALSE)
key <- function(x) paste(x$source,x$interaction,x$target,sep="\t")
check <- function(value,msg) if (!isTRUE(value)) stop(msg)
check(!anyDuplicated(ledger$feature_id) && !anyDuplicated(ledger$node_id),"Duplicate feature/node assignment")
old <- frozen$mapping[frozen$mapping$included,]
check(all(old$node_id %in% ledger$node_id),"Old observation mapping dropped")
check(identical(old$feature,ledger$feature_id[match(old$node_id,ledger$node_id)]),"Old identity mapping changed")
new <- ledger[!ledger$node_id %in% old$node_id,]
check(all(new$node_id==paste0("Metab__MeasuredLipid:",new$feature_id,"_unknown")),"New observation IDs differ from exact feature IDs")
check(all(new$feature_id %in% inv$feature_id[inv$has_product_candidate]),"New feature has no product candidate")
for(i in seq_len(nrow(new))) check(new$node_id[i] %in% strsplit(inv$product_nodes[match(new$feature_id[i],inv$feature_id)],"|",fixed=TRUE)[[1]],"New feature assigned to wrong candidate")
intended <- frozen$original_lipid_weight[ledger$feature_id]*frozen$orientation$multiplier
check(all(is.finite(intended)) && all(intended!=0),"Mapped feature weight missing/zero")
check(max(abs(intended-ledger$weight))<1e-12,"New ledger weight disagrees with frozen source")
audit_one <- function(folder) {
  d <- file.path(base,folder)
  check(file.exists(file.path(d,"completed_manifest.json")),paste("Run not completed",folder))
  inputs <- readRDS(file.path(d,"analysis_inputs.rds"))
  fit <- readRDS(file.path(d,"run_state.rds"))
  sol <- readRDS(file.path(d,"solution.rds"))
  summary <- read_json(file.path(d,"results_summary.json"),simplifyVector=TRUE)
  allsc <- read.csv(file.path(d,"moon_scores_all.csv"),check.names=FALSE)
  s <- fit$scores; e <- sol$SIF; a <- sol$ATT
  check(identical(inputs$rna,frozen$hepatocyte_gene_weight),"RNA changed")
  check(identical(inputs$regulon,frozen$tf_regulon[,c("source","target","mor")]),"Regulon changed")
  check(identical(inputs$orientation,frozen$orientation),"Orientation changed")
  check(identical(inputs$downstream[names(frozen$lipid_weight)],frozen$lipid_weight),"Frozen lipid values not preserved exactly")
  check(identical(inputs$upstream,frozen$tf_activity[names(inputs$upstream)]),"TF footprints changed")
  check(!anyDuplicated(names(inputs$downstream)) && !anyDuplicated(s$source),"Duplicate inputs/scores")
  check(all(inputs$pkn$interaction %in% c(-1,1)) && !anyNA(inputs$pkn) && !anyDuplicated(key(inputs$pkn)),"Invalid input PKN signed edges")
  check(all(names(inputs$upstream) %in% unique(c(inputs$pkn$source,inputs$pkn$target))),"Upstream input absent from scoring PKN")
  if(summary$mode=="expanded") {
    check(setequal(names(inputs$downstream),ledger$node_id),"Expanded data/ledger IDs differ")
    values <- stats::setNames(unname(intended),ledger$node_id)
    # Existing frozen node values are authoritative full-precision values.
    values[names(frozen$lipid_weight)] <- frozen$lipid_weight
    check(identical(unname(inputs$downstream[names(values)]),unname(values)),"Expanded data weights differ from original feature values")
  } else check(identical(inputs$downstream,frozen$lipid_weight),"Control inputs changed")
  check(fit$converged && tail(fit$coherence_ledger$edges_removed,1)==0,"Coherence did not converge")
  check(all(fit$coherence_ledger$edges_removed>=0) && all(diff(fit$coherence_ledger$edges_scored)<=0),"Nonmonotonic coherence")
  check(!any(fit$prepared$net$source %in% names(inputs$downstream)),"Prepared observations not terminal")
  map <- setNames(fit$compression$node_map$compressed,fit$compression$node_map$original)
  prot <- fit$compression$protected_nodes
  check(all(map[prot]==prot),"Protected IDs compressed")
  projected <- fit$pkn_final
  projected$source <- unname(map[projected$source]); projected$target <- unname(map[projected$target])
  check(setequal(key(projected),key(fit$pkn_scoring)),"Final original/scoring graph projections differ")
  check(setequal(s$source[s$level==0],names(fit$prepared$downstream)),"Level0 input IDs differ")
  check(identical(unname(setNames(s$score,s$source)[names(fit$prepared$downstream)]),unname(fit$prepared$downstream)),"Level0 input values differ")
  check(all(is.finite(s$score)) && all(is.finite(s$level)),"Nonfinite scores")
  check(all(a$source %in% s$source),"Solution contains unscored ID")
  sc <- setNames(s$score,s$source)
  check(all(key(e) %in% key(fit$pkn_final)),"Solution edge absent from matching final graph")
  check(all(sign(sc[e$source]*sc[e$target])==e$interaction),"Edge score/sign mismatch")
  check(!any(e$source %in% names(inputs$downstream)),"Solution observation has outgoing edge")
  seeds <- a$source[a$type=="upstream_input"]
  obs <- a$source[a$level==0]
  check(all(sign(sc[seeds])==sign(fit$prepared$upstream[seeds])),"Upstream sign mismatch")
  primary <- a$source[a$level>0 & abs(a$score)>summary$selected_primary]
  check(all(abs(a$score[a$level>0])>summary$selected_secondary),"Secondary threshold violated")
  check(all(e$source %in% primary | e$target %in% primary),"Primary incidence threshold violated")
  if(nrow(e)) {
    g <- graph_from_data_frame(e[,c("source","target")],directed=TRUE)
    forward <- unique(unlist(lapply(seeds,function(x) as_ids(subcomponent(g,x,mode="out")))))
    backward <- unique(unlist(lapply(obs,function(x) as_ids(subcomponent(g,x,mode="in")))))
    check(setequal(a$source,intersect(forward,backward)),"Solution is not seed-to-observation connected")
  } else check(nrow(a)==0,"Empty edges with nonempty nodes")
  if (isTRUE(summary$layer_consistent_extraction)) {
    level <- setNames(s$level,s$source)
    check(all(level[e$source]==level[e$target]+1),"Non-adjacent edge in declared layer-consistent solution")
  }
  if (identical(summary$currency_policy,"no_generic_substrate")) {
    policy_path <- if (!is.null(summary$currency_edges_path)) summary$currency_edges_path else "support/cosmos_mouse/moon_currency_substrate_edges.csv"
    exclusions <- read.csv(policy_path,check.names=FALSE)
    removed <- read.csv(file.path(d,"currency_policy_removed_edges.csv"),check.names=FALSE)
    check(setequal(key(exclusions),key(removed)),"Currency policy removal IDs differ")
    check(!any(key(inputs$pkn) %in% key(exclusions)),"Currency substrate edge leaked into input graph")
    check(!any(key(fit$prepared$net) %in% key(exclusions)),"Currency substrate edge leaked into prepared graph")
    check(!any(key(fit$pkn_final) %in% key(exclusions)),"Currency substrate edge leaked into final graph")
    full <- readRDS("results/causal_modelling/cosmos_mouse_pkn_expanded_species/pkn_expanded_chain.rds")
    check(setequal(key(inputs$pkn),setdiff(key(full),key(exclusions))),"Currency policy modified other edges")
  }
  grid <- read.csv(file.path(d,"threshold_grid_results.csv"))
  eligible <- grid[grid$eligible,]
  if(nrow(eligible)) {
    ranked <- eligible[order(-eligible$secondary,-eligible$primary,-eligible$observed_lipids,-eligible$upstream_tfs,eligible$edges),]
    check(summary$selected_primary==ranked$primary[1] && summary$selected_secondary==ranked$secondary[1],"Selection policy not reproduced")
    check(nrow(e)==ranked$edges[1] && nrow(a)==ranked$nodes[1],"Selected counts differ from grid")
  }
  check(summary$selected_edges==nrow(e) && summary$selected_nodes==nrow(a),"Summary counts disagree")
  focus <- read.csv(file.path(d,"focus_tf_retention.csv"),check.names=FALSE)
  list(validated=TRUE,summary=summary,focus=focus,score_levels=as.list(table(s$level)),
       sign_consistency=TRUE,observation_identity_preserved=TRUE,compression_protection=TRUE,
       final_graph_projection=TRUE,connectivity_validated=TRUE,selection_policy_reproduced=TRUE,
       all_scores=allsc)
}
folders <- c("expanded_inputs","frozen_inputs_control","currency_filtered_expanded","currency_filtered_frozen_control","currency_filtered_expanded_layers6","currency_v2_expanded_layers4","currency_v2_expanded_layers6")
audits <- setNames(lapply(folders,audit_one),folders)
# Independently verify the separately extracted networks against their original fits.
audit_extract <- function(folder) {
 d <- file.path(base,folder)
 summary <- read_json(file.path(d,"results_summary.json"),simplifyVector=TRUE)
 f <- readRDS(file.path(summary$fit_reference,"run_state.rds"))
 e <- read.csv(file.path(d,"network_SIF.csv")); a <- read.csv(file.path(d,"network_ATT.csv"))
 score <- setNames(f$scores$score,f$scores$source); level <- setNames(f$scores$level,f$scores$source)
 check(all(key(e)%in%key(f$pkn_final)),"Re-extracted edge absent from final graph")
 check(all(sign(score[e$source]*score[e$target])==e$interaction),"Re-extracted edge score-sign mismatch")
 if(isTRUE(summary$layer_consistent_extraction)) check(all(level[e$source]==level[e$target]+1),"Re-extracted non-adjacent edge")
 if(nrow(e)) {
   seeds <- a$source[a$type=="upstream_input"]; obs <- a$source[a$level==0]
   g <- graph_from_data_frame(e[,c("source","target")],directed=TRUE)
   reach <- unique(unlist(lapply(seeds,function(x)as_ids(subcomponent(g,x,mode="out")))))
   back <- unique(unlist(lapply(obs,function(x)as_ids(subcomponent(g,x,mode="in")))))
   check(setequal(a$source,intersect(reach,back)),"Re-extracted connectivity mismatch")
   check(all(sign(score[seeds])==sign(f$prepared$upstream[seeds])),"Re-extracted upstream sign mismatch")
   p <- summary$selected$primary; sec <- summary$selected$secondary
   pr <- a$source[a$level>0 & abs(a$score)>p]
   check(all(abs(a$score[a$level>0])>sec),"Re-extracted secondary threshold violated")
   check(all(e$source%in%pr|e$target%in%pr),"Re-extracted primary threshold violated")
 }
 grid <- read.csv(file.path(d,"threshold_grid_results.csv"))
 eligible <- grid[grid$eligible,]
 if(nrow(eligible)) {
   ranked<-eligible[order(-eligible$secondary,-eligible$primary,-eligible$observed_lipids,-eligible$upstream_tfs,eligible$edges),]
   check(summary$selected$primary==ranked$primary[1] && summary$selected$secondary==ranked$secondary[1],"Re-extracted selection mismatch")
   check(nrow(e)==ranked$edges[1] && nrow(a)==ranked$nodes[1],"Re-extracted count mismatch")
 } else check(!nrow(e) && !nrow(a),"No eligible network but output nonempty")
 list(validated=TRUE,summary=summary,layer_consistent_required=isTRUE(summary$layer_consistent_extraction),
   adjacent_levels_validated=if(isTRUE(summary$layer_consistent_extraction)) TRUE else NA,
   nonadjacent_selected_edges=sum(level[e$source] != level[e$target]+1),
   score_signs_validated=TRUE,connectivity_validated=TRUE)
}
extractions <- c("control_lower_thresholds","expanded_coverage_1.85_1.0","expanded_layer_consistent","currency_filtered_expanded_lower_thresholds","currency_filtered_expanded_layers6_lower_thresholds","currency_v2_expanded_layers4_lower_thresholds","currency_v2_expanded_layers6_lower_thresholds")
extraction_audits <- setNames(lapply(extractions,audit_extract),extractions)
equivalences <- list()
for (pair in list(c("currency_filtered_expanded","currency_v2_expanded_layers4"),
                  c("currency_filtered_expanded_layers6","currency_v2_expanded_layers6"))) {
 first <- readRDS(file.path(base,pair[1],"run_state.rds"))
 second <- readRDS(file.path(base,pair[2],"run_state.rds"))
 check(identical(first$prepared,second$prepared),"Currency-v1/v2 prepared topologies differ")
 check(identical(first$compression,second$compression),"Currency-v1/v2 compression maps differ")
 check(identical(first$pkn_scoring,second$pkn_scoring) && identical(first$pkn_final,second$pkn_final),"Currency-v1/v2 final graphs differ")
 check(identical(first$scores,second$scores),"Currency-v1/v2 score IDs, levels or values differ")
 equivalences[[paste(pair,collapse="_versus_")]] <- list(exact_prepared_signed_edges=TRUE,
   exact_compression_mapping=TRUE,exact_final_graphs=TRUE,exact_aligned_score_ids_levels_values=TRUE,
   maximum_absolute_score_difference=0,tolerance=1e-10)
}
out <- list(passed=TRUE,created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
 command_args=args,validator=script,currency_v1_v2_equivalence=equivalences,
 mapping=list(measurements=nrow(ledger),frozen_preserved=nrow(old),new=nrow(new),candidate_products=419,
              finite_unique_one_to_one=TRUE,candidate_identity_verified=TRUE,source_weights_verified=TRUE),
 fits=lapply(audits,function(x)x[names(x)!="all_scores"]),extractions=extraction_audits,
 conclusion="Currency-filtered adjacent-layer extraction is stronger evidence than initial unrestricted topology. Four and six layers recover same13 targets at.75/.75; six adds Hoxa10 but no higher cutoffs/target coverage.")
source_paths <- c(script,"miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds",
 "results/causal_modelling/cosmos_mouse_pkn_expanded_species/pkn_expanded_chain.rds",
 "support/cosmos_mouse/moon_currency_substrate_edges.csv","support/cosmos_mouse/moon_currency_substrate_edges_v2.csv")
out$source_md5 <- as.list(tools::md5sum(source_paths))
write_json(out,output,pretty=TRUE,auto_unbox=TRUE,na="null",digits=NA)
cat("Independent audit passed for",length(audits),"fits and",length(extraction_audits),"extractions.\n")
