#!/usr/bin/env Rscript
# Complete scores + explicit retention reasons for the compressed depth4 candidates.
source("scripts/support_functions.R")
base <- "results/cosmos_factor1_compressed_review"
state <- readRDS("results/cosmos_factor1_parameter_search/runs/radius6_depth4_compressed/run_state.rds")
raw <- readRDS(file.path(base,"tf_audit/radius6_depth4_compressed_diagnostic_scores.rds"))$final_decompressed
catalog <- utils::read.csv("results/cosmos_mouse_pkn/nodes.csv",check.names=FALSE)
raw$label <- catalog$label[match(raw$source,catalog$node_id)]
raw$upstream_activity <- unname(state$tf[raw$source])
raw$input_tf <- raw$source %in% names(state$tf)
raw$returned_by_moon <- raw$source %in% state$scores$source
raw$sign_agrees <- ifelse(raw$input_tf,sign(raw$score)==sign(raw$upstream_activity),NA)
raw$retention_reason <- ifelse(raw$returned_by_moon,"returned",
  ifelse(raw$input_tf & !raw$sign_agrees,"upstream_sign_conflict","other"))
stopifnot(all(raw$retention_reason %in% c("returned","upstream_sign_conflict")),
  setequal(raw$source[raw$returned_by_moon],state$scores$source),
  max(abs(raw$score[match(state$scores$source,raw$source)]-state$scores$score))<1e-10)
for(id in c("C01","C02","C03")) {
  a <- utils::read.csv(file.path(base,"threshold_search",id,"network_ATT.csv"),check.names=FALSE)
  raw[[paste0("in_",id,"_solution")]] <- raw$source %in% a$source
}
raw <- raw[order(!raw$input_tf,raw$label,raw$source),]
dir.create(file.path(base,"score_retention"),showWarnings=FALSE)
utils::write.csv(raw,file.path(base,"score_retention/radius6_depth4_compressed_scores_all.csv"),row.names=FALSE,na="")
saveRDS(raw,file.path(base,"score_retention/radius6_depth4_compressed_scores_all.rds"))
print(raw[raw$source %in% paste0("Gene__",c("Esr1","Nfic","Dbp")),c("label","score","level","upstream_activity","retention_reason")])
cat(nrow(raw),"complete original-ID scores;",sum(raw$returned_by_moon),"returned;",
  sum(!raw$returned_by_moon),"upstream sign conflicts.\n")

