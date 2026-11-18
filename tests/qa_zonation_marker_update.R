source("scripts/config/analysis_release.R")
# Independent saved-fit and manuscript-number QA for the16 November marker update.
# No pipeline, limma fitting, plot generation or bootstrap is rerun here.
# Run: Rscript tests/qa_zonation_marker_update.R [VERSIONED_ROOT]

qa_zonation_marker_update <- function(root = "results/zonation_marker_update_20261116") {
  if (dir.exists(".local/R-library")) .libPaths(c(".local/R-library", .libPaths()))
  checks <- list(); metrics <- list()
  same <- function(a, b, label, tolerance = 1e-8) {
    a <- as.numeric(a); b <- as.numeric(b)
    stopifnot(length(a) == length(b), identical(is.na(a), is.na(b)))
    valid <- is.finite(a) & is.finite(b)
    stopifnot(all(valid | (is.na(a) & is.na(b)) | (is.infinite(a) & a == b)))
    difference <- if (any(valid)) max(abs(a[valid] - b[valid])) else 0
    if (!is.finite(difference) || difference > tolerance) stop(label, ": ", difference)
    checks[[length(checks) + 1L]] <<- data.frame(check = label, max_abs_difference = difference)
  }
  metric <- function(key, old, new, source, note = "") {
    metrics[[length(metrics) + 1L]] <<- data.frame(metric_key = key, old_value = old,
      new_value = new, status = "INDEPENDENTLY_REPLAYED", source_file = source, note = note)
  }
  read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gls <- function(y, mouse, x, contrast, fit, label) {
    y <- as.matrix(y); x <- as.matrix(x); contrast <- as.matrix(contrast)
    stopifnot(identical(rownames(y), rownames(fit$coefficients)), length(mouse) == ncol(y))
    rho <- attr(x, "rho")
    stopifnot(length(rho) == 1L, abs(rho) < 1)
    v <- matrix(0, ncol(y), ncol(y)); v[outer(mouse, mouse, "==")] <- rho; diag(v) <- 1
    vi <- solve(v); cov <- solve(crossprod(x, vi %*% x))
    beta <- t(cov %*% crossprod(x, vi %*% t(y)))
    res <- y - beta %*% t(x); df <- ncol(y) - qr(x)$rank
    variance <- rowSums((res %*% vi) * res) / df
    effect <- beta %*% contrast
    prior_df <- rep_len(fit$df.prior, nrow(y)); prior <- rep_len(fit$s2.prior, nrow(y))
    post <- (prior_df * prior + df * variance) / (prior_df + df)
    post[is.infinite(prior_df)] <- prior[is.infinite(prior_df)]
    cc <- crossprod(contrast, cov %*% contrast)
    se <- sqrt(post) %o% sqrt(diag(cc)); tt <- effect / se
    pp <- matrix(2 * stats::pt(-abs(as.vector(tt)), rep(rep_len(fit$df.total, nrow(y)), ncol(tt))), nrow(y))
    same(fit$coefficients, effect, paste(label, "all coefficients from independent GLS"))
    same(fit$sigma^2, variance, paste(label, "all residual variances"))
    same(fit$cov.coefficients, cc, paste(label, "full coefficient covariance"))
    same(fit$s2.post, post, paste(label, "posterior variances"))
    same(fit$t, tt, paste(label, "all moderated statistics"))
    same(fit$p.value, pp, paste(label, "all moderated probabilities"))
    list(effect = effect, se = se, t = tt, p = pp, beta = beta, covariance = cov)
  }
  fcg <- readRDS(file.path(root, "fcg05/fcg05_checkpoint.rds"))
  e <- new.env(); load(file.path(root, "fcg05/zonation_limma_results.RData"), envir = e)
  fit <- e$limmaRes[[1L]]; x <- e$limmaRes[[3L]]$design
  attr(x, "rho") <- fcg$corfit$consensus.correlation
  y <- as.matrix(fcg$rna[, fcg$metadata$sample, drop = FALSE])
  f <- gls(y, fcg$metadata$sample_name, x, fcg$contrasts, fit, "FCG")
  for (j in seq_along(fcg$ttop)) {
    tab <- fcg$ttop[[j]]; ii <- match(tab$ID, rownames(y))
    same(tab$logFC, f$effect[ii, j], paste("FCG exported effect", names(fcg$ttop)[j]))
    same(tab$t, f$t[ii, j], paste("FCG exported t", names(fcg$ttop)[j]))
    same(tab$P.Value, f$p[ii, j], paste("FCG exported p", names(fcg$ttop)[j]))
    same(tab$adj.P.Val, stats::p.adjust(f$p[, j], "BH")[ii], paste("FCG per-contrast BH", names(fcg$ttop)[j]))
  }
  mt <- fcg$threshold_metrics
  same(mt$combined_score, mt$entropy_associated * mt$stability * mt$prop_associated * mt$specificity^3,
       "FCG entropy operating-point objective")
  same(fcg$threshold, mt$threshold[which.max(mt$combined_score)], "FCG operating point is objective maximum")
  t <- fcg$t_table; genes <- t[[if ("gene" %in% names(t)) "gene" else "ID"]]
  for (context in c("TPP", "TPC", "OPP", "OPC")) {
    male <- substr(context, 1L, 1L) == "T"; portal <- substr(context, 2L, 3L) == "PP"
    own <- if (male) t$portaltestes_vs_centraltestes else t$portalova_vs_centralova
    other <- if (male) t$portalova_vs_centralova else t$portaltestes_vs_centraltestes
    selected <- genes[abs(t$testes_effect_zonation) > fcg$threshold & abs(other) < 1 &
                        if (portal) own > 0 else own < 0]
    stopifnot(setequal(selected, fcg$candidates[[context]]))
    gonad <- if (male) "Testes" else "Ovaries"; zone <- if (portal) "Periportal" else "Pericentral"
    md <- fcg$metadata[fcg$metadata$gonadal_type == gonad & fcg$metadata$zonation_type == zone, , drop = FALSE]
    stem <- c(TPP="testes_portal",TPC="testes_central",OPP="Ovaries_portal",OPC="ova_central")[[context]]
    net <- read_csv(file.path(root, "fcg05/networks", paste0("cormat_",stem,"_full_SIF.csv")))
    net <- net[net$node2 %in% rownames(fcg$lipids), , drop = FALSE]
    stopifnot(nrow(net) == length(selected) * nrow(fcg$lipids), setequal(unique(net$node1), selected))
    r <- stats::cor(t(as.matrix(fcg$rna[selected, md$sample, drop = FALSE])),
                    t(as.matrix(fcg$lipids[, md$sample_name, drop = FALSE])))
    ii <- cbind(match(net$node1, rownames(r)), match(net$node2, colnames(r)))
    same(net$weight, r[ii], paste("FCG all selected gene-lipid coefficients including non-parenthetical lipids", context))
  }
  old_t <- read_csv(file.path(analysis_path("historical_fcg05"),"t_table_zonation.csv"))
  for (which in c("Fig4B", "Fig5A")) {
    columns <- if (which == "Fig4B") c("portal_vs_Central", "testes_vs_ova") else
      c("portalova_vs_centralova", "portaltestes_vs_centraltestes")
    a <- stats::cor.test(old_t[[columns[1L]]], old_t[[columns[2L]]])
    b <- stats::cor.test(t[[columns[1L]]], t[[columns[2L]]])
    pv <- function(v) if (v < 1e-300) "<1e-300" else format(v, scientific = TRUE, digits = 6)
    if (which == "Fig5A") {
      metric("female_male_zonation_r", sprintf("%.6f",a$estimate), sprintf("%.6f",b$estimate), file.path(root,"fcg05/t_table_zonation.csv"))
      metric("female_male_zonation_n_genes", nrow(old_t), nrow(t), file.path(root,"fcg05/t_table_zonation.csv"))
      metric("female_male_zonation_p", pv(a$p.value), pv(b$p.value), file.path(root,"fcg05/t_table_zonation.csv"))
    } else metric("figure4B_additional_correlation_audit", paste(a$estimate,pv(a$p.value),nrow(old_t)),
      paste("r=",signif(b$estimate,7),";p=",pv(b$p.value),";n=",nrow(t)), file.path(root,"fcg05/t_table_zonation.csv"),
      "Additional audit statistic; no r/p was printed in the embedded old Figure4B.")
  }
  cf <- readRDS(file.path(root,"cf453_models/models.rds")); attr(cf$design,"rho") <- cf$rho
  cg <- gls(cf$modeled_rna, cf$profile_metadata$mouse_id, cf$design, cf$contrasts, cf$fit,"CF453 model cache")
  for (family in unique(cf$rna_results$family)) {
    j <- cf$rna_results$family == family
    same(cf$rna_results$q[j],stats::p.adjust(cf$rna_results$p[j],"BH"),paste("CF453 model BH",family))
  }
  direct <- readRDS(file.path(root,"cf453_direct_diet/diet_analysis.rds"))
  stopifnot(identical(cf$modeled_rna,direct$modeled_rna))
  for (kind in c("primary","conditional")) {
    fit <- direct[[paste0("fit_",kind)]]; contrast <- direct[[paste0(kind,"_contrasts")]]
    gls(cf$modeled_rna,cf$profile_metadata$mouse_id,cf$design,contrast,fit,paste("CF453 direct",kind))
    result <- if(kind=="primary")direct$overview else direct$conditional
    same(result$q,stats::p.adjust(result$p,"BH"),paste("CF453 direct BH",kind))
  }
  rw <- readRDS(file.path(root,"cf453_rewiring/rewiring_analysis.rds"))
  chosen <- list()
  for (diet in c("SD","WD")) for(sex in c("female","male")) for(zone in c("PP","PC")) {
    rr <- cf$rna_results; inter <- rr[rr$contrast==paste0("sex_zone_",diet),]
    own <- rr[rr$contrast==paste("zone",sex,diet,sep="_"),];own<-own[match(inter$gene,own$gene),]
    other <- rr[rr$contrast==paste("zone",if(sex=="female")"male" else "female",diet,sep="_"),];other<-other[match(inter$gene,other$gene),]
    ii <- which(inter$q<=.05 & own$q<=.05 & abs(other$t)<1 & ifelse(own$effect>0,"PP","PC")==zone)
    ii <- head(ii[order(-abs(own$effect[ii]),inter$q[ii],own$q[ii],inter$gene[ii])],5L)
    chosen[[length(chosen)+1L]]<-paste(inter$gene[ii],diet,sex,zone)
  }
  stopifnot(identical(paste(rw$anchors$gene,rw$anchors$anchor_diet,rw$anchors$sex,rw$anchors$zone),unlist(chosen,use.names=FALSE)))
  three <- direct$overview[direct$overview$contrast=="diet_sex_zone",]
  same(rw$anchors$threeway_q,three$q[match(rw$anchors$gene,three$gene)],"CF453 anchors use new direct interaction annotations")
  for(context in unique(rw$edges$context)) {
    ee<-rw$edges[rw$edges$context==context,];s<-ee$sex[1];z<-ee$zone[1];gg<-unique(ee$gene)
    for(d in c("SD","WD")) {
      mice<-cf$metadata$mouse_id[cf$metadata$sex==s & cf$metadata$diet==d]
      pp<-cf$profile_metadata$profile_id[match(paste(mice,z),paste(cf$profile_metadata$mouse_id,cf$profile_metadata$zone))]
      r<-stats::cor(t(cf$modeled_rna[gg,pp,drop=FALSE]),t(cf$lipid_matrix[,mice,drop=FALSE]))
      same(ee[[paste0("r_",d)]],r[cbind(match(ee$gene,rownames(r)),match(ee$lipid,colnames(r)))],paste("CF453 saved programme correlations",context,d))
    }
  }
  val <- read_csv(file.path(root,"fcg07/Mtmr11_male_PP_prespecified_source.csv"))
  kd <- read_csv("results/validation/knockdown_t_table.csv")
  kd_id <- kd[[1L]]; kd_t <- kd[[which(tolower(names(kd))=="mtmr11")]]
  same(val$KD_t,kd_t[match(val$lipid,kd_id)],"Every plotted Mtmr11 experimental KD t remains unchanged")
  summary <- read_csv(file.path(root,"comparisons/Fig5F_prespecified_male_PP_old_vs_new.csv"))
  old_net <- read_csv(file.path(analysis_path("historical_fcg05"),"networks/cormat_testes_portal_full_SIF.csv"))
  old_net <- old_net[old_net$node1=="Mtmr11" & old_net$node2 %in% rownames(fcg$lipids),]
  new_net <- read_csv(file.path(root,"fcg05/networks/cormat_testes_portal_full_SIF.csv"))
  new_net <- new_net[new_net$node1=="Mtmr11" & new_net$node2 %in% rownames(fcg$lipids),]
  for(version in c("historical","revised")) {
    net<-if(version=="historical")old_net else new_net
    ii<-match(kd_id,net$node2);pred<-numeric(length(kd_id));ok<-!is.na(ii);pred[ok]<-net$weight[ii[ok]]
    r<-stats::cor(pred,kd_t);score<-r*sqrt((length(kd_t)-2)/(1-r^2));p<-2*stats::pt(-abs(score),length(kd_t)-2)
    saved<-summary[summary$version==version,]
    same(saved$ULM_score,score,paste(version,"ULM independently reconstructed on all661 experimental targets"))
    same(saved$ULM_p,p,paste(version,"all-background ULM probability"))
    same(saved$n_lipids,sum(ok),paste(version,"matched lipid point count"))
    same(saved$n_KD_background,length(kd_t),paste(version,"complete KD background count"))
    same(saved$lipid_vector_Pearson,stats::cor(pred[ok],kd_t[ok]),paste(version,"matched-only plotted Pearson r"))
  }
  old_kd_hash<-digest::digest("results/validation/knockdown_t_table.csv",file=TRUE,algo="sha256")
  same(nrow(val),345,"Figure5F has345 matched points, distinct from22 reference mice and661 ULM background features")
  out<-file.path(root,"report/qa");dir.create(out,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(do.call(rbind,checks),file.path(out,"independent_fit_checks.csv"),row.names=FALSE)
  utils::write.csv(do.call(rbind,metrics),file.path(out,"independent_numeric_summary.csv"),row.names=FALSE)
  saveRDS(list(checks=checks,KD_source_sha256=old_kd_hash,session=utils::sessionInfo()),file.path(out,"independent_fit_qa.rds"))
  cat("Versioned marker-update independent QA passed:",length(checks),"groups; no resampling or model fitting rerun.\n")
}
if(sys.nframe()==0L){args<-commandArgs(TRUE);qa_zonation_marker_update(if(length(args))args[1] else "results/zonation_marker_update_20261116")}
