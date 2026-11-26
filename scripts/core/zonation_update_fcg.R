source("scripts/config/analysis_release.R")
# Current FCG zonation models, panels and zonation-dependent validation.
zmu_require <- function() {
  if (dir.exists(".local/R-library")) .libPaths(c(".local/R-library", .libPaths()))
  pkgs <- c("limma", "dplyr", "tidyr", "ggplot2", "ggrepel", "cowplot",
    "gridExtra", "patchwork", "decoupleR", "reshape2", "pheatmap", "GSEABase",
    "readr", "readxl", "broom", "digest")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

zmu_validate <- function(input_dir) {
  paths <- file.path(input_dir, c("scaled_RNA_zonation.csv", "meta_data_df_zonation.csv"))
  if (!all(file.exists(paths)))
    stop("Missing FCG zonation input: ", paste(paths[!file.exists(paths)], collapse = ", "))
  x <- utils::read.csv(paths[1], row.names = 1, check.names = FALSE)
  m <- utils::read.csv(paths[2], check.names = FALSE)
  required <- c("sample", "sample_name", "zonation_type", "gonadal_type")
  stopifnot(all(required %in% names(m)), !anyDuplicated(m$sample),
    !anyDuplicated(rownames(x)), setequal(m$sample, colnames(x)),
    all(is.finite(as.matrix(x))), !anyNA(m[, required]),
    all(table(m$sample_name, m$zonation_type) == 1L),
    all(m$zonation_type %in% c("Pericentral", "Periportal")),
    all(m$gonadal_type %in% c("Ovaries", "Testes")))
  cat("FCG input dimensions:", nrow(x), "genes x", ncol(x), "profiles;",
    length(unique(m$sample_name)), "mice\n")
  print(table(m$gonadal_type, m$zonation_type))
  invisible(list(metadata = m, genes = rownames(x)))
}

zmu_runtime05 <- function(input_dir, output) {
  output <- file.path(normalizePath(dirname(output), mustWork = FALSE), basename(output))
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  text <- readLines(analysis_path("engine_fcg05",mustWork=TRUE), warn = FALSE)
  text <- gsub("results/zonation", output, text, fixed = TRUE)
  text <- gsub("results/lipidomic/", "results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/", text, fixed=TRUE)
  for (f in c("scaled_RNA_zonation.csv", "meta_data_df_zonation.csv")) {
    text <- gsub(paste0("data/RNA/", f), file.path(input_dir, f), text, fixed = TRUE)
  }
  # Positional legacy contrast vectors require this explicit condition order.
  at <- grep('^cat\\("Running PCA', text)[1]
  ordering <- c(
    'condition_order <- c("Pericentral_Ovaries", "Periportal_Ovaries", "Pericentral_Testes", "Periportal_Testes")',
    'stopifnot(setequal(meta_data_df$condition, condition_order))',
    'meta_data_df <- meta_data_df[order(match(meta_data_df$condition, condition_order), meta_data_df$sample_name), , drop=FALSE]',
    'scaled_RNA <- scaled_RNA[, meta_data_df$sample, drop=FALSE]',
    'stopifnot(identical(unique(meta_data_df$condition), condition_order))')
  text <- append(text, ordering, after = at - 1L)
  start <- grep('^# 7\\. PATHWAY ENRICHMENT', text)[1]
  end <- grep('^# 8\\. LOAD LIPIDOMICS', text)[1]
  stopifnot(length(start) == 1L, length(end) == 1L, start < end)
  fresh <- c('# 7. FRESH PATHWAY ENRICHMENT (no historical result cache)',
    'gene_sets <- import_gmt("support/m5.all.v2023.2.Mm.symbols.gmt")',
    'gene_sets <- dplyr::rename(gene_sets, target=gene, source=term)',
    'tval_df <- ttop_list_to_t_table(ttop_list)',
    'rownames(tval_df) <- tval_df$ID; tval_df <- tval_df[,setdiff(names(tval_df),"ID"),drop=FALSE]',
    'pathway_scores <- decoupleR::run_ulm(mat=tval_df,network=gene_sets,.source="source",.target="target",minsize=15)',
    'save(pathway_scores,gene_sets,file=file.path(zonation_dir,"pathway_enrichment.RData"))',
    'write.csv(pathway_scores,file.path(zonation_dir,"pathway_enrichment_fresh.csv"),row.names=FALSE)',
    'cat("Fresh pathway results:",nrow(pathway_scores),"rows\\n")', '')
  text <- c(text[seq_len(start - 1L)], fresh, text[end:length(text)])
  # A category absent at a threshold has zero mass. New dplyr rejects the
  # historical zero-length expression; sum() preserves a scalar denominator.
  stopifnot(sum(grepl('prop_associated = 1 - proportion[zonation == "no_association"]',text,fixed=TRUE))==1L)
  text <- gsub('prop_associated = 1 - proportion[zonation == "no_association"]',
    'prop_associated = 1 - sum(proportion[zonation == "no_association"])',text,fixed=TRUE)
  if (any(grepl('"results/zonation(/|")', text))) stop("Unredirected WF05 output path.")
  text
}

zmu_write_runtime <- function(text, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # Parse before saving so a partially valid runtime cannot start analysis.
  parse(text = text)
  writeLines(text, path)
  invisible(path)
}

zmu_run05 <- function(input_dir, output_root, current_panels = FALSE) {
  zmu_require(); zmu_validate(input_dir)
  output <- file.path(output_root, "fcg05")
  if (file.exists(file.path(output, "fcg05_checkpoint.rds")))
    stop("A completed FCG05 destination cannot be overwritten: ", output)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  runtime <- file.path(output_root, "runtime", "05_zonation_analysis_versioned.R")
  zmu_write_runtime(zmu_runtime05(input_dir, output), runtime)
  e <- new.env(parent = globalenv())
  warnings <- list()
  on.exit({
    warning_table <- if (length(warnings)) do.call(rbind, warnings) else
      data.frame(message = character(), call = character())
    utils::write.csv(warning_table, file.path(output, "warnings.csv"), row.names = FALSE)
  }, add = TRUE)
  withCallingHandlers(sys.source(runtime, envir = e), warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- data.frame(
      message = conditionMessage(w), call = paste(deparse(conditionCall(w)), collapse = " "))
    invokeRestart("muffleWarning")
  })
  warning_table <- if (length(warnings)) do.call(rbind, warnings) else
    data.frame(message = character(), call = character())
  utils::write.csv(warning_table, file.path(output, "warnings.csv"), row.names = FALSE)
  result <- list(rna = e$scaled_RNA, metadata = e$meta_data_df, ttop = e$ttop_list,
    t_table = e$t_table, lipids = e$lipidomic, TF_scores = e$TF_scores_df,
    threshold = e$optimal_threshold, threshold_metrics = e$optimal_metrics,
    corfit = e$corfit, contrasts = e$cont.matrix,
    candidates = list(TPP = e$testes_periportal_genes$gene,
      TPC = e$testes_pericentral_genes$gene, OPP = e$ovaries_periportal_genes$gene,
      OPC = e$ovaries_pericentral_genes$gene))
  if (current_panels) zmu_science_finish(output_root, result, include_validation = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
  saveRDS(result, file.path(output, "fcg05_checkpoint.rds"))
  invisible(result)
}

zmu_full_gene_cor <- function(rna, metadata, lipids, gene, context) {
  contexts <- list(TPP=c("Testes","Periportal"), TPC=c("Testes","Pericentral"),
    OPP=c("Ovaries","Periportal"), OPC=c("Ovaries","Pericentral"))
  z <- contexts[[context]]
  m <- metadata[metadata$gonadal_type == z[1] & metadata$zonation_type == z[2], , drop=FALSE]
  m <- m[m$sample_name %in% colnames(lipids), , drop=FALSE]
  stopifnot(gene %in% rownames(rna), nrow(m) >= 3L, !anyDuplicated(m$sample_name))
  r <- stats::cor(t(rna[gene,m$sample,drop=FALSE]), t(lipids[,m$sample_name,drop=FALSE]))
  data.frame(gene=gene, context=context, lipid=colnames(r), r=as.numeric(r),
    n_animals=nrow(m), stringsAsFactors=FALSE)
}

zmu_mean_cl_normal <- function(x,conf.int=.95,mult=NULL,na.rm=TRUE) {
  # Plot-only equivalent of Hmisc::smean.cl.normal / ggplot2 wrapper:
  # https://search.r-project.org/CRAN/refmans/Hmisc/html/smean.sd.html
  if(na.rm)x<-x[!is.na(x)]
  n<-length(x);m<-mean(x)
  if(n<2L)return(data.frame(y=m,ymin=NA_real_,ymax=NA_real_))
  if(is.null(mult))mult<-stats::qt((1+conf.int)/2,n-1L)
  delta<-mult*stats::sd(x)/sqrt(n)
  data.frame(y=m,ymin=m-delta,ymax=m+delta)
}

zmu_zonal_network_paths <- function(checkpoint_path) {
  file.path(dirname(checkpoint_path), "networks", paste0(
    c("cormat_testes_portal", "cormat_testes_central", "cormat_Ovaries_portal", "cormat_ova_central"),
    "_full_SIF.csv"))
}

zmu_run07 <- function(output_root, checkpoint = NULL,
                      checkpoint_path = file.path(output_root, "fcg05/fcg05_checkpoint.rds")) {
  zmu_require()
  networks <- zmu_zonal_network_paths(checkpoint_path)
  reference_files <- networks
  if (is.null(checkpoint)) {
    checkpoint <- readRDS(checkpoint_path)
    reference_files <- c(checkpoint_path, reference_files)
  }
  output <- file.path(output_root,"fcg07")
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  output <- normalizePath(output)
  text <- readLines(analysis_path("engine_validation07",mustWork=TRUE),warn=FALSE)
  start2 <- grep('^# 2\\. LOAD KNOCKDOWN',text)[1]
  start8 <- grep('^# 8\\. VALIDATION:',text)[1]
  start10 <- grep('^# 10\\. SUMMARY',text)[1]
  start11 <- grep('^# 11\\. STRUCTURAL',text)[1]
  end11 <- grep('^cat\\("Output directory: results/validation/structural',text)[1]
  stopifnot(all(is.finite(c(start2,start8,start10,start11,end11))))
  header <- text[seq_len(start2-1L)]
  header <- gsub('results/validation',output,header,fixed=TRUE)
  branch <- c(text[start8:(start10-1L)],text[start11:end11])
  branch <- gsub('results/zonation/networks',
    normalizePath(dirname(networks[1L]), mustWork = TRUE),branch,fixed=TRUE)
  branch <- gsub("results/lipidomic/", "results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/", branch, fixed=TRUE)
  # Only the zonation-dependent validation and structural branch is executed.
  setup <- c(
    't_table <- read.csv(analysis_resolve_recorded_path("results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv"),check.names=FALSE)',
    'rownames(t_table) <- t_table$ID',
    'comparisons <- setNames(as.list(rep(NA,ncol(t_table)-1L)),names(t_table)[-1L])',
    'available_global <- character(0)',
    'mean_cl_normal <- zmu_mean_cl_normal',
    'global_xlim <- c(-1.1,1.1); global_ylim <- c(-25,25)')
  runtime <- file.path(output_root,"runtime","07_zonation_validation_versioned.R")
  zmu_write_runtime(c(header,setup,branch),runtime)
  e <- new.env(parent=globalenv()); warnings <- list()
  on.exit({
    w <- if(length(warnings))do.call(rbind,warnings) else data.frame(message=character(),call=character())
    write.csv(w,file.path(output,"warnings.csv"),row.names=FALSE)
  },add=TRUE)
  withCallingHandlers(sys.source(runtime,envir=e),warning=function(w){
    warnings[[length(warnings)+1L]] <<- data.frame(message=conditionMessage(w),
      call=paste(deparse(conditionCall(w)),collapse=" "))
    invokeRestart("muffleWarning")
  })
  # Primary manuscript context is fixed before considering validation scores.
  kd_path <- analysis_resolve_recorded_path("results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv")
  kd <- utils::read.csv(kd_path, check.names = FALSE)
  fixed <- zmu_full_gene_cor(checkpoint$rna,checkpoint$metadata,checkpoint$lipids,"Mtmr11","TPP")
  fixed <- merge(fixed,kd[,c("ID","mtmr11")],by.x="lipid",by.y="ID")
  names(fixed)[names(fixed)=="mtmr11"] <- "KD_t"
  net <- data.frame(source="Mtmr11_TPP_prespecified",target=fixed$lipid,mor=fixed$r)
  # decoupleR's historical ULM background is every KD lipid; unmatched targets
  # receive zero prediction weights. Do not change this to matched-only ULM.
  mat <- matrix(kd$mtmr11,ncol=1L,dimnames=list(kd$ID,"mtmr11"))
  ulm <- decoupleR::run_ulm(mat,net,.source="source",.target="target",.mor="mor",minsize=2)
  ulm$n_lipids <- nrow(fixed); ulm$n_KD_background <- nrow(kd)
  ulm$n_reference_animals <- unique(fixed$n_animals)
  ulm$context <- "male periportal (prespecified)"
  ulm$selected_new <- "Mtmr11" %in% checkpoint$candidates$TPP
  write.csv(ulm,file.path(output,"Mtmr11_male_PP_prespecified_ULM.csv"),row.names=FALSE)
  write.csv(fixed,file.path(output,"Mtmr11_male_PP_prespecified_source.csv"),row.names=FALSE)
  background<-merge(kd[,c("ID","mtmr11")],fixed[,c("lipid","r")],by.x="ID",by.y="lipid",all.x=TRUE)
  background$matched_in_vivo<-!is.na(background$r)
  background$ULM_weight<-ifelse(background$matched_in_vivo,background$r,0)
  write.csv(background,file.path(output,"Mtmr11_male_PP_ULM_background_source.csv"),row.names=FALSE)
  selected <- c("BMP(22:6_22:6)","PC(18:0_18:1)","PI(36:1)","PS(18:0_18:1)","PI(18:0_18:1)")
  fixed$label <- ifelse(fixed$lipid %in% selected,fixed$lipid,NA)
  p <- ggplot2::ggplot(fixed,ggplot2::aes(r,KD_t))+
    ggplot2::geom_hline(yintercept=0,colour="grey85",linewidth=.3)+
    ggplot2::geom_vline(xintercept=0,colour="grey85",linewidth=.3)+
    ggplot2::geom_smooth(method="lm",formula=y~x,colour="grey35",fill="grey85",linewidth=.6)+
    ggplot2::geom_point(colour="#7ba7ff",size=1.7,alpha=.8)+
    ggrepel::geom_text_repel(ggplot2::aes(label=label),size=2.6,max.overlaps=Inf,
      min.segment.length=0,box.padding=.5,seed=42)+
    ggplot2::coord_cartesian(xlim=c(-1.05,1.05),ylim=c(-25,25))+
    ggplot2::labs(title="Mtmr11: male periportal network validation",
      subtitle=sprintf("%d reference mice; %d matched lipids; native ULM = %.2f, p = %.2g",
        unique(fixed$n_animals),nrow(fixed),ulm$score[1],ulm$p_value[1]),
      x="In-vivo gene-lipid correlation (Pearson r)",
      y="Mtmr11 knockdown lipid effect (t-statistic)",
      caption=sprintf("Prespecified male-periportal context. Native ULM uses all %d knockdown lipids; unmatched network weights are zero.",nrow(kd)))+
    ggplot2::theme_classic(base_size=11,base_family="Helvetica")+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=12),
      plot.subtitle=ggplot2::element_text(size=8.3),plot.caption=ggplot2::element_text(size=7.3,hjust=0))
  ggplot2::ggsave(file.path(output,"Fig5F_Mtmr11_male_PP_prespecified.pdf"),p,
    width=7,height=6,device=grDevices::pdf,family="Helvetica",useDingbats=FALSE)
  result <- list(scores=e$sszb_validation_scores,on_target=e$validation_set_only_sszb,
    structural=e$all_regression_stats,prespecified=ulm,prespecified_source=fixed)
  panel <- zmu_fig5f_panel(output_root)
  utils::write.csv(analysis_provenance(c(reference_files, kd_path)),
    file.path(output, "reference_provenance.csv"), row.names = FALSE)
  utils::write.csv(panel, file.path(output, "panel_manifest.csv"), row.names = FALSE)
  combined <- file.path(output_root, "fcg_panel_manifest.csv")
  same_root_parent <- file.path(output_root, "fcg05/fcg05_checkpoint.rds")
  if (file.exists(combined) && file.exists(same_root_parent) &&
      identical(normalizePath(checkpoint_path), normalizePath(same_root_parent))) {
    zmu_append_panel(combined, panel)
  }
  saveRDS(result,file.path(output,"fcg07_checkpoint.rds"))
  invisible(result)
}

zmu_candidates <- function(t, threshold) {
  rows <- list()
  for (context in c("OPP","OPC","TPP","TPC")) {
    own <- if(substr(context,1,1)=="O")"portalova_vs_centralova" else "portaltestes_vs_centraltestes"
    other <- if(substr(context,1,1)=="O")"portaltestes_vs_centraltestes" else "portalova_vs_centralova"
    selected <- abs(t$testes_effect_zonation)>threshold & abs(t[[other]])<1 &
      if(substr(context,2,3)=="PP")t[[own]]>0 else t[[own]]<0
    rows[[context]] <- data.frame(gene=t$gene[selected],context=context,
      own_t=t[[own]][selected],other_t=t[[other]][selected],
      interaction_t=t$testes_effect_zonation[selected],threshold=threshold)
  }
  do.call(rbind,rows)
}

zmu_old <- function() {
  old <- new.env(); load(file.path(analysis_path("historical_fcg05"),"zonation_ttop_list.RData"),envir=old)
  metrics <- utils::read.csv(file.path(analysis_path("historical_fcg05"),"optimal_threshold_metrics.csv"),check.names=FALSE)
  list(ttop=old$ttop_list,t_table=utils::read.csv(file.path(analysis_path("historical_fcg05"),"t_table_zonation.csv"),check.names=FALSE),
    threshold=metrics$threshold[which.max(metrics$combined_score)])
}

zmu_replay_old_model <- function(output_root) {
  zmu_require()
  suppressPackageStartupMessages({library(limma);library(dplyr);library(reshape2);library(readr)})
  source("scripts/support_functions.R")
  e<-new.env(parent=globalenv())
  e$scaled_RNA<-read.csv(analysis_path("historical_fcg_zonation_csv"),row.names=1,check.names=FALSE)
  e$meta_data_df<-read.csv(analysis_path("historical_fcg_zonation_metadata"),check.names=FALSE)
  e$meta_data_df$condition<-paste(e$meta_data_df$zonation_type,e$meta_data_df$gonadal_type,sep="_")
  ord<-c("Pericentral_Ovaries","Periportal_Ovaries","Pericentral_Testes","Periportal_Testes")
  e$meta_data_df<-e$meta_data_df[order(match(e$meta_data_df$condition,ord),e$meta_data_df$sample_name),]
  e$scaled_RNA<-e$scaled_RNA[,e$meta_data_df$sample,drop=FALSE]
  code<-readLines(analysis_path("engine_fcg05",mustWork=TRUE),warn=FALSE)
  start<-grep('^# 5\\. DIFFERENTIAL',code)[1];end<-grep('^# 6\\. VOLCANO',code)[1]
  eval(parse(text=code[start:(end-1L)]),envir=e)
  old<-zmu_old()
  replay<-do.call(rbind,lapply(names(old$ttop),function(n){
    a<-old$ttop[[n]];b<-e$ttop_list[[n]];b<-b[match(a$ID,b$ID),]
    data.frame(contrast=n,genes=nrow(a),max_abs_delta_t=max(abs(a$t-b$t)),
      max_abs_delta_effect=max(abs(a$logFC-b$logFC)),
      max_abs_delta_p=max(abs(a$P.Value-b$P.Value)),
      max_abs_delta_q=max(abs(a$adj.P.Val-b$adj.P.Val)))
  }))
  dir.create(file.path(output_root,"comparisons"),recursive=TRUE,showWarnings=FALSE)
  write.csv(replay,file.path(output_root,"comparisons","fcg_historical_model_exact_replay.csv"),row.names=FALSE)
  saveRDS(list(ttop=e$ttop_list,contrasts=e$cont.matrix,consensus=e$corfit$consensus.correlation),
    file.path(output_root,"comparisons","fcg_historical_model_replayed.rds"))
  stopifnot(all(replay$max_abs_delta_t<1e-7),all(replay$max_abs_delta_q<1e-7))
  print(replay)
  invisible(replay)
}

zmu_compare <- function(output_root, checkpoint=NULL) {
  if(is.null(checkpoint))checkpoint <- readRDS(file.path(output_root,"fcg05","fcg05_checkpoint.rds"))
  output <- file.path(output_root,"comparisons");dir.create(output,recursive=TRUE,showWarnings=FALSE)
  old <- zmu_old()
  long <- list();summary <- list()
  for(contrast in names(checkpoint$ttop)) {
    a <- old$ttop[[contrast]];b <- checkpoint$ttop[[contrast]]
    vars <- c("ID","logFC","t","P.Value","adj.P.Val")
    a <- a[,vars];b <- b[,vars]
    names(a) <- c("gene","effect_old","t_old","p_old","q_old")
    names(b) <- c("gene","effect_new","t_new","p_new","q_new")
    z <- merge(a,b,by="gene",all=TRUE)
    z$contrast <- contrast;z$delta_t<-z$t_new-z$t_old;z$delta_effect<-z$effect_new-z$effect_old
    z$FDR05_old <- z$q_old<.05;z$FDR05_new <- z$q_new<.05
    z$direction_changed <- sign(z$effect_old)!=sign(z$effect_new)
    long[[contrast]] <- z
    summary[[contrast]] <- data.frame(contrast=contrast,n_common=sum(complete.cases(z[,c("t_old","t_new")])),
      significant_old=sum(z$FDR05_old,na.rm=TRUE),significant_new=sum(z$FDR05_new,na.rm=TRUE),
      significant_gained=sum(!z$FDR05_old & z$FDR05_new,na.rm=TRUE),
      significant_lost=sum(z$FDR05_old & !z$FDR05_new,na.rm=TRUE),
      t_Pearson=stats::cor(z$t_old,z$t_new,use="complete.obs"),
      t_Spearman=stats::cor(z$t_old,z$t_new,use="complete.obs",method="spearman"),
      median_abs_delta_t=median(abs(z$delta_t),na.rm=TRUE),max_abs_delta_t=max(abs(z$delta_t),na.rm=TRUE),
      n_abs_delta_t_below_1e_8=sum(abs(z$delta_t)<1e-8,na.rm=TRUE),
      n_abs_delta_t_below_0_1=sum(abs(z$delta_t)<.1,na.rm=TRUE))
  }
  utils::write.csv(do.call(rbind,long),file.path(output,"fcg_gene_statistics_old_vs_new.csv"),row.names=FALSE)
  utils::write.csv(do.call(rbind,summary),file.path(output,"fcg_contrast_summary_old_vs_new.csv"),row.names=FALSE)
  oldcand <- zmu_candidates(old$t_table,old$threshold)
  newcand <- zmu_candidates(checkpoint$t_table,checkpoint$threshold)
  fixedcand <- zmu_candidates(checkpoint$t_table,old$threshold)
  cc <- list(historical=oldcand,revised_reoptimized=newcand,revised_fixed_historical_threshold=fixedcand)
  counts <- do.call(rbind,lapply(names(cc),function(v){
    z<-cc[[v]];data.frame(version=v,context=c("OPP","OPC","TPP","TPC"),
      n_genes=as.integer(table(factor(z$context,levels=c("OPP","OPC","TPP","TPC")))) ,
      threshold=if(v=="revised_reoptimized")checkpoint$threshold else old$threshold)
  }))
  write.csv(counts,file.path(output,"fcg_selected_gene_counts.csv"),row.names=FALSE)
  write.csv(do.call(rbind,lapply(names(cc),function(v){z<-cc[[v]];z$version<-v;z})),
    file.path(output,"fcg_selected_genes_all_selection_rules.csv"),row.names=FALSE)
  a<-oldcand[,c("gene","context")];a$selected_old<-TRUE
  b<-newcand[,c("gene","context")];b$selected_new<-TRUE
  membership<-merge(a,b,by=c("gene","context"),all=TRUE)
  membership$selected_old[is.na(membership$selected_old)]<-FALSE
  membership$selected_new[is.na(membership$selected_new)]<-FALSE
  membership$status<-ifelse(membership$selected_old & membership$selected_new,"retained",
    ifelse(membership$selected_new,"gained","lost"))
  write.csv(membership,file.path(output,"fcg_selected_gene_membership_changes.csv"),row.names=FALSE)
  # Figure 4B counts use FDR in the zonation-resolved model, not the WF02 fit.
  figure4 <- do.call(rbind,lapply(c("historical","revised"),function(v){
    t<-if(v=="historical")old$ttop else checkpoint$ttop
    s<-t$testes_vs_ova;z<-t$portal_vs_Central
    genes<-s$ID[s$adj.P.Val<.05];z<-z[match(genes,z$ID),]
    data.frame(version=v,sex_biased=nrow(z),dual_bias=sum(z$adj.P.Val<.05),
      periportal=sum(z$adj.P.Val<.05 & z$logFC>0),
      pericentral=sum(z$adj.P.Val<.05 & z$logFC<0))
  }))
  write.csv(figure4,file.path(output,"Fig4B_FDR_counts_old_vs_new.csv"),row.names=FALSE)
  wf02<-new.env();load("results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/ttop_full_list.RData",envir=wf02)
  whole<-wf02$ttop_full_list[["Testes vs Ovaries"]][["Hepatocytes"]]
  claims<-list();claim_genes<-list()
  for(v in c("historical","revised"))for(definition in c("WF05_zonation_model_sex_FDR","WF02_fixed_whole_hepatocyte_sex_FDR")){
    tt<-if(v=="historical")old$ttop else checkpoint$ttop
    sex<-if(definition=="WF05_zonation_model_sex_FDR")tt$testes_vs_ova else whole
    selected<-sex$ID[sex$adj.P.Val<=.05];zone<-tt$portal_vs_Central
    z<-zone[match(selected,zone$ID),]
    claims[[length(claims)+1L]]<-data.frame(version=v,sex_gene_set=definition,
      n_sex_biased=length(selected),n_measured_zone=sum(!is.na(z$ID)),
      n_also_zone_FDR05=sum(z$adj.P.Val<=.05,na.rm=TRUE),
      PC=sum(z$adj.P.Val<=.05 & z$logFC<0,na.rm=TRUE),
      PP=sum(z$adj.P.Val<=.05 & z$logFC>0,na.rm=TRUE),
      proportion_zone=sum(z$adj.P.Val<=.05,na.rm=TRUE)/length(selected))
    claim_genes[[length(claim_genes)+1L]]<-data.frame(version=v,sex_gene_set=definition,
      gene=selected,zone_effect=z$logFC,zone_t=z$t,zone_q=z$adj.P.Val,
      zone_FDR05=z$adj.P.Val<=.05,stringsAsFactors=FALSE)
  }
  write.csv(do.call(rbind,claims),file.path(output,"Fig4B_claim_denominator_definition_audit.csv"),row.names=FALSE)
  write.csv(do.call(rbind,claim_genes),file.path(output,"Fig4B_claim_gene_sets_old_vs_new.csv"),row.names=FALSE)
  stems<-c(TPP="cormat_testes_portal",TPC="cormat_testes_central",OPP="cormat_Ovaries_portal",OPC="cormat_ova_central")
  netcounts<-list();fixedpairs<-list();networkchanges<-list()
  for(context in names(stems)) {
    for(version in c("historical","revised"))for(type in c("displayed","full")) {
      dir<-if(version=="historical")file.path(analysis_path("historical_fcg05"),"networks") else file.path(output_root,"fcg05","networks")
      z<-utils::read.csv(file.path(dir,paste0(stems[[context]],if(type=="full")"_full" else "","_SIF.csv")),check.names=FALSE)
      lipid<-z$node2 %in% rownames(checkpoint$lipids)
      netcounts[[length(netcounts)+1L]]<-data.frame(version=version,context=context,type=type,
        all_edges=nrow(z),gene_lipid_edges=sum(lipid),TF_gene_edges=sum(!lipid),
        genes=length(unique(z$node1[lipid])),lipids=length(unique(z$node2[lipid])))
    }
    oldnet<-read.csv(file.path(analysis_path("historical_fcg05"),"networks",paste0(stems[[context]],"_full_SIF.csv")),check.names=FALSE)
    oldnet<-oldnet[oldnet$node2 %in% rownames(checkpoint$lipids),]
    newnet<-read.csv(file.path(output_root,"fcg05","networks",paste0(stems[[context]],"_full_SIF.csv")),check.names=FALSE)
    newnet<-newnet[newnet$node2 %in% rownames(checkpoint$lipids),]
    colnames(oldnet)<-c("gene","lipid","r_old");colnames(newnet)<-c("gene","lipid","r_new")
    network<-merge(oldnet,newnet,by=c("gene","lipid"),all=TRUE);network$context<-context
    networkchanges[[context]]<-network
    current<-do.call(rbind,lapply(unique(oldnet$gene),function(g)
      zmu_full_gene_cor(checkpoint$rna,checkpoint$metadata,checkpoint$lipids,g,context)))
    names(current)[names(current)=="r"]<-"r_new"
    pairs<-merge(oldnet,current,by=c("gene","lipid"),all.x=TRUE)
    pairs$selected_new<-pairs$gene %in% checkpoint$candidates[[context]]
    pairs$displayed_old<-abs(pairs$r_old)>=.65;pairs$displayed_new_fixed_genes<-abs(pairs$r_new)>=.65
    pairs$delta_r<-pairs$r_new-pairs$r_old
    fixedpairs[[context]]<-pairs
  }
  fixedpairs<-do.call(rbind,fixedpairs)
  write.csv(do.call(rbind,netcounts),file.path(output,"fcg_network_counts_old_vs_new.csv"),row.names=FALSE)
  write.csv(do.call(rbind,networkchanges),file.path(output,"fcg_reselected_network_edges_old_vs_new.csv"),row.names=FALSE)
  write.csv(fixedpairs,file.path(output,"fcg_fixed_historical_pairs_old_vs_new.csv"),row.names=FALSE)
  fixedsummary<-do.call(rbind,lapply(split(fixedpairs,fixedpairs$context),function(z)data.frame(
    context=z$context[1],n_pairs=nrow(z),genes=length(unique(z$gene)),
    coefficient_Pearson=cor(z$r_old,z$r_new,use="complete.obs"),
    median_abs_delta_r=median(abs(z$delta_r),na.rm=TRUE),max_abs_delta_r=max(abs(z$delta_r),na.rm=TRUE),
    displayed_old=sum(z$displayed_old),displayed_new=sum(z$displayed_new_fixed_genes),
    display_gained=sum(!z$displayed_old & z$displayed_new_fixed_genes),
    display_lost=sum(z$displayed_old & !z$displayed_new_fixed_genes))))
  write.csv(fixedsummary,file.path(output,"fcg_fixed_historical_pair_summary.csv"),row.names=FALSE)
  examples<-data.frame(gene=c("Fmo3","Sult3a1","Mtmr11","Kitl"),context=c("OPP","OPC","TPP","TPC"))
  audit<-merge(examples,checkpoint$t_table,by="gene",all.x=TRUE)
  audit$selected_old<-mapply(function(g,c)any(oldcand$gene==g & oldcand$context==c),audit$gene,audit$context)
  audit$selected_new<-mapply(function(g,c)any(newcand$gene==g & newcand$context==c),audit$gene,audit$context)
  audit$selected_new_fixed_threshold<-mapply(function(g,c)any(fixedcand$gene==g & fixedcand$context==c),audit$gene,audit$context)
  write.csv(audit,file.path(output,"fcg_historical_four_examples_audit.csv"),row.names=FALSE)
  key<-merge(fixedpairs,examples,by=c("gene","context"))
  key$top6_old<-FALSE;key$top6_new<-FALSE
  for(g in unique(key$gene)) {
    ii<-which(key$gene==g)
    key$top6_old[ii[order(-abs(key$r_old[ii]),key$lipid[ii])][seq_len(min(6,length(ii)))]]<-TRUE
    key$top6_new[ii[order(-abs(key$r_new[ii]),key$lipid[ii])][seq_len(min(6,length(ii)))]]<-TRUE
  }
  write.csv(key,file.path(output,"fcg_historical_four_examples_all_lipid_correlations.csv"),row.names=FALSE)
  invisible(list(contrasts=do.call(rbind,summary),counts=counts,membership=membership,
    fixed_pairs=fixedpairs,examples=key,example_audit=audit))
}

zmu_class_sources <- function(output_root,checkpoint,versions=c("historical","revised")) {
  stems<-c(TPP="cormat_testes_portal",TPC="cormat_testes_central",OPP="cormat_Ovaries_portal",OPC="cormat_ova_central")
  species_all<-list();class_all<-list()
  le<-new.env();load("results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData",envir=le)
  tt<-le$ttop_full_list[["Testes vs Ovaries"]];if(is.list(tt)&&!is.data.frame(tt))tt<-tt[[1]]
  for(v in versions){
    netdir<-if(v=="historical")file.path(analysis_path("historical_fcg05"),"networks") else file.path(output_root,"fcg05","networks")
    sets<-lapply(stems,function(stem){
      z<-read.csv(file.path(netdir,paste0(stem,"_SIF.csv")),check.names=FALSE)
      unique(z$node2[z$node2 %in% rownames(checkpoint$lipids)])
    })
    s<-data.frame(lipid=rownames(checkpoint$lipids),stringsAsFactors=FALSE)
    s$lipid_class<-sub("\\(.*","",s$lipid)
    for(c in names(sets))s[[c]]<-s$lipid %in% sets[[c]]
    s$n_contexts<-rowSums(s[,names(sets)])
    pp<-s$TPP|s$OPP;pc<-s$TPC|s$OPC
    s$zonation<-ifelse(pp&pc,"both",ifelse(pp,"periportal",ifelse(pc,"pericentral","no_association")))
    s$gonad_bias<-ifelse(tt$logFC[match(s$lipid,tt$ID)]>0,"Testes","Ovaries")
    s$version<-v
    s$unique_context<-ifelse(s$n_contexts==0,"unassigned",ifelse(s$n_contexts>1,"shared",
      names(sets)[max.col(s[,names(sets)],ties.method="first")]))
    species_all[[v]]<-s
    cl<-aggregate(list(count=rep(1L,nrow(s))),s[,c("lipid_class","gonad_bias","zonation")],sum)
    total<-aggregate(count~lipid_class+gonad_bias,cl,sum);names(total)[3]<-"total"
    cl<-merge(cl,total,by=c("lipid_class","gonad_bias"));cl$proportion<-cl$count/cl$total;cl$version<-v
    class_all[[v]]<-cl
  }
  source<-file.path(output_root,"fcg05","panel_source_data");dir.create(source,recursive=TRUE,showWarnings=FALSE)
  s<-do.call(rbind,species_all);cl<-do.call(rbind,class_all)
  write.csv(s,file.path(source,"Fig5C_lipid_species_membership_old_vs_new.csv"),row.names=FALSE)
  write.csv(cl,file.path(source,"Fig5C_faithful_butterfly_proportions_old_vs_new.csv"),row.names=FALSE)
  counts<-do.call(rbind,lapply(split(s,s$version),function(z){
    data.frame(version=z$version[1],lipid_universe=nrow(z),associated=sum(z$n_contexts>0),
      unique_context=sum(z$n_contexts==1),shared_context=sum(z$n_contexts>1),
      TPP=sum(z$TPP),TPC=sum(z$TPC),OPP=sum(z$OPP),OPC=sum(z$OPC))
  }))
  write.csv(counts,file.path(output_root,"comparisons","Fig5C_membership_counts_old_vs_new.csv"),row.names=FALSE)
  invisible(list(species=s,classes=cl,counts=counts))
}

zmu_validation_compare <- function(output_root) {
  output<-file.path(output_root,"comparisons")
  a<-read.csv(analysis_resolve_recorded_path("results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/zonation_network_validation_scores.csv"),check.names=FALSE)
  b<-read.csv(file.path(output_root,"fcg07","zonation_network_validation_scores.csv"),check.names=FALSE)
  a$version<-"historical";b$version<-"revised"
  both<-rbind(a,b);both$on_target<-tolower(both$source)==tolower(both$condition)
  both$q_on_target_BH<-NA_real_
  for(v in unique(both$version)){
    idx<-which(both$version==v & both$on_target)
    both$q_on_target_BH[idx]<-p.adjust(both$p_value[idx],method="BH")
  }
  write.csv(both,file.path(output,"fcg07_all_ULM_scores_old_vs_new.csv"),row.names=FALSE)
  write.csv(both[both$on_target,],file.path(output,"fcg07_on_target_ULM_scores_old_vs_new.csv"),row.names=FALSE)
  kd<-read.csv(analysis_resolve_recorded_path("results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv"),check.names=FALSE)
  old<-read.csv(file.path(analysis_path("historical_fcg05"),"networks/cormat_testes_portal_full_SIF.csv"),check.names=FALSE)
  old<-old[tolower(old$node1)=="mtmr11" & grepl("(",old$node2,fixed=TRUE),]
  names(old)<-c("gene","lipid","r_old")
  z<-read.csv(file.path(output_root,"fcg07","Mtmr11_male_PP_prespecified_source.csv"),check.names=FALSE)
  names(z)[names(z)=="r"]<-"r_new"
  z<-merge(z,old[,c("lipid","r_old")],by="lipid")
  rows<-lapply(c("historical","revised"),function(v){
    r<-z[[if(v=="historical")"r_old" else "r_new"]]
    net<-data.frame(source="Mtmr11_male_PP",target=z$lipid,mor=r)
    mat<-matrix(kd$mtmr11,ncol=1,dimnames=list(kd$ID,"mtmr11"))
    u<-decoupleR::run_ulm(mat,net,.source="source",.target="target",.mor="mor",minsize=2)
    data.frame(version=v,gene="Mtmr11",context="male periportal (prespecified)",
      n_reference_mice=unique(z$n_animals),n_lipids=nrow(z),n_KD_background=nrow(kd),ULM_score=u$score[1],ULM_p=u$p_value[1],
      lipid_vector_Pearson=cor(r,z$KD_t),lipid_vector_Spearman=cor(r,z$KD_t,method="spearman"))
  })
  scores<-do.call(rbind,rows)
  oldscore<-a$score[a$source=="mtmr11" & a$condition=="mtmr11" & a$gonade_zonation=="testes_portal"]
  stopifnot(length(oldscore)==1L,abs(scores$ULM_score[1]-oldscore)<1e-8)
  write.csv(scores,file.path(output,"Fig5F_prespecified_male_PP_old_vs_new.csv"),row.names=FALSE)
  write.csv(z,file.path(output,"Fig5F_matched_lipid_coefficients_old_vs_new.csv"),row.names=FALSE)
  oldstruct<-read.csv(analysis_resolve_recorded_path("results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/structural_consistency/consistency_structural_regression_summary.csv"),check.names=FALSE)
  oldstruct<-oldstruct[tolower(oldstruct$gene) %in% c("mtmr11","sult3a1"),];oldstruct$version<-"historical"
  newstruct<-read.csv(file.path(output_root,"fcg07","structural_consistency","consistency_structural_regression_summary.csv"),check.names=FALSE)
  newstruct$version<-"revised"
  write.csv(rbind(oldstruct,newstruct),file.path(output,"fcg07_structural_consistency_old_vs_new.csv"),row.names=FALSE)
  invisible(scores)
}

zmu_network_plot <- function(z,gene,context,selected) {
  z<-z[is.finite(z$r_new)&abs(z$r_new)>=.65,]
  z<-z[order(-abs(z$r_new),z$lipid),];z<-head(z,6L)
  colours<-c(OPP="#ff6969",OPC="#bf2e2e",TPP="#7ba7ff",TPC="#205cbc")
  labels<-c(OPP="Female periportal",OPC="Female pericentral",TPP="Male periportal",TPC="Male pericentral")
  z$x<-rep(c(-2,2),length.out=nrow(z));z$y<-rep(c(2,1,0),each=2,length.out=nrow(z))
  z$label_y<-4.25-(2-z$y)*.55
  z$lipid_label<-vapply(z$lipid,function(s)if(nchar(s)>30)sub("_","_\n",s,fixed=TRUE) else s,character(1))
  p<-ggplot2::ggplot()+
    ggplot2::geom_segment(data=z,ggplot2::aes(x=0,y=4.55,xend=x,yend=y,
      linewidth=abs(r_new),linetype=r_new>0),colour=colours[[context]],alpha=.7)+
    ggplot2::geom_label(data=z,ggplot2::aes(x=sign(x)*.75,y=label_y,label=sprintf("r = %.2f",r_new)),
      size=2.8,linewidth=0,fill="white",label.padding=grid::unit(.08,"lines"))+
    ggplot2::geom_label(data=z,ggplot2::aes(x=x,y=y,label=lipid_label),size=2.8,
      label.padding=grid::unit(.4,"lines"),colour=colours[[context]],fill="white",linewidth=.5)+
    ggplot2::annotate("label",x=0,y=4.7,label=gene,fontface="bold.italic",size=4,
      fill=colours[[context]],colour="white",linewidth=0,label.padding=grid::unit(.45,"lines"))+
    ggplot2::scale_linewidth_continuous(range=c(.5,1.6),limits=c(.65,1),guide="none")+
    ggplot2::scale_linetype_manual(values=c(`TRUE`="solid",`FALSE`="dashed"),guide="none")+
    ggplot2::coord_cartesian(xlim=c(-4.2,4.2),ylim=c(-.5,5.15),clip="off")+
    ggplot2::labs(title=paste0(gene,": gene-lipid associations"),
      subtitle=paste0(labels[[context]],"; ",if(selected)"retained in revised selected programme" else "historical example; revised selection not met"),
      caption="Up to six strongest associations with |r| >= 0.65. Solid: positive; dashed: negative.\nWhole-liver lipid measurements; associations do not establish causality.")+
    ggplot2::theme_void(base_family="Helvetica",base_size=11)+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=13),
      plot.subtitle=ggplot2::element_text(size=9,margin=ggplot2::margin(b=8)),
      plot.caption=ggplot2::element_text(size=8,hjust=0),plot.margin=ggplot2::margin(12,18,12,18))
  list(plot=p,source=z)
}

zmu_fig5f_panel <- function(output_root) {
  base <- file.path(output_root, "fcg07")
  data.frame(panel_id = "Fig5F",
    generated_file = normalizePath(file.path(base, "Fig5F_Mtmr11_male_PP_prespecified.pdf"), mustWork = TRUE),
    source_table = normalizePath(file.path(base, "Mtmr11_male_PP_prespecified_source.csv"), mustWork = TRUE),
    note = "Prespecified male-periportal reference; unchanged knockdown fit.")
}

zmu_append_panel <- function(path, entry) {
  manifest <- utils::read.csv(path, check.names = FALSE)
  utils::write.csv(rbind(manifest[manifest$panel_id != entry$panel_id, ], entry), path, row.names = FALSE)
}

zmu_main_panels <- function(output_root,checkpoint=NULL,comparison=NULL,include_historical=TRUE,
                            include_validation=TRUE) {
  if(is.null(checkpoint))checkpoint<-readRDS(file.path(output_root,"fcg05","fcg05_checkpoint.rds"))
  if(is.null(comparison)) {
    if(include_historical)comparison<-zmu_compare(output_root,checkpoint) else {
      definitions<-c(Fmo3="OPP",Sult3a1="OPC",Mtmr11="TPP",Kitl="TPC")
      examples<-do.call(rbind,lapply(names(definitions),function(g){
        x<-zmu_full_gene_cor(checkpoint$rna,checkpoint$metadata,checkpoint$lipids,g,definitions[[g]])
        names(x)[names(x)=="r"]<-"r_new";x
      }))
      comparison<-list(examples=examples,counts=data.frame(version="revised_reoptimized",
        context=names(checkpoint$candidates),n_genes=lengths(checkpoint$candidates),threshold=checkpoint$threshold))
    }
  }
  base<-file.path(output_root,"fcg05");source<-file.path(base,"panel_source_data")
  dir.create(source,recursive=TRUE,showWarnings=FALSE)
  panel<-file.path(base,"panels");dir.create(panel,recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(output_root,"comparisons"),recursive=TRUE,showWarnings=FALSE)
  classes<-zmu_class_sources(output_root,checkpoint,
    versions=if(include_historical)c("historical","revised") else "revised")
  write.csv(checkpoint$t_table,file.path(source,"Fig4B_Fig5A_all_gene_t_statistics.csv"),row.names=FALSE)
  write.csv(comparison$counts[comparison$counts$version=="revised_reoptimized",],
    file.path(source,"Fig5B_selected_gene_counts.csv"),row.names=FALSE)
  rows<-list()
  add<-function(id,file,src,note=""){
    rows[[length(rows)+1L]]<<-data.frame(panel_id=id,generated_file=normalizePath(file),
      source_table=normalizePath(src),note=note,stringsAsFactors=FALSE)
  }
  add("Fig4B",file.path(base,"zonation_vs_sex_bias_scatter.pdf"),file.path(source,"Fig4B_Fig5A_all_gene_t_statistics.csv"),
    "Faithful WF05 display uses |t|>3 colour bins; manuscript FDR claim is exported separately.")
  add("Fig5A",file.path(base,"sex_specific_zonation_quadrant.pdf"),file.path(source,"Fig4B_Fig5A_all_gene_t_statistics.csv"))
  add("Fig5B_counts",file.path(base,"sex_specific_zonation_summary_bar.pdf"),file.path(source,"Fig5B_selected_gene_counts.csv"),
    "Count component for manual Figure 5B schematic assembly.")
  add("Fig5C",file.path(base,"lipid_class_zonation_proportions_horizontal.pdf"),
    file.path(source,"Fig5C_faithful_butterfly_proportions_old_vs_new.csv"),
    "Faithful historical display: lipid DA gonad direction, with PP/PC associations pooled across gene-sex contexts.")
  add("FigS5A",file.path(base,"optimal_threshold_detection.pdf"),file.path(base,"optimal_threshold_metrics.csv"))
  for(g in c("Fmo3","Sult3a1","Mtmr11","Kitl")){
    z<-comparison$examples[comparison$examples$gene==g,];context<-unique(z$context)
    selected<-g %in% checkpoint$candidates[[context]]
    result<-zmu_network_plot(z,g,context,selected)
    stem<-paste0("Fig5D_",g,"_",context)
    csv<-file.path(source,paste0(stem,".csv"));write.csv(result$source,csv,row.names=FALSE)
    pdf<-file.path(panel,paste0(stem,".pdf"))
    ggplot2::ggsave(pdf,result$plot,width=7,height=5.2,device=grDevices::pdf,
      family="Helvetica",useDingbats=FALSE)
    add(paste0("Fig5D_",g),pdf,csv,if(selected)"Historical gene retained; partners re-ranked using revised coefficients." else "Historical example audit; revised candidate rule not met.")
  }
  if (include_validation) rows[[length(rows) + 1L]] <- zmu_fig5f_panel(output_root)
  manifest<-do.call(rbind,rows)
  write.csv(manifest,file.path(output_root,"fcg_panel_manifest.csv"),row.names=FALSE)
  invisible(manifest)
}

zmu_science_finish <- function(output_root, checkpoint = readRDS(file.path(output_root, "fcg05/fcg05_checkpoint.rds")),
                               include_validation = TRUE) {
  zmu_main_panels(output_root, checkpoint, include_historical = FALSE, include_validation = include_validation)
  zmu_context_class_panel(output_root, checkpoint)
  zmu_repair_auxiliary_heatmaps(output_root, checkpoint)
  zmu_patch_manifest(output_root, include_validation)
  invisible(TRUE)
}

zmu_patch_manifest <- function(output_root, include_validation = TRUE) {
  m<-data.frame(scope=c("WF05","WF05","WF05","WF05","WF05","WF07","WF07","WF07"),
    adaptation=c("Versioned output paths","Versioned input paths","Explicit condition order",
      "Fresh pathway ULM","Zero-frequency category scalar","Reuse knockdown fit","Zonation-only branch","Prespecified Mtmr11 context"),
    detail=c("All results/zonation paths redirected into fcg05, preserving original outputs.",
      "FCG zonation uses the supplied scaled pseudobulk profiles through the compatibility CSV pair.",
      "Pericentral_Ovaries, Periportal_Ovaries, Pericentral_Testes, Periportal_Testes; exact old-input model replay verified.",
      "Restore documented minsize15 ULM on current t statistics; do not load the historical pathway cache.",
      "Absent no_association has zero mass: sum(proportion[category]) avoids zero-length summarise without changing denominators or entropy formula.",
      "Read historical knockdown_t_table.csv; no knockdown model, raw data or global network fit is changed.",
      "Execute existing sections8/9/11, with available_global empty so structural outputs cover zonation genes only.",
      "Add separate male-periportal Mtmr11 panel/score regardless of adaptive best-context selection."),stringsAsFactors=FALSE)
  m<-rbind(m,data.frame(scope=c("Comparison","Presentation","WF07 rendering"),
    adaptation=c("Measured target identities","Expanded selected-gene heatmap layout","Optional Hmisc plotting fallback"),
    detail=c("Comparison and membership ledgers identify features from the measured1254-target universe, retaining three names without parentheses. None crosses0.65 or overlapsKD, so faithful displayed/validation statistics are unchanged.",
      "Auxiliary selected-gene heatmaps show the exact173-gene t-statistic matrices as labelled rows on taller pages; no model or selection changes. The four-context class alternative displays the same20classes as the original class panel; all raw classes remain in its source table.",
      "Use base R mean plus/minus t-quantile times standard error for the legacy stat_summary layer, matching Hmisc smean.cl.normal documented formula and stats::t.test confidence limits. No model, ULM, network or primary Figure5F changes.")))
  files <- c("scripts/workflows/05_zonation_analysis.R", "scripts/core/zonation_update_fcg.R",
    "scripts/core/zonation_inputs.R", "scripts/reports/run_zonation_update_fcg.R",
    file.path(output_root, "runtime/05_zonation_analysis_versioned.R"))
  if (include_validation) {
    files <- c(files, "scripts/workflows/07_validation_experiments.R",
      "results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv",
      "results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_ttop_list.RData",
      file.path(output_root, "runtime/07_zonation_validation_versioned.R"))
  } else {
    m <- m[!startsWith(m$scope, "WF07"), ]
  }
  write.csv(m,file.path(output_root,"fcg_runtime_adapter_manifest.csv"),row.names=FALSE)
  hash<-data.frame(path=normalizePath(files),sha256=vapply(files,digest::digest,character(1),file=TRUE,algo="sha256"))
  write.csv(hash,file.path(output_root,"fcg_source_fit_sha256.csv"),row.names=FALSE)
  invisible(m)
}

zmu_replay_old_threshold <- function(output_root) {
  zmu_require()
  suppressPackageStartupMessages({library(dplyr);library(tidyr)})
  e<-new.env(parent=globalenv())
  e$scaled_RNA<-read.csv(analysis_path("historical_fcg_zonation_csv"),row.names=1,check.names=FALSE)
  e$t_table<-zmu_old()$t_table
  e$lipidomic<-readRDS(file.path(output_root,"fcg05","fcg05_checkpoint.rds"))$lipids
  code<-readLines(analysis_path("engine_fcg05",mustWork=TRUE),warn=FALSE)
  start<-grep('^# 11\\. OPTIMAL',code)[1];end<-grep('^# 12\\. IDENTIFY',code)[1]
  code<-code[start:(end-1L)]
  code<-gsub('prop_associated = 1 - proportion[zonation == "no_association"]',
    'prop_associated = 1 - sum(proportion[zonation == "no_association"])',code,fixed=TRUE)
  warnings<-list()
  withCallingHandlers(eval(parse(text=code),envir=e),warning=function(w){
    warnings[[length(warnings)+1L]]<<-data.frame(message=conditionMessage(w),
      call=paste(deparse(conditionCall(w)),collapse=" "))
    invokeRestart("muffleWarning")
  })
  w<-if(length(warnings))do.call(rbind,warnings) else data.frame(message=character(),call=character())
  write.csv(w,file.path(output_root,"comparisons","fcg_historical_threshold_replay_warnings.csv"),row.names=FALSE)
  historical<-read.csv(file.path(analysis_path("historical_fcg05"),"optimal_threshold_metrics.csv"),check.names=FALSE)
  historical$version<-"historical_archived"
  replay<-e$optimal_metrics;replay$version<-"historical_replayed_zero_category_fix"
  current<-read.csv(file.path(output_root,"fcg05","optimal_threshold_metrics.csv"),check.names=FALSE)
  current$version<-"revised"
  write.csv(rbind(historical,replay,current),
    file.path(output_root,"comparisons","fcg_threshold_objective_old_replay_new.csv"),row.names=FALSE)
  keys<-c("historical_archived","historical_replayed_zero_category_fix","revised")
  thresholds<-c(zmu_old()$threshold,e$optimal_threshold,current$threshold[which.max(current$combined_score)])
  write.csv(data.frame(version=keys,threshold=thresholds),
    file.path(output_root,"comparisons","fcg_threshold_replay_summary.csv"),row.names=FALSE)
  stopifnot(abs(e$optimal_threshold-zmu_old()$threshold)<1e-12)
  cat("Historical threshold replay passed:",e$optimal_threshold,"\n")
  invisible(thresholds)
}

zmu_fixed_threshold_membership <- function(output_root,checkpoint) {
  old<-zmu_old();cand<-zmu_candidates(checkpoint$t_table,old$threshold)
  stems<-c(TPP="cormat_testes_portal",TPC="cormat_testes_central",OPP="cormat_Ovaries_portal",OPC="cormat_ova_central")
  s<-data.frame(lipid=rownames(checkpoint$lipids),stringsAsFactors=FALSE)
  s$lipid_class<-sub("\\(.*","",s$lipid)
  for(context in names(stems)){
    net<-read.csv(file.path(output_root,"fcg05","networks",paste0(stems[[context]],"_full_SIF.csv")),check.names=FALSE)
    net<-net[net$node1 %in% cand$gene[cand$context==context] & net$node2 %in% rownames(checkpoint$lipids) & abs(net$weight)>=.65,]
    s[[context]]<-s$lipid %in% unique(net$node2)
  }
  s$n_contexts<-rowSums(s[,names(stems)]);pp<-s$TPP|s$OPP;pc<-s$TPC|s$OPC
  s$zonation<-ifelse(pp&pc,"both",ifelse(pp,"periportal",ifelse(pc,"pericentral","no_association")))
  le<-new.env();load("results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData",envir=le)
  tt<-le$ttop_full_list[["Testes vs Ovaries"]];if(is.list(tt)&&!is.data.frame(tt))tt<-tt[[1]]
  s$gonad_bias<-ifelse(tt$logFC[match(s$lipid,tt$ID)]>0,"Testes","Ovaries")
  s$version<-"revised_fixed_historical_threshold"
  s$unique_context<-ifelse(s$n_contexts==0,"unassigned",ifelse(s$n_contexts>1,"shared",
    names(stems)[max.col(s[,names(stems)],ties.method="first")]))
  source<-file.path(output_root,"fcg05","panel_source_data")
  write.csv(s,file.path(source,"Fig5C_fixed_historical_threshold_membership.csv"),row.names=FALSE)
  cl<-aggregate(list(count=rep(1L,nrow(s))),s[,c("lipid_class","gonad_bias","zonation")],sum)
  total<-aggregate(count~lipid_class+gonad_bias,cl,sum);names(total)[3]<-"total"
  cl<-merge(cl,total,by=c("lipid_class","gonad_bias"));cl$proportion<-cl$count/cl$total;cl$version<-s$version[1]
  write.csv(cl,file.path(source,"Fig5C_fixed_historical_threshold_class_proportions.csv"),row.names=FALSE)
  count<-data.frame(version=s$version[1],lipid_universe=nrow(s),associated=sum(s$n_contexts>0),
    unique_context=sum(s$n_contexts==1),shared_context=sum(s$n_contexts>1),
    TPP=sum(s$TPP),TPC=sum(s$TPC),OPP=sum(s$OPP),OPC=sum(s$OPC))
  counts<-read.csv(file.path(output_root,"comparisons","Fig5C_membership_counts_old_vs_new.csv"),check.names=FALSE)
  write.csv(rbind(counts[counts$version!=count$version,],count),
    file.path(output_root,"comparisons","Fig5C_membership_counts_all_thresholds.csv"),row.names=FALSE)
  invisible(s)
}

zmu_context_class_panel <- function(output_root, checkpoint) {
  source<-file.path(output_root,"fcg05","panel_source_data")
  s<-read.csv(file.path(source,"Fig5C_lipid_species_membership_old_vs_new.csv"),check.names=FALSE)
  s<-s[s$version=="revised",]
  contexts<-c("OPP","OPC","TPP","TPC")
  z<-do.call(rbind,lapply(contexts,function(context){
    x<-aggregate(list(associated=s[[context]],total=rep(1L,nrow(s))),list(lipid_class=s$lipid_class),sum)
    x$context<-context;x$proportion<-x$associated/x$total;x
  }))
  labels<-c(OPP="Female\nperiportal",OPC="Female\npericentral",TPP="Male\nperiportal",TPC="Male\npericentral")
  z$context<-factor(z$context,levels=contexts)
  enrichment<-new.env();load("results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_scores.RData",envir=enrichment)
  displayed<-unique(enrichment$lipid_class_results$source[
    grepl("^Testes vs Ovaries",enrichment$lipid_class_results$condition)])
  z$displayed_class<-z$lipid_class %in% displayed
  classes<-sort(unique(z$lipid_class));z$lipid_class<-factor(z$lipid_class,levels=rev(classes))
  totals<-setNames(z$total[match(classes,z$lipid_class)],classes)
  threshold <- checkpoint$threshold
  p<-ggplot2::ggplot(z[z$displayed_class,],ggplot2::aes(context,lipid_class,fill=proportion))+
    ggplot2::geom_tile(colour="white",linewidth=.5)+
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.0f%%",100*proportion),
      colour=proportion>.55),size=2.7)+
    ggplot2::scale_colour_manual(values=c(`TRUE`="white",`FALSE`="#222222"),guide="none")+
    ggplot2::scale_fill_gradient(low="#F5F7FA",high="#315B6B",limits=c(0,1),
      labels=scales::percent,name="Class species\nassociated",guide=ggplot2::guide_colourbar(display="rectangles"))+
    ggplot2::scale_x_discrete(labels=labels,position="top")+
    ggplot2::scale_y_discrete(labels=function(x)paste0(x," (n=",totals[x],")"))+
    ggplot2::labs(title="Lipid associations with sex-specific zonation programmes",
      subtitle=paste0("Twenty original lipid classes; revised selection at |interaction t| > ",threshold),
      x=NULL,y=NULL,
      caption="Denominator: all measured species in each lipid class. |r| >= 0.65; species can enter multiple contexts.\nWhole-liver lipid measurements: these associations do not establish spatial lipid localization.")+
    ggplot2::theme_minimal(base_family="Helvetica",base_size=10)+
    ggplot2::theme(panel.grid=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold",size=12),
      plot.subtitle=ggplot2::element_text(size=9),plot.caption=ggplot2::element_text(size=8,hjust=0),
      axis.text.x=ggplot2::element_text(colour="#222222"),legend.position="right")
  csv<-file.path(source,"Fig5C_explicit_gene_context_class_fractions.csv")
  write.csv(z,csv,row.names=FALSE)
  pdf<-file.path(output_root,"fcg05","panels","Fig5C_explicit_gene_context_class_associations.pdf")
  ggplot2::ggsave(pdf,p,width=7.5,height=7.5,device=grDevices::pdf,family="Helvetica",useDingbats=FALSE)
  row<-data.frame(panel_id="Fig5C_gene_context_alternative",generated_file=normalizePath(pdf),source_table=normalizePath(csv),
    note="Unambiguous gene-sex/zone class association fractions; explicit alternative to historical lipid-DA-direction butterfly.")
  zmu_append_panel(file.path(output_root, "fcg_panel_manifest.csv"), row)
  invisible(z)
}

zmu_metric_ledger <- function(output_root) {
  dir<-file.path(output_root,"comparisons");rows<-list()
  add<-function(panel,metric,old,new,fixed=NA_real_,unit="count",note=""){
    rows[[length(rows)+1L]]<<-data.frame(panel_id=panel,metric=metric,historical=old,
      revised_reoptimized=new,revised_fixed_5_4=fixed,unit=unit,note=note,stringsAsFactors=FALSE)
  }
  cs<-read.csv(file.path(dir,"fcg_selected_gene_counts.csv"))
  for(context in unique(cs$context)){
    z<-cs[cs$context==context,];v<-setNames(z$n_genes,z$version)
    add("Fig5A",paste0("selected_genes_",context),v["historical"],v["revised_reoptimized"],v["revised_fixed_historical_threshold"])
  }
  add("Fig5A","selected_genes_total",sum(cs$n_genes[cs$version=="historical"]),
    sum(cs$n_genes[cs$version=="revised_reoptimized"]),sum(cs$n_genes[cs$version=="revised_fixed_historical_threshold"]))
  add("FigS5A","interaction_threshold",5.4,4.4,5.4,"absolute t")
  claim<-read.csv(file.path(dir,"Fig4B_claim_denominator_definition_audit.csv"))
  for(definition in unique(claim$sex_gene_set))for(metric in c("n_sex_biased","n_also_zone_FDR05","PC","PP","proportion_zone")){
    z<-claim[claim$sex_gene_set==definition,];v<-setNames(z[[metric]],z$version)
    add("Fig4B",paste(definition,metric,sep="__"),v["historical"],v["revised"],unit=if(metric=="proportion_zone")"fraction" else "count")
  }
  cl<-read.csv(file.path(dir,"Fig5C_membership_counts_all_thresholds.csv"))
  for(metric in setdiff(names(cl),"version")){
    v<-setNames(cl[[metric]],cl$version)
    add("Fig5C",metric,v["historical"],v["revised"],v["revised_fixed_historical_threshold"])
  }
  nc<-read.csv(file.path(dir,"fcg_network_counts_old_vs_new.csv"))
  for(context in unique(nc$context))for(metric in c("gene_lipid_edges","genes","lipids")){
    z<-nc[nc$context==context & nc$type=="displayed",];v<-setNames(z[[metric]],z$version)
    add("Fig5C",paste0("displayed_",metric,"_",context),v["historical"],v["revised"])
  }
  v<-read.csv(file.path(dir,"Fig5F_prespecified_male_PP_old_vs_new.csv"))
  for(metric in c("ULM_score","ULM_p","n_reference_mice","n_lipids","n_KD_background","lipid_vector_Pearson")){
    x<-setNames(v[[metric]],v$version)
    add("Fig5F",metric,x["historical"],x["revised"],unit=if(grepl("^n_",metric))"count" else "statistic",
      note="Native ULM background661; scatter and matched-vector r use345 lipids.")
  }
  old<-zmu_old()$t_table;new<-readRDS(file.path(output_root,"fcg05","fcg05_checkpoint.rds"))$t_table
  add("Fig5A","female_male_zonation_t_Pearson",cor(old$portalova_vs_centralova,old$portaltestes_vs_centraltestes),
    cor(new$portalova_vs_centralova,new$portaltestes_vs_centraltestes),unit="Pearson r")
  add("Fig4B","sex_zonation_t_Pearson",cor(old$portal_vs_Central,old$testes_vs_ova),
    cor(new$portal_vs_Central,new$testes_vs_ova),unit="Pearson r")
  ledger<-do.call(rbind,rows);rownames(ledger)<-NULL
  write.csv(ledger,file.path(output_root,"fcg_metric_ledger.csv"),row.names=FALSE)
  ex<-read.csv(file.path(dir,"fcg_historical_four_examples_all_lipid_correlations.csv"),check.names=FALSE)
  write.csv(ex[ex$top6_old|ex$top6_new,],file.path(dir,"Fig5D_top_partner_edge_changes.csv"),row.names=FALSE)
  e<-new.env();load(analysis_resolve_recorded_path("results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_ttop_list.RData"),envir=e)
  kd<-read.csv(analysis_resolve_recorded_path("results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv"),check.names=FALSE)
  audit<-do.call(rbind,lapply(names(e$ttop_list),function(g){
    z<-e$ttop_list[[g]];z<-z[match(kd$ID,z$ID),]
    data.frame(knockdown=g,n_lipids=nrow(kd),max_abs_delta_saved_t=max(abs(kd[[g]]-z$t)))
  }))
  stopifnot(all(audit$max_abs_delta_saved_t<1e-10))
  write.csv(audit,file.path(dir,"fcg07_unchanged_KD_fit_reuse_audit.csv"),row.names=FALSE)
  invisible(ledger)
}

zmu_repair_auxiliary_heatmaps <- function(output_root, z) {
  ids<-do.call(rbind,lapply(c("TPC","TPP","OPC","OPP"),function(context)
    data.frame(gene=z$candidates[[context]],context=context)))
  d<-merge(ids,z$t_table,by="gene")
  d<-d[order(match(d$context,c("TPC","TPP","OPC","OPP")),d$gene),]
  annotation<-data.frame(Programme=factor(d$context,levels=c("TPC","TPP","OPC","OPP")),row.names=d$gene)
  cols<-list(Programme=c(TPC="#205cbc",TPP="#7ba7ff",OPC="#bf2e2e",OPP="#ff6969"))
  height<-2.3+nrow(d)*9/72
  configurations<-list(
    zonation_biased_genes_heatmap=c("portalova_vs_centralova","portaltestes_vs_centraltestes"),
    sex_specific_zonation_heatmap_faceted=c("portalova_vs_centralova","portaltestes_vs_centraltestes","testes_effect_zonation","testes_vs_ova"))
  for(stem in names(configurations)){
    mat<-as.matrix(d[,configurations[[stem]],drop=FALSE]);rownames(mat)<-d$gene
    colnames(mat)<-c("Female zonation\n(PP - PC t)","Male zonation\n(PP - PC t)",
      "Sex x zonation\n(t)","Sex effect\n(t)")[seq_len(ncol(mat))]
    lim<-max(abs(mat));palette<-grDevices::colorRampPalette(c("#4173d8","white","#e54141"))(100)
    width<-if(ncol(mat)==2L)7.5 else 8.5
    path<-file.path(output_root,"fcg05",paste0(stem,".pdf"))
    grDevices::pdf(path,width=width,height=height,family="Helvetica",useDingbats=FALSE)
    h<-pheatmap::pheatmap(mat,annotation_row=annotation,annotation_colors=cols,
      cluster_rows=FALSE,cluster_cols=FALSE,color=palette,breaks=seq(-lim,lim,length.out=101),
        display_numbers=TRUE,number_format="%.1f",fontsize=8,fontsize_row=7.1,annotation_names_row=FALSE,
      fontsize_number=6.2,cellheight=9,cellwidth=if(ncol(mat)==2L)115 else 70,
      angle_col=0,main="Selected sex-specific zonation genes",border_color=NA,silent=TRUE)
    grid::grid.newpage();grid::pushViewport(grid::viewport(width=.96,height=.98))
    grid::grid.draw(h$gtable);grid::popViewport();grDevices::dev.off()
    write.csv(d[,c("gene","context",configurations[[stem]])],
      file.path(output_root,"fcg05","panel_source_data",paste0(stem,".csv")),row.names=FALSE)
  }
  invisible(TRUE)
}

zmu_finish <- function(output_root,rerun07=FALSE) {
  if(rerun07 || !file.exists(file.path(output_root,"fcg07","fcg07_checkpoint.rds")))zmu_run07(output_root)
  if(!file.exists(file.path(output_root,"comparisons","fcg_historical_model_exact_replay.csv")))zmu_replay_old_model(output_root)
  if(!file.exists(file.path(output_root,"comparisons","fcg_threshold_replay_summary.csv")))zmu_replay_old_threshold(output_root)
  checkpoint<-readRDS(file.path(output_root,"fcg05","fcg05_checkpoint.rds"))
  comparison<-zmu_compare(output_root,checkpoint)
  zmu_validation_compare(output_root)
  zmu_main_panels(output_root,checkpoint,comparison)
  zmu_fixed_threshold_membership(output_root,checkpoint)
  zmu_context_class_panel(output_root, checkpoint)
  zmu_repair_auxiliary_heatmaps(output_root, checkpoint)
  zmu_metric_ledger(output_root)
  zmu_patch_manifest(output_root)
  invisible(TRUE)
}
