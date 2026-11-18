# Scientific invariants for current expanded-marker results; no historical fits,
# private manuscript, resampling replay or external Python runtime is required.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/config/analysis_release.R")
qa_current_zonation_release <- function(root=analysis_current(),cf453_root=root,output=NULL) {
  rows<-list()
  check<-function(label,value){rows[[length(rows)+1L]]<<-data.frame(check=label,passed=isTRUE(value));if(!isTRUE(value))stop(label)}
  f<-readRDS(file.path(root,"fcg05/fcg05_checkpoint.rds"))
  check("FCG43 paired animals",length(unique(f$metadata$sample_name))==43L&&all(table(f$metadata$sample_name)==2L))
  check("FCG expression and metadata order",identical(colnames(f$rna),f$metadata$sample)&&all(is.finite(as.matrix(f$rna))))
  check("Seven FCG contrasts",length(f$ttop)==7L&&ncol(f$contrasts)==7L)
  for(n in names(f$ttop)){
    z<-f$ttop[[n]];check(paste("FCG BH family",n),max(abs(z$adj.P.Val-p.adjust(z$P.Value,"BH")))<1e-12)
  }
  check("FCG selected genes measured",all(unlist(f$candidates)%in%rownames(f$rna)))
  check("FCG threshold is objective maximum",f$threshold==f$threshold_metrics$threshold[which.max(f$threshold_metrics$combined_score)])
  bg<-read.csv(file.path(root,"fcg07/Mtmr11_male_PP_ULM_background_source.csv"),check.names=FALSE)
  u<-read.csv(file.path(root,"fcg07/Mtmr11_male_PP_prespecified_ULM.csv"),check.names=FALSE)
  r<-cor(bg$ULM_weight,bg$mtmr11);stat<-r*sqrt((nrow(bg)-2)/(1-r^2))
  check("Prespecified ULM independent t replay",abs(stat-u$score)<1e-10&&abs(2*pt(-abs(stat),nrow(bg)-2)-u$p_value)<1e-12)
  check("Matched and background lipids distinguished",sum(bg$matched_in_vivo)==u$n_lipids&&nrow(bg)==u$n_KD_background&&all(bg$ULM_weight[!bg$matched_in_vivo]==0))
  g<-readRDS(file.path(root,"gdx_zonation/gdx_paired_analysis.rds"))$new
  check("Current GDX22 paired animals",length(unique(g$metadata$mouse_id))==22L&&all(table(g$metadata$mouse_id)==2L))
  for(family in unique(g$results$family)){
    z<-g$results[g$results$family==family,];check(paste("GDX BH family",family),max(abs(z$q-p.adjust(z$p,"BH")))<1e-12)
  }
  a<-readRDS(file.path(cf453_root,"cf453_models/models.rds"))
  check("CF45328 independent animals",nrow(a$metadata)==28L&&all(table(a$metadata$sex,a$metadata$diet)==7L))
  source("scripts/core/pseudobulk_csv.R",local=TRUE)
  current_manifest<-jsonlite::read_json(analysis_path("current_wdsd_zonation_transport",mustWork=TRUE))
  recorded<-a$input_paths[["zonation"]]
  if(grepl("[.]rds$",recorded,ignore.case=TRUE)){
    v<-analysis_verify_inputs(recorded,current_manifest$source$sha256)
    same_source<-all(v$verified)
  }else{
    transported<-pb_components(recorded,verify=TRUE)
    same_source<-identical(transported$source$sha256,current_manifest$source$sha256)
  }
  check("CF453 uses current expanded-marker bytes",same_source)
  d<-readRDS(file.path(cf453_root,"cf453_direct_diet/diet_analysis.rds"))
  check("Direct diet overview contrasts",setequal(unique(d$overview$contrast),c("diet_sex","diet_zone","diet_sex_zone")))
  n<-readRDS(file.path(cf453_root,"cf453_rewiring/rewiring_analysis.rds"))
  check("RNA-only anchor bins preserve five-gene cap",all(n$anchor_bins$selected==pmin(n$anchor_bins$eligible,5L)))
  check("Class resampling count retained",n$settings$B==999L)
  result<-do.call(rbind,rows)
  if(!is.null(output)){dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE);write.csv(result,output,row.names=FALSE)}
  cat(nrow(result),"current-release scientific checks passed\n")
  invisible(result)
}
if(sys.nframe()==0L){
  args<-commandArgs(TRUE)
  if(any(!grepl("^--(root|cf453-root|output)=.+",args)))stop("Use --root=, --cf453-root= and/or --output=")
  value<-function(key,default){x<-args[startsWith(args,paste0(key,"="))];if(length(x))sub("^[^=]*=","",x[1]) else default}
  root<-value("--root",analysis_current())
  qa_current_zonation_release(root,value("--cf453-root",root),value("--output",file.path(root,"QA/current_release_scientific_checks.csv")))
}
