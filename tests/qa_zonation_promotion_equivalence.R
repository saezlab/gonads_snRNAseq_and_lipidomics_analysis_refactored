# Compare a fresh current-only FCG/GDX run with the approved release. Differences
# in file locations, retrospective columns and display row order are excluded.
source("scripts/config/analysis_release.R")
qa_zonation_promotion_equivalence <- function(fresh,approved=analysis_path("archived_zonation_bundle"),output=NULL){
  rows<-list()
  normalize<-function(x){
    keys<-intersect(c("version","gene","context","lipid","lipid_class","gonad_bias","zonation"),names(x))
    if(length(keys)&&nrow(x))x<-x[do.call(order,x[,keys,drop=FALSE]),,drop=FALSE]
    rownames(x)<-NULL;x
  }
  check<-function(name,x,y){
    pass<-isTRUE(all.equal(normalize(x),normalize(y),tolerance=1e-12,check.attributes=FALSE))
    rows[[length(rows)+1L]]<<-data.frame(check=name,passed=pass)
    if(!pass)stop("Fresh scientific run differs: ",name)
  }
  files<-c("fcg05/t_table_zonation.csv","fcg05/ttop_full_long.csv","fcg05/sex_specific_zonation_genes_summary.csv",
    "fcg05/optimal_threshold_metrics.csv","fcg05/pathway_enrichment_fresh.csv","fcg05/summary_statistics.csv",
    "fcg07/zonation_network_validation_scores.csv","fcg07/Mtmr11_male_PP_prespecified_ULM.csv",
    "fcg07/Mtmr11_male_PP_prespecified_source.csv","fcg07/Mtmr11_male_PP_ULM_background_source.csv",
    "gdx_zonation/new_results.csv","gdx_zonation/new_summary.csv","gdx_zonation/new_eligibility.csv")
  for(f in files)check(f,read.csv(file.path(approved,f),check.names=FALSE),read.csv(file.path(fresh,f),check.names=FALSE))
  for(f in list.files(file.path(approved,"fcg05/networks"),pattern="[.]csv$"))
    check(paste("Network",f),read.csv(file.path(approved,"fcg05/networks",f),check.names=FALSE),read.csv(file.path(fresh,"fcg05/networks",f),check.names=FALSE))
  for(f in list.files(file.path(fresh,"fcg05/panel_source_data"),pattern="[.]csv$")){
    old<-read.csv(file.path(approved,"fcg05/panel_source_data",f),check.names=FALSE)
    new<-read.csv(file.path(fresh,"fcg05/panel_source_data",f),check.names=FALSE)
    if("version"%in%names(old))old<-old[old$version%in%c("revised","revised_reoptimized"),,drop=FALSE]
    common<-intersect(names(old),names(new))
    check(paste("Panel source",f),old[,common,drop=FALSE],new[,common,drop=FALSE])
  }
  result<-do.call(rbind,rows)
  if(!is.null(output)){dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE);write.csv(result,output,row.names=FALSE)}
  cat(nrow(result),"fresh current-only scientific/source equivalence checks passed\n")
  invisible(result)
}
if(sys.nframe()==0L){
  args<-commandArgs(TRUE)
  if(length(args)<1L||length(args)>2L)stop("Usage: Rscript tests/qa_zonation_promotion_equivalence.R FRESH_OUTPUT [QA_CSV]")
  qa_zonation_promotion_equivalence(args[1],output=if(length(args)==2L)args[2] else NULL)
}
