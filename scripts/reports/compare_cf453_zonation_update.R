source("scripts/config/analysis_release.R")
# Observed old/new CF453 audit tables; no model or resampling is executed here.
cf453_compare_zonation_update <- function(root=analysis_current()) {
  output<-file.path(root,"comparisons");dir.create(output,recursive=TRUE,showWarnings=FALSE)
  paths<-c(old_models=file.path(analysis_path("historical_cf453_models"),"models.rds"),
    new_models=file.path(root,"cf453_models/models.rds"),
    old_direct=file.path(analysis_path("historical_cf453_direct"),"diet_analysis.rds"),
    new_direct=file.path(root,"cf453_direct_diet/diet_analysis.rds"),
    old_rewiring=file.path(analysis_path("historical_cf453_rewiring"),"rewiring_analysis.rds"),
    new_rewiring=file.path(root,"cf453_rewiring/rewiring_analysis.rds"))
  objects<-lapply(paths,readRDS);names(objects)<-names(paths)
  om<-objects$old_models;nm<-objects$new_models;od<-objects$old_direct;nd<-objects$new_direct
  old<-objects$old_rewiring;new<-objects$new_rewiring
  write<-function(x,name)utils::write.csv(x,file.path(output,paste0("cf453_",name,".csv")),row.names=FALSE,na="")
  joined<-function(x,y,by){x$evaluated_old<-TRUE;y$evaluated_new<-TRUE;merge(x,y,by=by,all=TRUE,suffixes=c("_old","_new"),sort=TRUE)}
  compare_tests<-function(x,y,stage) {
    result<-do.call(rbind,lapply(union(unique(x$contrast),unique(y$contrast)),function(k){
      a<-x[x$contrast==k,];b<-y[y$contrast==k,];shared<-intersect(a$gene,b$gene)
      ai<-match(shared,a$gene);bi<-match(shared,b$gene);sa<-a$gene[a$q<=.05];sb<-b$gene[b$q<=.05]
      data.frame(stage=stage,contrast=k,tested_old=nrow(a),tested_new=nrow(b),shared_tested=length(shared),
        q05_old=length(sa),q05_new=length(sb),q05_shared=length(intersect(sa,sb)),
        q05_old_only=length(setdiff(sa,sb)),q05_new_only=length(setdiff(sb,sa)),
        effect_pearson=stats::cor(a$effect[ai],b$effect[bi]),t_pearson=stats::cor(a$t[ai],b$t[bi]),
        max_abs_effect_change=max(abs(b$effect[bi]-a$effect[ai])))
    }));result
  }
  write(rbind(compare_tests(om$rna_results,nm$rna_results,"Workflow14"),
              compare_tests(od$overview,nd$overview,"Workflow15_direct")),"rna_contrast_comparison")
  write(joined(om$primary_interactions,nm$primary_interactions,c("gene","contrast","family")),"primary_gene_changes")
  write(joined(od$overview,nd$overview,c("gene","contrast","family")),"direct_diet_gene_changes")
  write(joined(om$gene_eligibility,nm$gene_eligibility,"gene"),"gene_eligibility")
  write(data.frame(version=c("old","updated_score"),eligible_genes=c(nrow(om$modeled_rna),nrow(nm$modeled_rna)),
    rho=c(om$rho,nm$rho),joint_direct_q05=c(sum(od$omnibus$q<=.05),sum(nd$omnibus$q<=.05)),
    lipid_interaction_q05=c(sum(od$lipid_scatter$q_interaction<=.05),sum(nd$lipid_scatter$q_interaction<=.05))),"model_summary")
  write(joined(old$anchors,new$anchors,c("gene","anchor_diet","sex","zone","context")),"anchor_changes")
  write(joined(old$programmes,new$programmes,c("context","sex","zone","gene")),"programme_gene_changes")
  cl<-joined(old$class_summary,new$class_summary,c("context","sex","zone","lipid_class"))
  cl$delta_E_change<-cl$delta_E_new-cl$delta_E_old
  cl$comparison_scope<-"Each version's selected programmes; gene denominators may differ"
  write(cl,"class_association_changes")
  write(joined(old$network_genes,new$network_genes,c("context","sex","zone")),"network_representatives")
  selected<-unique(rbind(old$network_edges[,c("context","gene","lipid")],new$network_edges[,c("context","gene","lipid")]))
  key<-function(x)paste(x$context,x$gene,x$lipid,sep="|")
  pairs<-joined(old$edges,new$edges,c("context","sex","zone","gene","lipid","lipid_label","lipid_class"))
  pairs$selected_old<-key(pairs)%in%key(old$network_edges);pairs$selected_new<-key(pairs)%in%key(new$network_edges)
  write(pairs[key(pairs)%in%key(selected),],"selected_network_edge_changes")
  pair_summary<-matched_classes<-list()
  for(context in union(old$programmes$context,new$programmes$context)){
    g<-intersect(old$programmes$gene[old$programmes$context==context],new$programmes$gene[new$programmes$context==context])
    z<-pairs[pairs$context==context & pairs$gene%in%g & !is.na(pairs$evaluated_old) & !is.na(pairs$evaluated_new),]
    for(diet in c("SD","WD")){
      x<-z[[paste0("r_",diet,"_old")]];y<-z[[paste0("r_",diet,"_new")]]
      pair_summary[[length(pair_summary)+1L]]<-data.frame(context=context,diet=diet,n_shared_genes=length(g),
        n_identical_pairs=nrow(z),coefficient_pearson=if(nrow(z)>1)stats::cor(x,y) else NA_real_,
        median_abs_change=if(nrow(z))stats::median(abs(y-x)) else NA_real_,max_abs_change=if(nrow(z))max(abs(y-x)) else NA_real_)
    }
    if(!nrow(z))next
    for(lc in unique(z$lipid_class)){
      inclass<-z$lipid_class==lc
      E<-function(v)mean(abs(v[inclass])>=.65)-mean(abs(v[!inclass])>=.65)
      esdo<-E(z$r_SD_old);ewdo<-E(z$r_WD_old);esdn<-E(z$r_SD_new);ewdn<-E(z$r_WD_new)
      matched_classes[[length(matched_classes)+1L]]<-data.frame(context=context,lipid_class=lc,n_genes=length(g),
        n_pairs_class=sum(inclass),n_pairs_rest=sum(!inclass),E_SD_old=esdo,E_WD_old=ewdo,E_SD_new=esdn,E_WD_new=ewdn,
        delta_E_old=ewdo-esdo,delta_E_new=ewdn-esdn,delta_E_change=(ewdn-esdn)-(ewdo-esdo),
        comparison_scope="Identical intersection programme genes and lipid pairs; observed only, no bootstrap")
    }
  }
  write(do.call(rbind,pair_summary),"matched_pair_correlations")
  write(do.call(rbind,matched_classes),"matched_programme_class_changes")
  genes<-unique(c("Mtmr11","Bmper","Cyp2c38","Cyp2b9","Cyp2b13","Cyp2d9","Cyp3a41a","Fmo3",old$network_genes$gene,new$network_genes$gene))
  write(joined(om$rna_results[om$rna_results$gene%in%genes,],nm$rna_results[nm$rna_results$gene%in%genes,],c("gene","contrast","family")),"named_gene_statistics")
  checks<-data.frame(check=c("Mouse metadata unchanged","Paired profile metadata unchanged","Primary lipid matrix unchanged",
    "Lipid model results unchanged","All anchor bins preserve five-gene cap","B uses999 bootstrap draws"),
    passed=c(identical(om$metadata,nm$metadata),identical(om$profile_metadata,nm$profile_metadata),
      identical(om$lipid_matrix,nm$lipid_matrix),identical(om$lipid_results,nm$lipid_results),
      all(new$anchor_bins$selected==pmin(new$anchor_bins$eligible,5L)),new$settings$B==999L))
  stopifnot(all(checks$passed));write(checks,"chain_checks")
  write(data.frame(role=names(paths),path=normalizePath(paths),sha256=vapply(paths,digest::digest,character(1L),file=TRUE,algo="sha256")),"comparison_provenance")
  cat("CF453 old/new comparisons written:",length(list.files(output,pattern="^cf453_.*[.]csv$")),"CSV files; chain checks pass\n")
  invisible(output)
}
if(sys.nframe()==0L){args<-commandArgs(TRUE);cf453_compare_zonation_update(if(length(args))args[1L] else analysis_current())}
