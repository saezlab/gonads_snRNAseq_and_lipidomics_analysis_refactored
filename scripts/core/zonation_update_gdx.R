source("scripts/config/analysis_release.R")
# Paired GDX zonation on the supplied pseudobulk profiles.
zu_gdx_read_input <- function(path) {
  source("scripts/core/pseudobulk_csv.R", local = TRUE)
  o <- pb_components(path)
  m <- o$meta.data
  required <- c("group", "sample_name", "zonation_type")
  stopifnot(all(required %in% names(m)), !anyNA(m[, required]),
    setequal(unique(m$group), c("f-s-v", "f-g-v", "m-s-v", "m-g-v")),
    setequal(unique(m$zonation_type), c("Periportal", "Pericentral")),
    length(o$cells) == 44L, all(table(m$sample_name, m$zonation_type) == 1L),
    identical(colnames(o$layers$scale.data), rownames(m)),
    identical(dimnames(o$layers$scale.data), dimnames(o$layers$counts)))
  o
}

zu_gdx_fit <- function(path,version) {
  o <- zu_gdx_read_input(path)
  m <- o$meta.data
  m$mouse_id<-m$sample_name;m$sex<-ifelse(substr(m$group,1,1)=="f","female","male")
  m$surgery<-ifelse(substr(m$group,3,3)=="g","gdx","sham")
  m$zone<-ifelse(m$zonation_type=="Periportal","PP","PC")
  m$cell<-paste(m$sex,m$surgery,m$zone,sep="_")
  x<-as.matrix(o$layers$scale.data)
  cts<-as.matrix(o$layers$counts)
  eligibility<-data.frame(gene=rownames(x),detected_profiles=rowSums(cts>0),variable=apply(x,1,stats::sd)>1e-8)
  eligibility$keep<-eligibility$detected_profiles>=5L & eligibility$variable
  x<-x[eligibility$keep,,drop=FALSE]
  cells<-as.vector(outer(as.vector(outer(c("female","male"),c("sham","gdx"),paste,sep="_")),c("PC","PP"),paste,sep="_"))
  design<-stats::model.matrix(~0+factor(m$cell,levels=cells));colnames(design)<-cells;rownames(design)<-colnames(x)
  v<-function(s,t,z)as.numeric(cells==paste(s,t,z,sep="_"))
  cs<-list();families<-character()
  add<-function(name,value,family){cs[[name]]<<-value;families[name]<<-family}
  for(s in c("female","male"))for(t in c("sham","gdx"))add(paste("zone",s,t,sep="_"),v(s,t,"PP")-v(s,t,"PC"),"simple_zonation")
  for(s in c("female","male"))for(z in c("PC","PP"))add(paste("gdx",s,z,sep="_"),v(s,"gdx",z)-v(s,"sham",z),"conditional_gdx")
  for(s in c("female","male"))add(paste0("gdx_zone_",s),v(s,"gdx","PP")-v(s,"gdx","PC")-v(s,"sham","PP")+v(s,"sham","PC"),"gdx_zonation_interactions")
  add("gdx_sex_zone",cs$gdx_zone_male-cs$gdx_zone_female,"gdx_zonation_interactions")
  cm<-do.call(cbind,cs);rownames(cm)<-cells
  rho<-limma::duplicateCorrelation(x,design,block=m$mouse_id)$consensus.correlation
  base<-limma::lmFit(x,design,block=m$mouse_id,correlation=rho)
  fit<-limma::eBayes(limma::contrasts.fit(base,cm),trend=FALSE,robust=FALSE)
  results<-do.call(rbind,lapply(seq_len(ncol(cm)),function(j)data.frame(gene=rownames(x),contrast=colnames(cm)[j],
    effect=fit$coefficients[,j],se=fit$stdev.unscaled[,j]*sqrt(fit$s2.post),t=fit$t[,j],p=fit$p.value[,j],
    family=families[j],version=version,stringsAsFactors=FALSE)))
  results$q<-ave(results$p,results$family,FUN=function(p)stats::p.adjust(p,"BH"))
  summary<-do.call(rbind,lapply(unique(results$contrast),function(co){z<-results[results$contrast==co,];data.frame(version=version,contrast=co,genes=nrow(z),FDR_005=sum(z$q<=.05),positive=sum(z$q<=.05 &z$effect>0),negative=sum(z$q<=.05 &z$effect<0),rho=rho)}))
  cat("GDX",version,":",nrow(x),"eligible genes;22mice/44profiles;rho",rho,"\n");print(summary)
  list(results=results,summary=summary,eligibility=eligibility,metadata=m,rna=x,design=design,contrasts=cm,rho=rho,fit=fit,
    settings=list(input=normalizePath(path),units="difference in gene-wise scaled normalized expression",min_detected_profiles=5L,
      model="limma mouse duplicateCorrelation; default eBayes on scaled expression",primary_FDR="BH genes x three GDX-zone interactions",
      historical_status="new paired supporting analysis of both old/new inputs; not a replay of the legacy hormone-replacement notebook"))
}

zu_gdx_run <- function(root=analysis_current(), science_only=FALSE, input=NULL) {
  analysis_assert_not_approved_output(root)
  output<-file.path(root,"gdx_zonation")
  if(file.exists(file.path(output,"gdx_paired_analysis.rds")))stop("Preserve completed GDX results: choose a fresh output root")
  if (is.null(input)) {
    staged <- file.path(root, c("inputs/gdx/manifest.json", "inputs/gdx/gdx_zonation_pseudobulks_score.rds"))
    available <- staged[file.exists(staged)]
    if (!length(available)) stop("The retrospective GDX driver needs a staged GDX transport")
    input <- available[1L]
  }
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  if(science_only){
    new<-zu_gdx_fit(input,"expanded_marker_score")
    for(nm in c("results","summary","eligibility","metadata"))
      utils::write.csv(new[[nm]],file.path(output,paste0("new_",nm,".csv")),row.names=FALSE)
    for(family in unique(new$results$family)){
      z<-new$results[new$results$family==family,]
      stopifnot(max(abs(z$q-stats::p.adjust(z$p,"BH")))<1e-15)
    }
    writeLines(c("Current expanded-marker GDX supporting analysis:22 mice/44 paired profiles.",
      "Historical comparisons and the manuscript audit are optional; no old input is required for this fit."),file.path(output,"README.txt"))
    saveRDS(list(new=new),file.path(output,"gdx_paired_analysis.rds"))
    return(invisible(list(new=new)))
  }
  old<-zu_gdx_fit(analysis_path("historical_gdx_zonation_rds"),"original_two_marker")
  new<-zu_gdx_fit(input,"expanded_marker_score")
  saveRDS(list(old=old,new=new),file.path(output,"gdx_paired_analysis.rds"))
  for(nm in c("results","summary","eligibility","metadata")) {
    utils::write.csv(old[[nm]],file.path(output,paste0("old_",nm,".csv")),row.names=FALSE)
    utils::write.csv(new[[nm]],file.path(output,paste0("new_",nm,".csv")),row.names=FALSE)
  }
  comp<-merge(old$results,new$results,by=c("gene","contrast","family"),suffixes=c("_old","_new"))
  comp$delta_effect<-comp$effect_new-comp$effect_old;comp$delta_t<-comp$t_new-comp$t_old
  comp$FDR_transition<-ifelse(comp$q_old<=.05 &comp$q_new<=.05,"retained",ifelse(comp$q_old<=.05,"lost",ifelse(comp$q_new<=.05,"gained","neither")))
  utils::write.csv(comp,file.path(root,"comparisons/gdx_old_new.csv"),row.names=FALSE)
  summary<-merge(old$summary,new$summary,by="contrast",suffixes=c("_old","_new"))
  utils::write.csv(summary,file.path(root,"comparisons/gdx_summary.csv"),row.names=FALSE)
  checks<-list()
  for(a in list(old,new)) {
    V<-matrix(0,nrow(a$metadata),nrow(a$metadata))
    for(id in unique(a$metadata$mouse_id)){i<-which(a$metadata$mouse_id==id);V[i,i]<-a$rho;diag(V)[i]<-1}
    W<-solve(V);bh<-solve(t(a$design)%*%W%*%a$design,t(a$design)%*%W%*%t(a$rna))
    expected<-t(a$contrasts)%*%bh
    err<-max(abs(t(expected)-a$fit$coefficients))
    stopifnot(err<1e-10,identical(a$metadata$mouse_id,if(a$settings$input==old$settings$input)old$metadata$mouse_id else new$metadata$mouse_id))
    checks[[length(checks)+1L]]<-data.frame(check=paste(a$settings$input,"baseGLScontrast_replay"),max_abs_difference=err)
    for(fam in unique(a$results$family)){z<-a$results[a$results$family==fam,];e<-max(abs(z$q-stats::p.adjust(z$p,"BH")));stopifnot(e<1e-15);checks[[length(checks)+1L]]<-data.frame(check=paste(a$settings$input,fam,"BH"),max_abs_difference=e)}
  }
  utils::write.csv(do.call(rbind,checks),file.path(output,"independent_checks.csv"),row.names=FALSE)
  chosen<-comp[comp$family=="gdx_zonation_interactions",]
  p<-ggplot2::ggplot(chosen,ggplot2::aes(effect_old,effect_new))+
    ggplot2::geom_abline(slope=1,intercept=0,color="#9C9C9C",linewidth=.4)+
    ggplot2::geom_point(color="#6E5A83",size=.35,alpha=.3)+ggplot2::facet_wrap(~contrast,nrow=1)+
    ggplot2::labs(title="Expanded-marker GDX zonation sensitivity",subtitle="Same paired model fitted separately to original and updated assignments",
      x="Original assignment effect (scaled RNA)",y="Updated assignment effect (scaled RNA)",
      caption="Supporting analysis; no current main-figure panel consumes this GDX branch. GDX effects are gonadectomy minus sham.")+
    ggplot2::theme_classic(base_size=10,base_family="Helvetica")
  ggplot2::ggsave(file.path(output,"GDX_zonation_interaction_comparison.pdf"),p,width=11,height=4.2,device=grDevices::pdf,useDingbats=FALSE)
  ggplot2::ggsave(file.path(output,"GDX_zonation_interaction_comparison.png"),p,width=11,height=4.2,dpi=300,device=ragg::agg_png)
  writeLines(c("Supporting paired GDX zonation analysis:22mice,44profiles,fourvehicle-only sex/surgery groups.",
    "Original and expanded-marker pseudobulks are fitted with the same scaled-RNA model and a minimum of five detected profiles.",
    "Mouse duplicateCorrelation handles the two zones per animal. Positive GDX contrast means gonadectomy minus sham; zone is PP minus PC.",
    "Primary BH family: all eligible genes x female/male GDX-zone and three-way GDX-sex-zone interactions. Simple zonation and conditional GDX effects have separate four-contrast families.",
    "The legacy hormone-replacement notebook is not applicable to these inputs. These results are newly generated support, not numerical replacements for whole-hepatocyte GDX manuscript results."),file.path(output,"README.txt"))
  invisible(list(old=old,new=new,comparison=comp))
}
