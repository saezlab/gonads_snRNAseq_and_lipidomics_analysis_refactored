source("scripts/config/analysis_release.R")
# Direct CF453 diet interactions. Read Workflow 15 plan before editing.
cfdi_fit <- function(parent=analysis_current("cf453_models/models.rds"), provenance_path=NULL) {
  a<-readRDS(parent)
  stopifnot(nrow(a$metadata)==28L,nrow(a$profile_metadata)==56L,
    all(table(a$metadata$sex,a$metadata$diet)==7L),nrow(a$modeled_rna)>0L,
    identical(rownames(a$modeled_rna),a$gene_eligibility$gene[a$gene_eligibility$keep]),
    identical(colnames(a$modeled_rna),a$profile_metadata$profile_id))
  if(is.null(provenance_path))provenance_path<-file.path(dirname(parent),"provenance_sha256.csv")
  prior<-utils::read.csv(provenance_path,stringsAsFactors=FALSE)
  paths<-unname(a$input_paths)
  expected<-prior$sha256[match(paths,prior$path)]
  input_verification<-analysis_verify_inputs(paths,expected)
  pm<-a$profile_metadata; cells<-colnames(a$design)
  cell<-function(s,d,z)as.numeric(cells==paste(s,d,z,sep="_"))
  ds<-function(z)cell("male","WD",z)-cell("female","WD",z)-cell("male","SD",z)+cell("female","SD",z)
  dz<-function(s)cell(s,"WD","PP")-cell(s,"WD","PC")-cell(s,"SD","PP")+cell(s,"SD","PC")
  primary<-cbind(diet_sex=(ds("PP")+ds("PC"))/2,
    diet_zone=(dz("female")+dz("male"))/2,diet_sex_zone=ds("PP")-ds("PC"))
  conditional<-cbind(diet_sex_PC=ds("PC"),diet_sex_PP=ds("PP"),
    diet_zone_female=dz("female"),diet_zone_male=dz("male"))
  rownames(primary)<-rownames(conditional)<-cells
  stopifnot(qr(primary)$rank==3L)
  base_fit<-limma::lmFit(a$modeled_rna,a$design,block=pm$mouse_id,correlation=a$rho)
  fit_primary<-limma::eBayes(limma::contrasts.fit(base_fit,primary),trend=TRUE,robust=TRUE)
  fit_conditional<-limma::eBayes(limma::contrasts.fit(base_fit,conditional),trend=TRUE,robust=TRUE)
  overview<-cf6_tests(fit_primary,rep("primary_direct_diet",3L))
  detail<-cf6_tests(fit_conditional,rep("conditional_direct_diet",4L))
  omnibus<-data.frame(gene=rownames(a$modeled_rna),F=fit_primary$F,p=fit_primary$F.p.value,
    q=stats::p.adjust(fit_primary$F.p.value,"BH"),df1=3L,df2=fit_primary$df.prior+fit_primary$df.residual,stringsAsFactors=FALSE)
  ranked<-omnibus[order(omnibus$q,omnibus$p,omnibus$gene),]
  heatmap_selection<-head(ranked[ranked$q<=.05,],20L)
  heatmap_selection$rank<-seq_len(nrow(heatmap_selection))
  heatmap<-detail[detail$gene %in% heatmap_selection$gene,]
  heatmap$rank<-heatmap_selection$rank[match(heatmap$gene,heatmap_selection$gene)]
  heatmap$joint_q<-omnibus$q[match(heatmap$gene,omnibus$gene)]
  three<-overview[overview$contrast=="diet_sex_zone",]
  heatmap$threeway_q<-three$q[match(heatmap$gene,three$gene)]
  selection<-list(); selected<-character()
  for(name in c("diet_sex","diet_zone","diet_sex_zone","diet_sex_zone")) {
    r<-overview[overview$contrast==name & overview$q<=.05 & !overview$gene %in% selected,]
    r<-r[order(-abs(r$effect),r$q,r$p,r$gene),]
    if(!nrow(r))next
    chosen<-r[1L,];selected<-c(selected,chosen$gene)
    selection[[length(selection)+1L]]<-data.frame(gene=chosen$gene,selected_contrast=name,
      q_primary=chosen$q,p_primary=chosen$p,effect_primary=chosen$effect,
      metric=if(name=="diet_sex")"Mean-zone expression" else "PP - PC zonation",
      example_rank=length(selection)+1L,stringsAsFactors=FALSE)
  }
  example_selection<-if(length(selection))do.call(rbind,selection) else
    data.frame(gene=character(),selected_contrast=character(),q_primary=numeric(),p_primary=numeric(),effect_primary=numeric(),metric=character(),example_rank=integer())
  points<-means<-list()
  for(j in seq_len(nrow(example_selection))) {
    s<-example_selection[j,]; g<-match(s$gene,rownames(a$modeled_rna))
    for(sex in c("female","male"))for(diet in c("SD","WD")) {
      mice<-a$metadata$mouse_id[a$metadata$sex==sex & a$metadata$diet==diet]
      pp<-pm$profile_id[match(paste(mice,"PP"),paste(pm$mouse_id,pm$zone))]
      pc<-pm$profile_id[match(paste(mice,"PC"),paste(pm$mouse_id,pm$zone))]
      v<-if(s$selected_contrast=="diet_sex")(cell(sex,diet,"PP")+cell(sex,diet,"PC"))/2 else cell(sex,diet,"PP")-cell(sex,diet,"PC")
      value<-if(s$selected_contrast=="diet_sex")as.numeric(a$modeled_rna[g,pp]+a$modeled_rna[g,pc])/2 else as.numeric(a$modeled_rna[g,pp]-a$modeled_rna[g,pc])
      estimate<-as.numeric(base_fit$coefficients[g,]%*%v)
      se<-sqrt(as.numeric(t(v)%*%base_fit$cov.coefficients%*%v)*fit_primary$s2.post[g])
      crit<-stats::qt(.975,fit_primary$df.total[g])
      points[[length(points)+1L]]<-data.frame(mouse_id=mice,gene=s$gene,sex=sex,diet=diet,value=value,
        metric=s$metric,example_rank=s$example_rank,selected_contrast=s$selected_contrast,q_primary=s$q_primary)
      means[[length(means)+1L]]<-data.frame(gene=s$gene,sex=sex,diet=diet,estimate=estimate,se=se,
        low=estimate-crit*se,high=estimate+crit*se,metric=s$metric,example_rank=s$example_rank,
        selected_contrast=s$selected_contrast,q_primary=s$q_primary)
    }
  }
  lr<-a$lipid_results
  li<-lr[lr$contrast=="sex_diet",];lf<-lr[lr$contrast=="diet_female",];lm<-lr[lr$contrast=="diet_male",]
  stopifnot(identical(li$lipid,lf$lipid),identical(li$lipid,lm$lipid))
  sxsd<-lr[lr$contrast=="sex_SD",];sxwd<-lr[lr$contrast=="sex_WD",]
  lipids<-data.frame(lipid=li$lipid,lipid_label=li$lipid_label,lipid_class=li$lipid_class,
    effect_female=lf$effect,effect_male=lm$effect,q_interaction=li$q,p_interaction=li$p,
    effect_interaction=li$effect,low_interaction=li$low,high_interaction=li$high,
    low_female=lf$low,high_female=lf$high,low_male=lm$low,high_male=lm$high,
    sex_SD=sxsd$effect[match(li$lipid,sxsd$lipid)],sex_WD=sxwd$effect[match(li$lipid,sxwd$lipid)])
  stopifnot(max(abs(lipids$effect_interaction-(lipids$effect_male-lipids$effect_female)))<1e-12)
  old<-a$rna_results[a$rna_results$contrast=="sex_diet_zone",]
  new<-three[match(old$gene,three$gene),]
  replay<-data.frame(metric=c("threeway_effect","threeway_p"),max_abs_difference=c(max(abs(old$effect-new$effect)),max(abs(old$p-new$p))))
  stopifnot(all(replay$max_abs_difference<1e-12))
  cat("Direct diet interaction RNA discoveries (new pooled 3-contrast FDR family):\n")
  print(stats::aggregate(I(q<=.05)~contrast,overview,sum))
  cat("Joint 3-df interaction genes:",sum(omnibus$q<=.05),"; lipid diet x sex interactions:",sum(lipids$q_interaction<=.05),"\n")
  print(example_selection)
  list(metadata=a$metadata,profile_metadata=pm,modeled_rna=a$modeled_rna,design=a$design,rho=a$rho,
    base_fit=base_fit,primary_contrasts=primary,conditional_contrasts=conditional,
    fit_primary=fit_primary,fit_conditional=fit_conditional,overview=overview,conditional=detail,
    omnibus=omnibus,heatmap_selection=heatmap_selection,heatmap=heatmap,example_selection=example_selection,
    example_points=if(length(points))do.call(rbind,points) else data.frame(),
    example_means=if(length(means))do.call(rbind,means) else data.frame(),lipid_scatter=lipids,
    lipid_matrix=a$lipid_matrix,lipid_fit=a$lipid_fit,lipid_design=a$lipid_design,
    lipid_contrasts=a$lipid_contrasts,lipid_results=lr,parent_threeway_replay=replay,
    source_paths=c(parent=normalizePath(parent),parent_provenance=normalizePath(provenance_path),a$input_paths),
    input_verification=input_verification,
    settings=list(date="2026-11-02",primary_FDR="BH across genes x 3 direct diet interactions",
      conditional_FDR="BH across genes x 4 conditional interactions",omnibus_FDR="BH across genes for joint 3-df diet interaction",
      RNA_units="Difference in natural-log normalized expression",lipid_units="Difference in log2 whole-liver lipid abundance",
      examples="Largest absolute estimated interaction among primary FDR<=0.05 genes; ties by q,p,gene; distinct genes for diet_sex, diet_zone and diet_sex_zone twice",
      heatmap="First up to20 joint-FDR genes ordered q,p,gene",rho_source="Validated Workflow14 mouse-block correlation"))
}

cfdi_save<-function(a,output) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  saveRDS(a,file.path(output,"diet_analysis.rds"))
  for(name in names(a)[vapply(a,is.data.frame,logical(1L))])cf6_write(a[[name]],file.path(output,"tables",paste0(name,".csv")))
  cf6_write(data.frame(cell=rownames(a$primary_contrasts),a$primary_contrasts,check.names=FALSE),file.path(output,"tables","primary_contrasts.csv"))
  cf6_write(data.frame(cell=rownames(a$conditional_contrasts),a$conditional_contrasts,check.names=FALSE),file.path(output,"tables","conditional_contrasts.csv"))
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  cfdi_provenance(a,output)
  invisible(a)
}

cfdi_provenance<-function(a,output) {
  code<-c("scripts/core/cf453_diet_impact.R","scripts/core/cf453_diet_impact_plots.R",
    "scripts/core/cf453_diet_impact_report.R","scripts/core/cf453_standalone.R",
    "scripts/workflows/15_cf453_diet_impact_figure6.R","tests/qa_cf453_diet_impact.R",
    "tests/qa_cf453_diet_impact_cli.R","scripts/reports/package_cf453_diet_impact.py",
    "docs/workflow_plans/15_cf453_diet_impact_figure6.md")
  files<-c(unname(a$source_paths),code[file.exists(code)])
  cf6_write(analysis_provenance(files),
    file.path(output,"provenance_sha256.csv"))
  invisible(TRUE)
}
