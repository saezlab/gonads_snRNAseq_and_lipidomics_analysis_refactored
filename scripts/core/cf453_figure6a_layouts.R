source("scripts/config/analysis_release.R")
# Figure 6A visual alternatives only. Read Workflow 16 layout-comparison section.
cf6a_contexts <- c("female_PP","female_PC","male_PP","male_PC")
cf6a_context_labels <- c(female_PP="Female periportal selection",female_PC="Female pericentral selection",
                         male_PP="Male periportal selection",male_PC="Male pericentral selection")

cf6a_prepare <- function(path=analysis_current("cf453_rewiring/rewiring_analysis.rds")) {
  a<-readRDS(path);x<-a$anchors
  required<-c("gene","context","sex","zone","anchor_diet","threeway_effect","threeway_p","threeway_q",
    "effect_female_SD","effect_female_WD","effect_male_SD","effect_male_WD")
  stopifnot(all(required %in% names(x)),nrow(x)>0L)
  invariant<-setdiff(required,"anchor_diet")
  rows<-lapply(unique(x$gene),function(g) {
    z<-x[x$gene==g,]
    if(nrow(unique(z[,invariant,drop=FALSE]))!=1L)stop("Inconsistent duplicate gene:",g)
    r<-z[1L,invariant,drop=FALSE]
    if(anyDuplicated(z$anchor_diet))stop("Repeated gene/anchor-diet membership:",g)
    r$selection<-if(setequal(z$anchor_diet,c("SD","WD")))"Both" else z$anchor_diet[1L]
    r$n_anchors<-nrow(z);r$direct_supported<-r$threeway_q<=.05
    r
  })
  d<-do.call(rbind,rows)
  d<-d[order(match(d$context,cf6a_contexts),-abs(d$threeway_effect),d$gene),]
  row.names(d)<-NULL
  d$display_rank<-as.integer(ave(seq_len(nrow(d)),d$context,FUN=seq_along))
  d$global_rank<-seq_len(nrow(d))
  stopifnot(nrow(d)==length(unique(x$gene)),!anyDuplicated(d$gene),
    sum(d$selection %in% c("SD","Both"))==sum(x$anchor_diet=="SD"),
    sum(d$selection %in% c("WD","Both"))==sum(x$anchor_diet=="WD"))
  cat("Figure6A layout ledger:",nrow(d),"unique genes;",sum(d$direct_supported),"direct three-way q<=0.05\n")
  print(table(d$context,d$selection))
  d
}

cf6a_long <- function(d) {
  rows<-list()
  for(s in c("female","male"))for(diet in c("SD","WD")) {
    z<-d;z$anchor_sex<-z$sex;z$sex<-s;z$diet<-diet
    z$effect<-z[[paste0("effect_",s,"_",diet)]]
    z$state<-paste(s,diet,sep="_")
    rows[[length(rows)+1L]]<-z
  }
  do.call(rbind,rows)
}

cf6a_theme <- function(size=11) {
  ggplot2::theme_minimal(base_size=size,base_family="Helvetica")+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=size+2,hjust=0),
      plot.subtitle=ggplot2::element_text(size=size-1,color="#555555",margin=ggplot2::margin(b=10)),
      plot.caption=ggplot2::element_text(size=size-2,color="#555555",hjust=0,lineheight=1.15),
      plot.title.position="plot",plot.caption.position="plot",
      panel.grid=ggplot2::element_blank(),axis.text=ggplot2::element_text(color="#303030"),
      legend.position="bottom",legend.title=ggplot2::element_text(face="bold",size=size-1),
      plot.margin=ggplot2::margin(12,16,12,12))
}

cf6a_export <- function(result,output,stem,width,height,option) {
  dir.create(file.path(output,"panels"),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(output,"source_data"),recursive=TRUE,showWarnings=FALSE)
  csv<-file.path(output,"source_data",paste0(stem,".csv"))
  utils::write.csv(result$source,csv,row.names=FALSE,na="")
  files<-file.path(output,"panels",paste0(stem,c(".pdf",".svg",".png")))
  ggplot2::ggsave(files[1L],result$plot,width=width,height=height,units="in",device=grDevices::pdf,
    family="Helvetica",useDingbats=FALSE,bg="white")
  ggplot2::ggsave(files[2L],result$plot,width=width,height=height,units="in",device=svglite::svglite,bg="white")
  ggplot2::ggsave(files[3L],result$plot,width=width,height=height,units="in",device=ragg::agg_png,dpi=300,bg="white")
  data.frame(option=option,panel=stem,file=normalizePath(files),format=c("PDF","SVG","PNG"),
    width_in=width,height_in=height,source_table=normalizePath(csv),stringsAsFactors=FALSE)
}

cf6a_original_files <- function(parent) {
  root<-dirname(parent)
  files<-c(parent,file.path(root,"CF453_Figure6_zonation_diet_panels.zip"),
    file.path(root,"panels",c("Fig6A_anchor_gene_shifts.pdf","Fig6B_lipid_class_shift.pdf","Fig6C_network_gain_loss.pdf")),
    file.path(root,"source_data",c("Fig6A_anchor_gene_shifts.csv","Fig6B_lipid_class_shift.csv","Fig6C_network_gain_loss.csv")))
  files[file.exists(files)]
}

cf6a_provenance <- function(parent,output,original) {
  files<-c(parent,"scripts/core/cf453_figure6a_layouts.R","scripts/core/cf453_figure6a_heatmap.R",
    "scripts/core/cf453_figure6a_dots.R","scripts/reports/cf453_figure6a_layouts.R",
    "scripts/reports/review_cf453_figure6a_layouts.py",
    "tests/qa_cf453_figure6a_layouts.R","docs/workflow_plans/16_cf453_diet_rewiring_figure6.md")
  files<-files[file.exists(files)]
  ledger<-data.frame(path=normalizePath(files),sha256=vapply(files,digest::digest,character(1L),file=TRUE,algo="sha256"))
  utils::write.csv(ledger,file.path(output,"provenance_sha256.csv"),row.names=FALSE)
  verification<-analysis_verify_inputs(original$path,original$sha256)
  utils::write.csv(verification,file.path(output,"retained_original_verification.csv"),row.names=FALSE)
  utils::write.csv(original,file.path(output,"retained_original_sha256.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
}
