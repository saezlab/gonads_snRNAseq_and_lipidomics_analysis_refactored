# Stage verified pseudobulk CSV transports; explicit local RDS imports remain possible.
source("scripts/support_functions.R")
source("scripts/config/analysis_release.R")
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))

zu_source_input <- function(cohort,source_root=analysis_path("pseudobulk_csv_root")) {
  stopifnot(cohort%in%c("fcg","gdx","wdsd"))
  # CSV always wins when both representations are present in an explicit source root.
  candidates<-c(file.path(source_root,paste0("current_",cohort,"_zonation"),"manifest.json"),
    file.path(source_root,cohort,"manifest.json"),file.path(source_root,paste0(cohort,"_zonation_pseudobulks_score.rds")))
  found<-candidates[file.exists(candidates)]
  if(!length(found))stop("Missing pseudobulk CSV transport for ",cohort," under ",source_root)
  normalizePath(found[1L],mustWork=TRUE)
}

zu_stage_inputs <- function(output=analysis_current(),source_root=analysis_path("pseudobulk_csv_root")) {
  source("scripts/core/pseudobulk_csv.R",local=TRUE)
  if(file.exists(file.path(output,"fcg05/fcg05_checkpoint.rds")))
    stop("Preserve completed release inputs: staging requires a fresh output directory")
  qc<-imports<-sources<-list()
  for(cohort in c("fcg","gdx","wdsd")) {
    s<-zu_source_input(cohort,source_root)
    files<-pb_source_files(s)
    source_dir<-dirname(s);destination<-file.path(output,"inputs",cohort)
    dir.create(destination,recursive=TRUE,showWarnings=FALSE)
    for(f in files){
      relative<-if(identical(normalizePath(f),normalizePath(s)))basename(s) else substring(normalizePath(f),nchar(normalizePath(source_dir))+2L)
      d<-file.path(destination,relative);dir.create(dirname(d),recursive=TRUE,showWarnings=FALSE)
      h<-digest::digest(f,file=TRUE,algo="sha256")
      if(file.exists(d)&&!identical(h,digest::digest(d,file=TRUE,algo="sha256")))stop("Refusing differing staged input: ",d)
      if(!file.exists(d))stopifnot(file.copy(f,d,overwrite=FALSE,copy.date=TRUE))
      stopifnot(identical(h,digest::digest(d,file=TRUE,algo="sha256")))
      imports[[length(imports)+1L]]<-data.frame(cohort=cohort,source=normalizePath(f),staged=normalizePath(d),sha256=h,bytes=file.info(d)$size)
    }
    d<-file.path(destination,basename(s));sources[[cohort]]<-normalizePath(d)
    o<-pb_components(d);m<-o$meta.data
    cts<-as.matrix(o$layers$counts);dat<-as.matrix(o$layers$data)
    err<-max(abs(dat-log1p(sweep(cts,2,colSums(cts),"/")*10000)))
    id<-if(cohort=="wdsd")m$cf453_id else m$sample_name
    stopifnot(all(is.finite(cts)),all(cts>=0),max(abs(cts-round(cts)))==0,all(is.finite(dat)),err<1e-12,
      all(table(id)==2L),all(table(m$zonation_type)==length(o$cells)/2))
    if(cohort %in% c("fcg","gdx")) {
      x<-as.matrix(o$layers$scale.data)
      m$sample<-gsub("-","_",as.character(m$orig.ident));m$sample_name<-gsub("-","_",m$sample_name)
      m<-m[,c("sample",setdiff(colnames(m),c("sample","orig.ident")))]
      colnames(x)<-m$sample
      tag<-if(cohort=="fcg")"" else "gdx_"
      utils::write.csv(data.frame(gene=rownames(x),x,check.names=FALSE),file.path(destination,paste0("scaled_RNA_",tag,"zonation.csv")),row.names=FALSE)
      utils::write.csv(m,file.path(destination,paste0("meta_data_df_",tag,"zonation.csv")),row.names=FALSE)
    }
    qc[[cohort]]<-data.frame(cohort=cohort,genes=length(o$features),profiles=length(o$cells),mice=length(unique(id)),paired_zones=TRUE,normalization_max_error=err)
  }
  utils::write.csv(do.call(rbind,imports),file.path(output,"inputs/import_manifest.csv"),row.names=FALSE)
  utils::write.csv(do.call(rbind,qc),file.path(output,"inputs/cohort_qc.csv"),row.names=FALSE)
  print(do.call(rbind,qc));invisible(list(qc=qc,sources=unlist(sources),imports=imports))
}
if(sys.nframe()==0L)zu_stage_inputs()
