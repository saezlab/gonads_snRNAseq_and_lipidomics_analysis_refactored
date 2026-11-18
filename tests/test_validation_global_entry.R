# Source guard, input scope and unchanged scientific results for WF07 global-only.
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/workflows/07_validation_experiments.R")
zonation_validation_run(c("--global-only","--validate-only"))
stopifnot(!any(grepl('"results/zonation/',validation_global_runtime("/private/tmp/validation_global_test"),fixed=TRUE)))
args<-commandArgs(TRUE)
if(length(args)>1L)stop("Usage: Rscript tests/test_validation_global_entry.R [FRESH_GLOBAL_OUTPUT]")
if(length(args)){
  output<-args[1L]
  read<-function(root,file)utils::read.csv(file.path(root,file),check.names=FALSE)
  old<-read("results/validation","knockdown_t_table.csv")
  new<-read(output,"knockdown_t_table.csv")
  stopifnot(setequal(old$ID,new$ID),isTRUE(all.equal(old,new[match(old$ID,new$ID),],check.attributes=FALSE,tolerance=1e-12)))
  a<-new.env();b<-new.env()
  load("results/validation/knockdown_ttop_list.RData",envir=a)
  load(file.path(output,"knockdown_ttop_list.RData"),envir=b)
  stopifnot(identical(names(a$ttop_list),names(b$ttop_list)))
  for(n in names(a$ttop_list)){
    x<-a$ttop_list[[n]];y<-b$ttop_list[[n]]
    stopifnot(setequal(rownames(x),rownames(y)),isTRUE(all.equal(x,y[match(rownames(x),rownames(y)),],check.attributes=FALSE,tolerance=1e-12)))
  }
  stopifnot(isTRUE(all.equal(read("results/validation","global_network_validation_scores.csv"),read(output,"global_network_validation_scores.csv"),tolerance=1e-12)))
  f<-"structural_consistency/consistency_structural_regression_summary.csv"
  x<-read("results/validation",f);y<-read(output,f);x<-x[x$gene%in%y$gene,]
  stopifnot(setequal(y$gene,c("CYP2B9","CYP2C38","CYP2D9")),
    isTRUE(all.equal(x,y[match(paste(x$gene,x$lipid_class),paste(y$gene,y$lipid_class)),],check.attributes=FALSE,tolerance=1e-12)))
  err<-tryCatch(zonation_validation_run(c("--global-only",paste0("--output=",output))),error=identity)
  stopifnot(inherits(err,"error"),grepl("fresh --output",conditionMessage(err),fixed=TRUE))
  cat("KD fits, global ULM, global structural regressions and overwrite guard match the independent preserved branch.\n")
}
cat("WF07 independent entry checks passed.\n")
