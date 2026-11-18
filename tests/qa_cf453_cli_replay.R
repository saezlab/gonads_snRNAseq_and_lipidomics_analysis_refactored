source("scripts/config/analysis_release.R")
# Verify fresh reduced-draw CLI reproduction without replacing production draws.
args<-commandArgs(TRUE)
production<-if(length(args)>=1L)args[1L] else analysis_path("historical_cf453_models")
smoke<-if(length(args)>=2L)args[2L] else "/private/tmp/cf453_full_cli_smoke_20261030"
a<-readRDS(file.path(production,"analysis.rds"))
b<-readRDS(file.path(smoke,"analysis.rds"))
rows<-list()
check<-function(name,x,y) {
  same<-isTRUE(all.equal(x,y,tolerance=1e-12,check.attributes=TRUE))
  if(!same)stop("Fresh CLI replay mismatch: ",name)
  rows[[length(rows)+1L]]<<-data.frame(check=name,passed=same)
}
for(n in c("rna_results","lipid_results","candidates","zonation","membership",
  "class_summary","class_pair_density","representatives","selection_ledger","examples","diet_edge_changes"))
  check(n,a[[n]],b[[n]])
fields<-setdiff(names(a$edges),grep("^bootstrap",names(a$edges),value=TRUE))
check("All observed edges and omission summaries",a$edges[,fields],b$edges[,fields])
check("All seven omission coefficients",a$edge_loo,b$edge_loo)
m<-utils::read.csv(file.path(smoke,"panel_manifest.csv"))
check("Fresh CLI all 30 panel groups",length(unique(m$panel)),30L)
stopifnot(nrow(m)==90L,a$settings$bootstraps==999L,b$settings$bootstraps==19L)
for(f in list.files(file.path(production,"source_data"),pattern="[.]csv$",full.names=TRUE)) {
  other<-file.path(smoke,"source_data",basename(f))
  x<-utils::read.csv(f,check.names=FALSE);y<-utils::read.csv(other,check.names=FALSE)
  # Network sidecars include uncertainty summaries, intentionally based on
  # different draw counts in the smoke and production runs.
  columns<-setdiff(names(x),grep("^bootstrap",names(x),value=TRUE))
  check(paste("Observed panel source data",basename(f)),x[,columns,drop=FALSE],y[,columns,drop=FALSE])
}
qa<-file.path(production,"QA","fresh_cli");dir.create(qa,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(do.call(rbind,rows),file.path(qa,"replay.csv"),row.names=FALSE)
writeLines(c("Fresh full CLI with 19 draws completed successfully.",
  "The saved-input marker fallback was added during that smoke and checked by --plot-only.",
  "Final observed statistics and all 30 source tables reproduce production, excluding bootstrap-summary columns with deliberately different draw counts.",
  "Production retains all 999 draws; temporary smoke uncertainty is not reported."),file.path(qa,"README.txt"))
cat(length(rows),"fresh CLI replay checks passed; production bootstrap count remains 999.\n")
