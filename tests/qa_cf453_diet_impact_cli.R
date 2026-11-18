source("scripts/config/analysis_release.R")
# Exact fresh-CLI replay for the compact direct-diet analysis.
args<-commandArgs(TRUE)
production<-if(length(args))args[1]else analysis_current("cf453_direct_diet")
smoke<-if(length(args)>=2)args[2]else"/private/tmp/cf453_diet_impact_smoke_20261102"
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
a<-readRDS(file.path(production,"diet_analysis.rds"));b<-readRDS(file.path(smoke,"diet_analysis.rds"))
checks<-list()
check<-function(label,x,y) {
  stopifnot(isTRUE(all.equal(x,y,tolerance=1e-12,check.attributes=TRUE)))
  checks[[length(checks)+1L]]<<-data.frame(check=label,passed=TRUE)
}
for(name in c("overview","conditional","omnibus","heatmap_selection","heatmap","example_selection",
  "example_points","example_means","lipid_scatter","primary_contrasts","conditional_contrasts","parent_threeway_replay"))
  check(name,a[[name]],b[[name]])
ma<-utils::read.csv(file.path(production,"panel_manifest.csv"))
mb<-utils::read.csv(file.path(smoke,"panel_manifest.csv"))
stopifnot(nrow(ma)==12L,nrow(mb)==12L,length(unique(ma$panel))==4L)
check("All four panel IDs",ma$panel,mb$panel)
for(file in unique(ma$source_table))check(paste("Panel source",basename(file)),
  utils::read.csv(file,check.names=FALSE),utils::read.csv(file.path(smoke,"source_data",basename(file)),check.names=FALSE))
qa<-file.path(production,"QA","fresh_cli");dir.create(qa,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(do.call(rbind,checks),file.path(qa,"replay.csv"),row.names=FALSE)
cat(length(checks),"fresh complete CLI replay checks passed.\n")
