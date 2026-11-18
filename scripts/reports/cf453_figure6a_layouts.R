source("scripts/config/analysis_release.R")
# Plot-only comparison of Figure6A candidates. Read Workflow16 plan first.
source("scripts/support_functions.R")
if(dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"),.libPaths()))
source("scripts/core/cf453_figure6a_layouts.R")
cf6a_run <- function(args=commandArgs(TRUE)) {
  if(any(!grepl("^--(output|parent)=.+",args))||anyDuplicated(sub("=.*$","",args)))stop("Use --output= and/or --parent= only")
  value<-function(key,default){x<-args[startsWith(args,paste0(key,"="))];if(length(x))sub("^[^=]*=","",x) else default}
  parent<-value("--parent",analysis_current("cf453_rewiring/rewiring_analysis.rds"))
  output<-value("--output",analysis_current("cf453_layouts"))
  analysis_assert_not_approved_output(output)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  original_path<-file.path(output,"retained_original_sha256.csv")
  original<-if(file.exists(original_path))utils::read.csv(original_path,stringsAsFactors=FALSE) else {
    files<-cf6a_original_files(parent)
    data.frame(path=normalizePath(files),sha256=unname(vapply(files,digest::digest,character(1L),file=TRUE,algo="sha256")))
  }
  d<-cf6a_prepare(parent)
  utils::write.csv(d,file.path(output,"common_gene_ledger.csv"),row.names=FALSE)
  source("scripts/core/cf453_figure6a_heatmap.R")
  source("scripts/core/cf453_figure6a_dots.R")
  manifest<-rbind(
    cf6a_export(cf6a_heatmap(d),output,"Fig6A_option1_heatmap",6.8,9.5,"1: paired heatmap"),
    cf6a_export(cf6a_dots(d),output,"Fig6A_option2_paired_dots",9.0,8.0,"2: paired dots"))
  utils::write.csv(manifest,file.path(output,"panel_manifest.csv"),row.names=FALSE)
  writeLines(c("# Figure6A: two layout candidates", "",
    paste("Both views use the same",nrow(d),"unique genes from",sum(d$n_anchors),"saved SD/WD anchor selections, the same model estimates, and the same direct-test annotations."), "",
    "- [Option1: four-state heatmap](panels/Fig6A_option1_heatmap.pdf)",
    "- [Option2: paired dots](panels/Fig6A_option2_paired_dots.pdf)", "",
    "Each gene is shown in female SD, female WD, male SD and male WD. SD means standard diet; WD means Western diet. PP means periportal and PC pericentral.", "",
    "The plotted effect is modelled PP-minus-PC natural-log normalized RNA expression difference. Positive values indicate PP bias and negative values PC bias. This is not log2 fold change. Both plots use untransformed saved effects; the heatmap has a single symmetric scale without row standardization.", "",
    paste("Selection badges SD, WD and Both retain the anchor-diet membership.",sum(d$selection=="Both"),"genes were selected in both diets. Programme headings describe the selection programme, not a claim that the gene retains that classification under both diets."), "",
    paste("The direct diet x sex x zone q<=0.05 mark identifies",sum(d$direct_supported),"genes with evidence that diet changes zonation differently between sexes. The other selected genes remain visible. This annotation uses the saved Workflow15 primary FDR family, with no new tests or reselection."), "",
    paste("Within each programme, genes are ordered by decreasing absolute direct three-way effect and then symbol. This display ordering is identical in both candidates and does not change gene eligibility. All",4*nrow(d),"model estimates are retained in each source CSV."), "",
    "Paired dots connect SD and WD group estimates for the same gene and sex, not repeated dietary measurements on the same animal. Female and male estimates have separate tracks. No confidence intervals are drawn in these overview candidates; statistical support is conveyed by the separately defined direct-test mark.", "",
    "PDF/SVG/300dpiPNG alternatives and source CSVs are indexed by panel_manifest.csv. The original Figure6A/B/C, statistical results and main figure bundle are preserved. These alternatives are for layout selection before replacing a manuscript panel."),file.path(output,"README.md"))
  cf6a_provenance(parent,output,original)
  invisible(list(data=d,manifest=manifest))
}
if(sys.nframe()==0L)cf6a_run()
