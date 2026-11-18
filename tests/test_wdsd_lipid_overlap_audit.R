# Independent checks for the lipid overlap audit, including collision protection.
source('scripts/support_functions.R')
source('scripts/core/cyp2c38_reference.R')
source('scripts/core/wdsd_lipid_overlap_audit.R')
stopifnot(identical(wp1_canonical_lipid('PC 34:1|PC 16:0_18:1'), 'PC(16:0_18:1)'))
stopifnot(identical(wp2_overlap_chain_key('PC(18:1_16:0)'), 'PC(16:0_18:1)'))
stopifnot(identical(wp2_overlap_chain_key('PC(18:1/16:0)'), 'PC(18:1/16:0)'))
stopifnot(identical(wp2_overlap_chain_key('Cer(18:1;O2/16:0)',TRUE), 'Cer(18:1;O2/16:0)'))
neutral <- c(600,600,600,600,600,600,600)
adduct <- c('[M-H]-','[M+H]+','[M+Na]+','[M+NH4]+','[M+CH3COO]-','[M+H-H2O]+','[M-2H]2-')
ion <- c(598.992724,601.007276,622.989218,618.033823,659.013851,582.996711,(600-2*1.007276)/2)
stopifnot(max(abs(wp2_overlap_neutral_mass(ion,adduct)-neutral))<1e-10)
stopifnot(is.na(wp2_overlap_neutral_mass(600,'[M+2H]+')))

# The same literal sum alias does not establish the same constituent chains.
make_annotations <- function(labels, mz, features=NULL) {
 z <- data.frame('Metabolite name'=labels,'Adduct type'=rep('[M+H]+',length(labels)),
                 'Average Mz'=mz,'Average Rt(min)'=rep(4,length(labels)),
                 Formula=rep('C34H66NO8P',length(labels)),INCHIKEY=rep('LIBRARY-REPRESENTATIVE',length(labels)),
                 SMILES=rep('C',length(labels)),check.names=FALSE,stringsAsFactors=FALSE)
 if(!is.null(features)) {
  z$feature_id<-features;z$canonical_lipid<-wp1_canonical_lipid(labels)
  z$lipid_class<-sub('\\(.*$','',z$canonical_lipid);z$keep_primary<-TRUE;z$exclusion_reason<-''
 }
 z
}
old <- make_annotations(c('PC 34:1|PC 16:0_18:1','PC 34:1|PC 14:0_20:1','PC 38:4|PC 18:0_20:4'), c(700,700,800))
new <- make_annotations(c('PC 34:1|PC 16:0_18:1','PC 38:4|PC 18:1_20:3'),c(700,800),c('N1','N2'))
z <- wp2_lipid_overlap_audit(wp1_canonical_lipid(old[['Metabolite name']]),new,old)
stopifnot(nrow(z$safe_map)==1L,nrow(z$conditional_compatibility_map)==2L,
          identical(z$conditional_additional_pairs$feature_id,'N2'))
stopifnot(!'PC(14:0_20:1)' %in% z$conditional_compatibility_map$reference_lipid)
stopifnot(any(z$candidate_pairs$reference_lipid=='PC(14:0_20:1)' & z$candidate_pairs$feature_id=='N1'))
# No subtotal alias substring matching is allowed.
old2<-make_annotations('PC 34:1|PC 16:0_18:1',700)
new2<-make_annotations('PC 134:1|PC 116:0_18:1',700,'X')
z2<-wp2_lipid_overlap_audit('PC(16:0_18:1)',new2,old2)
stopifnot(!any(z2$candidate_pairs$recorded_annotation_alias))
# Duplicate identities fail instead of arbitrarily selecting the first row.
bad<-rbind(new,new[1,,drop=FALSE])
stopifnot(inherits(try(wp2_lipid_overlap_audit(wp1_canonical_lipid(old[['Metabolite name']]),bad,old),silent=TRUE),'try-error'))

actual <- wp2_load_lipid_overlap_audit()
expected_local_paths <- c(
  original_annotation = "data/metabolomics/fcg/lipidomics_fcg.csv",
  new_annotation = "data/metabolomics/wdsd/lipidomics_cf453.csv",
  original_network_features = "data/metabolomics/Lipid_data_file_JPC_cleaned.csv"
)
stopifnot(identical(actual$input_paths,
                   vapply(expected_local_paths, normalizePath, character(1L), mustWork = TRUE)),
          !any(grepl("/results/", actual$input_paths, fixed = TRUE)))
stopifnot(nrow(actual$original_ledger)==1254L,nrow(actual$new_ledger)==1240L,
          sum(actual$new_ledger$keep_primary)==1177L,nrow(actual$safe_map)==752L,
          nrow(actual$conditional_compatibility_map)==854L,
          !anyDuplicated(actual$conditional_compatibility_map$feature_id),
          !anyDuplicated(actual$conditional_compatibility_map$reference_lipid))
stopifnot(all(actual$conditional_additional_pairs$recorded_annotation_alias),
          all(actual$conditional_additional_pairs$same_adduct_mass_within_5ppm),
          all(actual$conditional_additional_pairs$compatible_alias_mass_old_degree==1L),
          all(actual$conditional_additional_pairs$compatible_alias_mass_new_degree==1L),
          !any(actual$conditional_additional_pairs$feature_id %in% actual$safe_map$feature_id))
# Independent exact join uses plain final label equality against every supplied row.
old_plain <- actual$original_ledger$canonical_lipid
new_plain <- actual$new_ledger$canonical_lipid
stopifnot(sum(old_plain %in% new_plain)==752L,
          sum(old_plain %in% new_plain[actual$new_ledger$keep_primary])==752L)
# Independently enumerate exact literal alias sets; no helper parser is used.
canon_independent <- function(s) {
 vapply(s,function(t){t<-trimws(t);if(grepl(' ',t,fixed=TRUE))sub('^([^ ]+) (.*)$','\\1(\\2)',t) else t},character(1L))
}
raw_old_alias<-lapply(strsplit(actual$original_ledger[['Metabolite name']],'|',fixed=TRUE),canon_independent)
raw_new_alias<-lapply(strsplit(actual$new_ledger[['Metabolite name']],'|',fixed=TRUE),canon_independent)
all_new_alias<-unique(unlist(raw_new_alias[actual$new_ledger$keep_primary]))
independent_alias_covered<-vapply(raw_old_alias,function(t)any(t %in% all_new_alias),logical(1L))
stopifnot(sum(independent_alias_covered)==959L,
          sum(independent_alias_covered & !actual$original_ledger$exact_keep_primary)==207L)
category<-setNames(actual$category_summary$n_original_rows,actual$category_summary$overlap_category)
stopifnot(category[['unordered_chain_CL_candidate_unresolved']]==5L,
          category[['position_resolution_candidate_only']]==4L,
          category[['shared_coarse_annotation_alias_only']]==202L,
          category[['library_identifier_only_structure_mismatch']]==7L,
          category[['same_class_formula_only_possible_isomers']]==4L,
          category[['no_exported_annotation_or_same_class_formula_candidate']]==280L,
          sum(category)==1254L)
exact <- actual$candidate_pairs[actual$candidate_pairs$safe_for_primary_join,,drop=FALSE]
stopifnot(all(exact$same_formula),sum(exact$same_INCHIKEY)==740L,
          sum(!exact$same_adduct)==63L,max(exact$neutral_mass_difference_ppm)<8,
          stats::median(exact$neutral_mass_difference_ppm)>1,
          sum(actual$original_ledger$overlap_category=='unordered_chain_CL_candidate_unresolved' & actual$original_ledger$lipid_class!='CL')==0L)
cat('WDSD lipid overlap audit checks passed; primary exact752, conditional854, aliases959.\n')
