# WP2 follow-up: annotation-level FCG versus WDSD overlap audit.
# No abundance is merged and no ambiguous structural candidate is rescued.
# Source scripts/support_functions.R, then cyp2c38_reference.R first.

wp2_overlap_aliases <- function(x) {
  lapply(strsplit(as.character(x), "|", fixed = TRUE), function(z) unique(wp1_canonical_lipid(z)))
}

wp2_overlap_chain_key <- function(x, allow_position_loss = FALSE) {
  vapply(wp1_canonical_lipid(x), function(s) {
    cls <- sub("\\(.*$", "", s)
    if (!grepl("(", s, fixed = TRUE)) return(s)
    structure <- sub("^[^(]+\\(", "", s)
    structure <- sub("\\)$", "", structure)
    if (allow_position_loss) {
      # Diagnostic only: positional information is deliberately lost.
      if (!cls %in% c("DG", "TG", "PC", "PE", "PA", "PG", "PI", "PS", "BMP", "CL", "MLCL")) return(s)
      parts <- strsplit(structure, "[_/]", perl = TRUE)[[1L]]
    } else {
      if (!grepl("_", structure, fixed = TRUE) || grepl("/", structure, fixed = TRUE)) return(s)
      parts <- strsplit(structure, "_", fixed = TRUE)[[1L]]
    }
    paste0(cls, "(", paste(sort(parts), collapse = "_"), ")")
  }, character(1L), USE.NAMES = FALSE)
}

wp2_overlap_zero_chain_key <- function(x) {
  vapply(wp1_canonical_lipid(x), function(s) {
    cls <- sub("\\(.*$", "", s)
    if (!cls %in% c("LPC", "LPE", "LPG", "LPI", "LPS", "LPA")) return(s)
    # This removes position resolution and therefore remains a candidate only.
    s <- sub("/0:0\\)$", ")", s)
    sub("\\(0:0/", "(", s)
  }, character(1L), USE.NAMES = FALSE)
}

wp2_overlap_neutral_mass <- function(mz, adduct) {
  # Fiehn Lab ion-adduct calculator mass shifts (Da); unsupported labels remain NA.
  # [M+2H]+ has inconsistent charge notation and is deliberately not guessed.
  shifts <- c("[M-H]-" = -1.007276, "[M+H]+" = 1.007276,
              "[M+Na]+" = 22.989218, "[M+NH4]+" = 18.033823,
              "[M+CH3COO]-" = 59.013851, "[M+Hac-H]-" = 59.013851,
              "[M+H-H2O]+" = -17.003289)
  result <- as.numeric(mz) - shifts[as.character(adduct)]
  double <- !is.na(adduct) & adduct == "[M-2H]2-"
  result[double] <- 2 * mz[double] + 2 * 1.007276
  unname(result)
}

wp2_overlap_equal_nonempty <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  z <- outer(a, b, "==")
  z[is.na(z)] <- FALSE
  z[is.na(a) | !nzchar(a), ] <- FALSE
  z[, is.na(b) | !nzchar(b)] <- FALSE
  z
}

wp2_lipid_overlap_audit <- function(reference_lipids, wdsd_ledger, original_annotation) {
  if (!exists("wp1_canonical_lipid", mode = "function")) stop("Source the canonical lipid helper first.")
  old_names <- wp1_canonical_lipid(reference_lipids)
  required <- c("Metabolite name", "Adduct type", "Average Mz", "Average Rt(min)", "Formula", "INCHIKEY", "SMILES")
  if (!all(required %in% names(original_annotation)) ||
      !all(c(required, "feature_id", "canonical_lipid", "keep_primary", "exclusion_reason") %in% names(wdsd_ledger))) stop("Required original/new annotation fields are missing.")
  old_raw_names <- wp1_canonical_lipid(original_annotation[["Metabolite name"]])
  if (anyNA(old_names) || anyDuplicated(old_names) || anyDuplicated(old_raw_names) ||
      anyDuplicated(wdsd_ledger$feature_id) || anyDuplicated(wdsd_ledger$canonical_lipid)) stop("Unique old labels and new feature/annotation keys are required.")
  index <- match(old_names, old_raw_names)
  if (anyNA(index)) stop("Every original network label must have an exact raw FCG annotation row.")
  old <- original_annotation[index, , drop = FALSE]
  new <- wdsd_ledger
  if (!identical(wp1_canonical_lipid(new[["Metabolite name"]]), as.character(new$canonical_lipid))) stop("New canonical annotation ledger disagrees with raw names.")
  if (!is.logical(new$keep_primary) || anyNA(new$keep_primary)) stop("New eligibility flags must be complete logical values.")
  old$reference_lipid <- as.character(reference_lipids)
  old$canonical_lipid <- old_names
  old$lipid_class <- sub("\\(.*$", "", old_names)
  old$raw_annotation_row <- index
  exact <- wp2_overlap_equal_nonempty(old_names, new$canonical_lipid)
  chain_old <- wp2_overlap_chain_key(old_names)
  chain_new <- wp2_overlap_chain_key(new$canonical_lipid)
  chain <- wp2_overlap_equal_nonempty(chain_old, chain_new)
  position <- wp2_overlap_equal_nonempty(wp2_overlap_chain_key(old_names, TRUE), wp2_overlap_chain_key(new$canonical_lipid, TRUE))
  zero <- wp2_overlap_equal_nonempty(wp2_overlap_zero_chain_key(old_names), wp2_overlap_zero_chain_key(new$canonical_lipid))
  old_alias <- wp2_overlap_aliases(old[["Metabolite name"]])
  new_alias <- wp2_overlap_aliases(new[["Metabolite name"]])
  alias <- matrix(FALSE, nrow(old), nrow(new))
  # Alias inventory joins independently of final annotation resolution.
  key_to_new <- split(rep(seq_len(nrow(new)), lengths(new_alias)), unlist(new_alias, use.names = FALSE))
  for (i in seq_len(nrow(old))) {
    candidates <- unique(unlist(key_to_new[old_alias[[i]]], use.names = FALSE))
    if (length(candidates)) alias[i, candidates] <- TRUE
  }
  inchi <- wp2_overlap_equal_nonempty(old$INCHIKEY, new$INCHIKEY)
  smiles <- wp2_overlap_equal_nonempty(old$SMILES, new$SMILES)
  formula <- wp2_overlap_equal_nonempty(old$Formula, new$Formula)
  same_class <- outer(old$lipid_class, new$lipid_class, "==")
  same_adduct <- wp2_overlap_equal_nonempty(old[["Adduct type"]], new[["Adduct type"]])
  ppm <- outer(old[["Average Mz"]], new[["Average Mz"]], function(a, b) abs(a - b)/a*1e6)
  mass5 <- same_adduct & is.finite(ppm) & ppm <= 5
  old$neutral_mass_diagnostic <- wp2_overlap_neutral_mass(old[["Average Mz"]], old[["Adduct type"]])
  new$neutral_mass_diagnostic <- wp2_overlap_neutral_mass(new[["Average Mz"]], new[["Adduct type"]])
  neutral_ppm <- outer(old$neutral_mass_diagnostic, new$neutral_mass_diagnostic, function(a,b) abs(a-b)/a*1e6)
  # Nonempty supplied IDs can be representative library structures. Conflicting
  # reported acyl-chain/position/ether resolution blocks rescue by ID equality.
  identifier <- (inchi | smiles) & same_class
  formula_class <- formula & same_class
  retained <- matrix(rep(new$keep_primary, each = nrow(old)), nrow(old), nrow(new))
  annotation_candidate <- chain | position | zero | alias | identifier
  composition_candidate <- annotation_candidate | formula_class
  any_candidate <- composition_candidate | mass5
  pair <- which(any_candidate, arr.ind = TRUE)
  i <- pair[, 1L]; j <- pair[, 2L]
  flags <- list(exact_annotation = exact, unordered_chain_key = chain,
                positional_resolution_key = position, zero_chain_resolution_key = zero,
                recorded_annotation_alias = alias, same_INCHIKEY = inchi,
                same_SMILES = smiles, same_formula = formula, same_class = same_class,
                same_adduct = same_adduct, same_adduct_mass_within_5ppm = mass5)
  pairs <- data.frame(reference_lipid = old$reference_lipid[i], reference_canonical = old_names[i],
                      feature_id = new$feature_id[j], wdsd_canonical = new$canonical_lipid[j],
                      original_raw_annotation = old[["Metabolite name"]][i],
                      wdsd_raw_annotation = new[["Metabolite name"]][j],
                      new_keep_primary = new$keep_primary[j],
                      original_class = old$lipid_class[i], wdsd_class = new$lipid_class[j],
                      original_adduct = old[["Adduct type"]][i], wdsd_adduct = new[["Adduct type"]][j],
                      observed_mz_difference_ppm = ppm[cbind(i,j)],
                      neutral_mass_difference_ppm = neutral_ppm[cbind(i,j)],
                      retention_time_difference_min = new[["Average Rt(min)"]][j] - old[["Average Rt(min)"]][i],
                      stringsAsFactors = FALSE)
  for (nm in names(flags)) pairs[[nm]] <- flags[[nm]][cbind(i,j)]
  pairs$shared_recorded_aliases <- vapply(seq_along(i), function(k) paste(intersect(old_alias[[i[k]]], new_alias[[j[k]]]), collapse = ";"), character(1L))
  pairs$old_annotation_candidate_degree <- rowSums(annotation_candidate & retained)[i]
  pairs$new_annotation_candidate_degree <- colSums(annotation_candidate & retained)[j]
  pairs$unordered_old_multiplicity <- as.integer(table(chain_old)[chain_old[i]])
  pairs$unordered_new_multiplicity <- as.integer(table(chain_new)[chain_new[j]])
  pairs$original_has_exact_primary_match <- old_names[i] %in% new$canonical_lipid[new$keep_primary]
  pairs$safe_for_primary_join <- pairs$exact_annotation & pairs$new_keep_primary
  pairs$candidate_interpretation <- ifelse(pairs$exact_annotation, "same_reported_final_annotation",
     ifelse(pairs$unordered_chain_key, "chain_order_candidate_CL_pairing_or_key_collisions_unresolved",
     ifelse(pairs$positional_resolution_key | pairs$zero_chain_resolution_key, "position_resolution_compatibility_only",
     ifelse(pairs$recorded_annotation_alias, "shared_coarse_alias_not_same_molecular_feature",
     ifelse((pairs$same_INCHIKEY | pairs$same_SMILES) & pairs$same_class, "library_identifier_equal_but_reported_structure_differs",
     ifelse(pairs$same_formula & pairs$same_class, "same_class_formula_possible_isomers_not_identity",
            "mass_adduct_candidate_not_annotation_identity"))))))
  counts <- function(x) rowSums(x & retained)
  old$exact_retained_candidates <- counts(exact)
  old$unordered_retained_candidates <- counts(chain)
  old$position_resolution_candidates <- counts(position | zero)
  old$recorded_alias_candidates <- counts(alias)
  old$identifier_candidates <- counts(identifier)
  old$same_class_formula_candidates <- counts(formula_class)
  old$same_adduct_mass5_candidates <- counts(mass5)
  old$annotation_candidate_degree <- counts(annotation_candidate)
  old$exact_feature_id <- new$feature_id[match(old_names, new$canonical_lipid)]
  old$exact_keep_primary <- old$exact_retained_candidates > 0L
  old$overlap_category <- ifelse(old$exact_keep_primary, "exact_final_annotation",
     ifelse(old$unordered_retained_candidates > 0L, "unordered_chain_CL_candidate_unresolved",
     ifelse(old$position_resolution_candidates > 0L, "position_resolution_candidate_only",
     ifelse(old$recorded_alias_candidates > 0L, "shared_coarse_annotation_alias_only",
     ifelse(old$identifier_candidates > 0L, "library_identifier_only_structure_mismatch",
     ifelse(old$same_class_formula_candidates > 0L, "same_class_formula_only_possible_isomers",
            "no_exported_annotation_or_same_class_formula_candidate"))))))
  # Exact primary keys alone form the safe one-to-one join used by analysis.
  safe <- data.frame(reference_lipid = old$reference_lipid[old$exact_keep_primary],
                     canonical_lipid = old_names[old$exact_keep_primary],
                     feature_id = old$exact_feature_id[old$exact_keep_primary],
                     mapping_rule = rep("formatting_only_exact_final_annotation",sum(old$exact_keep_primary)), stringsAsFactors = FALSE)
  if (anyDuplicated(safe$reference_lipid) || anyDuplicated(safe$feature_id)) stop("Primary join must be one-to-one without merging features.")
  # Conditional compatibility sensitivity: exact map plus an entirely separate,
  # bidirectionally unique recorded-alias/same-adduct/5 ppm candidate graph.
  # Degrees are counted across ALL original rows and retained new rows before
  # excluding exact-assigned nodes, preventing reuse of a strict-matched feature.
  compatible <- alias & same_class & same_adduct & mass5 & retained
  unique_compatible <- compatible & rowSums(compatible)==1L &
    matrix(rep(colSums(compatible)==1L, each=nrow(old)), nrow(old))
  unique_compatible[old$exact_keep_primary, ] <- FALSE
  unique_compatible[, new$feature_id %in% safe$feature_id] <- FALSE
  sensitivity_pair <- which(unique_compatible, arr.ind=TRUE)
  additional <- data.frame(reference_lipid=old$reference_lipid[sensitivity_pair[,1L]],
                           canonical_lipid=old_names[sensitivity_pair[,1L]],
                           feature_id=new$feature_id[sensitivity_pair[,2L]],
                           mapping_rule=rep("conditional_unique_coarse_alias_same_adduct_5ppm_NOT_identity",nrow(sensitivity_pair)),
                           stringsAsFactors=FALSE)
  sensitivity_map <- rbind(safe, additional)
  if(anyDuplicated(sensitivity_map$reference_lipid) || anyDuplicated(sensitivity_map$feature_id)) stop("Conditional compatibility sensitivity reuses a feature.")
  pairs$in_conditional_compatibility_sensitivity <- paste(pairs$reference_lipid,pairs$feature_id,sep="__") %in%
    paste(additional$reference_lipid,additional$feature_id,sep="__")
  pairs$compatible_alias_mass_old_degree <- rowSums(compatible)[i]
  pairs$compatible_alias_mass_new_degree <- colSums(compatible)[j]
  identifier_examples <- pairs[!pairs$exact_annotation & pairs$new_keep_primary &
                                (pairs$same_INCHIKEY | pairs$same_SMILES), , drop=FALSE]
  new$exact_original_match <- new$canonical_lipid %in% old_names
  new$original_annotation_candidate_degree <- colSums(annotation_candidate)
  new$original_same_class_formula_degree <- colSums(formula_class)
  new$unordered_key_multiplicity <- as.integer(table(chain_new)[chain_new])
  stages <- list(formatting_only_exact = exact,
                 plus_unordered_chain_candidates = exact | chain,
                 plus_position_resolution_candidates = exact | chain | position | zero,
                 plus_recorded_alias_candidates = exact | chain | position | zero | alias,
                 plus_library_identifier_candidates = annotation_candidate,
                 plus_same_class_formula_candidates = composition_candidate,
                 plus_same_adduct_mass5_candidates = any_candidate)
  summary <- do.call(rbind, lapply(names(stages), function(nm) {
    x <- stages[[nm]]; p <- x & retained
    data.frame(rule = nm, original_rows_covered_all_new = sum(rowSums(x)>0L),
               original_rows_covered_retained_new = sum(rowSums(p)>0L),
               original_rows_with_unique_candidate = sum(rowSums(p)==1L),
               one_to_one_candidate_pairs = sum(p & rowSums(p)==1L & matrix(rep(colSums(p)==1L,each=nrow(old)),nrow(old))),
               primary_identity_safe_pairs = nrow(safe), stringsAsFactors = FALSE)
  }))
  class_summary <- as.data.frame(table(old$lipid_class, old$overlap_category), stringsAsFactors = FALSE)
  names(class_summary) <- c("lipid_class", "overlap_category", "n_original_rows")
  class_summary <- class_summary[class_summary$n_original_rows>0L, , drop = FALSE]
  category_summary <- as.data.frame(table(old$overlap_category), stringsAsFactors = FALSE)
  names(category_summary) <- c("overlap_category", "n_original_rows")
  classes <- sort(unique(c(old$lipid_class,new$lipid_class[new$keep_primary])))
  class_coverage <- do.call(rbind,lapply(classes,function(cl) {
    ix <- old$lipid_class==cl
    data.frame(lipid_class=cl,n_original=sum(ix),n_new_retained=sum(new$lipid_class==cl & new$keep_primary),
               n_exact=sum(ix & old$exact_keep_primary),n_original_unmatched=sum(ix & !old$exact_keep_primary),
               n_recorded_alias_covered=sum(ix & old$recorded_alias_candidates>0L),stringsAsFactors=FALSE)
  }))
  metrics <- data.frame(metric=c("original_supplied_rows","original_network_rows","new_supplied_rows","new_retained_rows",
                                  "new_excluded_rows","exact_against_all_supplied_new","exact_against_retained_new",
                                  "literal_recorded_alias_covered","conditional_compatibility_map_rows","conditional_additional_rows"),
                        value=c(nrow(original_annotation),nrow(old),nrow(new),sum(new$keep_primary),sum(!new$keep_primary),
                                sum(rowSums(exact)>0L),nrow(safe),sum(counts(alias)>0L),nrow(sensitivity_map),nrow(additional)),
                        stringsAsFactors=FALSE)
  limitations <- c(
    "No raw peak absence or limit-of-detection inference follows from absence of an exported annotation.",
    "Same sum composition, formula or precursor mass cannot distinguish acyl-chain isomers; candidates are not extra validated features.",
    "The supplied INCHIKEY/SMILES fields can be identical despite different reported chain/ether/positional resolution; these library identifiers are not independent chemical confirmation.",
    "Underscores normally omit sn-position, but both datasets have multiple CL rows with the same unordered chain quartet; CL group/pair assignment and feature identity are unresolved, so none is automatically rescued.",
    "Retention time and 5 ppm precursor checks are descriptive diagnostics; no cross-run retention-time calibration or peak/spectrum re-identification was performed.",
    "Zero-chain removal and slash/underscore conversion lose positional resolution and remain candidates only.")
  list(original_ledger = old, new_ledger = new, candidate_pairs = pairs,
       safe_map = safe, conditional_compatibility_map = sensitivity_map,
       conditional_additional_pairs = pairs[pairs$in_conditional_compatibility_sensitivity, , drop=FALSE],
       identifier_structure_mismatch_examples = identifier_examples, stage_summary = summary, class_summary = class_summary,
       category_summary = category_summary, class_coverage = class_coverage, metrics = metrics, limitations = limitations,
       sources = c("https://pmc.ncbi.nlm.nih.gov/articles/PMC7707175/",
                   "https://www.lipidmaps.org/lipid_nomenclature/rules/hierarchy",
                   "https://www.fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator"))
}

wp2_load_lipid_overlap_audit <- function(root = ".", data = NULL) {
  root <- normalizePath(root, mustWork = TRUE)
  # Derive the frozen original feature universe and new annotation ledger from
  # source inputs so the pooled phase can run before any zonation output exists.
  if (is.null(data)) {
    intake_helpers <- new.env(parent = baseenv())
    sys.source(file.path(root, "scripts/core/wdsd_intake.R"), envir = intake_helpers)
    data <- intake_helpers$wp2_load_inputs(root)
  }
  if (!is.list(data) || !is.data.frame(data$lipid_ledger)) {
    stop("data must be a current WDSD intake with a lipid annotation ledger.")
  }
  paths <- c(
    original_annotation = file.path(root, "data/metabolomics/fcg/lipidomics_fcg.csv"),
    new_annotation = file.path(root, "data/metabolomics/wdsd/lipidomics_cf453.csv"),
    original_network_features = file.path(root, "data/metabolomics/Lipid_data_file_JPC_cleaned.csv")
  )
  if (any(!file.exists(paths))) stop("Missing lipid overlap input(s): ",
    paste(paths[!file.exists(paths)], collapse = ", "))
  reference <- utils::read.csv(paths[["original_network_features"]], nrows = 0L,
                              check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(reference)[1L], "Unnamed: 0") || ncol(reference) != 1255L ||
      anyDuplicated(names(reference)[-1L])) {
    stop("Expected an original mouse-ID column and 1,254 unique lipid annotations.")
  }
  ledger <- data$lipid_ledger
  original <- utils::read.csv(paths[["original_annotation"]], check.names = FALSE, stringsAsFactors = FALSE)
  current_new <- utils::read.csv(paths[["new_annotation"]],check.names=FALSE,stringsAsFactors=FALSE)
  annotation_columns <- intersect(names(current_new),names(ledger))
  if(nrow(current_new)!=nrow(ledger) || !isTRUE(all.equal(current_new[,annotation_columns,drop=FALSE],ledger[,annotation_columns,drop=FALSE],check.attributes=FALSE))) stop("Saved WDSD annotation ledger does not replay current source annotations.")
  audit <- wp2_lipid_overlap_audit(names(reference)[-1L], ledger, original)
  audit$input_paths <- vapply(paths, normalizePath, character(1L), mustWork = TRUE)
  audit
}

wp2_write_lipid_overlap_audit <- function(audit, output) {
  if (dir.exists(output)) stop("Use a fresh overlap audit output directory.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  tables <- c("original_ledger", "new_ledger", "candidate_pairs", "safe_map", "conditional_compatibility_map", "conditional_additional_pairs",
              "identifier_structure_mismatch_examples", "stage_summary", "class_summary", "category_summary", "class_coverage", "metrics")
  for (nm in tables) utils::write.csv(audit[[nm]], file.path(output,paste0(nm,".csv")), row.names = FALSE, na = "")
  if(!is.null(audit$input_paths)) utils::write.csv(data.frame(input=names(audit$input_paths),path=unname(audit$input_paths)),file.path(output,"input_paths.csv"),row.names=FALSE)
  writeLines(audit$limitations, file.path(output,"limitations.txt"))
  writeLines(audit$sources, file.path(output,"nomenclature_sources.txt"))
  saveRDS(audit, file.path(output,"lipid_overlap_audit.rds"))
  invisible(output)
}
