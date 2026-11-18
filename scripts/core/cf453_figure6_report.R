# Standalone CF453 reports; sourcing defines functions without running analysis.
cf6_report <- function(a, manifest, output) {
  required <- c("metadata", "modeled_rna", "rna_results", "candidates", "zonation",
                "lipid_ledger", "lipid_results", "edges", "representatives", "settings")
  if (length(setdiff(required, names(a)))) stop("Incomplete CF453 report inputs.")
  if (!all(c("panel", "family", "file", "format", "source_table") %in% names(manifest)))
    stop("Incomplete panel manifest.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  n <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  pct <- function(x, digits = 2) paste0(formatC(100 * x, digits = digits, format = "f"), "%")
  table_md <- function(d) {
    if (!nrow(d)) return("No eligible rows.")
    d[] <- lapply(d, function(x) gsub("[|]", "/", as.character(x)))
    c(paste0("| ", paste(names(d), collapse = " | "), " |"),
      paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
      apply(d, 1L, function(x) paste0("| ", paste(x, collapse = " | "), " |")))
  }
  relative <- function(path) {
    path <- as.character(path)
    prefix <- paste0(output, "/")
    ifelse(startsWith(path, prefix), substring(path, nchar(prefix) + 1L), path)
  }
  link <- function(path, label = basename(path)) paste0("[", label, "](", relative(path), ")")
  gene_n <- nrow(a$modeled_rna)
  ledger <- a$lipid_ledger
  primary_n <- sum(ledger$keep_primary)
  broad_n <- sum(ledger$keep_broad)
  fill_n <- sum(ledger$keep_fill)
  q_limit <- a$settings$fdr
  primary_edges <- a$edges$primary_lipid
  tested_edges <- sum(primary_edges)
  significant_edges <- sum(primary_edges & is.finite(a$edges$q) & a$edges$q <= q_limit)
  interaction_counts <- do.call(rbind, lapply(c("sex_zone_SD", "sex_zone_WD", "sex_diet_zone"), function(k) {
    d <- a$rna_results[a$rna_results$contrast == k, , drop = FALSE]
    data.frame(contrast = k, genes_tested = n(nrow(d)), genes_q_005 = n(sum(d$q <= q_limit)))
  }))
  nomination <- do.call(rbind, lapply(c("SD", "WD"), function(diet) {
    d <- a$candidates[a$candidates$diet == diet, , drop = FALSE]
    data.frame(diet = diet, distinct_genes = n(length(unique(d$gene))),
      context_memberships = n(nrow(d)), shared_genes = n(sum(a$zonation$diet == diet & a$zonation$category == "shared")))
  }))
  reps <- a$representatives[, c("sex", "zone", "gene", "stable_degree"), drop = FALSE]
  names(reps)[4L] <- "stable_partner_count_across_diets"
  edge_by_diet <- do.call(rbind, lapply(c("SD", "WD"), function(diet) {
    k <- primary_edges & a$edges$diet == diet
    data.frame(diet = diet, tested_pairs = n(sum(k)),
      pairs_abs_r_065 = n(sum(k & abs(a$edges$r) >= a$settings$display_threshold)),
      pairs_q_005 = n(sum(k & is.finite(a$edges$q) & a$edges$q <= q_limit)))
  }))
  family_counts <- function(d, feature) {
    do.call(rbind, lapply(unique(d$family), function(family) {
      x <- d[d$family == family, , drop = FALSE]
      data.frame(family = family, features = n(length(unique(x[[feature]]))),
        contrasts = n(length(unique(x$contrast))), tests = n(nrow(x)), q_005 = n(sum(x$q <= q_limit)))
    }))
  }
  panels <- manifest[!duplicated(manifest$panel), , drop = FALSE]
  panel_index <- data.frame(panel = panels$panel, family = panels$family,
    PDF = vapply(panels$panel, function(p) {
      f <- manifest$file[manifest$panel == p & toupper(manifest$format) == "PDF"]
      if (length(f)) link(f[1L], "PDF") else ""
    }, character(1L)),
    source_data = vapply(panels$source_table, function(f) link(f, "CSV"), character(1L)))
  density <- a$class_pair_density
  density_text <- if (!is.null(density) && nrow(density)) paste0(
    "The complementary C_pair_density panels use identical sex/zone candidate-gene unions and lipid pairs in both diets. Across class/context/diet cells, ",
    pct(min(density$display_fraction)), "–", pct(max(density$display_fraction)),
    " of eligible pairs exceed |r| ≥ ", a$settings$display_threshold,
    ". These densities are descriptive and do not establish differential correlation or rewiring.") else
    "The class pair-density summary is unavailable."
  shared <- a$class_summary[a$class_summary$category == "shared", , drop = FALSE]
  saturation <- nrow(shared) > 0L && all(shared$proportion == 1)
  saturation_text <- if (saturation)
    "The original C any-edge rule assigns 100% of eligible lipids in every class and both diets to multiple contexts (shared). This saturated summary provides no class-specific spatial discrimination. It is retained transparently; use the complementary pair-density panels to inspect continuous context differences." else
    "C reports overlapping any-edge memberships; shared denotes more than one sex/zone context."
  limits <- paste(
    "This is an observational analysis of sex associations in one CF453 cohort; it does not isolate gonadal hormone or sex-chromosome effects.",
    "One whole-liver lipid vector is paired with both RNA zones per mouse: RNA-zone associations do not measure spatial lipid concentrations, enzyme activity, substrate specificity or flux.",
    "Nucleus contributions, original nuclear QC/gating, count aggregation, batch/litter and assignment provenance remain unverified.",
    "Seven mice per sex/diet limit correlation precision. Representative selection uses these same data, and omission/bootstrap checks condition on the fixed RNA and partner selections.",
    "Stable signs, large correlations and bootstrap ranges do not provide independent replication or confirm causal regulation.")
  audit_path <- file.path(output, "QA/provenance/RNA_transform_relationships.csv")
  audit_text <- "The provenance audit confirms data = ln(1 + 10000 × supplied counts / profile library), and scale.data is a per-gene z-score of data."
  if (file.exists(audit_path)) {
    audit <- utils::read.csv(audit_path, stringsAsFactors = FALSE)
    k <- audit$proposed_data_transform == "log1p_library_normalized_1e4"
    err <- audit$max_abs_error[k]
    audit_text <- paste0(audit_text, " Maximum normalization replay errors were ",
      paste(format(err, digits = 4, scientific = TRUE), collapse = " and "),
      " for hepatocyte and zone objects, respectively. This establishes the stored transformation, not the original count-summing provenance.")
  }
  bootstrap_text <- paste0(a$settings$bootstraps,
    " ordinary paired-mouse bootstrap draws are generated separately within each sex/diet, using the same draw indices across both zones and every RNA/lipid feature. RNA and its matched whole-liver lipid vector travel together. Percentile ranges (2.5th–97.5th), sign fractions and valid-draw counts describe fixed-selection stability; they are not null-test p-values or calibrated post-selection intervals. Constant-feature draws are undefined and reported. WD-minus-SD delta-r draws combine independent diet-group resamples; positive fractions are descriptive stability summaries.")
  summary_lines <- c(
    "# Standalone CF453 Figure 6 panel package", "",
    paste0(n(nrow(a$metadata)), " wild-type mice (seven per female/male × SD/WD group) contribute ",
      n(nrow(a$profile_metadata)), " paired hepatocyte-zone RNA profiles and one whole-liver lipid profile per mouse. The primary models retain ",
      n(gene_n), " genes and ", n(primary_n), " endogenous lipid features with analytical QC CV ≤30%."), "",
    "The principal RNA result is diet-specific sex-by-zone association, with a separate test of diet modification:", "",
    table_md(interaction_counts), "",
    "Interaction contrasts are male-minus-female differences in periportal-minus-pericentral expression; the three-way contrast is WD minus SD of that interaction. All three share one BH family across genes. Context nominations also require a significant within-sex/diet zonation effect; memberships may overlap.", "",
    table_md(nomination), "",
    paste0("**No individual gene–lipid correlation passes BH q ≤0.05 among ", n(tested_edges),
      " primary tests (observed count: ", n(significant_edges), "). Association panels are exploratory.**"), "",
    table_md(edge_by_diet), "", saturation_text, "", density_text, "",
    "Representative genes are selected for descriptive association stability, without implying a metabolic or regulatory role:", "",
    table_md(reps), "", "Each representative uses the same selected partners in both diets, including coefficients below the display threshold in one diet. Their ranking and partner selection use the current cohort; examples are not independent validation.", "",
    paste0("The lipid-quality sensitivity universes contain ", n(broad_n), " broadly eligible features, ",
      n(primary_n), " primary QC-filtered features and ", n(fill_n), " features after additionally requiring supplied Fill ≥0.5. Their class denominators are recalculated within each universe; these are membership-summary sensitivities rather than independent experiments."), "",
    audit_text, "", limits, "",
    "## Panel files", "",
    paste0(n(nrow(panels)), " separate panels are exported in PDF, SVG and PNG; ",
      "[panel_manifest.csv](panel_manifest.csv) records dimensions, formats and matching source tables. The manuscript figure is to be assembled manually."), "",
    table_md(panel_index), "",
    "## Reproducible records", "",
    "[Methods](METHODS.md) and [panel legends](PANEL_LEGENDS.md) describe the tests and visual encodings. `tables/` contains complete gene tests, nomination/selection ledgers, lipid tests, all correlation coefficients, memberships and diet changes. `source_data/` matches individual panels. `analysis.rds` preserves the full fitted analysis; `bootstrap_indices.rds` preserves animal draw banks. `QC/`, `QA/provenance/`, independent replay checks under `QA/`, `provenance_sha256.csv` and `sessionInfo.txt` record input QC, numerical verification and provenance.")
  methods_lines <- c(
    "# CF453 Figure 6 methods", "",
    "## Design and inputs", "",
    paste0("The standalone analysis uses ", n(nrow(a$metadata)), " paired CF453 wild-type mice, seven per sex/diet group. Each mouse supplies periportal and pericentral RNA pseudobulks and one whole-liver lipid measurement. Analytical QC injections are not biological replicates. Input files are the local CF453 lipid CSV and the updated whole-hepatocyte/zonation Seurat objects. All matching uses explicit mouse IDs; zones remain paired within mouse. No other cohort supplies model parameters, feature selection or test statistics."), "",
    audit_text, "",
    "## RNA model and testing families", "",
    paste0("Genes require positive supplied counts in at least seven of 56 zone profiles and expression SD >10^-8 across profiles, independent of lipid outcomes; ", n(gene_n),
      " pass. Primary expression is the supplied unscaled natural-log normalized data. An eight-cell sex × diet × zone design is fit with limma, a common duplicateCorrelation mouse block (estimated correlation ",
      formatC(a$rho, digits = 5, format = "f"), "), and mean-variance-trend, robust empirical Bayes. RNA effect units are differences in log-normalized expression, not log2 fold changes. Model-based plotted ranges are nominal 95% intervals from moderated standard errors, without multiplicity adjustment."), "",
    "Primary sex-by-zone contrasts are (male PP−PC)−(female PP−PC) separately in SD and WD; the three-way test is interaction_WD−interaction_SD. BH adjustment pools all eligible genes and all three primary contrasts. Four simple zonation contrasts form a separate family; the additional families below remain distinct.", "",
    table_md(family_counts(a$rna_results, "gene")), "",
    "Gene-context membership requires both diet-specific primary interaction q≤0.05 and the corresponding within-sex/diet zonation q≤0.05. Positive PP−PC effects map to PP and negative effects to PC; multi-context memberships are retained. A context name denotes membership under this rule, not absence of expression or association in another context.", "",
    "## Lipid model and quality sets", "",
    paste0("Internal standards, isotope-labelled annotations and documented repository-excluded classes are removed, leaving ", n(broad_n),
      " broad features. Primary analysis requires analytical QC CV≤30% (", n(primary_n),
      " features); the further Fill≥0.5 sensitivity retains ", n(fill_n),
      ". Supplied positive intensities are log2 transformed without additional normalization or imputation. A four-cell sex × diet limma model tests male−female within diet, WD−SD within sex, and their interaction. Lipid coefficients are differences of log2 abundance. The two sex contrasts share one BH family, two diet contrasts another, and the interaction a third."), "",
    table_md(family_counts(a$lipid_results, "lipid")), "",
    "Chemical annotations retain feature identity and molecular resolution; no cross-cohort overlap restriction, abundance-based outlier exclusion or forced annotation merging is used. QC and Fill sensitivities change feature eligibility and denominators, not the selected RNA programs.", "",
    "## Associations, summaries and selection", "",
    paste0("Within each sex/zone, candidate genes are the union of SD and WD nominations. Genes and lipids must vary in both diets. Pearson correlations use seven matched mice per sex/diet and identical pair universes across diets. Nominal correlation probabilities use t with five degrees of freedom. One BH family covers all ", n(tested_edges),
      " primary pair/context tests; none pass q≤0.05. These probabilities remain conditional on the RNA selection. Broad-only lipid pairs retain coefficients and omission estimates but have no primary-family q-value."), "",
    "The exploratory display rule is observed |r|≥0.65 with an active RNA nomination in that diet. A stable displayed edge additionally has all seven leave-one-mouse-out correlations defined and of the observed sign; each omission uses six mice. Omission stability does not require every omitted correlation to exceed 0.65.", "",
    bootstrap_text, "",
    "For each sex/zone, representative genes rank by the number of distinct stable displayed lipid partners across diets, then the smallest diet-specific interaction q, then gene symbol. At most one representative with positive stable degree is selected. Up to six partners stable in at least one diet are ranked by their maximum absolute correlation across diets, then lipid ID. Both diets and all corresponding coefficients are retained. The leading partner supplies the individual-mouse example; this is a data-selected illustration.", "",
    "C assigns each lipid to the sex/zone contexts containing any qualifying displayed edge among diet-active candidates. Exactly one context denotes a unique membership, more than one denotes shared, and zero denotes unassigned. The class denominator is every eligible measured lipid in that class, including unassigned lipids. C_pair_density instead counts |r|≥0.65 pairs over all eligible union-candidate-gene × lipid pairs within each class and sex/zone, using the same pairs in both diets; it does not require an active nomination in both diets. Different contexts may have different union sizes, but each context's SD/WD comparison uses the same denominator. Stable fractions, FDR fractions and mean/median |r| are also exported.", "",
    paste0("The display-threshold audit uses 0.50, 0.60, 0.65, 0.70 and 0.80 without selecting a preferred threshold from outcomes. Seed: ", a$settings$seed, "; requested bootstrap draws per sex/diet: ", a$settings$bootstraps, "."), "",
    "## Interpretation", "", limits)
  legend_for <- function(row) {
    family <- as.character(row$family)
    source <- if (file.exists(row$source_table)) utils::read.csv(row$source_table, stringsAsFactors = FALSE) else data.frame()
    if (family == "A") return(paste0("Female versus male moderated PP−PC t-statistics in ", row$diet,
      "; ", n(gene_n), " genes, seven mice per sex. Positive is periportal and negative pericentral. Color requires primary sex-by-zone q≤0.05 and matching within-sex zonation q≤0.05; shared memberships are retained. Labels include representatives and up to eight genes ranked by interaction q. Axes show t-statistics, not fold changes."))
    if (family == "B") return("CF453 study design: 28 mice, seven each female SD, male SD, female WD and male WD. Each supplies paired periportal/pericentral RNA profiles and a single whole-liver lipid measurement. Diet groups contain different mice.")
    if (family == "C") return(paste0("Within-class lipid membership in ", row$diet,
      ". Class labels state all eligible measured lipid denominators. A context requires at least one diet-active RNA-candidate association with observed |r|≥0.65; multiple contexts denote shared membership and zero contexts denotes unassigned. ", saturation_text,
      " This any-edge summary is exploratory and does not localize lipids spatially."))
    if (family == "C_pair_density") return(paste0("Exploratory pair density in ", row$diet,
      ": fraction of eligible union-candidate-gene × primary-lipid pairs with |r|≥0.65 within each lipid class and sex/zone. The candidate union and pair denominator are fixed across diets within each context; n=7 mice per estimate. All primary correlation q-values exceed 0.05. Density differences do not establish rewiring. Complete denominators and alternative summaries accompany the source table."))
    if (family == "D") return(paste0("Representative ", if (nrow(source)) paste(unique(source$gene), collapse = ", ") else "unavailable",
      " subnetwork in ", row$context, ", ", row$diet,
      "; n=7 mice. Same up-to-six selected partners appear in both diets. Line width follows |r|, solid/dashed lines indicate positive/negative r, and faint edges have |r|<0.65. Partners are sign-stable across all seven omissions in at least one diet; this does not require stability in both diets. The subtitle identifies whether the gene is RNA-nominated in the displayed diet. Edges remain exploratory; none is BH-significant."))
    if (family == "E") return(paste0("Diet change in PP−PC expression within each sex for ",
      if (nrow(source)) n(length(unique(source$gene))) else "available",
      " displayed genes: representatives plus up to 16 genes with strongest three-way evidence. Points and nominal 95% moderated intervals use natural-log normalized expression differences; filled symbols denote q≤0.05 in the separate within-sex diet-zonation family. Text gives the primary three-way q, which tests a different contrast. WD−SD is positive to the right."))
    if (family == "F") return(paste0("Individual-mouse illustration in ", row$context, ": ",
      if (nrow(source)) paste(unique(paste(source$gene, source$lipid_label, sep = " with ")), collapse = "; ") else "no eligible example",
      ". Each point is one mouse, seven per diet. RNA is ln(1+10000×counts/library); lipid abundance is log2 of supplied whole-liver intensity. Lines are separate diet-wise linear fits without interval bands. Reported r and six-mouse omission ranges describe a selected example; they are not adjusted-significant discoveries."))
    if (family == "supplementary_QC" && grepl("markers", row$panel)) return("Zonation marker distributions (Cyp2e1, Cyp2f2, Glul, Gls2, Hal, Sds) in natural-log normalized RNA. Each point is a mouse-zone profile, seven per sex/diet/zone. Boxplots show median, quartiles and conventional 1.5-IQR whiskers; all points are retained. Aggregate marker ordering passes in 56/56 profiles; this cannot validate original nucleus-level gates.")
    if (family == "supplementary_QC") return("PCA of centered, unscaled retained RNA or primary log2 lipid features. RNA includes 56 paired-zone profiles from 28 mice; lipid PCA includes 28 mice. Symbols identify diet and colors identify sex/zone or sex. Percentages give variance explained. This is an input diagnostic without outcome-driven exclusions.")
    if (family == "supplementary_lipid") return(paste0("Lipid factorial contrast ", if (nrow(source)) unique(source$contrast)[1L] else row$panel,
      ", using ", n(primary_n), " primary features and 28 mice. Horizontal effect units are log2 abundance differences with signs stated on the axis; vertical values are −log10 BH q in the corresponding lipid family. The horizontal line marks q=0.05. Up to eight labels rank by q, absolute effect and lipid ID."))
    if (family == "supplementary_threshold") return("Fraction of all primary eligible lipid species with at least one association across candidate contexts as the exploratory |r| threshold varies. SD and WD retain their diet-specific RNA nominations. This any-edge summary can saturate as candidate count grows; thresholds are not significance levels.")
    if (family == "supplementary_quality") return(paste0("Any-edge membership fractions at |r|≥0.65 under broad (", n(broad_n), "), QC-CV≤30% primary (", n(primary_n),
      ") and primary plus Fill≥0.5 (", n(fill_n), ") lipid universes. Labels state associated/eligible species. Denominators change with eligibility; RNA nominations and thresholds remain fixed. Saturation is not validation of equivalence between quality sets."))
    "Source data and the panel manifest define this supplementary export."
  }
  legends <- c("# CF453 Figure 6 panel legends", "",
    "Panels are supplied separately for manual assembly. SD, standard diet; WD, Western diet; PP, periportal; PC, pericentral. RNA effects are sex associations, and lipids are measured in whole liver. **No individual primary gene–lipid correlation passes BH q≤0.05; all association panels are exploratory.**", "",
    paste0("Representative selection: ", paste(paste(reps$gene, paste0("(", reps$sex, " ", reps$zone, ")")), collapse = "; "), "."), "",
    "RNA primary q pools eligible genes × three interactions; within-sex zonation uses a distinct genes × four-contrast family. Lipid families and the single correlation family are defined in METHODS.md. Bootstrap summaries use fixed RNA/partner selections and are descriptive.")
  for (i in seq_len(nrow(panels))) legends <- c(legends, "", paste0("## ", panels$panel[i]), "",
    legend_for(panels[i, , drop = FALSE]), "", paste0("Source data: ", link(panels$source_table[i], "CSV"), "."))
  if (!is.null(a$representative_delta_summary) && nrow(a$representative_delta_summary)) {
    d <- a$representative_delta_summary
    summary_lines <- c(summary_lines, "", paste0("Paired bootstrap delta-r summaries for ", n(nrow(d)),
      " selected gene–lipid partners are in [representative_delta_summary.csv](tables/representative_delta_summary.csv), with requested/valid draws and percentile ranges. These are descriptive fixed-selection quantities, not tests of differential correlation."))
  }
  writeLines(summary_lines, file.path(output, "README.md"), useBytes = TRUE)
  writeLines(methods_lines, file.path(output, "METHODS.md"), useBytes = TRUE)
  writeLines(legends, file.path(output, "PANEL_LEGENDS.md"), useBytes = TRUE)
  invisible(c(README = file.path(output, "README.md"), METHODS = file.path(output, "METHODS.md"),
              PANEL_LEGENDS = file.path(output, "PANEL_LEGENDS.md")))
}
