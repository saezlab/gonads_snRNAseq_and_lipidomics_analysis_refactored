# Reporting for selected CF453 programmes; no analysis runs when sourced.
cfdr_report <- function(a, manifest, output) {
  required <- c("anchors", "programmes", "class_summary", "edges", "network_genes",
                "network_edges", "metadata", "settings")
  if (length(setdiff(required, names(a)))) stop("Incomplete selected-programme report inputs.")
  stopifnot(all(c("gene", "anchor_diet", "context", "threeway_q") %in% names(a$anchors)),
            all(c("panel", "file") %in% names(manifest)))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  n <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  fixed <- function(x, digits = 2) formatC(x, digits = digits, format = "f")
  md_table <- function(d) {
    if (!nrow(d)) return("No eligible rows.")
    d[] <- lapply(d, function(x) gsub("[|]", "/", as.character(x)))
    c(paste0("| ", paste(names(d), collapse = " | "), " |"),
      paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
      apply(d, 1L, function(x) paste0("| ", paste(x, collapse = " | "), " |")))
  }
  relative <- function(path) {
    prefix <- paste0(output, "/")
    ifelse(startsWith(path, prefix), substring(path, nchar(prefix) + 1L), path)
  }
  anchors <- a$anchors
  programmes <- a$programmes
  classes <- a$class_summary
  anchor_counts <- do.call(rbind, lapply(unique(anchors$anchor_diet), function(diet) {
    d <- anchors[anchors$anchor_diet == diet, , drop = FALSE]
    data.frame(selection_diet = diet, arrows = n(nrow(d)),
      direct_threeway_q_005 = n(sum(d$threeway_q <= .05, na.rm = TRUE)))
  }))
  context_counts <- if (nrow(programmes)) stats::aggregate(gene ~ context, programmes,
    function(x) length(unique(x))) else data.frame(context = character(), gene = integer())
  names(context_counts)[2L] <- "fixed_union_genes"
  networks <- a$network_genes[, intersect(c("context", "gene", "n_edges", "status"),
                                         names(a$network_genes)), drop = FALSE]
  lipid_count <- if ("lipid" %in% names(a$edges)) length(unique(a$edges$lipid)) else NA_integer_
  q_columns <- intersect(c("q_SD", "q_WD"), names(a$edges))
  q_values <- unlist(a$edges[q_columns], use.names = FALSE)
  q_summary <- if (length(q_columns)) paste0("The explicitly selected pair/context/diet family contains ",
    n(sum(is.finite(q_values))), " correlation tests; ", n(sum(q_values <= .05, na.rm = TRUE)),
    " pass its BH q≤0.05 threshold. This is conditional on the selected RNA programmes and is not a test of differential correlation.") else
    "Nominal and adjusted correlation evidence is exported for the explicitly selected programmes; it does not test differential correlation."
  finite_class <- is.finite(classes$delta_E)
  class_summary_text <- if (any(finite_class)) paste0("Across the class/context cells, ΔE ranges from ",
    fixed(100 * min(classes$delta_E[finite_class])), " to ",
    fixed(100 * max(classes$delta_E[finite_class])), " percentage points. ",
    "Defined bootstrap draws range from ", n(min(classes$bootstrap_valid)), " to ",
    n(max(classes$bootstrap_valid)), " of ", n(max(classes$bootstrap_requested)),
    " per context; intervals and directional fractions condition on those defined draws.") else
    "No finite class/context change is available; undefined summaries remain explicit."
  class_uncertainty_text <- character()
  if (any(finite_class)) {
    defined_bounds <- is.finite(classes$bootstrap_low) & is.finite(classes$bootstrap_high)
    exclude_zero <- defined_bounds & (classes$bootstrap_low > 0 | classes$bootstrap_high < 0)
    largest <- classes[which.max(abs(classes$delta_E)), , drop = FALSE]
    class_uncertainty_text <- paste0(n(sum(exclude_zero)), "/", n(sum(defined_bounds)),
      " conditional 95% percentile ranges exclude zero; ",
      n(sum(classes$bootstrap_direction >= .90, na.rm = TRUE)),
      " cells meet the descriptive ≥0.90 directional marker. The largest absolute displayed change is ",
      largest$lipid_class, " in ", gsub("_", " ", largest$context), " (", n(largest$n_lipids),
      " measured lipid", if (largest$n_lipids == 1L) "" else "s",
      "). Small class denominators can generate extreme rates; marker stability is not class-level significance.")
  }
  panels <- manifest[!duplicated(manifest$panel), , drop = FALSE]
  file_map <- data.frame(panel = panels$panel,
    PDF = vapply(panels$panel, function(panel) {
      f <- manifest$file[manifest$panel == panel & grepl("[.]pdf$", manifest$file, ignore.case = TRUE)]
      if (length(f)) paste0("[PDF](", relative(f[1L]), ")") else ""
    }, character(1L)))
  if ("source_table" %in% names(panels)) file_map$source <- if (nrow(panels))
    paste0("[CSV](", relative(panels$source_table), ")") else character()
  limits <- paste(
    "This is a standalone observational analysis of approximately 10-month-old CF453 wild-type mice, seven per sex/diet.",
    "The specificity gate and correlation thresholds define displayed programmes, not biological exclusivity, causal regulation or established rewiring.",
    "Whole-liver lipids are not spatial lipid measurements. Nucleus contributions, raw count aggregation, original nuclear QC/gating, batch/litter and assignment provenance remain unverified.",
    "Selection and bootstrap summaries use these same animals; fixed-selection stability is not independent validation or a calibrated post-selection significance test.")
  readme <- c("# Liver zonation and gene-lipid associations under standard and Western diets", "",
    paste0("The CF453 cohort contains ", n(nrow(a$metadata)), " approximately 10-month-old wild-type mice: seven females and seven males on each of standard diet (SD) and Western diet (WD). Each mouse contributes paired periportal/pericentral liver RNA and one whole-liver lipid measurement. Three panel groups connect selected zonation genes, relative lipid-class associations and exploratory correlation changes."), "",
    "## RNA anchors and direct diet evidence", "",
    paste0(n(nrow(anchors)), " anchor selections represent ", n(length(unique(anchors$gene))),
      " distinct genes. Selection uses anchor-diet RNA specificity and zonation magnitude; neither the SD–WD shift nor lipid correlations ranks the anchors."), "",
    md_table(anchor_counts), "",
    "Arrows track the same genes across diets even if their destination classification changes. Solid arrows have direct three-way diet×sex×zone q≤0.05; other arrows remain visible and dashed. Changes in t-statistic coordinates combine effect and precision, and anchor selection can produce regression toward weaker destination estimates.", "",
    "## Fixed selected programmes and class changes", "", md_table(context_counts), "",
    paste0("Both diets use the same per-context union of anchor genes", if (is.finite(lipid_count))
      paste0(" and ", n(lipid_count), " eligible lipid features") else "", ". ", class_summary_text), "",
    class_uncertainty_text, if (length(class_uncertainty_text)) "" else character(),
    "For class C, E is its fraction of gene–lipid pairs with |r|≥0.65 minus the corresponding fraction across all other lipid classes. Panel B shows ΔE=E_WD−E_SD. Positive ΔE means increased relative association-display density, not increased lipid abundance or spatial localization.", "",
    q_summary, "", "## Selected gain/loss examples", "", md_table(networks), "",
    "Network examples satisfy explicit strength, omission and bootstrap-stability gates. Gain/loss describes a threshold-defined change in observed |r|. The selected illustrations do not establish gene regulation or a significant rewiring event.", "",
    "## Panel files", "", md_table(file_map), "",
    "PDF/SVG/PNG versions and matching CSVs are listed in `panel_manifest.csv`. No extra or assembled figure is produced. [Methods](METHODS.md) and [legends](PANEL_LEGENDS.md) define all selection and uncertainty rules. Tables retain exclusions, exact denominators, correlation evidence, bootstrap valid counts and observed sensitivity summaries.", "",
    limits)
  methods <- c("# Methods: selected-programme diet analysis", "",
    "The validated eight-cell sex×diet×zone limma model uses supplied ln(1+10000×counts/profile-library) RNA, mouse duplicateCorrelation, mean-variance trend and robust empirical Bayes. RNA effects are natural-log normalized-expression differences; arrow axes are moderated t-statistics. Lipids are log2 supplied intensities, filtered to endogenous supported classes and analytical QC CV≤30%. No new normalization or other-cohort reference is used.", "",
    "Anchor eligibility requires anchor-diet sex×zone q≤0.05, own-sex PP−PC q≤0.05 and other-sex |t|<1. The last is a descriptive gate, not an equivalence test. Within each sex/zone/diet, select up to five by decreasing |own zonation effect|, then interaction q, own q and gene symbol. Positive PP−PC defines PP. Anchor q uses the earlier genes×three family (SD sex×zone, WD sex×zone, three-way); own zonation uses its separate four-test family. Arrow diet-modification q uses the genes×three direct diet-interaction family, independent of anchor ranking.", "",
    "For each sex/zone, the union of selected anchors is fixed across diets. Class and rest hit rates use all eligible gene×lipid pairs in the respective classes, including pairs below threshold. Rates are not fractions of independent animals. Nominal Pearson p uses five degrees of freedom; one BH family covers the selected context/diet/gene/lipid tests. Class changes are descriptive, rather than a genome-wide enrichment test.", "",
    "999 independent-arm bootstrap draws resample matched RNA/lipid mouse vectors within sex/diet, sharing indices across both zones and all features. For class summaries, any constant modeled feature invalidates its whole context/draw, preserving fixed denominators. Edge-specific strength-change summaries require that individual pair to be defined in both diet draws and may therefore have more valid draws. Percentile bounds and directional fractions condition on defined draws and fixed programme selection; they are not null-test p-values. Thresholds 0.60/0.65/0.70/0.80 and the broader all-strict-eligible RNA union are table-only sensitivities.", "",
    "A gained edge has |r_SD|≤0.4, |r_WD|≥0.8 and Δ|r|≥0.5; loss reverses the diets. All seven high-diet omissions must retain |r|≥0.65 and the observed sign; all seven low-diet omissions require |r|≤0.65. Minimum high minus maximum low strength equals the minimum over all 49 cross-arm omission combinations. At least 90% of bootstrap draws must be defined and at least 90% of defined draws support the observed strength-change direction. Genes rank by qualifying-edge count, maximum strength change, best anchor interaction q and symbol. Up to three gained and three lost partners are retained, filling up to six from remaining qualifiers. Unsupported contexts remain empty; partner selection is exploratory.")
  legends <- c("# Figure 6 legends", "",
    "SD, standard diet; WD, Western diet; PP, periportal; PC, pericentral. Seven independent mice per sex/diet contribute paired RNA zones and one whole-liver lipid vector. All association summaries concern selected RNA programmes.", "",
    "**A. Anchor-gene zonation shifts.** Female PP−PC t is x and male PP−PC t is y. Each facet selects five PP and five PC genes per sex in its anchor diet and tracks their other-diet coordinates. Arrows start in the selection diet; open points denote SD and filled points WD. Color identifies anchor sex/zone. Solid arrows have directly tested diet×sex×zone q≤0.05; dashed arrows do not. Arrow displacement alone is not a test of diet modification. Selection favors large anchor effects and can produce weaker destination estimates.", "",
    "**B. Relative lipid-class association change.** Fixed anchor-gene unions are used in each sex/zone and both diets. For each class, E is class pair hit rate minus all-other-class pair hit rate at |r|≥0.65; tiles show 100×(E_WD−E_SD), in percentage points. Denominators include every eligible pair in class/rest. Bootstrap directional markers, where present, indicate ≥0.90 defined-draw stability, not significance; invalid-draw counts and conditional percentile bounds are exported. This is neither lipid abundance enrichment nor spatial lipid localization.", "",
    "**C. Selected threshold-defined gain/loss networks.** Up to one gene per sex/zone is chosen from the same anchor union, using qualifying association changes. Each partner retains both SD and WD r values. Gain requires high WD |r|≥0.8, low SD |r|≤0.4 and strength change≥0.5; loss reverses these. Seven omissions per diet must pass the stated strength/sign gates; ≥90% of bootstrap draws must be defined and ≥90% of those support the observed change direction. Up to six partners are shown; empty contexts retain their unsupported status. These data-selected examples are exploratory and do not demonstrate causal regulation or calibrated differential correlation.")
  files <- c("README.md", "METHODS.md", "PANEL_LEGENDS.md")
  for (i in seq_along(files)) writeLines(list(readme, methods, legends)[[i]],
                                        file.path(output, files[i]), useBytes = TRUE)
  invisible(file.path(output, files))
}
