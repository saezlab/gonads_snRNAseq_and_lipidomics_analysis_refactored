# Compact reporting for direct CF453 diet interactions; no analysis on source().
cfdi_report <- function(a, manifest, output) {
  required <- c("overview", "conditional", "omnibus", "heatmap_selection",
                "example_selection", "lipid_scatter", "metadata")
  if (length(setdiff(required, names(a)))) stop("Incomplete direct-diet report input.")
  stopifnot(all(c("gene", "contrast", "effect", "q") %in% names(a$overview)),
            all(c("panel", "file") %in% names(manifest)))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  number <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  fixed <- function(x, digits = 3) formatC(x, digits = digits, format = "f")
  qtext <- function(x) formatC(x, digits = 3, format = "g")
  markdown_table <- function(d) {
    if (!nrow(d)) return("No selected features.")
    d[] <- lapply(d, function(x) gsub("[|]", "/", as.character(x)))
    c(paste0("| ", paste(names(d), collapse = " | "), " |"),
      paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
      apply(d, 1L, function(x) paste0("| ", paste(x, collapse = " | "), " |")))
  }
  relative <- function(path) {
    prefix <- paste0(output, "/")
    ifelse(startsWith(path, prefix), substring(path, nchar(prefix) + 1L), path)
  }
  labels <- c(diet_sex = "Diet × sex, averaged over zones",
              diet_zone = "Diet × zone, averaged over sexes",
              diet_sex_zone = "Diet × sex × zone")
  counts <- do.call(rbind, lapply(names(labels), function(contrast) {
    d <- a$overview[a$overview$contrast == contrast, , drop = FALSE]
    data.frame(test = unname(labels[contrast]), genes_tested = number(nrow(d)),
      q_005 = number(sum(d$q <= .05, na.rm = TRUE)))
  }))
  gene_n <- length(unique(a$overview$gene))
  replay_text <- character()
  if (!is.null(a$parent_threeway_replay) &&
      all(a$parent_threeway_replay$max_abs_difference == 0)) {
    count_threeway <- sum(a$overview$contrast == "diet_sex_zone" & a$overview$q <= .05)
    replay_text <- paste0("Three-way estimates and unadjusted probabilities exactly replay the earlier paired model. The current ",
      number(count_threeway), " BH discoveries replace the earlier 194-count reporting family: this analysis pools three-way tests with diet×sex and diet×zone, whereas the earlier family pooled them with the separate SD/WD sex×zone tests. The change in adjusted count reflects the testing family, not changed raw three-way evidence.")
  }
  lipid <- a$lipid_scatter
  lipid_sig <- is.finite(lipid$q_interaction) & lipid$q_interaction <= .05
  examples <- a$example_selection
  selection_rule <- if (!is.null(a$settings$examples)) as.character(a$settings$examples) else
    "Distinct primary-significant genes, selected separately for diet×sex, diet×zone and two three-way illustrations"
  examples_text <- if (nrow(examples)) paste0(paste(paste0(examples$gene, " (",
    unname(labels[examples$selected_contrast]), "; q=", qtext(examples$q_primary), ")"), collapse = "; "), ".") else
    "No genes pass the prespecified example rule."
  example_table <- data.frame()
  if (nrow(examples)) {
    at <- match(paste(examples$gene, examples$selected_contrast),
                paste(a$overview$gene, a$overview$contrast))
    stopifnot(!anyNA(at))
    d <- a$overview[at, , drop = FALSE]
    example_table <- data.frame(gene = d$gene, contrast = d$contrast,
      effect = fixed(d$effect), nominal_95_interval = paste0("[", fixed(d$low), ", ", fixed(d$high), "]"),
      primary_q = qtext(d$q))
  }
  panels <- manifest[!duplicated(manifest$panel), , drop = FALSE]
  file_map <- data.frame(panel = panels$panel,
    PDF = vapply(panels$panel, function(panel) {
      f <- manifest$file[manifest$panel == panel & grepl("[.]pdf$", manifest$file, ignore.case = TRUE)]
      if (length(f)) paste0("[PDF](", relative(f[1L]), ")") else ""
    }, character(1L)))
  if ("source_table" %in% names(panels)) file_map$source <- if (nrow(panels))
    paste0("[CSV](", relative(panels$source_table), ")") else character()
  definitions <- c(
    "Let Δ(s,z)=mean(WD,s,z)−mean(SD,s,z), with z=PP or PC. The direct contrasts are:", "",
    "- Diet × sex: average over zones of Δ(male,z)−Δ(female,z).",
    "- Diet × zone: average over sexes of Δ(s,PP)−Δ(s,PC).",
    "- Three-way: [Δ(male,PP)−Δ(male,PC)]−[Δ(female,PP)−Δ(female,PC)].", "",
    "Positive values follow these signed differences; they need not imply a positive diet effect in every group. The interaction is the test of a difference in response; significance in only one diet or sex is not such a test.")
  limits <- paste(
    "These results describe diet-associated differences within CF453, with seven mice per sex/diet and paired RNA zones.",
    "They do not isolate hormonal or chromosome effects, and whole-liver lipids are not spatial lipid measurements.",
    "Source count aggregation, nucleus contributions/QC, batch/litter and assignment provenance remain unverified.",
    "Displayed genes are selected from the same data; nominal intervals do not account for that selection.")
  readme <- c("# CF453: direct diet-impact panels", "",
    paste0(number(nrow(a$metadata)), " mice (seven per female/male × SD/WD group), ",
      number(gene_n), " modeled genes and ", number(nrow(lipid)),
      " analytical-QC-filtered whole-liver lipid features. Four separate panels focus on direct diet interactions."), "",
    "The fitted model compares WD with SD and tests whether that difference depends on sex, zone, or both. This tests diet impact directly; comparing separate SD and WD significance labels would not.", "",
    "## Results", "", markdown_table(counts), "",
    paste0("The separate three-degree-of-freedom RNA omnibus test identifies ",
      number(sum(a$omnibus$q <= .05, na.rm = TRUE)), " genes at BH q≤0.05. Its ranking supplies the ",
      number(nrow(a$heatmap_selection)), " genes decomposed in panel B."), "",
    paste0("Diet × sex is supported for ", number(sum(lipid_sig)), "/", number(nrow(lipid)),
      " lipids (", fixed(100 * mean(lipid_sig), 2), "%); ",
      number(sum(lipid_sig & lipid$effect_interaction > 0)), " interactions are positive and ",
      number(sum(lipid_sig & lipid$effect_interaction < 0)),
      " negative. Positive means the male WD−SD response is more positive, or less negative, than the female response."), "",
    "Panel C examples and their direct interaction estimates (natural-log normalized-expression units; intervals are nominal 95%):", "",
    markdown_table(example_table), "",
    "## Direct comparisons", "", definitions, "",
    "## Files and interpretation", "", markdown_table(file_map), "",
    "PDF, SVG and PNG versions and matching source CSVs are listed in `panel_manifest.csv`. Assemble the final figure manually. [Methods](METHODS.md) and [legends](PANEL_LEGENDS.md) define the models, contrasts, uncertainty and selection.", "",
    limits)
  methods <- c("# Methods: direct CF453 diet effects", "",
    paste0("RNA uses ", number(gene_n),
      " genes detected in ≥7/56 zone profiles with nonzero expression variation. Supplied data numerically replays ln(1+10000×counts/profile-library) to <10^-15; this does not establish count-aggregation provenance. The eight sex/diet/zone means are fit using limma, mouse duplicateCorrelation, mean-variance trend and robust empirical Bayes. RNA effects are natural-log normalized-expression differences, not log2 fold changes."), "",
    paste0("BH adjusts the three direct RNA contrasts jointly across ", number(3 * gene_n),
      " tests. A separate genes × four-contrast family covers diet×sex within PP/PC and diet×zone within female/male groups. The omnibus F jointly tests the three direct coefficients and has its own across-gene BH adjustment; omnibus significance need not imply that an individual contrast passes its different testing family."), "",
    replay_text, if (length(replay_text)) "" else character(),
    paste0("Panel B ranks omnibus q, p and gene symbol, retaining at most 20 q≤0.05 genes. Panel C uses distinct primary-significant genes: one diet×sex, one diet×zone and two three-way examples. Its saved ranking rule is: ",
      selection_rule, ". Selection and full results are exported. Points are individual mice; fitted means and effect ranges use nominal 95% model intervals, not simultaneous or selection-adjusted intervals."), "",
    paste0("Lipid intensities are supplied positive values on a log2 scale, modeled by sex×diet limma. Primary eligibility excludes standards/isotopes and unsupported classes and requires analytical QC CV≤30% (", number(nrow(lipid)),
      " features). Diet×sex is (male WD−SD)−(female WD−SD), adjusted across all primary lipids. Displayed female/male simple effects retain their declared shared two-contrast BH family. No extra normalization or imputation is applied. This figure uses direct mean-effect tests, without correlation-network inference or network-selected examples."))
  legends <- c("# Figure 6: direct diet impact in CF453", "",
    "SD, standard diet; WD, Western diet; PP, periportal; PC, pericentral. There are seven independent mice per sex/diet; PP and PC measurements are paired within mouse. RNA units are differences in natural-log normalized expression; lipids use log2 abundance. Δ denotes WD−SD.", "",
    "**A. Direct RNA interactions.** Genome-wide estimates and adjusted evidence for diet×sex (averaged over zones), diet×zone (averaged over sexes) and diet×sex×zone. BH q uses one genes×three-test family. The contrasts follow Δmale−Δfemale, ΔPP−ΔPC and their difference-in-differences, respectively; the three-way contrast is [Δmale,PP−Δmale,PC]−[Δfemale,PP−Δfemale,PC].", "",
    paste0("**B. Interaction decomposition.** The top ", number(nrow(a$heatmap_selection)),
      " omnibus-selected genes are shown across four conditional interactions: diet×sex within PC/PP and diet×zone within female/male groups. Effect colors use natural-log normalized-expression differences. Conditional q-values use a separate genes×four-test BH family. The omnibus is a joint test; individual conditional tests may remain nonsignificant."), "",
    paste0("**C. Mouse-level examples.** ", examples_text,
      " Each point represents one of 28 mice. The diet×sex example averages that mouse's RNA over zones; diet×zone/three-way examples use its PP−PC expression difference. Group means and nominal 95% model intervals accompany the raw points. Reported q-values test the selected direct contrast, rather than a difference between separate significance decisions. Examples are data-selected."), "",
    "**D. Sex difference in lipid response.** Female WD−SD effect on the x-axis and male WD−SD effect on the y-axis for every primary lipid. The identity line represents equal estimated diet responses. Color marks the directly tested lipid diet×sex q≤0.05; selected labels identify the strongest interaction evidence. Above the line, the male response is more positive; below it, the female response is more positive. Neither axis measures zonal lipid abundance.")
  files <- c(README = "README.md", METHODS = "METHODS.md", PANEL_LEGENDS = "PANEL_LEGENDS.md")
  writeLines(readme, file.path(output, files["README"]), useBytes = TRUE)
  writeLines(methods, file.path(output, files["METHODS"]), useBytes = TRUE)
  writeLines(legends, file.path(output, files["PANEL_LEGENDS"]), useBytes = TRUE)
  invisible(file.path(output, unname(files)))
}
