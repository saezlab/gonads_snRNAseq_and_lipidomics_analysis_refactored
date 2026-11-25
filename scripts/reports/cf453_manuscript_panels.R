# Render the current manuscript Figure 6B and 6E from saved CF453 tables.
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/config/analysis_release.R")
source("scripts/config/plotting_themes.R")
source("scripts/core/cf453_diet_impact_plots.R")
source("scripts/core/cf453_diet_rewiring_plots.R")
source("scripts/core/cf453_manuscript_panels.R")

cf453_manuscript_export <- function(panel, output) {
  source_table <- file.path("source_data", paste0(panel$stem, "_source.csv"))
  pdf <- file.path("panels", paste0(panel$stem, "_A4.pdf"))
  png <- file.path("panels", paste0(panel$stem, "_A4.png"))
  utils::write.csv(panel$source, file.path(output, source_table), row.names = FALSE)
  ggplot2::ggsave(
    file.path(output, pdf), panel$plot,
    width = panel$width_mm, height = panel$height_mm, units = "mm",
    device = grDevices::pdf, family = "Helvetica", encoding = "WinAnsi.enc",
    useDingbats = FALSE, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output, png), panel$plot,
    width = panel$width_mm, height = panel$height_mm, units = "mm",
    device = ragg::agg_png, dpi = panel$dpi, bg = "white"
  )
  data.frame(
    panel = sub("_.*$", "", panel$stem),
    file = c(pdf, png), format = c("PDF", "PNG"),
    width_mm = panel$width_mm, height_mm = panel$height_mm,
    source_table = source_table
  )
}

cf453_manuscript_panels_run <- function(args = commandArgs(TRUE)) {
  keys <- sub("=.*$", "", args)
  if (any(!grepl("^--(output|lipids|edges)=.+$", args)) || anyDuplicated(keys) ||
      !"--output" %in% keys || sum(c("--lipids", "--edges") %in% keys) == 1L) {
    stop("Supply --output=FRESH_DIR and optionally both --lipids=WF15_CSV and --edges=WF16_CSV.",
         call. = FALSE)
  }
  value <- function(key, default = NULL) {
    match <- args[startsWith(args, paste0(key, "="))]
    if (length(match)) sub("^[^=]*=", "", match) else default
  }
  output <- analysis_assert_not_approved_output(value("--output"))
  if (file.exists(output) || dir.exists(output)) stop("Choose a fresh output directory.")
  lipid_path <- value("--lipids", analysis_current("cf453_direct_diet/tables/lipid_results.csv"))
  edge_path <- value("--edges", analysis_current("cf453_rewiring/source_data/Fig6C_network_gain_loss.csv"))
  input_paths <- c(lipid_path, edge_path)
  lipids <- cf453_manuscript_lipid_contrasts(utils::read.csv(lipid_path))
  b <- cf453_manuscript_figure6b(lipids)
  e <- cf453_manuscript_figure6e(lipids, utils::read.csv(edge_path))

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  staged <- tempfile("cf453_manuscript_", tmpdir = dirname(output))
  dir.create(staged)
  on.exit(unlink(staged, recursive = TRUE), add = TRUE)
  dir.create(file.path(staged, "panels"))
  dir.create(file.path(staged, "source_data"))
  manifest <- rbind(cf453_manuscript_export(b, staged), cf453_manuscript_export(e, staged))
  utils::write.csv(manifest, file.path(staged, "panel_manifest.csv"), row.names = FALSE)
  inputs <- data.frame(
    role = c("WF15 lipid contrasts", "WF16 selected gene-lipid pairs"),
    file = basename(input_paths), bytes = file.info(input_paths)$size,
    sha256 = unname(vapply(input_paths, digest::digest, character(1L), file = TRUE, algo = "sha256"))
  )
  utils::write.csv(inputs, file.path(staged, "input_sha256.csv"), row.names = FALSE)
  writeLines(c(
    "# CF453 manuscript Figure 6B and 6E", "",
    "- [Figure 6B vector panel](panels/Fig6B_diet_specificity_A4.pdf)",
    "- [Figure 6E vector panel](panels/Fig6E_four_examples_A4.pdf)", "",
    "**B. Sex-polarised whole-liver lipid responses to Western diet.** Each dot shows the saved limma moderated t statistics for WD minus SD in females (x) and males (y). Shading marks absolute t below 1 in the opposing sex. Coloured points also have absolute diet-by-sex t above 4.4. The sign of the within-sex t distinguishes higher and lower abundance under WD. The fixed gates carry over the Figure 5A display rule. The axes use asinh(t/2) with ticks labelled in t units. These are descriptive selections.", "",
    "**E. Gene-lipid correlations by diet.** Each selected pair is placed at the within-sex whole-liver lipid log2(WD/SD) estimate. Arrows connect the saved Pearson correlations between zonal gene expression and whole-liver lipid abundance in SD (open) and WD (filled). The independent diet groups contain seven mice per sex. The selected pairs illustrate exploratory association changes.", "",
    "[Panel manifest](panel_manifest.csv) records the A4 placement dimensions and the plotted CSV for each export. [Input hashes](input_sha256.csv) identify the two saved tables read by the renderer."
  ), file.path(staged, "README.md"))
  if (!file.rename(staged, output)) stop("Could not publish the completed render to ", output)
  invisible(manifest)
}

if (sys.nframe() == 0L) cf453_manuscript_panels_run()
