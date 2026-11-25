source("scripts/reports/cf453_manuscript_panels.R")

test_cf453_manuscript_panels <- function() {
  lipid_path <- analysis_current("cf453_direct_diet/tables/lipid_results.csv", mustWork = TRUE)
  edge_path <- analysis_current("cf453_rewiring/source_data/Fig6C_network_gain_loss.csv", mustWork = TRUE)
  lipids <- cf453_manuscript_lipid_contrasts(utils::read.csv(lipid_path))
  edges <- utils::read.csv(edge_path)
  b <- cf453_manuscript_figure6b(lipids)
  e <- cf453_manuscript_figure6e(lipids, edges)

  stopifnot(
    identical(as.integer(table(b$source$category)), c(9L, 2L, 0L, 5L, 1113L)),
    nrow(e$source) == 24L, all(table(e$source$context) == 6L)
  )
  e_layers <- ggplot2::ggplot_build(e$plot)$data
  segments <- e_layers[[4L]]
  stopifnot(
    identical(segments$x, e$source$log2_WD_over_SD),
    identical(segments$xend, e$source$log2_WD_over_SD),
    identical(segments$y, e$source$r_SD),
    identical(sign(segments$yend - segments$y), sign(e$source$r_WD - e$source$r_SD))
  )
  incomplete <- cf453_manuscript_figure6e(lipids, subset(edges, context != "male_PC"))
  stopifnot(nrow(incomplete$source) == 18L,
            nrow(ggplot2::ggplot_build(incomplete$plot)$layout$layout) == 4L)
  empty <- cf453_manuscript_figure6e(lipids, edges[FALSE, ])
  stopifnot(nrow(empty$source) == 0L,
            nrow(ggplot2::ggplot_build(empty$plot)$layout$layout) == 4L)

  output <- tempfile("cf453_manuscript_test_")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  cf453_manuscript_panels_run(c(paste0("--output=", output),
    paste0("--lipids=", lipid_path), paste0("--edges=", edge_path)))
  release <- analysis_current("cf453_manuscript_panels", mustWork = TRUE)
  manifest <- utils::read.csv(file.path(output, "panel_manifest.csv"))
  sources <- unique(manifest$source_table)
  stopifnot(
    nrow(manifest) == 4L, all(file.exists(file.path(output, manifest$file))),
    identical(unname(tools::md5sum(file.path(output, sources))),
              unname(tools::md5sum(file.path(release, sources))))
  )
  cat("Figure 6B and 6E source replay, selection and arrow geometry passed.\n")
}

test_cf453_manuscript_panels()
