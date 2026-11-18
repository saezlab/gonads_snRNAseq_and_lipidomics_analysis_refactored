#!/usr/bin/env Rscript
# Read-only expression/feasibility evidence for the TF-to-lipid KO proposal.
# Counts, Factor 1 loadings and inferred activities are separate measurements.
options(stringsAsFactors = FALSE, digits = 17)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (is.na(script)) stop("Run this entry point with Rscript")
root <- normalizePath(file.path(dirname(script), "../../.."), mustWork = TRUE)
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R [--output=NEW_DIR] [--validate-only]\n",
      "Audits actual FCG pseudobulk counts and four historical neg2 cultured-hepatocyte RNA controls.\n",
      "Exports gene, celltype/sex and culture-library evidence plus source hashes. No MOON fitting.\n",
      "Default output is a fresh directory under R's temporary directory. Existing nonempty output is refused.\n",
      "--validate-only performs all input and count checks without creating output. Run from any directory.\n")
  quit(save = "no")
}
if (any(!grepl("^--output=.+$|^--validate-only$", args))) stop("Unknown or empty option; use --help")
if (anyDuplicated(sub("=.*$", "", args))) stop("Repeated option")
output_arg <- args[startsWith(args, "--output=")]
output <- if (length(output_arg)) sub("^--output=", "", output_arg) else tempfile("factor1_ko_expression_")
validate_only <- "--validate-only" %in% args
source("scripts/support_functions.R")
for (pkg in c("Matrix", "SeuratObject", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
}
check <- function(x, message) if (!isTRUE(x)) stop(message)
paths <- c(
  FCG_counts = "data/RNA/fcg_pseudobulks.rds",
  culture_counts = "miscellaneous/older_inputs/legacy_rna/merged_counts_matrix_250508_VH00211_433_AACVFM3HV_with_gene_symbols_and_groups.csv",
  selected_network = "results/causal_modelling/factor1_tf_to_lipid_expanded_species/currency_v2_expanded_layers4_lower_thresholds/network_ATT.csv",
  Factor1_inputs = "miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds"
)
check(all(file.exists(paths)), "Required input missing; see the paths declared in this script")
obj <- readRDS(paths[["FCG_counts"]])
assay <- obj@assays$RNA
check("counts" %in% names(assay@layers), "FCG RNA counts layer missing")
counts <- assay@layers$counts
features <- rownames(assay@features)[as.matrix(assay@features)[, "counts"]]
samples <- rownames(assay@cells)[as.matrix(assay@cells)[, "counts"]]
check(inherits(counts, "sparseMatrix"), "Expected actual sparse FCG counts, not scaled RNA")
check(identical(dim(counts), c(length(features), length(samples))), "Counts dimensions disagree with layer maps")
check(!anyDuplicated(features) && !anyDuplicated(samples), "Duplicate FCG feature/sample IDs")
check(all(is.finite(counts@x)) && all(counts@x >= 0) && all(counts@x == round(counts@x)), "FCG counts are not finite nonnegative integers")
dimnames(counts) <- list(features, samples)
md <- obj@meta.data[samples, , drop = FALSE]
check(identical(rownames(md), samples), "FCG metadata/sample mismatch")
check(all(c("sample_name", "cell_type", "gonadal_type") %in% names(md)), "FCG metadata columns missing")
hep <- which(md$cell_type == "Hepatocytes")
check(length(hep) == 43L && !anyDuplicated(md$sample_name[hep]), "Expected 43 unique FCG hepatocyte samples, without the five WT additions")
check(sum(md$gonadal_type[hep] == "Ovaries") == 21L && sum(md$gonadal_type[hep] == "Testes") == 22L, "FCG gonadal group counts changed")
libraries <- Matrix::colSums(counts)
check(all(is.finite(libraries)) && all(libraries > 0), "Invalid FCG pseudobulk library size")
cat("FCG actual counts:", class(counts)[1], nrow(counts), "genes x", ncol(counts), "pseudobulks; 43 hepatocyte mice\n")

att <- utils::read.csv(paths[["selected_network"]], check.names = FALSE)
primary <- sub("^Gene__", "", att$source[startsWith(att$source, "Gene__")])
check(length(primary) == 21L && !anyDuplicated(primary), "Expected 21 distinct primary-network genes")
extra <- c("Pld2", "Hoxa10", "Esr1", "Nfic", "Dbp", "Cyp2b9", "Cyp2b13", "Cyp2c38", "Cyp2d9", "Sult3a1", "Mtmr11", "Abcb1a", "Alb", "Hnf4a", "Pparg", "Cebpb", "Vdr")
genes <- unique(c(primary, extra))
inp <- readRDS(paths[["Factor1_inputs"]])
rows <- list()
for (gene in genes) {
  for (cell_type in unique(md$cell_type)) {
    j <- which(md$cell_type == cell_type)
    present <- gene %in% features
    vals <- if (present) as.numeric(counts[gene, j]) else rep(NA_real_, length(j))
    cpm <- vals / libraries[j] * 1e6
    for (gonads in c("all", "Ovaries", "Testes")) {
      k <- if (gonads == "all") seq_along(j) else which(md$gonadal_type[j] == gonads)
      rows[[length(rows) + 1L]] <- data.frame(gene = gene, cell_type = cell_type, gonads = gonads,
        n = length(k), source_present = present, mean_count = mean(vals[k]), median_count = stats::median(vals[k]),
        min_count = min(vals[k]), max_count = max(vals[k]), n_nonzero = sum(vals[k] > 0),
        mean_CPM = mean(cpm[k]), median_CPM = stats::median(cpm[k]), min_CPM = min(cpm[k]), max_CPM = max(cpm[k]),
        factor1_RNA_loading = unname(inp$hepatocyte_gene_weight[paste0("Gene__", gene)]),
        factor1_TF_activity = unname(inp$tf_activity[paste0("Gene__", gene)]),
        selected_inferred_score = att$score[match(paste0("Gene__", gene), att$source)])
    }
  }
}
celltype <- do.call(rbind, rows)

# This source has two metadata rows, not a normal one-row CSV header.
raw <- utils::read.csv(paths[["culture_counts"]], header = FALSE, check.names = FALSE,
  colClasses = "character", stringsAsFactors = FALSE)
groups <- unlist(raw[1L, -1L], use.names = FALSE)
culture_samples <- unlist(raw[2L, -1L], use.names = FALSE)
culture_genes <- raw[-c(1L, 2L), 1L]
culture_counts <- matrix(as.numeric(unlist(raw[-c(1L, 2L), -1L], use.names = FALSE)),
  nrow = length(culture_genes), ncol = length(culture_samples))
check(!anyNA(culture_counts) && all(is.finite(culture_counts)) && all(culture_counts >= 0) &&
  all(culture_counts == round(culture_counts)), "Culture counts are not finite nonnegative integers")
check(!anyDuplicated(culture_samples) && all(nzchar(culture_samples)), "Duplicate/empty culture library IDs")
ctrl <- which(groups == "neg2")
check(length(ctrl) == 4L, "Expected four historical neg2 control RNA libraries")
culture_libs <- colSums(culture_counts)
check(all(culture_libs > 0), "Empty historical RNA library")
cat("Historical culture raw counts:", nrow(culture_counts), "gene rows x", ncol(culture_counts),
    "libraries; four neg2 controls\n")
culture_rows <- list()
summary_rows <- list()
for (gene in genes) {
  found <- which(culture_genes == gene)
  vals <- if (length(found)) colSums(culture_counts[found, ctrl, drop = FALSE]) else rep(NA_real_, length(ctrl))
  cpm <- vals / culture_libs[ctrl] * 1e6
  for (k in seq_along(ctrl)) {
    culture_rows[[length(culture_rows) + 1L]] <- data.frame(gene = gene, sample_id = culture_samples[ctrl[k]],
      group = "neg2", source_gene_rows = length(found), raw_count = vals[k], library_size = culture_libs[ctrl[k]], CPM = cpm[k])
  }
  h <- celltype[celltype$gene == gene & celltype$cell_type == "Hepatocytes" & celltype$gonads == "all", ]
  other <- celltype[celltype$gene == gene & celltype$gonads == "all", ]
  strongest <- other[which.max(other$mean_CPM), ]
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(gene = gene, in_selected_primary_network = gene %in% primary,
    hepatocyte_pseudobulk_mean_CPM = h$mean_CPM, hepatocyte_nonzero_samples = h$n_nonzero,
    hepatocyte_n_samples = h$n, historical_culture_control_mean_CPM = mean(cpm),
    historical_culture_control_min_CPM = min(cpm), historical_culture_control_max_CPM = max(cpm),
    historical_culture_control_positive_samples = sum(vals > 0), historical_culture_control_n = length(ctrl),
    historical_culture_source_gene_rows = length(found), historical_culture_control_counts = paste(vals, collapse = "|"),
    historical_culture_control_CPM = paste(format(cpm, digits = 17, trim = TRUE), collapse = "|"),
    highest_celltype = strongest$cell_type, highest_celltype_mean_CPM = strongest$mean_CPM,
    factor1_RNA_loading = h$factor1_RNA_loading, factor1_TF_activity = h$factor1_TF_activity,
    selected_inferred_score = h$selected_inferred_score)
}
summary <- do.call(rbind, summary_rows)
culture <- do.call(rbind, culture_rows)
check(nrow(summary) == length(genes) && !anyDuplicated(summary$gene), "Gene-summary duplication")
check(all(summary$historical_culture_source_gene_rows == 1L), "Audit gene missing/duplicated in historical RNA source")
print(summary[, c("gene", "in_selected_primary_network", "hepatocyte_pseudobulk_mean_CPM", "historical_culture_control_mean_CPM")], row.names = FALSE)
source_paths <- c(paths, entry_point = normalizePath(script), support_functions = "scripts/config/support_functions.R")
manifest <- data.frame(role = names(source_paths), path = unname(source_paths),
  sha256 = vapply(source_paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)))
limitations <- c(
  "FCG samples are 43 unique mice (21 ovarian/22 testicular), excluding the separate five WT additions. Counts are finite nonnegative integers.",
  "FCG counts are per-mouse/per-celltype pseudobulks; nonzero across mice is not single-nucleus prevalence and does not exclude ambient or contaminating RNA.",
  "CPM uses total supplied counts per library. FCG and cultured bulk RNA are distinct assays/states; their CPM ratio is not a measured culture-induced fold change.",
  "Four historical neg2 RNA libraries establish detection in a prior culture system, not biological independence, protein activity, functional isoform, future expression or knockdown efficiency.",
  "Factor 1 RNA loadings, footprint TF activities and MOON inferred scores are separate fields; none is baseline abundance.",
  "No RNA expression filter was applied during MOON fitting; frozen RNA TF-target sign coherence does not establish mediator expression or cell autonomy.",
  "Pld2 and Pparg/Cebpb/Vdr context genes outside the selected graph do not acquire model-supported paths merely by inclusion in this expression audit.",
  "The upcoming experimental system is unspecified. Historical culture sex/timing should not be imposed on the new experiment.")
if (validate_only) {
  cat("Validation passed; no output written.\n")
  quit(save = "no")
}
if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) stop("Nonempty output directory: ", output)
if (file.exists(output) && !dir.exists(output)) stop("Output exists as a file: ", output)
if (!dir.exists(output)) check(dir.create(output, recursive = TRUE), "Cannot create output directory")
utils::write.csv(summary, file.path(output, "gene_expression.csv"), row.names = FALSE, na = "")
utils::write.csv(celltype, file.path(output, "gene_expression_by_celltype.csv"), row.names = FALSE, na = "")
utils::write.csv(culture, file.path(output, "culture_control_expression.csv"), row.names = FALSE, na = "")
utils::write.csv(manifest, file.path(output, "source_manifest.csv"), row.names = FALSE)
report <- list(schema_version = "1.0", sources = manifest, selected_primary_gene_count = length(primary),
  gene_summary = summary, FCG_celltype_sex_summaries = celltype,
  culture_control_samples = culture_samples[ctrl], culture_library_sizes = culture_libs[ctrl],
  integer_count_checks_passed = TRUE, FCG_hepatocyte_mice = length(hep), FCG_gonadal_counts = c(Ovaries = 21L, Testes = 22L),
  limitations = limitations)
jsonlite::write_json(report, file.path(output, "expression_audit.json"), pretty = TRUE, auto_unbox = TRUE,
  na = "null", digits = 17)
writeLines(c("# Expression evidence for the Factor 1 TF-to-lipid KO proposal", "",
  "Reproduce from the repository root:", "",
  "```sh", "Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R --validate-only",
  "Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R --output=NEW_DIR", "```", "",
  "`gene_expression.csv` keeps baseline count evidence distinct from RNA loadings and inferred activity. Full celltype/sex and control-library tables, source hashes and session details accompany it. No new fit or perturbation prediction is calculated.",
  "", paste0("- ", limitations)), file.path(output, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
cat("Audit written to", normalizePath(output), "\n")
