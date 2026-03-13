# Export ttop_full_list.RData to CSV for MCP server use
# Run from project root: Rscript mcp/export_ttop_csv.R

load("results/pseudobulk/ttop_full_list.RData")

cat("Contrasts:", names(ttop_full_list), "\n")
cat("Cell types:", names(ttop_full_list[[1]]), "\n")

# Build a single long-format CSV with all contrasts and cell types
all_rows <- list()
for (contrast in names(ttop_full_list)) {
  for (cell_type in names(ttop_full_list[[contrast]])) {
    tt <- ttop_full_list[[contrast]][[cell_type]]
    if (is.null(tt) || nrow(tt) == 0) next

    # Add gene column from rownames if not present
    if (!"gene" %in% colnames(tt) && !"ID" %in% colnames(tt)) {
      tt$gene <- rownames(tt)
    }
    if ("ID" %in% colnames(tt)) {
      colnames(tt)[colnames(tt) == "ID"] <- "gene"
    }

    tt$contrast <- contrast
    tt$cell_type <- cell_type

    # Keep only standard columns
    std_cols <- c("gene", "contrast", "cell_type", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
    std_cols <- std_cols[std_cols %in% colnames(tt)]
    tt <- tt[, std_cols, drop = FALSE]

    all_rows[[length(all_rows) + 1]] <- tt
  }
}

combined <- do.call(rbind, all_rows)

# Standardize gene column name
if ("ID" %in% colnames(combined) && !"gene" %in% colnames(combined)) {
  colnames(combined)[colnames(combined) == "ID"] <- "gene"
}

# Keep only key columns
keep_cols <- c("gene", "contrast", "cell_type", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
keep_cols <- keep_cols[keep_cols %in% colnames(combined)]
combined <- combined[, keep_cols]

cat("Total rows:", nrow(combined), "\n")
cat("Columns:", paste(colnames(combined), collapse=", "), "\n")
cat("Unique contrasts:", length(unique(combined$contrast)), "\n")
cat("Unique cell types:", length(unique(combined$cell_type)), "\n")
cat("Unique genes:", length(unique(combined$gene)), "\n")

# Save
outpath <- "results/pseudobulk/ttop_full_long.csv"
write.csv(combined, outpath, row.names = FALSE)
cat("Saved to:", outpath, "\n")
cat("File size:", file.size(outpath) / 1e6, "MB\n")
