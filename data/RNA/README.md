# RNA input files

Only aggregated pseudobulk CSV inputs are distributed in the Git tree. Every RNA CSV must be strictly smaller than 50,000,000 bytes. RNA/Seurat `.rds`, `.RData` and other single-cell object formats remain local and are excluded from the commit.

- `scaled_RNA_{1,2,3}.csv`, metadata: existing main-study pseudobulk interfaces.
- `scaled_RNA_gdxhr_{1,2,3}.csv`, metadata: existing whole-cell-type validation interfaces.
- `pseudobulk_csv/`: seven manifest-verified exports of supplied counts, normalized/scaled values and aggregate-profile metadata, including the current expanded-marker inputs.
- `zonation/current/`: local original expanded-marker objects, public scoring provenance and small compatibility CSVs.
- `wd/`, `wdsd/`: local original comparison-cohort objects; corresponding public exports are in `pseudobulk_csv/`.

CSV columns represent mouse-by-zone or mouse-by-cell-type profiles, not individual nuclei. Export does not reassign zones or renormalize expression. See `docs/PSEUDOBULK_CSV_INPUTS.md` for ordering, hashes, original-object provenance and the distinction between current and historical compartments.
