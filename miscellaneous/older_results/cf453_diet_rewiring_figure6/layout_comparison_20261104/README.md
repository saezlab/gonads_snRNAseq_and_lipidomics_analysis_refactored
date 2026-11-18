# Figure6A: two layout candidates

Both views use the same 34 unique genes from the original 40 SD/WD anchor selections, the same model estimates, and the same direct-test annotations.

- [Option1: four-state heatmap](panels/Fig6A_option1_heatmap.pdf)
- [Option2: paired dots](panels/Fig6A_option2_paired_dots.pdf)

Each gene is shown in female SD, female WD, male SD and male WD. SD means standard diet; WD means Western diet. PP means periportal and PC pericentral.

The plotted effect is modelled PP-minus-PC natural-log normalized RNA expression difference. Positive values indicate PP bias and negative values PC bias. This is not log2 fold change. Both plots use untransformed saved effects; the heatmap has a single symmetric scale without row standardization.

Selection badges SD, WD and Both retain the anchor-diet membership. Six genes were selected in both diets. Programme headings describe the selection programme, not a claim that the gene retains that classification under both diets.

The direct diet x sex x zone q<=0.05 mark identifies nine genes with evidence that diet changes zonation differently between sexes. The other selected genes remain visible. This annotation uses the existing Workflow15 primary FDR family, with no new tests or reselection.

Within each programme, genes are ordered by decreasing absolute direct three-way effect and then symbol. This display ordering is identical in both candidates and does not change gene eligibility. All 136 model estimates are retained in each source CSV.

Paired dots connect SD and WD group estimates for the same gene and sex, not repeated dietary measurements on the same animal. Female and male estimates have separate tracks. No confidence intervals are drawn in these overview candidates; statistical support is conveyed by the separately defined direct-test mark.

PDF/SVG/300dpiPNG alternatives and source CSVs are indexed by panel_manifest.csv. The original Figure6A/B/C, statistical results and main figure bundle are preserved. These alternatives are for layout selection before replacing a manuscript panel.
