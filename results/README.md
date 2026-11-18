# Main analysis results

**Current zonation and standalone diet release:** [expanded markers, 16 November 2026](zonation_marker_update_20261116/README.md).

- [Current panel and source-data navigation](current_panels.md)
- [Current numerical replacement ledger](zonation_marker_update_20261116/report/new_replacements.csv)
- [Production report](zonation_marker_update_20261116/report/zonation_update_report.pdf)
- [Release registry](../scripts/config/analysis_release.json)

| Root | Role |
| --- | --- |
| `qc/` | Main input preparation and QC |
| `pseudobulk/` | Main cell-type RNA models, pathways, TFs and validations |
| `lipidomic/` | Lipid models, enrichment, structural and variance analyses |
| `mofa/` | Current trained model set, weights and global networks |
| `zonation_marker_update_20261116/` | **Current** FCG/GDX zonation, Figure5 and standalone Figure6 models/panels/source tables |
| `validation/` | Unchanged experimental KD/global validation and Cyp2c38 in-vivo analysis |
| `xenium/`, `human/`, `in_vitro_validation/` | Independent validation/collaborator outputs |
| `causal_modelling/` | Active optional causal resources, selected models and sensitivities; see causal map |

Earlier `zonation/`, previous CF453 Figure6 roots and cross-cohort comparison roots have moved into [miscellaneous](../miscellaneous/README.md). Main scripts resolve frozen references explicitly; they no longer infer that an old generic `results/zonation` folder is the current default.

The approved release is intact, including its provenance, comparisons and audit records. Git includes curated outputs; large reproducible caches and duplicate runtime exports remain local. See [the commit guide](../docs/GITHUB_COMMIT_GUIDE.md) before staging.
