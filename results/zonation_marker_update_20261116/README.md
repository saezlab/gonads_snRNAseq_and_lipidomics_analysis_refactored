# Current expanded-marker zonation and diet release

Approved 16 November 2026; this is the default scientific release selected by `scripts/config/analysis_release.json`. Original production files, numerical values and provenance are retained. This index was added during repository organization; it does not replace the original production manifests.

| Folder | Main content |
| --- | --- |
| `fcg05/` | FCG zonation models, gene selection, lipid associations and Figure4B/5/S5A sources |
| `fcg07/` | Updated zonation-network KD validation; prespecified male-PP Mtmr11 Figure5F |
| `gdx_zonation/` | Paired gonadectomy supporting analysis; separate from the legacy hormone-replacement notebook |
| `cf453_models/` | Standalone diet-study upstream model cache |
| `cf453_direct_diet/` | Direct diet interaction models |
| `cf453_rewiring/` | Current Figure6A–C and source data |
| `cf453_layouts/` | Both alternative Figure6A layouts |
| `tables/` | Updated source Tables S10/S13/S14/S15 |
| `report/` | Completed panel trace, numerical replacement ledger and report |
| `inputs/`, `runtime/`, `QA/`, `comparisons/`, `logs/` | Original input/source snapshots, numerical evidence and execution history |

Start with [the panel/source navigation](../current_panels.md) and [the replacement ledger](report/new_replacements.csv). The [report](report/zonation_update_report.pdf) documents selection-threshold effects, denominator/sign/sample-size corrections and exploratory diet associations.

Figure4A and S4B–C and their nucleus-level statistics remain deferred to James. Public inputs are split pseudobulk CSVs under `data/RNA/pseudobulk_csv/`, each under 50 MB; the original Seurat objects remain local-only; the private manuscript is optional and excluded from scientific default runs.

Use a fresh output directory for reruns. Do not overwrite this approved release or silently rebase older cross-cohort analyses to its new compartment definition.
