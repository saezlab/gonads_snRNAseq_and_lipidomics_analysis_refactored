# Gonadal hormone coordination of sex-specific adult liver physiology

Analysis code and results for the Four Core Genotypes mouse liver study, with independent validation and a standalone sex × diet × liver-zonation analysis.

## Current analysis

**The expanded-marker zonation analysis approved on 16 November 2026 is the default.** Its immutable production run is [results/zonation_marker_update_20261116](results/zonation_marker_update_20261116/README.md). Current FCG, gonadectomy and standard/Western-diet pseudobulk inputs are distributed under `data/RNA/pseudobulk_csv/` as numbered CSV parts, each strictly under 50 MB. Seurat/RNA RDS objects stay local and are excluded from Git. See [the CSV input guide](docs/PSEUDOBULK_CSV_INPUTS.md) and [the public input manifest](docs/provenance/current_input_manifest.csv).

Start with:

- [Pipeline and dependency map](PIPELINE.md): execution order, current outputs and optional branches.
- [Current results](results/README.md): figure-panel and table locations.
- [Result lookup guide](docs/AGENT_GUIDE.md): concise query map.
- [Provenance and archive map](docs/REPOSITORY_ORGANIZATION.md): original paths, archive locations, hashes and restoration.
- [GitHub commit guide](docs/GITHUB_COMMIT_GUIDE.md): what is included, what is external and the pre-commit check.

The executable release registry is [scripts/config/analysis_release.json](scripts/config/analysis_release.json). Current scientific outputs and historical reference inputs are different entries; they must not be substituted for one another.

## Reproduce the current zonation and diet analysis

Run commands from the repository root in an R environment with the required packages installed.

```sh
# Check current inputs, historical reference prerequisites and setup; no models run.
Rscript scripts/workflows/17_expanded_marker_zonation_update.R --validate-only

# Recompute into a fresh directory; the approved release is preserved.
Rscript scripts/workflows/17_expanded_marker_zonation_update.R --output=results/runs/expanded_markers_reproduction
```

Workflow 17 is the main zonation entry. It runs the FCG zonation models and gene–lipid networks, zonation-dependent knockdown validation, paired gonadectomy supporting analysis, diet-study models, selected-programme Figure 6, and both Figure 6A layouts. Workflow 14 and 15 model stages are upstream components of the current three-panel Workflow 16 presentation. Earlier full Figure 6 presentations are archived.

Scientific reruns do not require the private manuscript. Manuscript auditing/report regeneration is an explicit optional mode; the completed report and its replacement ledger are already provided with the approved release. Read the [Workflow 17 plan](docs/workflow_plans/17_expanded_marker_zonation_update.md) for CLI details and prerequisites.

## Main study pipeline

| Branch | Entry | Main results |
| --- | --- | --- |
| Input preparation and QC | `01_data_preparation_qc.R` | `data/`, `results/qc/` |
| Cell-type RNA differential analysis and validations | `02_fcg_pseudobulk_analysis.R` | `results/pseudobulk/` |
| Lipid differential analysis | `03_fcg_lipidomic_analysis.R` | `results/lipidomic/` |
| MOFA integration and global networks | `04_mofa_integration.R` | `results/mofa/` |
| **Current zonation / standalone diet figures** | **`17_expanded_marker_zonation_update.R`** | **`results/zonation_marker_update_20261116/`** |
| Lipid structural analysis | `06_structural_analysis.R` | `results/lipidomic/structural_model/` |
| Experimental validation | `07_validation_experiments.R --global-only`, current zonal branch through 17 | `results/validation/`, current release `fcg07/` |
| Lipid variance | `08_variance_analysis.R` | `results/lipidomic/variance_analysis/` |
| In-vivo Cyp2c38 validation | `11_cyp2c38_in_vivo_validation.R` | `results/validation/cyp2c38_in_vivo/` |

Entries are under `scripts/workflows/`. The [workflow plans](docs/workflow_plans/) describe inputs, contrasts, units, outputs and limitations. Raw sci-RNA-seq preprocessing (00), human comparison (09) and collaborator Xenium scripts (10) need external raw data and are optional. MOFA defaults use the retained trained models rather than training new models. The committed whole-cell CSVs support starting downstream analyses at 02; regenerating them with 01 also requires the original whole-cell FCG/GDX Seurat objects specified in its plan.

## Figures and interpretation

The current release supplies revised Figure 4B; Figure 5A–D/F; S5A; source Tables S10/S13/S14/S15; and standalone Figure 6A–C with both 6A alternatives. Panels are provided separately for manual assembly. Independent spatial assays and knockdown measurements are preserved.

Figure 4A and S4B–C, their linked nucleus-level counts/composition statistics, and final confirmation of executed retained-marker lists for Methods remain James's tasks. Pseudobulks cannot reconstruct these quantities.

The diet cohort is analysed independently of the younger original cohort. Gene–lipid associations and illustrative gain/loss networks remain exploratory: no tested correlation passes the declared BH family, and all 120 class bootstrap intervals include zero. See the current report for the exact numerical replacement ledger and selection-threshold sensitivity.

## Optional analyses and historical references

[COSMOS/MOON and PKN analyses](docs/CAUSAL_ANALYSIS_MAP.md) are retained as an active, separate extension, including their frozen reference dependencies. They are not quick tests. Their pinned model-fitting environment is separate from the R environment used for the zonation update.

Earlier cross-cohort analyses (12/13) retain their original two-marker compartment definition and frozen references. They are not rebased to the expanded-marker inputs. Their saved results, previous Figure 6 versions, old model generations, exploratory outputs and refactoring records are under [miscellaneous/](miscellaneous/README.md). Archived source and provenance remain discoverable; large historical payloads are excluded from ordinary Git commits and can be located by their recorded path and SHA256.

## GitHub preparation

```sh
python3 scripts/maintenance/audit_github_readiness.py --strict
```

The check inventories the files Git would include, oversized files, local metadata and data/provenance availability. It does not stage, commit or push. See [the commit guide](docs/GITHUB_COMMIT_GUIDE.md) for large-data restoration and review steps. This repository already contains substantial scientific data; generated caches, duplicated bundles, bulk old results, local settings and private manuscript files are excluded.

## Contributing and reuse

Read [AGENTS.md](AGENTS.md) before changing analysis code. Use explicit package namespaces, test modified R scripts immediately, and update the corresponding workflow plan. Preserve frozen scientific artifacts and their recorded hashes; use the release resolver when reading a relocated historical path.

Licensed under [LICENSE](LICENSE).
