# Gonadal hormone coordination of sex-specific adult liver physiology

Analysis code and results for the Four Core Genotypes mouse liver study, with independent validation and a standalone sex × diet × liver-zonation analysis.

## Current analysis

**The expanded-marker zonation analysis approved on 16 November 2026 is the default.** Its verified outputs are in the [manuscript-organized Results tree](results/README.md). Current FCG, gonadectomy and standard/Western-diet pseudobulk inputs are distributed under `data/RNA/pseudobulk_csv/` as numbered CSV parts, each strictly under 50 MB. Seurat/RNA RDS objects stay local and are excluded from Git. See [the CSV input guide](docs/PSEUDOBULK_CSV_INPUTS.md) and [the public input manifest](docs/provenance/current_input_manifest.csv).

Start with:

- [Pipeline and dependency map](PIPELINE.md): execution order, current outputs and optional branches.
- [Complete workflow reference](docs/WORKFLOW_REFERENCE.md): each workflow’s inputs, producers, modes, outputs and restoration requirements.
- [Current results](results/README.md): figure-panel and table locations.
- [Result lookup guide](docs/AGENT_GUIDE.md): concise query map.
- [Provenance and archive map](docs/REPOSITORY_ORGANIZATION.md): original paths, archive locations, hashes and restoration.
- [GitHub commit guide](docs/GITHUB_COMMIT_GUIDE.md): what is included, what is external and the pre-commit check.
- [Parallel manuscript results replay](docs/MANUSCRIPT_RESULTS_REPLAY.md): subsection layout, handoffs and comparison procedure.

The executable release registry is [scripts/config/analysis_release.json](scripts/config/analysis_release.json). Current scientific outputs and historical reference inputs are different entries; they must not be substituted for one another.

## Reproduce the current zonation and diet analysis

Run commands from the repository root in an R environment with the required packages installed.

```sh
# FCG zonation, gene–lipid networks and Figure 5 panels
Rscript scripts/workflows/05_zonation_analysis.R --output=results/runs/fcg_zonation
Rscript scripts/workflows/07_validation_experiments.R --output=results/runs/fcg_zonation

# Paired sham/gonadectomy zonation
Rscript scripts/workflows/05_zonation_analysis.R --cohort=gdx --output=results/runs/gdx_zonation
```

Workflow 05 reads the supplied expanded-marker profiles directly. Workflow 07 uses
the checkpoint and networks from the FCG run, together with the saved knockdown
effects. Use `--parent=PATH` to use an FCG checkpoint from a different output root.

For standalone diet analysis, Workflow 14 reads the expanded-marker diet profiles
and supplies the model cache to Workflows 15 and 16. The saved lipid contrasts and
representative pairs feed the Figure 6B and 6E renderer. The full command sequence
is in the [pipeline guide](PIPELINE.md#standalone-diet-execution).

Workflow 17 is an optional batch entry that runs these branches together with
release QA and Figure 6 rendering. Scientific reruns use the shipped derived
prerequisites. The optional manuscript audit uses additional local references.
See the [Workflow 17 plan](docs/workflow_plans/17_expanded_marker_zonation_update.md).

## Main study pipeline

| Branch | Entry | Main results |
| --- | --- | --- |
| Input preparation and QC | `01_data_preparation_qc.R` | `data/`, `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc/` |
| Cell-type RNA differential analysis and validations | `02_fcg_pseudobulk_analysis.R` | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/` |
| Lipid differential analysis | `03_fcg_lipidomic_analysis.R` | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/` |
| MOFA integration and global networks | `04_mofa_integration.R` | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/` |
| FCG and paired gonadectomy zonation | `05_zonation_analysis.R` (FCG default, `--cohort=gdx` for gonadectomy) | Sections 04–05 under `results/current/` |
| Lipid structural analysis | `06_structural_analysis.R` | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/` |
| Experimental validation | `07_validation_experiments.R` (zonal default, `--global-only` for the independent fit) | `results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/`, current release `fcg07/` |
| Lipid variance | `08_variance_analysis.R` | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/variance_analysis/` |
| In-vivo Cyp2c38 validation | `11_cyp2c38_in_vivo_validation.R` | `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/` |
| Standalone diet figures | Workflow 14 `--models-only` → 15 `--models-only` → 16 | Section 06 under `results/current/` |

Entries are under `scripts/workflows/`. The [workflow plans](docs/workflow_plans/) describe inputs, contrasts, units, outputs and limitations. Raw sci-RNA-seq preprocessing (00), human comparison (09) and collaborator Xenium scripts (10) need external inputs. Workflow04 requires MOFA2 and supports `--reuse-models-only` to require the retained model set before any writes. Its default mode can train when caches are missing. Workflows01–04,06 and08 have fixed output paths. Use the parallel replay to run them without modifying the release. Workflow10 has known execution blockers documented in its plan. The committed whole-cell CSVs allow downstream analyses to start at 02. Workflow01 re-export requires the local whole-cell FCG/GDX pseudobulks.

## Figures and interpretation

The current release supplies Figure 4C, Figure 5A–D/F, S5A, source Tables S10/S13/S14/S15, and standalone [Figure 6 panels](results/current_panels.md). The approved A4 exports of Figure 6B and 6E are indexed with their source tables and legends. The panels are provided separately for manual assembly. Independent spatial assays and knockdown measurements are preserved.

Figure 4A and S4B–C, their linked nucleus-level counts/composition statistics, and final confirmation of executed retained-marker lists for Methods remain James's tasks. Pseudobulks cannot reconstruct these quantities.

The diet cohort is analysed independently of the younger original cohort. Gene–lipid associations and illustrative gain/loss networks remain exploratory: no tested correlation passes the declared BH family, and all 120 class bootstrap intervals include zero. The current panel source tables document the thresholds and the estimates.

## Optional analyses and historical references

[COSMOS/MOON and PKN analyses](docs/CAUSAL_ANALYSIS_MAP.md) are retained as an active, separate extension, including their frozen reference dependencies. They are not quick tests. Their pinned model-fitting environment is separate from the R environment used for the zonation update.

Earlier cross-cohort analyses (12/13) retain their original two-marker compartment definition and frozen references. They are not rebased to the expanded-marker inputs. Their saved results, previous Figure 6 versions, old model generations, exploratory outputs and refactoring records are under [miscellaneous/](miscellaneous/README.md). Archived source and provenance remain discoverable; large historical payloads are excluded from ordinary Git commits and can be located by their recorded path and SHA256.

## GitHub preparation

```sh
python3 scripts/maintenance/build_workflow_documentation.py --check
python3 scripts/maintenance/verify_repository_organization.py
python3 scripts/maintenance/audit_github_readiness.py --strict
```

The check inventories the files Git would include, oversized files, local metadata and data/provenance availability. It does not stage, commit or push. See [the commit guide](docs/GITHUB_COMMIT_GUIDE.md) for large-data restoration and review steps. This repository already contains substantial scientific data; generated caches, duplicated bundles, bulk old results, local settings and private manuscript files are excluded.

## Contributing and reuse

Read [AGENTS.md](AGENTS.md) before changing analysis code. Use explicit package namespaces, test modified R scripts immediately, and update the corresponding workflow plan. Preserve frozen scientific artifacts and their recorded hashes; use the release resolver when reading a relocated historical path.

Licensed under [LICENSE](LICENSE).
