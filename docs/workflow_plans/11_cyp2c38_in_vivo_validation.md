# Workflow 11 — Cyp2c38 in vivo validation (WP1)

## Release routing - 16 November 2026

This independent whole-hepatocyte/KO workflow retains its inputs and statistical
scope. The expanded-marker change does not rebase its frozen MOFA reference.
The optional original global-siRNA report reference is opened through the release
path resolver so recorded provenance identities survive archive relocation.

---

**Entry point:** `scripts/workflows/11_cyp2c38_in_vivo_validation.R`.
**Orchestration:** `scripts/core/wp1_workflow.R`.
**Consolidated:** 9 October 2026. Active complete WP1 workflow; former 07b
and 07c scripts remain thin phase compatibility wrappers. Their plans retain
full scientific specifications and original production results.

Read this plan, both phase plans and `_meta_cross_workflow_overview.md` before
changes. Consolidation changes paths and execution organization; statistical
helper implementations and primary numerical results are preserved. Workflow
07's primary-hepatocyte siRNA experiments remain separate.

## Sections and execution order

1. **Configuration:** source `scripts/support_functions.R`, enable the existing
   `.local/R-library` when present, and load `wp1_workflow.R`.
2. **Intake/reference replay:** `wp1_run()` reads local KO data, separates QC
   injections from animals, creates feature/sample ledgers and reconstructs
   the historical paired FCG reference.
3. **Analyses:** base weighted differential models and historical concordance,
   then `wp1_ulm_run(previous=out)` reading the newly created base results.
   Rebuild full/historical gonadal networks; perform native ULM, explicit null
   tests and calibration/injection diagnostics.
4. **Exports:** phase CSV/RDS/plots/Markdown, input/code hashes and session
   provenance; package `workflow_manifest.json` records phase arguments/order.

`source()` defines functions without executing analysis; `Rscript` calls
`wp1_package_run()`. Run from the repository root. Required packages include
limma, statmod, digest, jsonlite, decoupleR and the existing plot dependencies.

## Inputs and invariants

| Repository input | Role |
|---|---|
| `data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv` | Integrated KO export, including nine analytical QC injections |
| `data/RNA/scaled_RNA_{1,2,3}.csv`, `data/RNA/meta_data_df.csv` | Original whole-hepatocyte FCG profiles/metadata |
| `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | Original 1,254-feature lipid matrix |
| `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv`, `data/Suppl_Table_1_OverviewData.xlsx` | Independently checked original FCG groups |
| `results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | Historical 828-edge Cyp2c38 universe |

Required KO intake: 1,237 primary endogenous features, 1,215 QC-CV<=30%
features, 32 mice in F_WT/F_KO/M_WT/M_KO counts 8/8/7/9, and 552 strict
historical matches. Regenerated full-network matching has 856 features.
FCG reference retains 43 paired mice; XXO/XYO/XXT/XYT counts 10/11/11/11,
ovarian/testicular counts 21/22. Preserve unique keys, finite positive supplied
values, strict one-to-one molecular matching, contrast signs and reference replay.

Imported data are byte-identical to shared archive sources; no processing was
changed. See `data/WP1_WP2_INPUTS.md` and
`data/wp1_wp2_input_migration_manifest.csv`. No active WP1 input path uses the
sibling `original_datasets` directory; original archives remain available.

## Statistical specification

Base model: log2 supplied intensities, four group means, sex precision weights,
robust/trend limma moderation. Female/male effects are KO-minus-WT; pooled
is their equal-weight average. Preserve unweighted, separate-sex, QC and
supplied-scale sensitivities and historical whole-animal resampling.

Current requested reporting is in `ulm_followup/`: native intercept ULM scores,
nominal ULM p-values and explicit bootstrap p-values. Negative ULM indicates
alignment with reduced Cyp2c38. Female/male/pooled KO contrasts match
ovarian/testicular/pooled reference networks; historical and full universes
remain distinct and feature selection is frozen.

Primary HC2 wild bootstrap constrains the tested KO contrast to zero for every
lipid while preserving nuisance mean dimensions. One shared Rademacher sign
acts on each mouse vector; weights/moderation are refitted per draw. Reference
resampling and restricted-residual tests remain separate sensitivities. Retain
within-sex genotype permutations, conditional null calibration and aligned-effect
injections. See phase 07c plan for assumptions and the original findings.
Dependent-edge nominal p-values and fitted injection curves differ from
animal-level null evidence. Inconclusive KO concordance does not establish
absence of Cyp2c38 regulation or negate independent experimental evidence.

## CLI and defaults

```bash
Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --validate-only
Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --output=results/validation/cyp2c38_in_vivo_reproduction
```

| Setting | Base phase | ULM/calibration phase |
|---|---:|---:|
| Seed | 20261006 | 20261007 |
| Genotype permutations | 1,999 | 1,999 |
| Main bootstrap draws | 999 | 4,999 |
| Restricted-residual draws | — | min(main draws, 1,999) |
| Workers | 3 | 3 |
| Calibration datasets / inner draws | — | 120 / 199 |
| Injection datasets per strength/contrast | — | 399 |

Options: `--output`, `--permutations`, `--bootstraps`, `--workers`, `--seed`,
`--base-permutations`, `--base-bootstraps`, `--base-seed`,
`--calibration-datasets`, `--calibration-inner`, `--injection-datasets`,
`--validate-only`. Unknown/duplicate options are rejected. Explicit common
permutation/bootstrap options affect both phases unless the corresponding
`--base-*` override is supplied. Omitting them preserves distinct defaults;
explicit production settings include `--bootstraps=4999 --base-bootstraps=999`.
`--seed` controls ULM; `--base-seed` controls the base phase.

`--validate-only` passes fresh raw intake directly to ULM validation; it requires
no saved WP1 results and creates no outputs. Full runs reject nonempty output
directories. Default root is `results/validation/cyp2c38_in_vivo/`, with
ULM results in `ulm_followup/`.

## Outputs and helper ownership

Base: differential/sensitivity tables, fits and t-values, feature/QC/sample and
matching ledgers, historical concordance/influence diagnostics, plots and README.
ULM: `ULM_summary.csv`, `ULM_matched_lipid_edges.csv`, six network SIFs,
null draw banks, calibration/injection diagnostics and current README.
Helpers retain `cyp2c38_*` names in `scripts/core/`. Existing PDF builders under
`scripts/reports/` remain available as separate steps from R analysis.
Original historical run hashes remain historical records; fresh runs record
new local paths/code hashes and the package manifest.

## Verified consolidation and limitations

Actual validation of both phases passed from a nonexistent output directory;
source loading, CLI routing/defaults and overwrite protection passed independent
review. A fresh complete reduced-draw smoke run executed all analysis/report
phases. Four differential/t-value/QC/matching tables reproduced within 1e-10;
all ten observed native ULM scores/p-values matched exactly. Temporary reduced
Monte Carlo outputs were not promoted to primary results. Original full
production draws and numerical outputs were preserved; their verification is
recorded in phase plans.

Tests: `tests/test_cyp2c38_in_vivo.R`, `test_cyp2c38_sex_networks.R`,
`test_cyp2c38_ulm_bootstrap.R`, `test_wp1_wp2_local_inputs.R` and
`test_work_package_workflows.R`. Actual replay is available with
`Rscript tests/test_work_package_workflows.R --wp1-output=/path/to/fresh/output`.
Consolidation evidence is under `tests/outputs/wp1_wp2_migration/`.
Run every modified R script immediately and update this plan after changes.

Inherited scientific limits: heterozygous perturbation, unknown batch/litter and
upstream normalization provenance, fixed reference selection, correlated lipids,
and conditional bootstrap/simulation assumptions. Local import resolves external
paths and fresh validation resolves saved-output dependence; neither resolves
those scientific limitations.
