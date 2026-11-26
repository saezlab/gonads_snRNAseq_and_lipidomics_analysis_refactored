# Workflow 05: Current FCG and gonadectomy zonation

## Authoritative execution contract (25 November 2026)

The source-safe public entry is `scripts/workflows/05_zonation_analysis.R`.
The default cohort is FCG. Use `--cohort=gdx` for the paired gonadectomy model.
Both modes read their current expanded-marker transport directly from the release
registry. `--input=PATH` accepts a transport directory, its manifest, or an explicit
local pseudobulk RDS. FCG also accepts an explicit directory containing
`scaled_RNA_zonation.csv` and `meta_data_df_zonation.csv`, taking that pair if
both representations are present. Supply a fresh `--output=ROOT`. The approved
registry root and the legacy archive are protected from writes.

The FCG mode adapts the immutable engine in `scripts/core/engines/` through
`scripts/core/zonation_update_fcg.R`. It fits the model, builds the networks, and
generates the current Figure 5D and class panels in the same run. The resulting
checkpoint and networks feed zonal Workflow 07. Workflow 17 can run both entries
together with the diet workflows.

`--validate-only` checks the selected cohort’s RNA and metadata and its required
packages. FCG materializes any engine inputs in a temporary directory, checks the
engine syntax, and removes the temporary files. Full FCG runs also require the
GO, TF and Workflow 03 resources below.

## Actual inputs and producers

| Input | Producer/role |
| --- | --- |
| `current_fcg_zonation_transport` | Default verified FCG manifest and payloads. The engine reads scaled expression and paired mouse/zone/gonad metadata materialized under the run’s `inputs/fcg/`. |
| `current_fcg_zonation_csv`, `current_fcg_zonation_metadata` | Shipped compatibility pair for explicit `--input=DIR` use and historical report compatibility. |
| `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | Curated whole-liver lipid data, read through `core/data_loading.R::load_lipidomic_data()`. |
| `support/m5.all.v2023.2.Mm.symbols.gmt` | Curated GO-BP sets for fresh pathway ULM. |
| `support/TF_regulons.RData`, `support/Gonadal_Hormone_Genes.csv` | Curated regulons/custom ESR1/AR target replacement. |
| `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData` | Workflow03 lipid DE, used by the association presentation. |
| `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_scores.RData` | Workflow03 lipid-class enrichment, used for class ordering/presentation. |

These scientific prerequisites are distributed in the public tree with the R
packages; no source Seurat object, Workflow02 DE, Workflow04/MOFA output, historical
zonation model or historical pathway cache is required. Recompute03 first if its
lipid DE/enrichment products are being rebuilt. Paired RNA zones remain within
mouse; the lipid vector is one whole-liver measurement per mouse.

The wrapper sources the registry/adapter; the generated engine runtime sources
`analysis_config.R`, `plotting_themes.R`, `support_functions.R` and `core/data_loading.R`.
Their generic config lists do not replace the explicit input/output overrides.

## Runtime adaptations and order

`zmu_runtime05()` redirects all engine outputs and the input CSV pair, orders the
four gonad-zone conditions explicitly before positional limma contrasts, replaces
engine section7 with fresh `decoupleR::run_ulm(minsize=15)`, and repairs an absent
`no_association` category using a scalar zero. It preserves the seven paired-mouse
limma contrasts, threshold objective, candidate gates and correlation calculations.
The historical pathway cache is an output of the fresh calculation, not an input.

| Engine section | Current runtime behavior |
| --- | --- |
| 1-4 | Configuration, paired zonal CSVs, PCA and metadata enrichment. |
| 5-6 | Seven limma contrasts with mouse `duplicateCorrelation`, then volcanoes. |
| 7 | Replaced by fresh current-input pathway ULM; writes `pathway_enrichment.RData` and `pathway_enrichment_fresh.csv`. |
| 8-10b | Lipid import, TF scores and overall/conditional zonation/sex scatter plots. |
| 11-12 | Optimize interaction threshold and select gonad-specific PP/PC candidates. |
| 13 | Four RNA-lipid/TF networks; full, displayed and class-compressed exports. |
| 14-16 | Lipid association/class/threshold summaries, fitted objects and result tables. |

Engine line references refer to the frozen file, approximately lines 1–2720.
The adapter replaces its `results/zonation/` prefix in memory.

## Outputs under the selected root

- `fcg05/fcg05_checkpoint.rds`, `zonation_ttop_list.RData`,
  `zonation_limma_results.RData`, `t_table_zonation.csv`, `ttop_full_long.csv`,
  `metadata.csv`, `summary_statistics.csv`, `warnings.csv`, `sessionInfo.txt`.
- `fcg05/pathway_enrichment.RData`, `pathway_enrichment_fresh.csv`.
- `fcg05/optimal_threshold_metrics.csv`, `sex_specific_zonation_genes_summary.csv`,
  `lipid_species_zonation_association.csv`, `lipid_category_binary_matrix.csv`,
  `lipid_class_association_summary.csv`, and class/total threshold-sensitivity CSVs.
- `fcg05/networks/cormat_{testes_portal,testes_central,Ovaries_portal,ova_central}_{SIF,ATT}.csv`,
  the same four `_full_{SIF,ATT}.csv` pairs and `compressed_cormat_*_{SIF,ATT}.csv`.
- Engine PCA, volcano, scatter, candidate, class and threshold PDFs under `fcg05/`.
- Figure 5D and the explicit gene-context class panel under `fcg05/panels/`,
  their plotted data under `fcg05/panel_source_data/`, and
  `fcg_panel_manifest.csv` at the run root.
- `runtime/05_zonation_analysis_versioned.R` records the exact executed adaptation.
- `inputs/fcg/input_provenance.csv` records the hashes of the manifest, payloads
  and generated engine CSVs when reading the transport directly. The manifest
  retains the separate original RDS identity.

A completed `fcg05/fcg05_checkpoint.rds` marks the end of the fit and panel export
and blocks a rerun into the same root. Files from a partial run can be rewritten.
The retrospective report driver uses explicit historical registry roles and also
writes old/new comparisons and metric/replay ledgers.

## Paired gonadectomy mode

`--cohort=gdx` reads `current_gdx_zonation_transport` and executes
`scripts/core/zonation_update_gdx.R`. Each of the 22 mice supplies two zonal RNA
profiles. The model fits the eight sex × surgery × zone means on the supplied
scaled expression using limma with mouse blocking. Gene eligibility requires
counts above zero in at least five profiles and a scaled standard deviation above
1e−8. The gene × contrast BH families are separate for simple zonation effects,
conditional gonadectomy effects and the three gonadectomy–zonation interactions.

The mode writes `gdx_zonation/new_{results,summary,eligibility,metadata}.csv`
and `gdx_paired_analysis.rds`. The historical comparison is available through
the Workflow 17 audit. The published GDX directory is a supporting result.

## Contrasts and candidate semantics

Contrasts are `portal_vs_Central`, `testes_vs_ova`,
`portalova_vs_centralova`, `portaltestes_vs_centraltestes`,
`centraltestes_vs_centralova`, `portaltestes_vs_portalova`, and
`testes_effect_zonation`. Own-zone t signs assign PP positive and PC negative.
Candidates require interaction absolute-t above the optimized threshold and weak
other-gonad zonation absolute-t below1. These are the original selection rules,
not a direct FDR-defined differential-correlation test.

The faithful class butterfly groups by lipid DE direction and pooled PP/PC gene
association; it must not be relabelled as four gene-sex/zone programme fractions.
The science finalizer supplies a separate explicit-context class panel for that
interpretation. Networks contain measured whole-liver associations, not spatial
lipid localization or causal effects.

## Frozen consumers and vulnerabilities

Workflow12/13 use `historical_fcg05` and historical input roles to replay the old
two-marker Figure5 branches. They do not automatically consume this current
expanded-marker output. Current zonal07 consumes the same-root FCG checkpoint
and its adjacent networks by default. Its `--parent=` option selects a completed
FCG run in another root. Preserve the mixed `Ovaries`/`ova` network filenames.

Runtime adaptation is text/section based: engine header markers and exact literal
substitutions form a compatibility contract. Metadata sample suffixes and the
four-condition order matter. Threshold state created in section11 is consumed
later in section14; the locally defined correlation helper is reused there.
Fresh/current models must not be combined with old cached pathways or network
coefficients. The approved root and historical engine/reference bytes are frozen.

## Historical production and replay record

## Expanded-marker FCG rerun, 16 November 2026

The isolated adapter `scripts/core/zonation_update_fcg.R`, invoked through
`scripts/reports/run_zonation_update_fcg.R`, reruns this workflow against the
staged expanded-marker FCG zonation export under
`results/legacy/pre_reorganization_baseline/zonation_marker_update_20261116/inputs/fcg/`. This paragraph records the isolated production run before the later default
promotion. Historical inputs/results remain archived; current defaults are the
maintained wrapper and registry paths documented above. Its versioned
runtime copy and outputs are in `runtime/` and `fcg05/` under the update root.
See Workflow 17 for the coordinated multi-cohort revision scope.

The adapter preserves the seven blocked limma contrasts and original
selection/correlation calculations. It explicitly orders the four conditions
before using the legacy positional contrast vectors; an old-input replay
reproduces every archived contrast (maximum absolute t difference below
`5e-12`). Pathway ULM is computed fresh (`minsize=15`) instead of loading the
old cache. One compatibility repair changes absent-category lookup to a scalar
zero via `sum(proportion[zonation == "no_association"])`; lipid denominators
and the entropy formula are preserved. Reoptimization on the historical inputs
still selects 5.4, whereas revised inputs select 4.4.

Verified selected-gene counts are 78 historical, 173 at the revised optimum,
and 82 for revised inputs at the fixed historical cutoff 5.4. Both operating
points, membership changes and fixed historical gene–lipid pair coefficients
are exported under `comparisons/`; this separates input effects from changes
caused by threshold reoptimization. Measured target identity, rather than
parentheses in target names, identifies lipid features in comparison tables,
retaining CoQ8, CoQ9 and Uridine-5-diphosphoacetylglucosamine in the full
1,254-feature universe.

Panel paths and source tables are listed in `fcg_panel_manifest.csv`; numerical
keys are in `fcg_metric_ledger.csv`. The adapter supplies Figure 4B, Figure
5A/B-counts/C/D and the S5A threshold panel, with Figure 5F from the validation
branch. It does not regenerate nucleus composition or marker-QC panels.
The historical Figure 5C butterfly groups by lipid differential-abundance sex
direction and pooled PP/PC gene association. That calculation is retained and
explicitly documented; an additional class-fraction heatmap uses the actual
four gene-sex/zone contexts and has a separate source table/manifest entry.

The Figure 4B manuscript claim is reproduced from the WF05 sex-FDR set:
historical 2,783/4,276 (1,441 PC; 1,342 PP), revised 2,901/4,449
(1,539 PC; 1,362 PP). A separate audit holds the unchanged WF02 whole-hepatocyte
sex-FDR set fixed at 4,103 genes. Display colouring in the faithful WF05
scatter remains the original absolute-t threshold, rather than these FDR
claim sets. Existing legend sign and sample-size discrepancies are recorded
separately from the marker update.

The documented report driver also rejects approved-root full fitting before either
core fit; explicit validate/postprocess-only modes retain their stated behavior.
The report driver `scripts/reports/run_zonation_update_fcg.R` accepts
`--validate-only`, `--input=...`, `--output=...`, `--science-only`, and
`--postprocess-only`; these report-driver modes are not extra flags of the
numbered05 wrapper. `--postprocess-only --rerun-zonal-validation` is an explicit
audit refresh against saved models and archived references. It can rewrite
presentation/audit outputs and must not be confused with read-only validation.
Runtime adaptations, source/fit hashes, execution logs and warnings are retained.
