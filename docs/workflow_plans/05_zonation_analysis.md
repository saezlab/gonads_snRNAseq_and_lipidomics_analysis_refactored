# Workflow 05: Current FCG Zonation Analysis

## Authoritative execution contract (audited 18 November 2026)

- Public entry: `scripts/workflows/05_zonation_analysis.R` (19-line source-safe wrapper).
- Implementation: `scripts/core/zonation_update_fcg.R`, which reads the immutable
  `engine_fcg05` body at `scripts/core/engines/05_zonation_analysis_legacy.R`.
- Default input directory: the parent of registry `current_fcg_zonation_csv`.
  Default output root: registry `current_root` (`results/legacy/pre_reorganization_baseline/zonation_marker_update_20261116`).
- `source()` defines functions; `Rscript` calls `fcg_zonation_run()`. Supported
  arguments are only `--input=DIR`, `--output=ROOT`, and `--validate-only`.
- A fit writes `<ROOT>/fcg05/` plus `<ROOT>/runtime/`. A completed
  `fcg05/fcg05_checkpoint.rds` blocks rerunning into that root. Use a fresh root;
  the guard does not promise to preserve every file in a partial destination.
  `analysis_assert_not_approved_output()` additionally refuses the approved
  registry root and every descendant before scientific work, regardless of
  whether ignored checkpoints exist. This also protects a Git-only clone.
  Always supply a fresh `--output` outside the approved root.

Workflow17 is the main coordinated current zonation/diet entry. It stages CSV
transports and invokes the report driver in `--science-only` mode, running05,
zonal07 and current panel finalization. Running05 alone runs the model/network
engine and saves its checkpoint; it does not run07 or the report driver's extra
Figure5D/context-class panel finalization.

`--validate-only` verifies RNA matrix/metadata integrity, required R packages and
prepared runtime syntax without writing. It does not open/check every later
GO/regulon/Workflow03 resource; file presence there remains a full-run prerequisite.
The frozen engine must not be sourced directly: its old paths, condition ordering
and cache behavior are implementation history rather than maintained defaults.

## Actual inputs and producers

| Input | Producer/role |
| --- | --- |
| `data/RNA/zonation/current/fcg/scaled_RNA_zonation.csv` (`current_fcg_zonation_csv`) | Supplied current FCG compatibility matrix; Workflow01 section11 can regenerate it from `current_fcg_zonation_transport`. Workflow17 instead stages the same values into its fresh run's `inputs/fcg/`. |
| `data/RNA/zonation/current/fcg/meta_data_df_zonation.csv` (`current_fcg_zonation_metadata`) | Matching paired mouse/zone/gonad metadata. |
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

Engine line references are approximately1-2720 in the frozen file, not the19-line
public wrapper. Its old literal `results/zonation/` prefix is replaced in memory.

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
- `runtime/05_zonation_analysis_versioned.R` records the exact executed adaptation.

The report-driver/17 science finalizer additionally writes current standalone
panels/source data under `fcg05/panels/`, `fcg05/panel_source_data/`, the FCG panel
manifest and patch manifest. Its retrospective mode additionally needs explicit
historical registry roles and writes `comparisons/` and metric/replay ledgers.
It is not an upstream requirement of a current05 fit.

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
expanded-marker output. Current zonal07 consumes the selected root's current
checkpoint/networks. Preserve the mixed `Ovaries`/`ova` network filenames.

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
