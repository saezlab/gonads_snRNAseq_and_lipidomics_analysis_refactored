# Workflow execution and data reference

Code audit: **2026-11-18**, working tree based on commit `f45e2a9702dae4274d8581d647f219b548525bdb`.

Generated from [workflow_dependencies.json](workflow_dependencies.json). Use [PIPELINE.md](../PIPELINE.md) for the short execution map and the linked per-workflow plans for methods. This reference distinguishes reading an existing result from producing it, and orchestration from data flow.

**Availability is a checked-out-file inventory, not a full execution guarantee.** A shipped derived result can be consumed without rerunning its producer. Local model caches and private/external inputs must be regenerated or restored as stated. Validation flags have the exact scope shown below; none is a substitute for a scientific reproduction run.

**Numbering is not execution order.** Workflow17 orchestrates the current zonation/diet branch; compatibility wrappers and external/reference scripts are not extra sequential main-pipeline steps.

## Workflow index

| ID | Role | Entry point | Detailed plan |
| --- | --- | --- | --- |
| [00](#workflow-00) | optional_external_inactive_reference | [00_snrnaseq_preprocessing.R](../scripts/workflows/00_snrnaseq_preprocessing.R) | [00_snrnaseq_preprocessing.md](workflow_plans/00_snrnaseq_preprocessing.md) |
| [01](#workflow-01) | local_preparation_not_public_prerequisite | [01_data_preparation_qc.R](../scripts/workflows/01_data_preparation_qc.R) | [01_data_preparation_qc.md](workflow_plans/01_data_preparation_qc.md) |
| [02](#workflow-02) | public_science_fixed_output | [02_fcg_pseudobulk_analysis.R](../scripts/workflows/02_fcg_pseudobulk_analysis.R) | [02_fcg_pseudobulk_analysis.md](workflow_plans/02_fcg_pseudobulk_analysis.md) |
| [03](#workflow-03) | public_science_fixed_output | [03_fcg_lipidomic_analysis.R](../scripts/workflows/03_fcg_lipidomic_analysis.R) | [03_fcg_lipidomic_analysis.md](workflow_plans/03_fcg_lipidomic_analysis.md) |
| [04](#workflow-04) | public_science_with_retained_model_cache | [04_mofa_integration.R](../scripts/workflows/04_mofa_integration.R) | [04_mofa_integration.md](workflow_plans/04_mofa_integration.md) |
| [05](#workflow-05) | current_component_public_CSV_inputs | [05_zonation_analysis.R](../scripts/workflows/05_zonation_analysis.R)<br>[run_zonation_update_fcg.R](../scripts/reports/run_zonation_update_fcg.R) | [05_zonation_analysis.md](workflow_plans/05_zonation_analysis.md) |
| [06](#workflow-06) | public_science_fixed_output | [06_structural_analysis.R](../scripts/workflows/06_structural_analysis.R) | [06_structural_analysis.md](workflow_plans/06_structural_analysis.md) |
| [07](#workflow-07) | source_safe_two_mode_component | [07_validation_experiments.R](../scripts/workflows/07_validation_experiments.R) | [07_validation_experiments.md](workflow_plans/07_validation_experiments.md) |
| [08](#workflow-08) | independent_public_science_fixed_output | [08_variance_analysis.R](../scripts/workflows/08_variance_analysis.R) | [08_variance_analysis.md](workflow_plans/08_variance_analysis.md) |
| [09](#workflow-09) | Optional external analysis; no public-checkout end-to-end validation. | [09_human_analysis.py](../scripts/workflows/09_human_analysis.py) | [09_human_analysis.md](workflow_plans/09_human_analysis.md) |
| [10](#workflow-10) | Imported provenance; existing fresh-process blockers; intended order00→01→02/03 is not a validated standalone run. | [00_infer_doublets.py](../scripts/workflows/10_xenium_analysis/00_infer_doublets.py)<br>[01_merge_doublets.py](../scripts/workflows/10_xenium_analysis/01_merge_doublets.py)<br>[02_compositional_analysis.py](../scripts/workflows/10_xenium_analysis/02_compositional_analysis.py)<br>[03_DGE_analysis.py](../scripts/workflows/10_xenium_analysis/03_DGE_analysis.py) | [10_xenium_analysis.md](workflow_plans/10_xenium_analysis.md) |
| [11](#workflow-11) | active independent extension | [11_cyp2c38_in_vivo_validation.R](../scripts/workflows/11_cyp2c38_in_vivo_validation.R) | [11_cyp2c38_in_vivo_validation.md](workflow_plans/11_cyp2c38_in_vivo_validation.md) |
| [12](#workflow-12) | active historical-definition comparison | [12_wdsd_network_comparison.R](../scripts/workflows/12_wdsd_network_comparison.R) | [12_wdsd_network_comparison.md](workflow_plans/12_wdsd_network_comparison.md) |
| [13](#workflow-13) | active historical-definition comparison | [13_wd_crosscohort_robustness.R](../scripts/workflows/13_wd_crosscohort_robustness.R) | [13_wd_crosscohort_robustness.md](workflow_plans/13_wd_crosscohort_robustness.md) |
| [14](#workflow-14) | current upstream model stage; full presentation alternative | [14_cf453_standalone_figure6.R](../scripts/workflows/14_cf453_standalone_figure6.R) | [14_cf453_standalone_figure6.md](workflow_plans/14_cf453_standalone_figure6.md) |
| [15](#workflow-15) | current upstream direct-evidence stage; four-panel alternative | [15_cf453_diet_impact_figure6.R](../scripts/workflows/15_cf453_diet_impact_figure6.R) | [15_cf453_diet_impact_figure6.md](workflow_plans/15_cf453_diet_impact_figure6.md) |
| [16](#workflow-16) | current Figure6 presentation | [16_cf453_diet_rewiring_figure6.R](../scripts/workflows/16_cf453_diet_rewiring_figure6.R) | [16_cf453_diet_rewiring_figure6.md](workflow_plans/16_cf453_diet_rewiring_figure6.md) |
| [17](#workflow-17) | current main zonation/diet entry | [17_expanded_marker_zonation_update.R](../scripts/workflows/17_expanded_marker_zonation_update.R) | [17_expanded_marker_zonation_update.md](workflow_plans/17_expanded_marker_zonation_update.md) |
| [05b](#workflow-05b) | compatibility wrapper; historical scientific definition | [05b_wdsd_network_comparison.R](../scripts/workflows/05b_wdsd_network_comparison.R) | [05b_wdsd_network_comparison.md](workflow_plans/05b_wdsd_network_comparison.md) |
| [05c](#workflow-05c) | compatibility wrapper; historical scientific definition | [05c_wdsd_lipid_overlap_and_pooled_network.R](../scripts/workflows/05c_wdsd_lipid_overlap_and_pooled_network.R) | [05c_wdsd_lipid_overlap_and_pooled_network.md](workflow_plans/05c_wdsd_lipid_overlap_and_pooled_network.md) |
| [07b](#workflow-07b) | compatibility wrapper; historical scientific definition | [07b_cyp2c38_in_vivo_validation.R](../scripts/workflows/07b_cyp2c38_in_vivo_validation.R) | [07b_cyp2c38_in_vivo_validation.md](workflow_plans/07b_cyp2c38_in_vivo_validation.md) |
| [07c](#workflow-07c) | compatibility wrapper; historical scientific definition | [07c_cyp2c38_ulm_resampling.R](../scripts/workflows/07c_cyp2c38_ulm_resampling.R) | [07c_cyp2c38_ulm_resampling.md](workflow_plans/07c_cyp2c38_ulm_resampling.md) |
| [6A_layouts](#workflow-6a_layouts) | current optional display entry | [cf453_figure6a_layouts.R](../scripts/reports/cf453_figure6a_layouts.R) | [16_cf453_diet_rewiring_figure6.md](workflow_plans/16_cf453_diet_rewiring_figure6.md) |
| [MOON_EXPANDED](#workflow-moon_expanded) | Active optional sensitivity; retained preferred outputs differ from runner defaults. | [audit_measured_lipid_reaction_targets.py](../scripts/pkn/audit_measured_lipid_reaction_targets.py)<br>[build_expanded_lipid_rules.py](../scripts/pkn/build_expanded_lipid_rules.py)<br>[prepare_expanded_factor1_mapping.py](../scripts/moon/ancillary/prepare_expanded_factor1_mapping.py)<br>[build_currency_traversal_policy.py](../scripts/moon/ancillary/build_currency_traversal_policy.py)<br>[run_expanded_factor1_tf_to_lipid.R](../scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R)<br>[reextract_expanded_factor1_thresholds.R](../scripts/moon/ancillary/reextract_expanded_factor1_thresholds.R)<br>[export_expanded_factor1_network.py](../scripts/moon/ancillary/export_expanded_factor1_network.py) | [moon_factor1_tf_to_lipid_expanded_species.md](workflow_plans/moon_factor1_tf_to_lipid_expanded_species.md) |
| [MOON_KO](#workflow-moon_ko) | Optional no-fit evidence/proposal branch; raw-count regeneration needs local inputs. | [audit_factor1_ko_expression.R](../scripts/moon/ancillary/audit_factor1_ko_expression.R)<br>[propose_factor1_ko_validation.py](../scripts/moon/ancillary/propose_factor1_ko_validation.py) | [moon_factor1_ko_validation_prioritization.md](workflow_plans/moon_factor1_ko_validation_prioritization.md) |
| [MOON_LIPID_TO_TF](#workflow-moon_lipid_to_tf) | Active optional branch; main entry produces coarsegrid, not every retained refinement/export. | [run_factor1_lipid_to_tf.R](../scripts/moon/run_factor1_lipid_to_tf.R)<br>[refine_lipid_tf_thresholds.R](../scripts/moon/ancillary/refine_lipid_tf_thresholds.R)<br>[audit_lipid_tf_scores.R](../scripts/moon/ancillary/audit_lipid_tf_scores.R) | [moon_factor1_lipid_to_tf.md](workflow_plans/moon_factor1_lipid_to_tf.md) |
| [MOON_RELAXED](#workflow-moon_relaxed) | Active optional sensitivity; distinct from original-prior lipid→TF and selected baseline. | [run_factor1_lipid_to_tf_relaxed.R](../scripts/moon/run_factor1_lipid_to_tf_relaxed.R)<br>[search_factor1_lipid_to_tf_relaxed_parameters.R](../scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R) | [moon_factor1_lipid_to_tf_relaxed.md](workflow_plans/moon_factor1_lipid_to_tf_relaxed.md) |
| [MOON_TF_TO_LIPID](#workflow-moon_tf_to_lipid) | Active optional branch; saved result/reference checks available; current full-fit software gate fails. | [run_factor1_tf_to_lipid.R](../scripts/moon/run_factor1_tf_to_lipid.R)<br>[validate_frozen_archives.R](../scripts/moon/ancillary/validate_frozen_archives.R) | [moon_factor1_tf_to_lipid.md](workflow_plans/moon_factor1_tf_to_lipid.md) |
| [PKN](#workflow-pkn) | Optional resource generation; retained priors/evidence ship; full source rebuild has additional packages/raw prerequisites. | [run_build.py](../scripts/pkn/run_build.py)<br>[build_stitch_directional.py](../scripts/pkn/build_stitch_directional.py)<br>[validate_export.R](../scripts/pkn/validate_export.R) | [CAUSAL_ANALYSIS_MAP.md](CAUSAL_ANALYSIS_MAP.md) |

## Workflow 00

**Raw sci-RNA preprocessing reference — optional_external_inactive_reference**

Document raw import/QC/clustering/doublets, cell-type aggregation and historical two-marker zonation.

**Entry:** `scripts/workflows/00_snrnaseq_preprocessing.R`

**Validation:** No CLI/validate-only/source guard. Both data-processing blocks are if(FALSE); loading still loads packages/defines helpers/prints.

**Writes / overwrite behavior:** If blocks are manually enabled, fixed HPC filenames can be overwritten and setwd changes session state.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Raw alignments and source sample workbooks<br>Nucleus-level original sequencing/preprocessing inputs. | `/omics/odcf/analysis/OE0538_projects/DO-0012/processed_data/20240524_36047_sx42b/36047/alignment/`<br>`/omics/odcf/analysis/OE0538_projects/DO-0012/processed_data/41064_gdxhr/41064/alignment`<br>`/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/sci_metadata.xlsx`<br>`/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/sci_metadata.xlsx` | external<br>manually enabled processing | External / runtime-supplied<br>External DKFZ/HPC inputs, not a public-clone prerequisite; importer is an intentional stop placeholder. |
| Common annotation GTF<br>config$fcg$gtf_file is used by BOTH cohort imports. | `/omics/groups/OE0538/internal/users/j083d/sci-fcg/data/annotations/gencode.vM34.basic.annotation.gtf` | external<br>manually enabled processing | External / runtime-supplied<br>External/HPC. |

| Output | Location | Mode |
| --- | --- | --- |
| Potential FCG intermediates/pseudobulks | `/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/fcg_sci/*.rds` | manually enabled processing |
| Potential GDX intermediates/pseudobulks | `/omics/odcf/analysis/OE0538_projects/DO-0012/fcgl_paper/gdxhr_sci/*.rds` | manually enabled processing |

**Command examples** (not executed by this documentation audit):

Loads packages and defines/prints reference workflow; disabled processing blocks do not run.

```sh
Rscript scripts/workflows/00_snrnaseq_preprocessing.R
```

**Important distinctions:**

- Does not produce the current expanded-marker transports.
- HPC pseudobulks_scc.rds is not automatically handed to01 as data/RNA/fcg_pseudobulks.rds.
- Requires reviewed importer, external data and cluster mapping before enabling.

**Code checked:** [scripts/workflows/00_snrnaseq_preprocessing.R](../scripts/workflows/00_snrnaseq_preprocessing.R).


## Workflow 01

**Local pseudobulk preparation and QC — local_preparation_not_public_prerequisite**

Convert original whole-cell-type pseudobulks to main split CSVs; optional current FCG zonation and legacy GDX conversions.

**Entry:** `scripts/workflows/01_data_preparation_qc.R`

**Validation:** No CLI/source guard/validate-only. Whole body executes; optional sections cannot bypass earlier mandatory main RDS. DATASETS source validation is a separate helper, not entry dry-run.

**Writes / overwrite behavior:** Writes fixed data/RNA and results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc paths and can overwrite curated inputs; use isolated checkout/copy.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Main FCG pseudobulk object<br>Mandatory first input of complete script. | `data/RNA/fcg_pseudobulks.rds` | external<br>default full body | Local-only; restore or regenerate<br>Ignored local Seurat object; public scientific reruns use already supplied CSVs instead. |
| Overview workbook<br>Metadata/PCA merge. | `data/Suppl_Table_1_OverviewData.xlsx` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Current FCG transport<br>Optional section11 current scale.data/typed profile metadata source. | `data/RNA/pseudobulk_csv/current_fcg_zonation/manifest.json`<br>`registry:current_fcg_zonation_transport → data/RNA/pseudobulk_csv/current_fcg_zonation/manifest.json` | curated<br>section11 if manifest exists | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Legacy whole-cell-type GDX pseudobulk object<br>Optional section12 source for02 GDX RNA validation. | `data/RNA/gdx_pseudobulks.rds` | external<br>section12 if local RDS exists | Local-only; restore or regenerate<br>Ignored local object; distinct from17 paired-zone GDX source. |

| Output | Location | Mode |
| --- | --- | --- |
| Main RNA inputs | `data/RNA/scaled_RNA_1.csv`<br>`data/RNA/scaled_RNA_2.csv`<br>`data/RNA/scaled_RNA_3.csv`<br>`data/RNA/meta_data_df.csv` | default |
| Preparation QC | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc/data_quality_report.txt`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc/pseudobulk_per_biosample.pdf`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc/pseudobulk_per_celltype.pdf`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc/sample_celltype_design.csv`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc/samples_per_group_celltype.csv`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/qc/PCA_*.pdf` | default |
| Current FCG compatibility CSVs | `data/RNA/zonation/current/fcg/scaled_RNA_zonation.csv`<br>`data/RNA/zonation/current/fcg/meta_data_df_zonation.csv`<br>`registry:current_fcg_zonation_csv → data/RNA/zonation/current/fcg/scaled_RNA_zonation.csv`<br>`registry:current_fcg_zonation_metadata → data/RNA/zonation/current/fcg/meta_data_df_zonation.csv` | optional section11 |
| Legacy GDX RNA validation CSVs | `data/RNA/scaled_RNA_gdxhr_1.csv`<br>`data/RNA/scaled_RNA_gdxhr_2.csv`<br>`data/RNA/scaled_RNA_gdxhr_3.csv`<br>`data/RNA/meta_data_df_gdxhr.csv` | optional section12 |

**Command examples** (not executed by this documentation audit):

Executes full local preparation with fixed-path writes; needs ignored main RDS.

```sh
Rscript scripts/workflows/01_data_preparation_qc.R
```

**Important distinctions:**

- Main CSVs feed02/04; main RNA metadata also filters03 cohort.
- Only02 consumes the GDX RNA CSVs;03 GDX lipid analysis reads its own experimental CSV.
- Section11 transport can regenerate05 compatibility input, while17 stages current transports independently.
- DATASETS does not redirect hardcoded whole-FCG/GDX RDS paths.

**Code checked:** [scripts/workflows/01_data_preparation_qc.R](../scripts/workflows/01_data_preparation_qc.R), [scripts/config/datasets.R](../scripts/config/datasets.R), [scripts/core/pseudobulk_csv.R](../scripts/core/pseudobulk_csv.R).


## Workflow 02

**FCG RNA differential analysis — public_science_fixed_output**

Fit main cell-type RNA contrasts, pathway/TF activities and optional validation displays.

**Entry:** `scripts/workflows/02_fcg_pseudobulk_analysis.R`

**Validation:** No CLI/source guard/validate-only; source or Rscript executes all sections. vsn is loaded even if proteomics input is absent.

**Writes / overwrite behavior:** Overwrites fixed results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk paths; use isolated checkout/copy.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Main FCG RNA CSVs and profile metadata<br>Supplied main pseudobulk expression and annotation; no need to rerun01. | `data/RNA/scaled_RNA_1.csv`<br>`data/RNA/scaled_RNA_2.csv`<br>`data/RNA/scaled_RNA_3.csv`<br>`data/RNA/meta_data_df.csv` | 01<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Mouse overview<br>RNA sample-name-to-group/allosome merge. | `data/Suppl_Table_1_OverviewData.xlsx` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Pathway and TF resources<br>GO ULM minsize5; TF ULM minsize20. | `support/m5.all.v2023.2.Mm.symbols.gmt`<br>`support/TF_regulons.RData`<br>`support/Gonadal_Hormone_Genes.csv` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Legacy GDX RNA validation<br>Independent legacy GDX RNA comparisons, not17 paired-zone GDX. | `data/RNA/scaled_RNA_gdxhr_1.csv`<br>`data/RNA/scaled_RNA_gdxhr_2.csv`<br>`data/RNA/scaled_RNA_gdxhr_3.csv`<br>`data/RNA/meta_data_df_gdxhr.csv` | 01<br>optional when all four files exist | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Xenium gene panel<br>Select genes from FCG DE for a panel-restricted volcano; no Xenium measurement read. | `support/FCGL_Xenium_Panel.csv` | curated<br>optional if file exists | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| WT bulk proteomics<br>Optional bulk proteomics normalization/DE validation. | `data/proteomics/P3225_fractionated_protein.tsv`<br>`data/proteomics/metadata.csv` | curated<br>optional when both files exist | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |

| Output | Location | Mode |
| --- | --- | --- |
| Main RNA DE/statistics | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/ttop_full_list.RData`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/tval_full_df.RData`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/tval_full_df.csv`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | default |
| RNA enrichment | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/pathways/pathway_scores_df_all.RData`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/pathways/tf_scores_df_all.RData`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/pathways/*.pdf` | default |
| RNA figures | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/pca/*.pdf`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/volcanos/*.pdf` | default |
| Legacy GDX results | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/gdxrh/gdhxr_groups_ttop_list.RData`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/gdx_validation_comparison.pdf` | optional GDX |
| Optional validation displays and proteomics source tables | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/xenium_validation_volcano.pdf`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/proteo_validation_comparison.pdf`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/proteomics_differential_results.csv`<br>`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/proteomics_rna_concordance.csv` | optional panel/proteomics |

**Command examples** (not executed by this documentation audit):

Runs RNA analyses and available optional branches; writes fixed result paths.

```sh
Rscript scripts/workflows/02_fcg_pseudobulk_analysis.R
```

**Important distinctions:**

- Current04/05 do not consume02 DE.
- PCA checked group metadata is used by independent11/frozen-reference consumers; optional09 uses its separately curated hepatocyte signature.
- Xenium panel restriction is not independent spatial validation or a10 output dependency.

**Code checked:** [scripts/workflows/02_fcg_pseudobulk_analysis.R](../scripts/workflows/02_fcg_pseudobulk_analysis.R), [scripts/core/data_loading.R](../scripts/core/data_loading.R).


## Workflow 03

**FCG lipid differential analysis — public_science_fixed_output**

Fit lipid mean-abundance contrasts, class enrichment, optional GDX lipid comparison, network and bubble displays.

**Entry:** `scripts/workflows/03_fcg_lipidomic_analysis.R`

**Validation:** No CLI/source guard/validate-only; source/Rscript executes body.

**Writes / overwrite behavior:** Fixed results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic outputs can be overwritten; isolated checkout/copy required to preserve curated outputs.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Curated main lipid data and mouse overview<br>Log2 lipid data and experimental groups. | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv`<br>`data/Suppl_Table_1_OverviewData.xlsx` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Main RNA mouse roster<br>Restrict main lipid models to RNA-shared mice; no02 DE needed. | `data/RNA/meta_data_df.csv` | 01<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| GDX lipid validation<br>Read original lipid experimental measurements directly. | `data/metabolomics/lipidomics_gdxhr.csv` | curated<br>optional if file exists | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| LINEX2 adjacency<br>Significant lipid-reaction subnetwork. | `support/linex2_adjacency_matrix_OT.csv` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Bubble annotation tables<br>Class-specific structural bubble displays. | `data/metabolomics/for_bubbles/*.csv` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |

| Output | Location | Mode |
| --- | --- | --- |
| Main lipid DE | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData` | default |
| Class enrichment | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_scores.RData`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_plot.pdf` | default |
| Lipid figures/network | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/pca/*.pdf`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/volcanos/*.pdf`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/network/significant_lipids_network_SIF.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/network/significant_lipids_network_ATT.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/bubbles_DBC/*.pdf` | default |
| GDX lipid results | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gdxhr_t_table.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gdxhr_ttop_list.RData`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/fcg_gdx_correlation.pdf`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/fcg_gdx_correlation_by_sex.pdf`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/fcg_gdx_correlation_source.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/fcg_gdx_correlation_statistics.csv` | optional GDX |

**Command examples** (not executed by this documentation audit):

Runs main and available GDX lipid analysis, class/network/bubble displays.

```sh
Rscript scripts/workflows/03_fcg_lipidomic_analysis.R
```

**Important distinctions:**

- 03 intersects to43 RNA-shared FCG mice in supplied data;08 independently retains all48 lipid mice.
- Feeds05 lipid DE/enrichment and06 structural input.
- No dependency on02 DE or02 GDX RNA results.

**Code checked:** [scripts/workflows/03_fcg_lipidomic_analysis.R](../scripts/workflows/03_fcg_lipidomic_analysis.R).


## Workflow 04

**MOFA integration and global networks — public_science_with_retained_model_cache**

Integrate main cell-type RNA/lipids, reuse/train factor models, extract weights/activity and global correlation networks.

**Entry:** `scripts/workflows/04_mofa_integration.R`

**Validation:** Reuse-only mode verifies the eleven exact saved models before any writes. No standalone validate-only mode. MOFA2 is required.

**Writes / overwrite behavior:** Overwrites prepared input, weights, factor/network outputs; skips existing individual HDF5 training files. Use isolated checkout/copy.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Main RNA inputs<br>Supplied main cell-type pseudobulk CSVs; no need to rerun01. | `data/RNA/scaled_RNA_1.csv`<br>`data/RNA/scaled_RNA_2.csv`<br>`data/RNA/scaled_RNA_3.csv` | 01<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Lipid input<br>Whole-liver lipid vector, sample-intersected with RNA for six-view MOFA input. | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Overview and regulatory resources<br>Factor group association and lipid-metabolism/TF annotations. | `data/Suppl_Table_1_OverviewData.xlsx`<br>`support/m5.all.v2023.2.Mm.symbols.gmt`<br>`support/TF_regulons.RData`<br>`support/Gonadal_Hormone_Genes.csv` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Pretrained model set<br>At least11 regex matches skip training; selects12factor model. | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/RNA/mofa_res_*factor.hdf5` | frozen<br>default retained-cache mode | Shipped in Git<br>All11 factor5-15 models exist and are unignored in supplied public tree; MOFA2 R package required. |
| Fallback weights<br>Later fallback if model unavailable; does not bypass earlier training gate or unconditional MOFA2 helper import. | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/main_model_RNA_lipid_weight.RData` | frozen<br>model unavailable | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Python trainer<br>Python mofapy2 training command with GPU then CPU retry. | `scripts/mofa/mofa2.py` | curated<br>automatic if fewer than11 matching HDF5 files | Shipped in Git<br>Python/mofapy2 needed only in training contingency. |

| Output | Location | Mode |
| --- | --- | --- |
| Prepared integration input | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/MOFA_data.RData` | default |
| Training input/models | `data/mofa/MOFA_RNA_data.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/RNA/mofa_res_*factor.hdf5` | automatic missing-cache training |
| Weight/factor products | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/main_model_RNA_lipid_weight.RData`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/factor_analysis/factor_enrichment_gonadal.RData`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/factor_analysis/*.pdf` | default |
| Global network products | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/gene_lipid_correlations_factor1.RData`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/gene_lipid_correlations_factor1.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_full_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_full_ATT.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_ATT.csv` | default |

**Command examples** (not executed by this documentation audit):

Full integration analysis; normally reuse retained11models, otherwise automatically train missing factors.

```sh
Rscript scripts/workflows/04_mofa_integration.R
```

Requires exactly the named models for factor counts 5–15 before writing and skips training. Use the parallel replay for a separate output tree.

```sh
Rscript scripts/workflows/04_mofa_integration.R --reuse-models-only --model-dir=results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/RNA
```

**Important distinctions:**

- Global07 uses permissive gene-lipid full network;12 pooled uses separately frozen sparse TF-gene-lipid branch/precursor.
- Current05/06 do not depend on04.
- Active MOON uses frozen precision inputs historically derived from04, not an automatic04 rerun.
- Generic MOFA_SETTINGS does not control this entry hardcoded factor loop/trainer options.
- Read-only --model-dir is available in --reuse-models-only mode. A missing model stops the run before result writes.

**Code checked:** [scripts/workflows/04_mofa_integration.R](../scripts/workflows/04_mofa_integration.R), [scripts/core/multiomics.R](../scripts/core/multiomics.R).


## Workflow 05

**Current expanded-marker FCG zonation — current_component_public_CSV_inputs**

Run seven blocked zonation/gonad contrasts, fresh pathways, TFs, threshold selection, networks and class associations.

**Entry:** `scripts/workflows/05_zonation_analysis.R`; `scripts/reports/run_zonation_update_fcg.R`

**Validation:** --validate-only checks current RNA integrity, packages and runtime parse without writes; later GO/TF/03 resources are not all opened by validation.

**Writes / overwrite behavior:** Shared approved-root guard rejects current_root and every descendant before full fitting, independent of cache presence (also in documented FCG report driver). Validation and explicit postprocess refresh remain allowed. Outside approved root, completed fcg05 checkpoint blocks rerun; partial destinations may be rewritten. Use explicit fresh --output.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Current FCG compatibility CSV pair<br>Supplied paired-zone matrix/metadata;01section11 or17staging can regenerate from current transport. | `data/RNA/zonation/current/fcg/scaled_RNA_zonation.csv`<br>`data/RNA/zonation/current/fcg/meta_data_df_zonation.csv`<br>`registry:current_fcg_zonation_csv → data/RNA/zonation/current/fcg/scaled_RNA_zonation.csv`<br>`registry:current_fcg_zonation_metadata → data/RNA/zonation/current/fcg/meta_data_df_zonation.csv` | curated<br>default --input | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Current FCG transport<br>Upstream source copied/extracted by17; not opened by bare05 wrapper. | `data/RNA/pseudobulk_csv/current_fcg_zonation/manifest.json`<br>`registry:current_fcg_zonation_transport → data/RNA/pseudobulk_csv/current_fcg_zonation/manifest.json` | curated<br>coordinated17 staging | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Lipid and support data<br>Whole-liver associations; GO/TF annotation. | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv`<br>`support/m5.all.v2023.2.Mm.symbols.gmt`<br>`support/TF_regulons.RData`<br>`support/Gonadal_Hormone_Genes.csv` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Lipid DE and class enrichment<br>Class direction/order and lipid association presentation. | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_scores.RData` | 03<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Frozen implementation engine<br>Read as text and adapted, not directly sourced. | `scripts/core/engines/05_zonation_analysis_legacy.R`<br>`registry:engine_fcg05 → scripts/core/engines/05_zonation_analysis_legacy.R` | frozen<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |

| Output | Location | Mode |
| --- | --- | --- |
| Current model/network results | `results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/fcg05_checkpoint.rds`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/zonation_ttop_list.RData`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/t_table_zonation.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/ttop_full_long.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/optimal_threshold_metrics.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/networks/*.csv`<br>`registry:current_root → results/current` | default: generated below selected fresh --output; current_root paths identify the saved approved release only |
| Fresh pathway/runtime records | `results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/pathway_enrichment.RData`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/pathway_enrichment_fresh.csv`<br>`results/legacy/pre_reorganization_baseline/zonation_marker_update_20261116/runtime/05_zonation_analysis_versioned.R`<br>`registry:current_root → results/current` | default: generated below selected fresh --output; current_root paths identify the saved approved release only |

**Command examples** (not executed by this documentation audit):

Read-only current input/runtime validation.

```sh
Rscript scripts/workflows/05_zonation_analysis.R --validate-only
```

Fit current05 engine/checkpoint into explicit fresh root; does not call07 or science-panel finalizer.

```sh
Rscript scripts/workflows/05_zonation_analysis.R --output=results/runs/fcg_zonation_reproduction
```

**Important distinctions:**

- Source-safe wrapper; calls zmu_run05, not old direct body.
- No02/04/historical model/pathway cache needed for current fit.
- Fresh ULM minsize15 replaces frozen cache section; condition order and absent-category scalar fixes are preserved.
- 17/report-driver --science-only adds07 and current panel finalization; retrospective report-driver modes require archive references.
- 12/13 consume historical_fcg05, not current05 outputs.
- Guard-only hardening was tested with missing-cache Git-like fixtures and model/render stubs; engines/statistical functions were unchanged.

**Code checked:** [scripts/workflows/05_zonation_analysis.R](../scripts/workflows/05_zonation_analysis.R), [scripts/core/zonation_update_fcg.R](../scripts/core/zonation_update_fcg.R), [scripts/reports/run_zonation_update_fcg.R](../scripts/reports/run_zonation_update_fcg.R).


## Workflow 06

**Lipid-only structural regression — public_science_fixed_output**

Parse main lipid DE structures and fit class-specific logFC ~ chain_length + db_position + interaction.

**Entry:** `scripts/workflows/06_structural_analysis.R`

**Validation:** No CLI/source guard/validate-only. Parses nested03 result with structural helper.

**Writes / overwrite behavior:** Overwrites fixed structural root; conditional empty results do not remove stale previous CSVs. Use isolated checkout/copy.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Main lipid DE archive<br>Nested Testes vs Ovaries result table; active workflow has no RNA/MOFA input. | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData` | 03<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |

| Output | Location | Mode |
| --- | --- | --- |
| Structural features and regressions | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipids_with_structure.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipid_only_class_coefficients.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipid_only_class_summaries.csv` | default |
| Conditional significant terms | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipid_only_chain_length_effects.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipid_only_db_position_effects.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipid_only_interaction_effects.csv` | only when qualifying terms exist |
| Structural figures | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipid_only_summary_heatmap.pdf`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipid_only_plots/*_structural_effects.pdf` | default |

**Command examples** (not executed by this documentation audit):

Runs lipid-only structural models/figures on current03 DE.

```sh
Rscript scripts/workflows/06_structural_analysis.R
```

**Important distinctions:**

- Both07modes use lipids_with_structure.csv for shared structural axes.
- Older gene-lipid/MOFA structural concept in header is not active.
- Focus classes SM/BMP/PG/TG and paths are hardcoded despite sourced configs.

**Code checked:** [scripts/workflows/06_structural_analysis.R](../scripts/workflows/06_structural_analysis.R), [scripts/core/structural_analysis.R](../scripts/core/structural_analysis.R).


## Workflow 07

**siRNA global/current-zonal validation — source_safe_two_mode_component**

Validate original gene-lipid predictions with unchanged KD data, with separate global-refit and current-zonal modes.

**Entry:** `scripts/workflows/07_validation_experiments.R`

**Validation:** Default --validate-only checks05checkpoint, KD table/fit and packages (not later structural reference). --global-only --validate-only checks raw/global/structural inputs and runtime parse. Source-safe.

**Writes / overwrite behavior:** Shared approved-root guard rejects full zonal fitting at current_root/descendants independently of caches; validate-only remains allowed. Outside it, completed fcg07 checkpoint blocks rerun; partial destinations may be rewritten. Global-only requires explicit empty/fresh output.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Same-root current FCG reference<br>Paired current FCG reference and candidate networks. | `results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/fcg05_checkpoint.rds`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/networks/cormat_testes_portal_full_SIF.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/networks/cormat_testes_central_full_SIF.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/networks/cormat_Ovaries_portal_full_SIF.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/networks/cormat_ova_central_full_SIF.csv`<br>`registry:current_root → results/current` | 05<br>default zonal | Mixed: some files are local-only<br>Checkpoint is an ignored generated model cache; clone must first fit05/17 in the same fresh output root or restore it explicitly. |
| Unchanged canonical KD fit/table<br>Table is scored; fit RData is additionally checked by the public wrapper. | `results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_ttop_list.RData` | frozen<br>default zonal | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Structural reference<br>Axis reference for structural consistency; required by both modes. | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipids_with_structure.csv` | 06<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Original siRNA lipid measurements<br>Refit original KD effects. | `data/metabolomics/lipidomics_20250423transfections.csv` | curated<br>--global-only | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Global full correlation network<br>Permissive original global reference. | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | 04<br>--global-only | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |
| Frozen validation engine<br>Extract sections8/9/11 for zonal,1-7/10-11 forglobal. | `scripts/core/engines/07_validation_experiments_legacy.R`<br>`registry:engine_validation07 → scripts/core/engines/07_validation_experiments_legacy.R` | frozen<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |

| Output | Location | Mode |
| --- | --- | --- |
| Current zonal validation | `results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg07/fcg07_checkpoint.rds`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg07/zonation_network_validation_scores.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg07/Mtmr11_male_PP_prespecified_ULM.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg07/Mtmr11_male_PP_prespecified_source.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg07/Mtmr11_male_PP_ULM_background_source.csv`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg07/Fig5F_Mtmr11_male_PP_prespecified.pdf`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg07/structural_consistency/*`<br>`registry:current_root → results/current` | default zonal: generated below selected fresh --output; current_root paths identify the saved approved release only |
| Independent global refit and validation | `results/runs/global_validation_reproduction/knockdown_t_table.csv`<br>`results/runs/global_validation_reproduction/knockdown_ttop_list.RData`<br>`results/runs/global_validation_reproduction/global_network_validation_scores.csv`<br>`results/runs/global_validation_reproduction/validation_summary.csv`<br>`results/runs/global_validation_reproduction/structural_consistency/*` | --global-only explicit fresh output |

**Command examples** (not executed by this documentation audit):

Read-only selected-root zonal setup check; requires generated05checkpoint.

```sh
Rscript scripts/workflows/07_validation_experiments.R --validate-only
```

Run zonal validation after05 in the same fresh root.

```sh
Rscript scripts/workflows/07_validation_experiments.R --output=results/runs/fcg_zonation_reproduction
```

Read-only independent global setup check.

```sh
Rscript scripts/workflows/07_validation_experiments.R --global-only --validate-only
```

Refit KD/global/structural independent branch into fresh output; does not overwrite canonical KD inputs.

```sh
Rscript scripts/workflows/07_validation_experiments.R --global-only --output=results/runs/global_validation_reproduction
```

**Important distinctions:**

- Default mode does not read rawKD orglobalMOFA network; globalonly does not read zonation network.
- Prespecified Mtmr11 malePP is separate from adaptive best context; ULM fullKDbackground differs from matchedscatter universe.
- Workflow11 in-vivoKO is independent;07b/07c are compatibility phase entries for11.

**Code checked:** [scripts/workflows/07_validation_experiments.R](../scripts/workflows/07_validation_experiments.R), [scripts/core/zonation_update_fcg.R](../scripts/core/zonation_update_fcg.R), [scripts/core/validation_global.R](../scripts/core/validation_global.R), [scripts/core/engines/07_validation_experiments_legacy.R](../scripts/core/engines/07_validation_experiments_legacy.R).


## Workflow 08

**Whole-liver lipid variability — independent_public_science_fixed_output**

Compare log2-lipid CV/Brown-Forsythe variability between gonad groups and class means.

**Entry:** `scripts/workflows/08_variance_analysis.R`

**Validation:** No CLI/source guard/validate-only; direct full body.

**Writes / overwrite behavior:** Fixed output files can be overwritten; use isolated checkout/copy.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Whole-liver lipid matrix and mouse overview<br>All measured lipid mice; load_metadata default exclude_wt=FALSE. No RNA input. | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv`<br>`data/Suppl_Table_1_OverviewData.xlsx` | curated<br>default | Shipped in Git<br>Public tracked/nonignored files with the required software environment. |

| Output | Location | Mode |
| --- | --- | --- |
| Variance tables | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/variance_analysis/variance_analysis_results.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/variance_analysis/significant_differential_variance_lipids.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/variance_analysis/variance_analysis_by_class.csv`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/variance_analysis/variance_analysis_summary.csv` | default |
| Variance figures | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/variance_analysis/cv_difference_volcano.pdf`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/variance_analysis/cv_difference_volcano_by_class.pdf` | default |

**Command examples** (not executed by this documentation audit):

Run independent lipid-variability models/plots with fixed output paths.

```sh
Rscript scripts/workflows/08_variance_analysis.R
```

**Important distinctions:**

- Supplied roster48mice includes21ovaries,22FCGtestes and5WTtestes;03 intersects to43 instead.
- No03result, RNAzonation, MOFA or17 requirement.
- Global pairedtests pairfeature-level CVmetrics, notmatchedmice; classvalues average supplied log2 species values.

**Code checked:** [scripts/workflows/08_variance_analysis.R](../scripts/workflows/08_variance_analysis.R), [scripts/core/data_loading.R](../scripts/core/data_loading.R).


## Workflow 09

**Optional human bulk-liver comparison — Optional external analysis; no public-checkout end-to-end validation.**

Compare a whole-hepatocyte mouse gonad signature with human bulk-liver sex/age expression.

**Entry:** `scripts/workflows/09_human_analysis.py`

**Validation:** CLI help after base package imports. The explicit DE-only mode bypasses pypath. Full scoring checks for a nonempty mouse-to-human map. No validate-only mode.

**Writes / overwrite behavior:** Explicit --output requires an empty destination. The historical default HPC path can overwrite its named files.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Human read/TPM/SOFT files<br>Counts, normalized comparison layer and sample sex/age/title metadata | `/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/data/bulk/zhou/GSE179250_gene_reads.txt`<br>`/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/data/bulk/zhou/GSE179250_gene_tpm_normalized.txt`<br>`/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/data/bulk/zhou/GSE179250_family.soft.gz` | external<br>full run | External / runtime-supplied<br>External absolute HPC DATA_DIR; not shipped. |
| Mouse signature seed and provisioned copy<br>ID,t_Hepatocytes signature | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/hepa_gonad_signature.csv`<br>`/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/data/hepa_gonad_signature.csv` | frozen<br>full run | Shipped in Git; also runtime/external paths<br>Repository seed ships, but script reads a manually provisioned external copy. Its biological definition is the02gonadcontrast, but02doesnotwrite/copythisnamedseed. |
| Full mouse DE table<br>ID,t_Hepatocytes,adj.P.Val,logFC table | `/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/data/male_vs_female_DA_res.csv` | external<br>ranked concordance stage | External / runtime-supplied<br>Separate required external copy; not produced by09. |
| Orthology resources<br>Live package-resource/cache orthology lookup | `resource:Live package-resource/cache orthology lookup` | external<br>full run | External package/resource<br>Compatible pypath resource/cache or network access; frozen support/ortholog_resource.RData is not read. |

| Output | Location | Mode |
| --- | --- | --- |
| Analysis output | `/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/results/bulk_manuscript/adata_processed.h5ad`<br>`/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/results/bulk_manuscript/*.csv`<br>`/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/results/bulk_manuscript/*.pdf` | external configured output |
| Retained study artifacts | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/human/` | saved outputs; not the script default destination |

**Command examples** (not executed by this documentation audit):

Runs full analysis after external paths/data/packages are provisioned; not a preflight.

```sh
python3 scripts/workflows/09_human_analysis.py
```

Recompute all three human DE contrasts without mouse orthology or signature scoring.

```sh
python3 scripts/workflows/09_human_analysis.py --data-dir HUMAN_INPUTS --output results/runs/human_replay --de-only
```

**Important distinctions:**

- Independent of17 expanded-marker and diet inputs.
- Use --data-dir and --output to override the historical external HPC paths.
- Male−female DE uses female reference. Female age contrast is ≥50 versus12–49.

**Code checked:** [scripts/workflows/09_human_analysis.py](../scripts/workflows/09_human_analysis.py), [docs/workflow_plans/09_human_analysis.md](../docs/workflow_plans/09_human_analysis.md).


## Workflow 10

**Optional collaborator Xenium analysis — Imported provenance; existing fresh-process blockers; intended order00→01→02/03 is not a validated standalone run.**

Imported nucleus-boundary doublet inference, fused annotations, composition plots and cell-level DGE.

**Entry:** `scripts/workflows/10_xenium_analysis/00_infer_doublets.py`; `scripts/workflows/10_xenium_analysis/01_merge_doublets.py`; `scripts/workflows/10_xenium_analysis/02_compositional_analysis.py`; `scripts/workflows/10_xenium_analysis/03_DGE_analysis.py`

**Validation:** No CLI/no-write validation; do not use --help as preflight.03 line14 uses adata before definition;00 lines139–146 mix cell_id keys and sequentialid. Static findings only.

**Writes / overwrite behavior:** Fixed output names;02 overwrites doublet_freq_core.pdf sex plot with TMA plot;03 repeats some filenames. No fresh-output guard.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Raw Xenium and panel workbook<br>Raw boundaries/counts/cell positions and marker-selection columns | `scripts/workflows/10_xenium_analysis/data_raw/output-XETG00078__0018188__Region_1__20240301__143248/nucleus_boundaries.csv.gz`<br>`scripts/workflows/10_xenium_analysis/data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cell_feature_matrix.h5`<br>`scripts/workflows/10_xenium_analysis/data_raw/output-XETG00078__0018188__Region_1__20240301__143248/cells.csv.gz`<br>`scripts/workflows/10_xenium_analysis/data_raw/XeniumTargetList_Mouse-Panel-1.xlsx` | external<br>00/01/03 | Absent locally; restore or generate<br>Not shipped; requires folder-local cwd and data_raw/data_tidy/figures directories. |
| Fused state<br>00 components and01 fused AnnData | `scripts/workflows/10_xenium_analysis/inferred_components.json`<br>`scripts/workflows/10_xenium_analysis/fused_anndata_annotated.h5ad` | curated<br>01 then02/03 | Absent locally; restore or generate<br>Generated local products; existing raw-code blockers must be reviewed before reproduction. |

| Output | Location | Mode |
| --- | --- | --- |
| Components/fused object | `scripts/workflows/10_xenium_analysis/inferred_components.json`<br>`scripts/workflows/10_xenium_analysis/fused_anndata_annotated.h5ad`<br>`scripts/workflows/10_xenium_analysis/data_tidy/` | intended folder-local generated outputs |
| DGE/plots | `scripts/workflows/10_xenium_analysis/figures/`<br>`scripts/workflows/10_xenium_analysis/*.xlsx`<br>`scripts/workflows/10_xenium_analysis/*.pdf` | intended generated outputs |

**Command examples** (not executed by this documentation audit):

Required working directory for relative paths; does not validate data.

```sh
cd scripts/workflows/10_xenium_analysis
```

Intended first stage; review existing component-key blocker before running.

```sh
python3 00_infer_doublets.py
```

Consumes raw files +00JSON, generates fused object.

```sh
python3 01_merge_doublets.py
```

Consumes fused object; generates compositions, with known overwritten output.

```sh
python3 02_compositional_analysis.py
```

Not a supported fresh-process run until existing undefined-adata blocker is addressed.

```sh
python3 03_DGE_analysis.py
```

**Important distinctions:**

- No dependency on17 pseudobulks; not a replacement for James deferred4A/S4B-C nucleus QC.
- 01 workbook needs reason and gene name columns; spatial assignment and cluster/cell-type labels are hardcoded.
- Keep raw/object products outside ordinary Git commits.

**Code checked:** [scripts/workflows/10_xenium_analysis/00_infer_doublets.py](../scripts/workflows/10_xenium_analysis/00_infer_doublets.py), [scripts/workflows/10_xenium_analysis/01_merge_doublets.py](../scripts/workflows/10_xenium_analysis/01_merge_doublets.py), [scripts/workflows/10_xenium_analysis/02_compositional_analysis.py](../scripts/workflows/10_xenium_analysis/02_compositional_analysis.py), [scripts/workflows/10_xenium_analysis/03_DGE_analysis.py](../scripts/workflows/10_xenium_analysis/03_DGE_analysis.py).


## Workflow 11

**Cyp2c38 in-vivo validation (WP1) — active independent extension**

Independent KO differential models, fixed-reference concordance, native ULM/null bootstrap and calibration

**Entry:** `scripts/workflows/11_cyp2c38_in_vivo_validation.R`

**Validation:** --validate-only reconstructs intake/reference and fits observed weighted KO model; no saved WP1 cache, resampling, figures or output writes required.

**Writes / overwrite behavior:** Full phases reject nonempty destinations; default is completed lookup root, use fresh --output.

**Invokes other entries/phases:** `07b`, `07c`. These are control-flow calls, not additional cohorts.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Cyp2c38 in-vivo KO lipid input<br>required | `data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original whole-hepatocyte RNA and metadata<br>required | `data/RNA/scaled_RNA_[123].csv`<br>`data/RNA/meta_data_df.csv` | 01<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original lipid measurements<br>required | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original FCG group/PCA roster<br>required | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | 02<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original full-network Cyp2c38 selection<br>required | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | 04<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Overview workbook recorded in production hashes<br>production provenance dependency | `data/Suppl_Table_1_OverviewData.xlsx` | curated<br>full | Shipped in Git<br>Shipped current artifact or source input |

| Output | Location | Mode |
| --- | --- | --- |
| Complete WP1 package | `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/ulm_followup/` | full |

**Command examples** (not executed by this documentation audit):

Observed-data/reference validation without output

```sh
Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --validate-only
```

Both phases in fresh root

```sh
Rscript scripts/workflows/11_cyp2c38_in_vivo_validation.R --output=results/runs/wp1
```

**Important distinctions:**

- No dependency on zonation05/17 or source Seurat objects.
- Default base B999/seed20261006; followup B4999/seed20261007. Shared B/P options affect both unless base override; --seed controls followup, --base-seed base.
- No report-only mode; PDF builder separate.

**Code checked:** [scripts/core/wp1_workflow.R](../scripts/core/wp1_workflow.R), [scripts/core/cyp2c38_reference.R](../scripts/core/cyp2c38_reference.R), [scripts/core/cyp2c38_sex_networks.R](../scripts/core/cyp2c38_sex_networks.R).


## Workflow 12

**Frozen pooled and zonation network comparison (WP2) — active historical-definition comparison**

Original Figure3/5 conservation versus CF453 SD and descriptive WD-SD changes

**Entry:** `scripts/workflows/12_wdsd_network_comparison.R`

**Validation:** --validate-only runs both actual intake/reference phases without saved WP2 result caches, resampling or writes; historical FCG reference payloads are still required.

**Writes / overwrite behavior:** Complete full run rejects any nonempty root; explicit report-only rewrites saved combined Markdown.

**Invokes other entries/phases:** `05c`, `05b`. These are control-flow calls, not additional cohorts.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| CF453 whole-hepatocyte pseudobulk transport<br>required | `data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json`<br>`data/RNA/pseudobulk_csv/wdsd_hepatocytes/*.csv`<br>`registry:wdsd_hepatocytes_transport → data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF453 lipid measurements<br>required | `data/metabolomics/wdsd/lipidomics_cf453.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original CF453 zone transport<br>required | `data/RNA/pseudobulk_csv/historical_wdsd_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/historical_wdsd_zonation/*.csv`<br>`registry:historical_wdsd_zonation_transport → data/RNA/pseudobulk_csv/historical_wdsd_zonation/manifest.json` | frozen<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen whole-hepatocyte Figure3 networks and precursor<br>required | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_full_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_full_ATT.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/gene_lipid_correlations_factor1.csv` | 04<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original raw lipid annotations<br>required | `data/metabolomics/fcg/lipidomics_fcg.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original whole-hepatocyte RNA and metadata<br>required | `data/RNA/scaled_RNA_[123].csv`<br>`data/RNA/meta_data_df.csv` | 01<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original lipid measurements<br>required | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original FCG group/PCA roster<br>required | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | 02<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original FCG zone matrices<br>required | `miscellaneous/older_inputs/two_marker_zonation/RNA/scaled_RNA_zonation.csv`<br>`miscellaneous/older_inputs/two_marker_zonation/RNA/meta_data_df_zonation.csv`<br>`registry:historical_fcg_zonation_csv → miscellaneous/older_inputs/two_marker_zonation/RNA/scaled_RNA_zonation.csv`<br>`registry:historical_fcg_zonation_metadata → miscellaneous/older_inputs/two_marker_zonation/RNA/meta_data_df_zonation.csv` | frozen<br>all | Local-only; restore or regenerate<br>Retained local archive; payloads excluded from ordinary Git; restore through migration/recovery ledger |
| Four frozen Figure5 context references<br>required | `miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_full_ATT.csv`<br>`registry:historical_fcg05 → miscellaneous/older_results/zonation_two_marker` | frozen<br>all | Local-only; restore or regenerate<br>Retained local archive; numerical payloads mostly ignored; restore original selected SIF/full-SIF/full-ATT, not current expanded-marker networks |

| Output | Location | Mode |
| --- | --- | --- |
| Pooled-first package | `miscellaneous/comparisons/wdsd_network_comparison/`<br>`miscellaneous/comparisons/wdsd_network_comparison/followup/`<br>`registry:historical_wf12 → miscellaneous/comparisons/wdsd_network_comparison` | full |
| Report refresh | `miscellaneous/comparisons/wdsd_network_comparison/README.md`<br>`registry:historical_wf12 → miscellaneous/comparisons/wdsd_network_comparison` | --report-only |

**Command examples** (not executed by this documentation audit):

Validate frozen pooled and zonation prerequisites

```sh
Rscript scripts/workflows/12_wdsd_network_comparison.R --validate-only
```

Pooled then zonal phases into fresh root

```sh
Rscript scripts/workflows/12_wdsd_network_comparison.R --output=results/runs/wp2
```

Saved report refresh at historical default root

```sh
Rscript scripts/workflows/12_wdsd_network_comparison.R --report-only
```

**Important distinctions:**

- No MOFA refit. Pooled phase still loads CF453 historical zones through shared intake even though pooled coefficients use whole hepatocytes.
- Report-only needs root network_analysis.rds and followup/pooled_analysis.rds plus companion CSVs/linked artifacts; these caches are ignored.
- Public transports do not provide original FCG two-marker selection/matrices; restore archive for complete12.

**Code checked:** [scripts/core/wp2_workflow.R](../scripts/core/wp2_workflow.R), [scripts/core/wdsd_reference.R](../scripts/core/wdsd_reference.R), [scripts/core/wdsd_pooled_reference.R](../scripts/core/wdsd_pooled_reference.R), [scripts/core/wdsd_lipid_overlap_audit.R](../scripts/core/wdsd_lipid_overlap_audit.R), [scripts/core/wp2_report.R](../scripts/core/wp2_report.R).


## Workflow 13

**Male WD cross-cohort robustness and original-FCG anchor — active historical-definition comparison**

Compare CF459 WD with CF453 WD and original male FCG on fixed three-cohort pairs

**Entry:** `scripts/workflows/13_wd_crosscohort_robustness.R`

**Validation:** Fresh --validate-only checks input/reference/matching and native ULM without saved13 outputs or resampling; still needs ignored original FCG references.

**Writes / overwrite behavior:** Full run rejects nonempty output. Extension-only requires fresh extension destination; --comparison-output only valid with extension-only. Report-only rewrites Markdown/extension figures.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| CF453 whole-hepatocyte pseudobulk transport<br>required | `data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json`<br>`data/RNA/pseudobulk_csv/wdsd_hepatocytes/*.csv`<br>`registry:wdsd_hepatocytes_transport → data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF453 lipid measurements<br>required | `data/metabolomics/wdsd/lipidomics_cf453.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original CF453 zone transport<br>required | `data/RNA/pseudobulk_csv/historical_wdsd_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/historical_wdsd_zonation/*.csv`<br>`registry:historical_wdsd_zonation_transport → data/RNA/pseudobulk_csv/historical_wdsd_zonation/manifest.json` | frozen<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF459 WD cell-type and zone transports<br>required | `data/RNA/pseudobulk_csv/wd_cell_types/manifest.json`<br>`data/RNA/pseudobulk_csv/wd_cell_types/*.csv`<br>`data/RNA/pseudobulk_csv/wd_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/wd_zonation/*.csv`<br>`registry:wd_cell_types_transport → data/RNA/pseudobulk_csv/wd_cell_types/manifest.json`<br>`registry:wd_zonation_transport → data/RNA/pseudobulk_csv/wd_zonation/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF459 WD lipid export<br>required | `data/metabolomics/wd/lipidomics_wd.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original FCG zone matrices<br>required | `miscellaneous/older_inputs/two_marker_zonation/RNA/scaled_RNA_zonation.csv`<br>`miscellaneous/older_inputs/two_marker_zonation/RNA/meta_data_df_zonation.csv`<br>`registry:historical_fcg_zonation_csv → miscellaneous/older_inputs/two_marker_zonation/RNA/scaled_RNA_zonation.csv`<br>`registry:historical_fcg_zonation_metadata → miscellaneous/older_inputs/two_marker_zonation/RNA/meta_data_df_zonation.csv` | frozen<br>all | Local-only; restore or regenerate<br>Retained local archive; payloads excluded from ordinary Git; restore through migration/recovery ledger |
| Four frozen Figure5 context references<br>required | `miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_full_ATT.csv`<br>`registry:historical_fcg05 → miscellaneous/older_results/zonation_two_marker` | frozen<br>all | Local-only; restore or regenerate<br>Retained local archive; numerical payloads mostly ignored; restore original selected SIF/full-SIF/full-ATT, not current expanded-marker networks |
| Original lipid measurements<br>required | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original FCG group/PCA roster<br>required | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | 02<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Prior WP2 coefficient replay<br>optional replay only | `miscellaneous/comparisons/wdsd_network_comparison/network_edges.csv`<br>`registry:historical_wf12 → miscellaneous/comparisons/wdsd_network_comparison` | 12<br>all | Mixed: some files are local-only<br>If absent, comparison proceeds and records replay unavailable |
| Saved base draw banks and contexts<br>required | `miscellaneous/comparisons/wd_crosscohort_robustness/contexts.rds`<br>`miscellaneous/comparisons/wd_crosscohort_robustness/analysis.rds`<br>`registry:historical_wf13 → miscellaneous/comparisons/wd_crosscohort_robustness` | 13<br>--fcg-comparison-only | Mixed: some files are local-only<br>Ignored generated caches; restore or regenerate full13 |
| Saved report bundle<br>required | `miscellaneous/comparisons/wd_crosscohort_robustness/report_bundle.rds`<br>`miscellaneous/comparisons/wd_crosscohort_robustness/fcg_comparison/analysis.rds`<br>`registry:historical_wf13 → miscellaneous/comparisons/wd_crosscohort_robustness` | 13<br>--report-only | Mixed: some files are local-only<br>Bundle required; extension cache used when present, otherwise bundle extension |

| Output | Location | Mode |
| --- | --- | --- |
| Full comparison | `miscellaneous/comparisons/wd_crosscohort_robustness/`<br>`registry:historical_wf13 → miscellaneous/comparisons/wd_crosscohort_robustness` | full |
| Original-FCG extension | `miscellaneous/comparisons/wd_crosscohort_robustness/fcg_comparison/`<br>`registry:historical_wf13 → miscellaneous/comparisons/wd_crosscohort_robustness` | full or --fcg-comparison-only |
| Saved report/extension figures | `miscellaneous/comparisons/wd_crosscohort_robustness/README.md`<br>`miscellaneous/comparisons/wd_crosscohort_robustness/plots/`<br>`registry:historical_wf13 → miscellaneous/comparisons/wd_crosscohort_robustness` | --report-only |

**Command examples** (not executed by this documentation audit):

Observed fixed-reference validation

```sh
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --validate-only
```

Full comparison, sensitivities and original-FCG extension

```sh
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --output=results/runs/wp3
```

Reuse saved WD draws; compute extension/FCG draws

```sh
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --fcg-comparison-only --output=PATH --comparison-output=FRESH_PATH
```

Saved report and extension-plot refresh

```sh
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --report-only --output=PATH
```

**Important distinctions:**

- No CF459 SD arm and no expanded-marker CF459 transport.
- WF12 execution is not mandatory; optional replay is not feature selection.
- Extension-only rejects bootstrap-count overrides and reuses saved counts.
- Retain user-corrected ages and conditional, descriptive interpretation.

**Code checked:** [scripts/core/wd_crosscohort_workflow.R](../scripts/core/wd_crosscohort_workflow.R), [scripts/core/wd_crosscohort_setup.R](../scripts/core/wd_crosscohort_setup.R), [scripts/core/wd_crosscohort_presentation.R](../scripts/core/wd_crosscohort_presentation.R), [scripts/core/wd_fcg_comparison.R](../scripts/core/wd_fcg_comparison.R), [scripts/core/wd_intake.R](../scripts/core/wd_intake.R).


## Workflow 14

**CF453 standalone models and legacy complete presentation — current upstream model stage; full presentation alternative**

Supply current RNA/lipid models; optionally run full legacy30-panel association presentation

**Entry:** `scripts/workflows/14_cf453_standalone_figure6.R`

**Validation:** No --validate-only. Source guard defines functions; use intake tests or17 validation for no-fit checks.

**Writes / overwrite behavior:** Full/models-only call approved-root/descendant guard before mkdir/fitting, independent of ignored caches. Outside that root existing model/analysis cache guards remain. Explicit finish/plot-only skip new-fit guard and can mutate saved results; use copied/fresh run destinations.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| CF453 whole-hepatocyte pseudobulk transport<br>required | `data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json`<br>`data/RNA/pseudobulk_csv/wdsd_hepatocytes/*.csv`<br>`registry:wdsd_hepatocytes_transport → data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF453 lipid measurements<br>required | `data/metabolomics/wdsd/lipidomics_cf453.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Current CF453 zone transport<br>required | `data/RNA/pseudobulk_csv/current_wdsd_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/current_wdsd_zonation/*.csv`<br>`registry:current_wdsd_zonation_transport → data/RNA/pseudobulk_csv/current_wdsd_zonation/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Saved model/full analysis<br>required | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/models.rds`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/analysis.rds`<br>`registry:current_root → results/current` | 14<br>--finish/--plot-only | Mixed: some files are local-only; incomplete: declared paths absent<br>Ignored caches; finish needs models and no analysis, plot-only needs analysis; approved current model-only stage has no analysis.rds |

| Output | Location | Mode |
| --- | --- | --- |
| Model stage | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/models.rds`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/tables/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/provenance_sha256.csv`<br>`registry:current_root → results/current` | all fresh modes; --models-only stops here |
| Legacy full association/panel set | `<output>/analysis.rds`<br>`<output>/panels/`<br>`<output>/source_data/`<br>`<output>/panel_manifest.csv` | default/full, --finish or saved --plot-only exports |

**Command examples** (not executed by this documentation audit):

Fit complete parent cache only

```sh
Rscript scripts/workflows/14_cf453_standalone_figure6.R --models-only --output=results/runs/cf14
```

Run legacy network/panel remainder from saved models

```sh
Rscript scripts/workflows/14_cf453_standalone_figure6.R --finish --output=FRESH_MODEL_STAGE
```

**Important distinctions:**

- Whole-hepatocyte28-profile transport is a hard shared-intake dependency, even models-only.
- Primary RNA uses supplied normalized zone layer; no FCG, CF459 or MOFA reference.
- Options --input-zonation for fresh fit only, --bootstraps=999, --seed=20261030.
- Source-loading does not run pipeline; saved caches excluded from public Git.

**Code checked:** [scripts/core/cf453_standalone.R](../scripts/core/cf453_standalone.R), [scripts/core/wdsd_intake.R](../scripts/core/wdsd_intake.R), [scripts/core/pseudobulk_csv.R](../scripts/core/pseudobulk_csv.R).


## Workflow 15

**CF453 direct diet interactions — current upstream direct-evidence stage; four-panel alternative**

Extend14 cache with direct diet interactions; optional four-panel presentation

**Entry:** `scripts/workflows/15_cf453_diet_impact_figure6.R`

**Validation:** No --validate-only; full or models-only refits direct contrasts.

**Writes / overwrite behavior:** New fit requires empty output except QA/. Plot-only reads cache from output and rewrites exports.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Selected Workflow14 model cache and adjacent provenance<br>required | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/models.rds`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/provenance_sha256.csv`<br>`registry:current_root → results/current` | 14<br>fresh fit | Mixed: some files are local-only<br>Generated model RDS ignored in Git; regenerate14 models-only into fresh directory or restore exact cache; explicit --parent does not infer any other parent |
| Saved direct analysis<br>required | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_direct_diet/diet_analysis.rds`<br>`registry:current_root → results/current` | 15<br>--plot-only | Mixed: some files are local-only<br>Ignored generated cache; must restore/regenerate |

| Output | Location | Mode |
| --- | --- | --- |
| Direct model evidence | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_direct_diet/diet_analysis.rds`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_direct_diet/tables/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_direct_diet/provenance_sha256.csv`<br>`registry:current_root → results/current` | full or --models-only |
| Four-panel alternative | `<output>/panels/`<br>`<output>/source_data/`<br>`<output>/panel_manifest.csv` | full or saved --plot-only |

**Command examples** (not executed by this documentation audit):

Fit direct evidence from explicit14 parent

```sh
Rscript scripts/workflows/15_cf453_diet_impact_figure6.R --models-only --parent=results/runs/cf14/models.rds --output=results/runs/cf15
```

**Important distinctions:**

- --parent-provenance overrides adjacent14 provenance; raw-input identities are checked, not all historical code/plan rows.
- No automatic14 execution and no need for14 legacy all-edge analysis.
- Public clone must regenerate ignored14/15 caches.
- --plot-only ignores parent choice for loading: reads output/diet_analysis.rds.

**Code checked:** [scripts/core/cf453_diet_impact.R](../scripts/core/cf453_diet_impact.R), [scripts/config/analysis_release.R](../scripts/config/analysis_release.R).


## Workflow 16

**CF453 selected programmes and exploratory diet association shifts — current Figure6 presentation**

Select zonation anchors and compute fixed-programme class/network displays

**Entry:** `scripts/workflows/16_cf453_diet_rewiring_figure6.R`

**Validation:** No --validate-only. --models-only computes complete999-bootstrap/LOO statistics, not just models or validation.

**Writes / overwrite behavior:** New analysis requires empty output except QA/. Plot-only reads output cache and overwrites panels/reports.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Selected Workflow14 model cache and adjacent provenance<br>required | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/models.rds`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/provenance_sha256.csv`<br>`registry:current_root → results/current` | 14<br>fresh fit | Mixed: some files are local-only<br>Generated model RDS ignored in Git; regenerate14 models-only into fresh directory or restore exact cache; explicit --parent does not infer any other parent |
| Selected Workflow15 direct-diet evidence cache<br>required | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_direct_diet/diet_analysis.rds`<br>`registry:current_root → results/current` | 15<br>fresh analysis | Mixed: some files are local-only<br>Generated RDS ignored in Git; regenerate15 models-only using the same14 parent or restore |
| Saved selected-programme analysis<br>required | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/rewiring_analysis.rds`<br>`registry:current_root → results/current` | 16<br>--plot-only | Mixed: some files are local-only<br>Ignored cache; restore/regenerate16 |

| Output | Location | Mode |
| --- | --- | --- |
| Selected-programme statistics | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/rewiring_analysis.rds`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/tables/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/provenance_sha256.csv`<br>`registry:current_root → results/current` | full or --models-only |
| Three main panel groups | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/panels/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/source_data/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/panel_manifest.csv`<br>`registry:current_root → results/current` | full or --plot-only |

**Command examples** (not executed by this documentation audit):

Fresh complete selected-programme analysis and three panels

```sh
Rscript scripts/workflows/16_cf453_diet_rewiring_figure6.R --parent=results/runs/cf14/models.rds --direct=results/runs/cf15/diet_analysis.rds --output=results/runs/cf16
```

**Important distinctions:**

- Both explicit parent paths required when reproducing fresh chain; changing --parent does not infer sibling --direct.
- Metadata and modeled RNA must agree across parents.
- No CLI bootstrap/seed override; implementation B999/seed20261102.
- Association thresholds/bootstrap ranges are descriptive and exploratory.
- Normal parent verification checks data identities; strict historical source-snapshot exportQA separately checks mutable code/plan hashes.

**Code checked:** [scripts/core/cf453_diet_rewiring.R](../scripts/core/cf453_diet_rewiring.R), [scripts/core/cf453_diet_rewiring_plots.R](../scripts/core/cf453_diet_rewiring_plots.R), [scripts/config/analysis_release.R](../scripts/config/analysis_release.R).


## Workflow 17

**Current expanded-marker zonation and standalone diet pipeline — current main zonation/diet entry**

Reproduce current FCG zonal/KD, paired GDX and CF45314/15/16/layout stages; optional retrospective manuscript audit

**Entry:** `scripts/workflows/17_expanded_marker_zonation_update.R`

**Validation:** Entry --validate-only verifies3zone transports, balanced zones, engine paths and KD CSV/RData; no models or writes. It does NOT check all03/06/resources/whole-hepatocyte14 prerequisites. With --with-manuscript-audit it also checks Python and selected historical cache files, not all private/report dependencies.

**Writes / overwrite behavior:** Before any writes, full science/audit and report-only modes reject published current_root, legacy_results_root and all descendants including cache-absent/symlink aliases via analysis_assert_not_approved_output; ordinary completed-cache guards remain. Validation skips guard and writes nothing. Explicit report-only skips guard and rewrites saved audit/report outputs.

**Invokes other entries/phases:** `05`, `07`, `14`, `15`, `16`, `6A_layouts`. These are control-flow calls, not additional cohorts.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Three current expanded-marker zone transports<br>required | `data/RNA/pseudobulk_csv/current_fcg_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/current_gdx_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/current_wdsd_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/current_*_zonation/*.csv`<br>`registry:pseudobulk_csv_root → data/RNA/pseudobulk_csv`<br>`registry:current_fcg_zonation_transport → data/RNA/pseudobulk_csv/current_fcg_zonation/manifest.json`<br>`registry:current_gdx_zonation_transport → data/RNA/pseudobulk_csv/current_gdx_zonation/manifest.json`<br>`registry:current_wdsd_zonation_transport → data/RNA/pseudobulk_csv/current_wdsd_zonation/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF453 whole-hepatocyte pseudobulk transport<br>required | `data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json`<br>`data/RNA/pseudobulk_csv/wdsd_hepatocytes/*.csv`<br>`registry:wdsd_hepatocytes_transport → data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF453 lipid measurements<br>required | `data/metabolomics/wdsd/lipidomics_cf453.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| FCG lipid measurements<br>required | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Lipid DA and enrichment outputs<br>required | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData`<br>`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_scores.RData` | 03<br>full | Shipped in Git<br>Shipped result artifacts; actual FCG engine requires them, but17 validate-only does not preflight them |
| Lipid structural reference<br>required | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipids_with_structure.csv` | 06<br>full zonal07 structural branch | Shipped in Git<br>Shipped; needed when matched structural data exist;17 validate-only does not check it |
| Independent KD fitted effects<br>reused KD fit, not raw-experiment refitting | `results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_ttop_list.RData` | 07<br>all | Shipped in Git<br>Shipped current artifact or source input |
| TF/pathway resources<br>required | `support/TF_regulons.RData`<br>`support/Gonadal_Hormone_Genes.csv`<br>`support/m5.all.v2023.2.Mm.symbols.gmt` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen implementation engines<br>implementation dependency | `scripts/core/engines/05_zonation_analysis_legacy.R`<br>`scripts/core/engines/07_validation_experiments_legacy.R`<br>`registry:engine_fcg05 → scripts/core/engines/05_zonation_analysis_legacy.R`<br>`registry:engine_validation07 → scripts/core/engines/07_validation_experiments_legacy.R` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Retrospective FCG/CF453/GDX references<br>required | `miscellaneous/older_results/zonation_two_marker/*`<br>`miscellaneous/older_inputs/two_marker_zonation/RNA/*`<br>`miscellaneous/older_results/cf453_standalone_figure6/models.rds`<br>`miscellaneous/older_results/cf453_diet_impact_figure6/diet_analysis.rds`<br>`miscellaneous/older_results/cf453_diet_rewiring_figure6/rewiring_analysis.rds`<br>`registry:historical_fcg05 → miscellaneous/older_results/zonation_two_marker`<br>`registry:historical_fcg_zonation_csv → miscellaneous/older_inputs/two_marker_zonation/RNA/scaled_RNA_zonation.csv`<br>`registry:historical_fcg_zonation_metadata → miscellaneous/older_inputs/two_marker_zonation/RNA/meta_data_df_zonation.csv`<br>`registry:historical_gdx_zonation_rds → miscellaneous/older_inputs/two_marker_zonation/RNA/gdx_zonation_pseudobulks.rds`<br>`registry:historical_cf453_models → miscellaneous/older_results/cf453_standalone_figure6`<br>`registry:historical_cf453_direct → miscellaneous/older_results/cf453_diet_impact_figure6`<br>`registry:historical_cf453_rewiring → miscellaneous/older_results/cf453_diet_rewiring_figure6` | frozen<br>--with-manuscript-audit or --report-only | Mixed: some files are local-only<br>Mostly ignored local archive; restore exact original baselines. Old GDX RDS required for old/new supporting fit; current science does not need it. |
| Whole-hepatocyte DE comparison audit<br>denominator-definition comparison only | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/ttop_full_list.RData` | 02<br>--with-manuscript-audit | Shipped in Git<br>Shipped current artifact or source input |
| Private manuscript and supplied preprocessing evidence<br>required | `manuscript/Cleland and Dugourd et al.docx`<br>`<source-root>/james_scripts_20261104/preprocessing_{02_fcg,03_gdx,04_wdsd}.R`<br>`<source-root>/Zonation_Genes.csv`<br>`results/legacy/pre_reorganization_baseline/zonation_marker_update_20261116/inputs/preprocessing_reference/*` | external<br>audit/report modes | Local-only; restore or regenerate; also runtime/external paths<br>Private DOCX ignored; exact saved reference snapshots may satisfy completed report replay. Fresh audit source-root must also contain preprocessing evidence; default public CSV root alone does not. |

| Output | Location | Mode |
| --- | --- | --- |
| Current staged inputs and scientific result stages | `results/legacy/pre_reorganization_baseline/zonation_marker_update_20261116/inputs/`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg05/`<br>`results/current/05_sex_polarised_zonation_gene_lipid_interactions/fcg07/`<br>`results/current/04_pericentral_and_periportal_hepatocyte_bias/gdx_zonation/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_models/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_direct_diet/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_layouts/`<br>`results/runs/NAME/fcg_panel_manifest.csv`<br>`registry:current_root → results/current` | default science/full |
| Formal scientific table copies | `results/current/05_sex_polarised_zonation_gene_lipid_interactions/supplementary_tables/`<br>`registry:current_root → results/current` | science when system Python available; audit/report |
| Old/new comparisons and private revision report | `results/legacy/pre_reorganization_baseline/zonation_marker_update_20261116/comparisons/`<br>`results/legacy/pre_reorganization_baseline/zonation_marker_update_20261116/report/`<br>`registry:current_root → results/current` | --with-manuscript-audit or --report-only |

**Command examples** (not executed by this documentation audit):

No-fit limited input validation

```sh
Rscript scripts/workflows/17_expanded_marker_zonation_update.R --validate-only
```

Current science; does not run03/06 or refit rawKD

```sh
Rscript scripts/workflows/17_expanded_marker_zonation_update.R --output=results/runs/expanded_markers_reproduction
```

Saved-cache QA/table/inventory/report/bundle refresh; private prerequisites

```sh
Rscript scripts/workflows/17_expanded_marker_zonation_update.R --report-only
```

**Important distinctions:**

- Normal source-root is public pseudobulk_csv_root; staging copies zones only.14 whole-hepatocyte/lipid intake still comes from repository registry/defaults.
- Normal scientific run does not require private manuscript, original datasets or Codex runtime.
- Optional Python tables-only copies13tableCSVs; statistical source CSVs do not depend on that copy.
- Audit source-root is reused for preprocessing reference lookup: fresh audit needs a combined source folder or exact pre-staged references.
- Report-only implies audit and still requires saved model/checkpoint/comparison caches; it does not rebuild missing model stages.
- James retains nucleus-composition panels and exact executed marker-list confirmation.
- Guard-only public-clone protection was added and tested with cache-absent stub fixtures; no estimates or approved artifacts changed.
- Internal staging/GDX/core helpers retain cache guards and are not public approved-output-safe entries; maintained17 protects the chain before invoking them. Standalone maintained FCG driver also has the shared full guard.

**Code checked:** [scripts/workflows/17_expanded_marker_zonation_update.R](../scripts/workflows/17_expanded_marker_zonation_update.R), [scripts/reports/stage_zonation_marker_update.R](../scripts/reports/stage_zonation_marker_update.R), [scripts/reports/run_zonation_update_fcg.R](../scripts/reports/run_zonation_update_fcg.R), [scripts/reports/run_zonation_update_gdx.R](../scripts/reports/run_zonation_update_gdx.R), [scripts/core/zonation_update_fcg.R](../scripts/core/zonation_update_fcg.R), [scripts/core/zonation_update_gdx.R](../scripts/core/zonation_update_gdx.R), [scripts/core/engines/05_zonation_analysis_legacy.R](../scripts/core/engines/05_zonation_analysis_legacy.R), [scripts/core/engines/07_validation_experiments_legacy.R](../scripts/core/engines/07_validation_experiments_legacy.R), [scripts/reports/prepare_zonation_update_report_sources.py](../scripts/reports/prepare_zonation_update_report_sources.py), [scripts/config/analysis_release.R](../scripts/config/analysis_release.R).


## Workflow 05b

**WP2 zonation/diet phase — compatibility wrapper; historical scientific definition**

Compatibility phase of complete workflow

**Entry:** `scripts/workflows/05b_wdsd_network_comparison.R`

**Validation:** Observed intake/reference checks only; no savedWP2analysis or resampling.

**Writes / overwrite behavior:** Standalone phase rejects any existing root contents.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| CF453 whole-hepatocyte pseudobulk transport<br>required | `data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json`<br>`data/RNA/pseudobulk_csv/wdsd_hepatocytes/*.csv`<br>`registry:wdsd_hepatocytes_transport → data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF453 lipid measurements<br>required | `data/metabolomics/wdsd/lipidomics_cf453.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original CF453 zone transport<br>required | `data/RNA/pseudobulk_csv/historical_wdsd_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/historical_wdsd_zonation/*.csv`<br>`registry:historical_wdsd_zonation_transport → data/RNA/pseudobulk_csv/historical_wdsd_zonation/manifest.json` | frozen<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original FCG zone matrices<br>required | `miscellaneous/older_inputs/two_marker_zonation/RNA/scaled_RNA_zonation.csv`<br>`miscellaneous/older_inputs/two_marker_zonation/RNA/meta_data_df_zonation.csv`<br>`registry:historical_fcg_zonation_csv → miscellaneous/older_inputs/two_marker_zonation/RNA/scaled_RNA_zonation.csv`<br>`registry:historical_fcg_zonation_metadata → miscellaneous/older_inputs/two_marker_zonation/RNA/meta_data_df_zonation.csv` | frozen<br>all | Local-only; restore or regenerate<br>Retained local archive; payloads excluded from ordinary Git; restore through migration/recovery ledger |
| Four frozen Figure5 context references<br>required | `miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_portal_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_testes_central_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_Ovaries_portal_full_ATT.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_full_SIF.csv`<br>`miscellaneous/older_results/zonation_two_marker/networks/cormat_ova_central_full_ATT.csv`<br>`registry:historical_fcg05 → miscellaneous/older_results/zonation_two_marker` | frozen<br>all | Local-only; restore or regenerate<br>Retained local archive; numerical payloads mostly ignored; restore original selected SIF/full-SIF/full-ATT, not current expanded-marker networks |
| Original lipid measurements<br>required | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original FCG group/PCA roster<br>required | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | 02<br>all | Shipped in Git<br>Shipped current artifact or source input |

| Output | Location | Mode |
| --- | --- | --- |
| Zonation/diet phase | `miscellaneous/comparisons/wdsd_network_comparison/`<br>`registry:historical_wf12 → miscellaneous/comparisons/wdsd_network_comparison` | full |

**Command examples** (not executed by this documentation audit):

Phase-specific validation as documented

```sh
Rscript scripts/workflows/05b_wdsd_network_comparison.R --validate-only
```

**Important distinctions:**

- Only bootstraps/paired-bootstraps/seed/output/validate-only; no permutations or report-only option.

**Code checked:** [scripts/core/wp2_workflow.R](../scripts/core/wp2_workflow.R).


## Workflow 05c

**WP2 pooled/annotation phase — compatibility wrapper; historical scientific definition**

Compatibility phase of complete workflow

**Entry:** `scripts/workflows/05c_wdsd_lipid_overlap_and_pooled_network.R`

**Validation:** Observed intake/pooledreference/mapping checks only; no savedWP2analysis or resampling.

**Writes / overwrite behavior:** Rejects existing pooled_analysis.rds or README; partial folder is not generally fresh-guarded.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| CF453 whole-hepatocyte pseudobulk transport<br>required | `data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json`<br>`data/RNA/pseudobulk_csv/wdsd_hepatocytes/*.csv`<br>`registry:wdsd_hepatocytes_transport → data/RNA/pseudobulk_csv/wdsd_hepatocytes/manifest.json` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| CF453 lipid measurements<br>required | `data/metabolomics/wdsd/lipidomics_cf453.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original CF453 zone transport<br>required | `data/RNA/pseudobulk_csv/historical_wdsd_zonation/manifest.json`<br>`data/RNA/pseudobulk_csv/historical_wdsd_zonation/*.csv`<br>`registry:historical_wdsd_zonation_transport → data/RNA/pseudobulk_csv/historical_wdsd_zonation/manifest.json` | frozen<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen whole-hepatocyte Figure3 networks and precursor<br>required | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_full_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_full_ATT.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/gene_lipid_correlations_factor1.csv` | 04<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original raw lipid annotations<br>required | `data/metabolomics/fcg/lipidomics_fcg.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original whole-hepatocyte RNA and metadata<br>required | `data/RNA/scaled_RNA_[123].csv`<br>`data/RNA/meta_data_df.csv` | 01<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original lipid measurements<br>required | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original FCG group/PCA roster<br>required | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | 02<br>all | Shipped in Git<br>Shipped current artifact or source input |

| Output | Location | Mode |
| --- | --- | --- |
| Pooled/annotation phase | `miscellaneous/comparisons/wdsd_network_comparison/followup/`<br>`registry:historical_wf12 → miscellaneous/comparisons/wdsd_network_comparison` | full |

**Command examples** (not executed by this documentation audit):

Phase-specific validation as documented

```sh
Rscript scripts/workflows/05c_wdsd_lipid_overlap_and_pooled_network.R --validate-only
```

**Important distinctions:**

- Permutations/bootstraps/paired-bootstraps/seed/output/validate-only; no report-only.
- No FCGzone statistical reference, but sharedCF453intake still loads historicalzone transport.

**Code checked:** [scripts/core/wp2_workflow.R](../scripts/core/wp2_workflow.R).


## Workflow 07b

**WP1 base KO phase — compatibility wrapper; historical scientific definition**

Compatibility phase of complete workflow

**Entry:** `scripts/workflows/07b_cyp2c38_in_vivo_validation.R`

**Validation:** Intake and frozenreference checks; no savedWP1cache needed.

**Writes / overwrite behavior:** Rejects nonempty output.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Cyp2c38 in-vivo KO lipid input<br>required | `data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original whole-hepatocyte RNA and metadata<br>required | `data/RNA/scaled_RNA_[123].csv`<br>`data/RNA/meta_data_df.csv` | 01<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original lipid measurements<br>required | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original FCG group/PCA roster<br>required | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | 02<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original full-network Cyp2c38 selection<br>required | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | 04<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Overview workbook recorded in production hashes<br>production provenance dependency | `data/Suppl_Table_1_OverviewData.xlsx` | curated<br>full | Shipped in Git<br>Shipped current artifact or source input |

| Output | Location | Mode |
| --- | --- | --- |
| Base KO results | `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/` | full |

**Command examples** (not executed by this documentation audit):

Phase-specific validation as documented

```sh
Rscript scripts/workflows/07b_cyp2c38_in_vivo_validation.R --validate-only
```

**Important distinctions:**

- Permutations1999,bootstrap999,seed20261006,workers3; no ULMfollowup/report-only.

**Code checked:** [scripts/core/wp1_workflow.R](../scripts/core/wp1_workflow.R).


## Workflow 07c

**WP1 ULM/resampling phase — compatibility wrapper; historical scientific definition**

Compatibility phase of complete workflow

**Entry:** `scripts/workflows/07c_cyp2c38_ulm_resampling.R`

**Validation:** Reads fixed savedbase and refits/replays observedKOmodel; no resampling/writes.

**Writes / overwrite behavior:** Rejects nonempty followup output.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Cyp2c38 in-vivo KO lipid input<br>required | `data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original whole-hepatocyte RNA and metadata<br>required | `data/RNA/scaled_RNA_[123].csv`<br>`data/RNA/meta_data_df.csv` | 01<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original lipid measurements<br>required | `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | curated<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Original FCG group/PCA roster<br>required | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv` | 02<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Frozen original full-network Cyp2c38 selection<br>required | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv` | 04<br>all | Shipped in Git<br>Shipped current artifact or source input |
| Overview workbook recorded in production hashes<br>production provenance dependency | `data/Suppl_Table_1_OverviewData.xlsx` | curated<br>full | Shipped in Git<br>Shipped current artifact or source input |
| Fixed base parent outputs<br>required | `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/log2_analysis_matrix.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/sample_manifest.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/feature_annotation_and_QC_ledger.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/limma_fits.rds` | 07b<br>all | Mixed: some files are local-only<br>Required even validate-only; fit RDS ignored, must restore/regenerate |
| Previous resampling banks<br>required | `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/combined_bootstrap_draws.csv`<br>`results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/permutation_draws.csv` | 07b<br>full | Shipped in Git<br>Required full followup; base generated outputs |

| Output | Location | Mode |
| --- | --- | --- |
| ULM followup | `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/ulm_followup/` | full |

**Command examples** (not executed by this documentation audit):

Phase-specific validation as documented

```sh
Rscript scripts/workflows/07c_cyp2c38_ulm_resampling.R --validate-only
```

**Important distinctions:**

- --output changes destination only; no --parent/--previous CLI.
- Complete11 avoids standalone saved-parent validation dependency by passing raw intake.

**Code checked:** [scripts/core/wp1_workflow.R](../scripts/core/wp1_workflow.R).


## Workflow 6A_layouts

**Figure6A saved-anchor visual alternatives — current optional display entry**

Replot saved anchor effects as four-state heatmap and paired dots without reselection

**Entry:** `scripts/reports/cf453_figure6a_layouts.R`

**Validation:** No validate-only option; source-loading defines runner. Reads16anchors/effects, no models or resampling.

**Writes / overwrite behavior:** Approved-root/descendant guard runs before any writes, independent of ignored caches and symlink aliases. Default approved output is intentionally refused; supply fresh --output. Outside protected root existing review exports can be rewritten; retained-original identities checked at end.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Saved16 analysis<br>required | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/rewiring_analysis.rds`<br>`registry:current_root → results/current` | 16<br>all | Mixed: some files are local-only<br>Ignored generated cache; restore or run16 first |
| Retained-original artifacts recorded when present<br>optional retained-artifact checks | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/panels/Fig6A_anchor_gene_shifts.pdf`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/panels/Fig6B_lipid_class_shift.pdf`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/panels/Fig6C_network_gain_loss.pdf`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/source_data/*.csv`<br>`results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/CF453_Figure6_zonation_diet_panels.zip` | 16<br>all | Shipped in Git; incomplete: declared paths absent<br>Only existing files enter a new ledger; reusing oldoutput ledger can require its recorded original files |

| Output | Location | Mode |
| --- | --- | --- |
| Two layout groups, ledger and provenance | `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_layouts/`<br>`registry:current_root → results/current` | plot-only by design |

**Command examples** (not executed by this documentation audit):

Render two layouts from explicit saved16analysis

```sh
Rscript scripts/reports/cf453_figure6a_layouts.R --parent=results/runs/cf16/rewiring_analysis.rds --output=results/runs/cf6a_layouts
```

**Important distinctions:**

- Only --parent/--output options.
- No automatic16 execution; do not treat dots as repeated dietary measurements on one mouse.
- Plan/source hashes are recorded anew, not compared as retained-original data identities.

**Code checked:** [scripts/reports/cf453_figure6a_layouts.R](../scripts/reports/cf453_figure6a_layouts.R), [scripts/core/cf453_figure6a_layouts.R](../scripts/core/cf453_figure6a_layouts.R).


## Workflow MOON_EXPANDED

**Expanded-prior TF→lipid sensitivity — Active optional sensitivity; retained preferred outputs differ from runner defaults.**

Prepare chain-compatibleprior/mapping and fit explicit graph/currency/depth sensitivities; extract lower thresholds separately.

**Entry:** `scripts/pkn/audit_measured_lipid_reaction_targets.py`; `scripts/pkn/build_expanded_lipid_rules.py`; `scripts/moon/ancillary/prepare_expanded_factor1_mapping.py`; `scripts/moon/ancillary/build_currency_traversal_policy.py`; `scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R`; `scripts/moon/ancillary/reextract_expanded_factor1_thresholds.R`; `scripts/moon/ancillary/export_expanded_factor1_network.py`

**Validation:** Mapping --validate-only checks prepared files/identities without writes or scoring. Full-fit --validate-only enforces baseline software/input/graph/mapping/policy gates. Reextract --validate-only checks completed fit,source commit and grid without scoring/writes.

**Writes / overwrite behavior:** Builders/fits/mapping/exports require new empty destinations; currency policy refuses existing three output paths. Baseline PKN builder is a separate overwrite-risk command.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Base prior/resources<br>Frozen reactions,nodes,evidence,mappings | `results/causal_modelling/cosmos_mouse_pkn/`<br>`support/cosmos_mouse/raw/uniprot_mouse.tsv.gz` | PKN<br>expanded prior build | Shipped in Git<br>Retained artifacts ship; base builder does not run expandedbuilder. |
| Selected baseline seed/export<br>Selected TF seed/context table | `results/causal_modelling/factor1_tf_to_lipid/upstream_input.csv` | MOON_TF_TO_LIPID<br>expanded prior audit/build | Shipped in Git<br>Shipped savedCSV can satisfy prerequisite without refitting baseline. |
| Curated templates/inventory<br>Reaction templates and measured-target audit | `support/cosmos_mouse/expanded_species_templates.json`<br>`results/causal_modelling/lipid_reaction_target_mapping_audit/feature_reaction_target_inventory.csv` | curated<br>expanded prior build | Shipped in Git<br>Retained catalogue/inventory ship; target audit is separate preparation. |
| Frozen factor/mapping companions<br>Full precision values and one-to-one identities | `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_lipid_node_mapping.csv`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_lipid_weights_all.csv`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_model_audit.json` | frozen<br>mapping/fit | Shipped in Git<br>Included minimal archived payloads. |
| Retained expanded mapping/policy<br>Prepared observations and explicit V2 currency exclusions | `results/causal_modelling/factor1_tf_to_lipid_expanded_species/inputs/input_mapping.csv`<br>`support/cosmos_mouse/moon_currency_substrate_edges_v2.csv` | curated<br>conservative sensitivity | Shipped in Git<br>Shipped; default runner instead uses full currency and no adjacent-layer restriction. |

| Output | Location | Mode |
| --- | --- | --- |
| Expanded prior | `results/causal_modelling/cosmos_mouse_pkn_expanded_species/` | separate builder; new empty output |
| Fits/extractions | `results/causal_modelling/factor1_tf_to_lipid_expanded_species/` | default mode dirs expanded/ or frozen_control/; saved study names were explicit |
| Preferred saved sensitivity | `results/causal_modelling/factor1_tf_to_lipid_expanded_species/currency_v2_expanded_layers4_lower_thresholds/` | saved separate lower-threshold extraction |

**Command examples** (not executed by this documentation audit):

No-write mapping check using retained expanded graph/resources.

```sh
python3 scripts/moon/ancillary/prepare_expanded_factor1_mapping.py --validate-only
```

Full-fit gate; currentdecoupleR mismatch blocks.

```sh
Rscript scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R --mode=expanded --validate-only
```

Check saved fit/commit and lower-cutoff grid; no scores or writes.

```sh
Rscript scripts/moon/ancillary/reextract_expanded_factor1_thresholds.R --fit=COMPLETED_FIT_DIR --output=NEW_DIRECTORY --layer-consistent --validate-only
```

**Important distinctions:**

- Order:baseprior+selectedseeds→targetaudit/catalogue→expandedprior→mapping→optionalV2policy→fit→lower-thresholdreextraction→evidence/routeexport.
- Preferredconservativebehavior requires explicit no_generic_substrate, V2edgeCSV and --layer-consistent; not a default invocation.
- --alias-complete alone still targets existingV1policy filenames; supply three fresh outputpaths if regeneratingV2.
- Default fit outputsubdirs expanded/frozen_control differ from saved expanded_inputs/frozen_inputs_control and namedcurrency runs.

**Code checked:** [scripts/pkn/build_expanded_lipid_rules.py](../scripts/pkn/build_expanded_lipid_rules.py), [scripts/moon/ancillary/prepare_expanded_factor1_mapping.py](../scripts/moon/ancillary/prepare_expanded_factor1_mapping.py), [scripts/moon/ancillary/build_currency_traversal_policy.py](../scripts/moon/ancillary/build_currency_traversal_policy.py), [scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R](../scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R), [scripts/moon/ancillary/reextract_expanded_factor1_thresholds.R](../scripts/moon/ancillary/reextract_expanded_factor1_thresholds.R), [scripts/moon/config/factor1_tf_to_lipid_expanded_species.json](../scripts/moon/config/factor1_tf_to_lipid_expanded_species.json).


## Workflow MOON_KO

**Reviewed KO validation proposal and optional expression audit — Optional no-fit evidence/proposal branch; raw-count regeneration needs local inputs.**

Export reviewed target/path priorities separately from actual-count feasibility evidence.

**Entry:** `scripts/moon/ancillary/audit_factor1_ko_expression.R`; `scripts/moon/ancillary/propose_factor1_ko_validation.py`

**Validation:** Expression --validate-only performs its actual count/metadata checks and still needs local raw count inputs; proposal --validate-only checks reviewed paths/graph/source/expression without exporting. Neither fitsMOON.

**Writes / overwrite behavior:** Fresh proposal/expression output required; existing nonempty destinations refused.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Expanded selected networks/mapping<br>Four/six-layer frozen networks and mapped endpoint identities | `results/causal_modelling/factor1_tf_to_lipid_expanded_species/` | MOON_EXPANDED<br>proposal/expression context | Shipped in Git<br>Retained selected outputs ship. |
| Review and source metadata<br>Reviewed exclusions,hypotheses,ranking and source identities | `support/cosmos_mouse/factor1_ko_validation_review.json`<br>`support/cosmos_mouse/ko_validation_literature_sources.json` | curated<br>proposal | Shipped in Git<br>Shipped; no graph mutation or automatic causal intervention inference. |
| Original count evidence<br>Actual whole-cell FCG counts and historical culture controls; factor values remain separate | `data/RNA/fcg_pseudobulks.rds`<br>`miscellaneous/older_inputs/legacy_rna/merged_counts_matrix_250508_VH00211_433_AACVFM3HV_with_gene_symbols_and_groups.csv` | external<br>expression-audit regeneration | Local-only; restore or regenerate<br>FCGSeuratobject andcultureCSV arelocal-only; seven publictransports do not replace whole-cell FCG counts. |
| Retained expression audit<br>Reusable count-derived feasibility table | `results/causal_modelling/factor1_ko_validation_targets/expression_audit/gene_expression.csv` | frozen<br>proposal only | Shipped in Git<br>Shipped; proposal need not rerun raw-count audit. |
| Frozen factor annotations<br>FullprecisionFactor1values shown separately from actualcounts | `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds` | frozen<br>expression audit | Shipped in Git<br>Included compact active-model dependency. |

| Output | Location | Mode |
| --- | --- | --- |
| Prioritization/path exports | `results/causal_modelling/factor1_ko_validation_targets/alternating_prioritization/` | default existing; fresh output |
| Expression audit | See mode notes | default temporary or explicit fresh output |

**Command examples** (not executed by this documentation audit):

No-write checks; requires local original FCG/culture counts.

```sh
Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R --validate-only
```

Validate retained proposal inputs without model fit/output.

```sh
python3 scripts/moon/ancillary/propose_factor1_ko_validation.py --validate-only
```

**Important distinctions:**

- This proposal consumes expanded sensitivities, not solely the baseline388edge graph.
- Archived proposalPDF is local-only; frozen READMEoldlink is not a public download.
- Raw count audit and reviewedproposal exports are different rerun scopes.

**Code checked:** [scripts/moon/ancillary/audit_factor1_ko_expression.R](../scripts/moon/ancillary/audit_factor1_ko_expression.R), [scripts/moon/ancillary/propose_factor1_ko_validation.py](../scripts/moon/ancillary/propose_factor1_ko_validation.py), [support/cosmos_mouse/factor1_ko_validation_review.json](../support/cosmos_mouse/factor1_ko_validation_review.json).


## Workflow MOON_LIPID_TO_TF

**Exact-mapped lipid→TF hypotheses — Active optional branch; main entry produces coarsegrid, not every retained refinement/export.**

Fit original-direction regulatory/precursor feasibility policies; optionalph73 remains separate.

**Entry:** `scripts/moon/run_factor1_lipid_to_tf.R`; `scripts/moon/ancillary/refine_lipid_tf_thresholds.R`; `scripts/moon/ancillary/audit_lipid_tf_scores.R`

**Validation:** Main --validate-only checks frozen software/observations before coverage fitting; it does not run every later ancillary/export check. Refinement reuses scores; score audit computes diagnostic scores on saved graphs (not a new coherence fit).

**Writes / overwrite behavior:** Main/refinement/selected export require fresh directories; no automatic overwrite.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Frozen observations/contract<br>Exact lipid inputs,TF observations,RNA coherence,softwarepin | `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_lipid_node_mapping.csv`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_hepatocyte_tf_activities_mapped.csv`<br>`scripts/moon/config/factor1_tf_to_lipid.json` | frozen<br>fit/coverage | Shipped in Git<br>Included; forward config used as dependency contract, not a forward-run call. |
| Original and optional normalized priors<br>Original PKN plus --include-ph73 curated sensitivity | `results/causal_modelling/cosmos_mouse_pkn/` | PKN<br>fit/coverage/export | Shipped in Git<br>Retained resources/evidence ship. |

| Output | Location | Mode |
| --- | --- | --- |
| Coarse fitting/extraction | `results/causal_modelling/factor1_lipid_to_tf/` | default existing; fresh --output |
| Refinement/export | `results/causal_modelling/factor1_lipid_to_tf/refinement/` | separate ancillary extraction from saved fits |
| Score diagnostics | `results/causal_modelling/factor1_lipid_to_tf/score_audit/` | separate recomputation on final graphs |

**Command examples** (not executed by this documentation audit):

Input/software check; current pin blocks.

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf.R --validate-only --include-ph73
```

Fit explicit original+normalized sensitivities after setup.

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf.R --include-ph73 --output=results/runs/moon_lipid_to_tf
```

Refine saved coarsefits and materialize reviewed candidate; no refit.

```sh
Rscript scripts/moon/ancillary/refine_lipid_tf_thresholds.R --results=results/runs/moon_lipid_to_tf
```

**Important distinctions:**

- Config: scripts/moon/config/factor1_lipid_to_tf.json.
- No PKN edge reversal. Main does not automatically perform refinement or final-score audit.
- Omitting --include-ph73 preserves original-prior negative policies; later refinement is not needed for those no-route results.

**Code checked:** [scripts/moon/run_factor1_lipid_to_tf.R](../scripts/moon/run_factor1_lipid_to_tf.R), [scripts/moon/config/factor1_lipid_to_tf.json](../scripts/moon/config/factor1_lipid_to_tf.json), [scripts/moon/ancillary/refine_lipid_tf_thresholds.R](../scripts/moon/ancillary/refine_lipid_tf_thresholds.R), [scripts/moon/ancillary/audit_lipid_tf_scores.R](../scripts/moon/ancillary/audit_lipid_tf_scores.R).


## Workflow MOON_RELAXED

**Relaxed regulatory lipid→TF sensitivity — Active optional sensitivity; distinct from original-prior lipid→TF and selected baseline.**

Separate direct-regulatory conditional mapping/scoring analysis using directional STITCH.

**Entry:** `scripts/moon/run_factor1_lipid_to_tf_relaxed.R`; `scripts/moon/ancillary/search_factor1_lipid_to_tf_relaxed_parameters.R`

**Validation:** --validate-only checks exact pinned packages/source hashes,selected reference and direct-graph inputs without fit/writes; currentdecoupleR pin fails. Candidate/evidence audits are separate commands.

**Writes / overwrite behavior:** Selected and search outputs refuse existing nonempty directories.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Directional prior<br>Direct chemical→protein +protein prior | `results/causal_modelling/cosmos_mouse_pkn_stitch_directional/pkn_exploratory_stitch_directional.rds` | PKN<br>selected fit/search | Shipped in Git<br>Retained overlay ships; large rawSTITCH files are rebuild-only prerequisites. |
| Frozen factor/mapping/reference<br>Conditional feature/chemical assignments and selected saved comparator | `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds`<br>`results/causal_modelling/factor1_relaxed_lipid_regulatory_mapping/`<br>`miscellaneous/older_results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory/runs/direct_strong_tfs_pruning4_layers4/run_state.rds`<br>`miscellaneous/older_results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory/score_layer_constrained/solutions/direct_strong_tfs_pruning4_layers4__p0.5__s0.5/` | frozen<br>fit/export/reference | Mixed: some files are local-only<br>Required compact states,ledgers and current evidence ship. |

| Output | Location | Mode |
| --- | --- | --- |
| Selected sensitivity | `results/causal_modelling/factor1_lipid_to_tf_relaxed/` | default existing; fresh output |
| Parameter search | `results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search/` | separate search entry/config |

**Command examples** (not executed by this documentation audit):

No-write prerequisites; currentdecoupleR pin blocks.

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R --validate-only
```

New selected sensitivity fit after environment provisioning.

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf_relaxed.R --output=results/runs/moon_relaxed
```

**Important distinctions:**

- Config: scripts/moon/config/factor1_lipid_to_tf_relaxed.json; searchconfigseparate.
- Unlike the stale older README claim, this active configuration DOES use the STITCH overlay.
- Candidate chemical topology seeds and post-score abundance-sign annotations are not receptor perturbations or independent TF validation.

**Code checked:** [scripts/moon/run_factor1_lipid_to_tf_relaxed.R](../scripts/moon/run_factor1_lipid_to_tf_relaxed.R), [scripts/moon/lib/relaxed_lipid_tf.R](../scripts/moon/lib/relaxed_lipid_tf.R), [scripts/moon/config/factor1_lipid_to_tf_relaxed.json](../scripts/moon/config/factor1_lipid_to_tf_relaxed.json), [scripts/moon/config/factor1_lipid_to_tf_relaxed_parameter_search.json](../scripts/moon/config/factor1_lipid_to_tf_relaxed_parameter_search.json).


## Workflow MOON_TF_TO_LIPID

**Selected TF→lipid replay — Active optional branch; saved result/reference checks available; current full-fit software gate fails.**

Refit selected baseline from frozen mapped factor inputs and compare it with the original selected result.

**Entry:** `scripts/moon/run_factor1_tf_to_lipid.R`; `scripts/moon/ancillary/validate_frozen_archives.R`

**Validation:** Main --validate-only enforces exact package/commit/input contracts, reads prerequisites and stops before scoring/writes. Saved-artifact validator instead re-extracts from retained main+relaxed fits; default prints only, optional --output=NEW_JSON writes report.

**Writes / overwrite behavior:** Nonempty fit output refused; saved-artifact validator refuses existing report filename.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Full-precision factor observations<br>TF activities,lipid weights,RNA coherence and mapping | `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_lipid_node_mapping.csv` | frozen<br>fit/reference replay | Shipped in Git<br>Mandatory compact archived states are included via active_frozen_dependencies.json; not regenerated by17. |
| Base PKN and provenance<br>Signed prior plus node/edge/reaction evidence | `results/causal_modelling/cosmos_mouse_pkn/` | PKN<br>fit/export | Shipped in Git<br>Retained resource ships; no builder call by main runner. |
| Selected reference states<br>Exact saved scores,compression,graph and selectedsolution | `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth4_compressed/`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/topology_radius_6/compressed/compression.rds`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03/network_SIF.csv`<br>`miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/score_retention/radius6_depth4_compressed_scores_all.csv` | frozen<br>reference comparison | Mixed: some files are local-only<br>Minimal required files ship; remaining exploratory archive payloads can be ignored/local. |

| Output | Location | Mode |
| --- | --- | --- |
| Selected fit,exports | `results/causal_modelling/factor1_tf_to_lipid/` | default already populated; fresh --output required for rerun |

**Command examples** (not executed by this documentation audit):

No-write full-fit prerequisites; currently blocked at decoupleR pin.

```sh
Rscript scripts/moon/run_factor1_tf_to_lipid.R --validate-only
```

Saved hashes/reference/extraction checks; no fit or new scoring.

```sh
Rscript scripts/moon/ancillary/validate_frozen_archives.R
```

New fit/export only after pinned environment is provisioned.

```sh
Rscript scripts/moon/run_factor1_tf_to_lipid.R --output=results/runs/moon_selected_replay
```

**Important distinctions:**

- Config: scripts/moon/config/factor1_tf_to_lipid.json.
- Requires cosmosR1.21.2 commitb2ed195738479c4ad204cace37fd11ae2fd02bb8,decoupleR2.17.0,igraph2.3.3,jsonlite2.0.0.
- Installed metadata18Nov2026 hasdecoupleR2.9.7; this blocks full-fit validation without invalidating saved artifact checks.
- Does not run04,PKN builders or17; no expanded-marker RNA dependency.

**Code checked:** [scripts/moon/run_factor1_tf_to_lipid.R](../scripts/moon/run_factor1_tf_to_lipid.R), [scripts/moon/lib/workflow.R](../scripts/moon/lib/workflow.R), [scripts/moon/config/factor1_tf_to_lipid.json](../scripts/moon/config/factor1_tf_to_lipid.json), [scripts/moon/ancillary/validate_frozen_archives.R](../scripts/moon/ancillary/validate_frozen_archives.R), [docs/provenance/active_frozen_dependencies.json](../docs/provenance/active_frozen_dependencies.json).


## Workflow PKN

**Base PKN and separate directional overlay — Optional resource generation; retained priors/evidence ship; full source rebuild has additional packages/raw prerequisites.**

Build/audit signed prior resources separately from omic scoring.

**Entry:** `scripts/pkn/run_build.py`; `scripts/pkn/build_stitch_directional.py`; `scripts/pkn/validate_export.R`

**Validation:** validate_export.R --validate-only checks existing signed CSV graphs against installed COSMOS cleaner without conversion/writes. It is not a MOON/software-version full-fit gate.

**Writes / overwrite behavior:** run_build.py can overwrite its fixed base directory; no --validate-only/help parser. Overlay builder requires empty/new output. Default validate_export writes RDS and validation/session records.

| Required input / role | File or registry key | Producer / mode | Availability and action |
| --- | --- | --- | --- |
| Frozen source snapshots<br>Mouse-GEM,Rhea,regulation,SwissLipids,LINEX identity/reaction resources | `support/cosmos_mouse/raw/`<br>`results/causal_modelling/cosmos_pkn_audit/raw/` | frozen<br>base build | Mixed: some files are local-only<br>Retained base snapshots ship; do not infer that every optional upstream download ships. |
| Study feature inputs<br>Measured lipid DA/catalogue | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gonade_type_ttop.csv` | 03<br>base mapping | Shipped in Git<br>Shipped saved output. |
| MOFA selected-feature context<br>Factor1 feature selection context | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/gene_lipid_correlations_factor1.csv` | 04<br>base mapping | Shipped in Git<br>Shipped saved output; builder does not invoke04. |
| Large directional STITCH downloads<br>Raw acting-partner/action/link evidence | `support/cosmos_mouse/raw/stitch_actions_10090_v5.0.tsv.gz`<br>`support/cosmos_mouse/raw/stitch_links_detailed_10090_v5.0.tsv.gz` | external<br>overlay rebuild only | Local-only; restore or regenerate<br>Local-only/ignored; source URLs+hashes retained. Existing overlay consumers do not need rawdownloads. |

| Output | Location | Mode |
| --- | --- | --- |
| Base resource | `results/causal_modelling/cosmos_mouse_pkn/` | base generation; fixed destination |
| Audit cache | `results/causal_modelling/cosmos_pkn_audit/` | retained source evidence |
| Directional overlay | `results/causal_modelling/cosmos_mouse_pkn_stitch_directional/` | separate builder; new empty output |

**Command examples** (not executed by this documentation audit):

No-fit/no-write validation of retained base CSV exports.

```sh
Rscript scripts/pkn/validate_export.R --validate-only
```

Full resource rebuild in fixed directory; never use as a no-write preflight.

```sh
python3 scripts/pkn/run_build.py
```

Separate overlay rebuild; requires local rawSTITCH files.

```sh
python3 scripts/pkn/build_stitch_directional.py --output=NEW_EMPTY_DIRECTORY
```

**Important distinctions:**

- Base order: build_mouse_gem→build_rhea_lipid_reactions→build_lipid_regulation→map_lipids→build_lipid_rules→assemble_mouse_pkn→audit_mouse_pkn→validate_export.R→write_build_manifest.
- Base runner does not call overlay/expanded builders or fitMOON.
- Python pygoslin/pinned dependencies and RcosmosR/jsonlite are separate setup requirements.

**Code checked:** [scripts/pkn/run_build.py](../scripts/pkn/run_build.py), [scripts/pkn/map_lipids.py](../scripts/pkn/map_lipids.py), [scripts/pkn/build_stitch_directional.py](../scripts/pkn/build_stitch_directional.py), [scripts/pkn/validate_export.R](../scripts/pkn/validate_export.R).
