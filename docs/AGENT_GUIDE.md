# Result lookup guide

Use [the release registry](../scripts/config/analysis_release.json) as the source of current defaults. This guide replaces the previous static map, preserved under `miscellaneous/documentation/pre_cleanup_20261116/docs/AGENT_GUIDE.md`.

For input dependencies, producers, execution modes and fresh-checkout requirements, use the [complete workflow reference](WORKFLOW_REFERENCE.md). Its [CSV dependency inventory](workflow_dependencies.csv) can be filtered by consumer or producer; the [JSON contract](workflow_dependencies.json) is its editable source. This result guide alone is not an execution recipe.

## Current main result map

| Question | File or folder | Producer / interpretation |
| --- | --- | --- |
| Which release is current? | `scripts/config/analysis_release.json` | `current_root`; distinct historical keys |
| Which updated inputs were used? | `docs/provenance/current_input_manifest.csv` | Public CSV/JSON payload hashes; original source-object hashes remain in export manifests |
| Where are the current figure panels? | [Panel navigation](../results/current_panels.md) | Current 3J–K, 4C, 5, S5A and standalone Figure 6 panels |
| Which manuscript numbers need changing? | `results/manuscript/manuscript_updates.md` | Local revision audit, excluded from Git |
| Main cell-type RNA DE | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/ttop_full_list.RData` | Workflow02; separate from zonation-model sex effects |
| Main lipid DE | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gonade_type_ttop.csv` | Workflow03; use actual tables rather than old prose counts |
| MOFA weights / global networks | `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/main_model_RNA_lipid_weight.RData`, `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/` | Workflow04, selected trained model set |
| Current FCG zonation DE | `analysis_current()` · `fcg05/ttop_full_long.csv`, `t_table_zonation.csv` | Workflow17/FCG adapter; all seven contrasts |
| Selected spatial programmes | `analysis_current()` · `fcg05/sex_specific_zonation_genes_summary.csv` | Reoptimized cutoff4.4; fixed5.4 sensitivity is separately audited |
| Gene–lipid networks | `analysis_current()` · `fcg05/networks/` | Four gonad × zone contexts; SIF/ATT sources |
| Lipid-class associations | `analysis_current()` · `fcg05/panel_source_data/` | Faithful historical butterfly and explicit gene-context view have different sex encodings |
| Mtmr11 zonal validation | `analysis_current()` · `fcg07/Fig5F_Mtmr11_male_PP_prespecified.pdf`, `Mtmr11_male_PP_prespecified_source.csv` | MalePP,22reference animals; native ULM661KD targets;345matched scatter points |
| Updated formal source tables | [Sections 04–05](../results/README.md) | S10/S13/S14/S15 in the supplementary table directories |
| Standalone diet models | `analysis_current()` · `cf453_models/`, `cf453_direct_diet/` | Workflow14/15 model caches; seven mice per sex/diet |
| Current Figure6 | `analysis_current()` · `cf453_rewiring/panels/`, `source_data/` | Workflow16; class and edge associations exploratory |
| Figure6B alternatives | `analysis_current()` · `cf453_layouts/panels/`, `source_data/` | Heatmap and paired dots from the same new gene set |
| GDX supporting result | `analysis_current()` · `gdx_zonation/` | Paired sham/GDX model; not a replay of the mismatched hormone-replacement notebook |
| KD experimental measurements | `results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv`, `knockdown_ttop_list.RData` | Preserved experimental data |
| In-vivo Cyp2c38 validation | `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/` | Workflow11; independent whole-hepatocyte reference |
| Lipid-only structural / variance models | `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/`, `variance_analysis/` | Workflows06/08 |
| Independent spatial / human results | `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/xenium/`, `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/human/` | Separate assays and optional external-data analyses |
| Active causal results | [Causal analysis map](CAUSAL_ANALYSIS_MAP.md) | Selected MOON models, optional sensitivities and frozen priors |

Here **current root** is `results/current/`. Use `analysis_current()` with a component string to resolve its section path. Large fitted caches and raw archival payloads may be local-only; see [the commit guide](GITHUB_COMMIT_GUIDE.md).

## Historical lookup

Do not search the generic old `results/zonation` path for current findings. Original two-marker FCG networks now live at `miscellaneous/older_results/zonation_two_marker/`. Earlier CF453 Figure6 analyses are sibling archived folders. WF12/13 comparison outputs are under `miscellaneous/comparisons/`, retaining their original input definition.

[Organization/provenance](REPOSITORY_ORGANIZATION.md) maps every move. Use the release resolver when reading frozen absolute or old relative paths; do not replace their scientific meaning with the current release.

## Minimal current lookup

```r
source("scripts/config/analysis_release.R")
zone <- read.csv(analysis_current("fcg05/ttop_full_long.csv"))
subset(zone, gene == "Cyp2e1")

selected <- read.csv(analysis_current("fcg05/sex_specific_zonation_genes_summary.csv"))
head(selected)
```

Commands, prerequisites and fit definitions are in [PIPELINE.md](../PIPELINE.md) and the corresponding [workflow plans](workflow_plans/). A saved-result lookup does not require model refitting.

## Validation and provenance scope

A saved table or PDF can be inspected without its ignored model cache. Re-rendering or chaining analyses can require that cache and its adjacent provenance files; the complete workflow reference identifies these dependencies. A shipped upstream result can be consumed without rerunning its producer.

Validation flags differ by workflow:17 checks transport/engine/KD prerequisites without fitting, while11’s `--validate-only` fits the observed weighted model in memory. Several older entries have no validation CLI and run on source. Consult the documented mode before executing any entry.

Some historical export QA checks compare every recorded code/plan hash with the current working tree. After a documented source or plan edit, that strict source-snapshot comparison can fail even while the frozen numerical results and ordinary cache-parent checks remain valid. Preserve the original manifests and inspect the stated mismatch; do not rewrite their hashes to conceal it. The repository organization check separately verifies preserved approved results and inputs.
