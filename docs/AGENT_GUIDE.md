# Result lookup guide

Use [the release registry](../scripts/config/analysis_release.json) as the source of current defaults. This guide replaces the previous static map, preserved under `miscellaneous/documentation/pre_cleanup_20261116/docs/AGENT_GUIDE.md`.

## Current main result map

| Question | File or folder | Producer / interpretation |
| --- | --- | --- |
| Which release is current? | `scripts/config/analysis_release.json` | `current_root`; distinct historical keys |
| Which updated inputs were used? | `docs/provenance/current_input_manifest.csv` | Public CSV/JSON payload hashes; original source-object hashes remain in export manifests |
| Where are the current figure panels? | [Panel navigation](../results/current_panels.md) | Current 4B/5/S5A/6 and source tables |
| Which manuscript numbers need changing? | `results/zonation_marker_update_20261116/report/new_replacements.csv` | 89 claims plus24 pictured-edge coefficients; source manuscript unchanged |
| Main cell-type RNA DE | `results/pseudobulk/ttop_full_list.RData` | Workflow02; separate from zonation-model sex effects |
| Main lipid DE | `results/lipidomic/gonade_type_ttop.csv` | Workflow03; use actual tables rather than old prose counts |
| MOFA weights / global networks | `results/mofa/main_model_RNA_lipid_weight.RData`, `results/mofa/networks/` | Workflow04, selected trained model set |
| Current FCG zonation DE | current root `fcg05/ttop_full_long.csv`, `t_table_zonation.csv` | Workflow17/FCG adapter; all seven contrasts |
| Selected spatial programmes | current root `fcg05/sex_specific_zonation_genes_summary.csv` | Reoptimized cutoff4.4; fixed5.4 sensitivity is separately audited |
| Gene–lipid networks | current root `fcg05/networks/` | Four gonad × zone contexts; SIF/ATT sources |
| Lipid-class associations | current root `fcg05/panel_source_data/` | Faithful historical butterfly and explicit gene-context view have different sex encodings |
| Mtmr11 zonal validation | current root `fcg07/Fig5F_Mtmr11_male_PP_prespecified.pdf`, `Mtmr11_male_PP_prespecified_source.csv` | MalePP,22reference animals; native ULM661KD targets;345matched scatter points |
| Updated formal source tables | current root `tables/table_manifest.csv` | S10/S13/S14/S15; independent S11/S12 unchanged |
| Standalone diet models | current root `cf453_models/`, `cf453_direct_diet/` | Workflow14/15 model caches; seven mice per sex/diet |
| Current Figure6 | current root `cf453_rewiring/panels/`, `source_data/` | Workflow16; class and edge associations exploratory |
| Figure6A alternatives | current root `cf453_layouts/panels/`, `source_data/` | Heatmap and paired dots from the same new gene set |
| GDX supporting result | current root `gdx_zonation/` | Paired sham/GDX model; not a replay of the mismatched hormone-replacement notebook |
| KD experimental measurements | `results/validation/knockdown_t_table.csv`, `knockdown_ttop_list.RData` | Preserved experimental data |
| In-vivo Cyp2c38 validation | `results/validation/cyp2c38_in_vivo/` | Workflow11; independent whole-hepatocyte reference |
| Lipid-only structural / variance models | `results/lipidomic/structural_model/`, `variance_analysis/` | Workflows06/08 |
| Independent spatial / human results | `results/xenium/`, `results/human/` | Separate assays and optional external-data analyses |
| Active causal results | [Causal analysis map](CAUSAL_ANALYSIS_MAP.md) | Selected MOON models, optional sensitivities and frozen priors |

Here **current root** is `results/zonation_marker_update_20261116/`. Large fitted caches and raw archival payloads may be local-only; see [the commit guide](GITHUB_COMMIT_GUIDE.md).

## Historical lookup

Do not search the generic old `results/zonation` path for current findings. Original two-marker FCG networks now live at `miscellaneous/older_results/zonation_two_marker/`. Earlier CF453 Figure6 analyses are sibling archived folders. WF12/13 comparison outputs are under `miscellaneous/comparisons/`, retaining their original input definition.

[Organization/provenance](REPOSITORY_ORGANIZATION.md) maps every move. Use the release resolver when reading frozen absolute or old relative paths; do not replace their scientific meaning with the current release.

## Minimal current lookup

```r
source("scripts/config/analysis_release.R")
current <- analysis_path("current_root")
zone <- read.csv(file.path(current, "fcg05", "ttop_full_long.csv"))
subset(zone, gene == "Cyp2e1")

selected <- read.csv(file.path(current, "fcg05", "sex_specific_zonation_genes_summary.csv"))
head(selected)
```

Commands, prerequisites and fit definitions are in [PIPELINE.md](../PIPELINE.md) and the corresponding [workflow plans](workflow_plans/). A saved-result lookup does not require model refitting.
