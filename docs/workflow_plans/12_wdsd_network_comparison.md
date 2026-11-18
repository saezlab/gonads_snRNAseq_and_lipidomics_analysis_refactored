# Workflow 12 — WDSD pooled and zonation network comparison (WP2)

## CSV transport update - 18 November 2026

Both phases now explicitly request `historical_wdsd_zonation_transport`; the shared
intake defaults its whole-hepatocyte input to `wdsd_hepatocytes_transport`. These
lossless exports preserve the historical source RDS identities without opening
Seurat objects. All3 RNA layers, metadata, pairing, QC/marker summaries and stored
command/version provenance matched the previous RDS intake exactly. Explicit local
RDS overrides remain supported but public defaults choose CSV even while originals
exist. Each manifest and every payload file enters input provenance.

---

## Frozen comparison routing - 16 November 2026

Workflow12 remains a historical Figure3/Figure5 comparison and explicitly passes
the registry's `historical_wdsd_zonation_rds` into both analysis phases. Original
FCG zonation matrices and networks use `historical_fcg_zonation_*` / `historical_fcg05`.
The default archived output root is `historical_wf12`; fresh runs still require a
fresh `--output`. No expanded-marker reference is substituted into this comparison.
The exact original78-gene,708-displayed-edge reconstruction and coefficient replay
pass after physical archive relocation. Report-only input links also resolve to
the frozen input role. Current CF453 zonation science runs through Workflow17.

---

**Entry point:** `scripts/workflows/12_wdsd_network_comparison.R`.
**Orchestration:** `scripts/core/wp2_workflow.R`.
**Report:** `scripts/core/wp2_report.R`.
**Consolidated:** 9 October 2026. Active complete WP2 workflow; former 05c
(pooled/annotation) and 05b (zonation/diet) entries remain thin compatibility
wrappers. Their plans retain detailed scientific specifications and original
production results/QA. Read both, this plan and the meta overview before changes.

Organization and data paths changed; statistical implementations and primary
numerical outputs were preserved. CF459/WP3, nucleus-level processing and
TF-gene replication remain outside WP2.

## Sections and execution order

1. **Configuration:** source `scripts/support_functions.R`, enable existing
   `.local/R-library`, load `wp2_workflow.R`.
2. **Pooled Figure 3/annotation:** `wp2_pooled_run()` reads fresh WDSD intake,
   reconstructs frozen whole-hepatocyte networks and calls
   `wp2_load_lipid_overlap_audit()` directly on local raw/cleaned annotations.
   Pooled comparisons and sensitivity/candidate ledgers go to `followup/`.
3. **Sex/zonation and diet:** `wp2_run(allow_pooled=TRUE)` reconstructs four
   Figure 5 contexts and runs SD conservation/descriptive WD-SD comparisons.
   At this stage only the just-created `followup/` may already exist.
4. **Combined report:** `wp2_write_package_report()` composes pooled-first
   Markdown from CSV/RDS, archives the original zonation README once and validates
   local links. `workflow_manifest.json` records phase order/arguments/completion.

`source()` defines functions without executing analysis; `Rscript` calls
`wp2_package_run()`. Run from repository root. Pooled analysis does not depend
on earlier WP2 outputs or the zonation phase.

## Local inputs and frozen references

| Repository input | Role |
|---|---|
| `data/metabolomics/wdsd/lipidomics_cf453.csv` | CF453 lipidomics, controls and raw annotations |
| `data/RNA/wdsd/wdsd_hep_pseudobulks.rds` | Updated whole-hepatocyte Seurat pseudobulks |
| `data/RNA/wdsd/wdsd_zonation_pseudobulks.rds` | Updated paired zone pseudobulks |
| `data/metabolomics/fcg/lipidomics_fcg.csv` | Original raw annotation audit |
| `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | Original cleaned 1,254-feature matrix |
| `data/RNA/scaled_RNA_{1,2,3}.csv`, `meta_data_df.csv` | Original whole-hepatocyte FCG profiles |
| `data/RNA/scaled_RNA_zonation.csv`, `meta_data_df_zonation.csv` | Original zonation profiles |
| `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv`, overview workbook | Validated original FCG identities |
| Workflow 04 sparse species Figure 3 SIF and precursor | Frozen displayed/precursor pairs |
| `results/zonation/networks/cormat_*_{SIF,full_SIF}.csv` | Four frozen Figure 5 branches |
| `docs/reference_figures/FCGL_Fig3.ai` | Optional selection provenance; not required for numerical replay |

See `data/WP1_WP2_INPUTS.md` and `data/wp1_wp2_input_migration_manifest.csv`
for original source paths/verified hashes. Input bytes/processing are unchanged;
shared archives preserved. No active helper reads sibling `original_datasets`
or `Revised_Figures`. Only two frozen input RDSs receive gitignore exceptions.

Require 28 distinct cross-mapped CF453 mice, seven per sex/diet group, 56
profiles with two zones per mouse and 14 pooled SD animals. Original paired FCG
reference has 43 mice, 21 ovarian/22 testicular. Analytical QC is not an animal;
per-mouse bulk lipid abundance is shared across its RNA zones. Preserve feature
IDs, raw annotations, molecular resolution and full pairing/exclusion ledgers.

## Scientific specification and coverage

Original coefficients replay from scaled RNA and log2 lipids. No MOFA refit or
selection on SD/WD responses. Primary molecular matching remains strict,
one-to-one: 752/1,254 original lipids. Recorded-alias coverage (959) and qualified
conditional compatibility (854) are distinct evidence levels. Fresh annotation
audit reads source inputs directly, saving candidates, ambiguity and exclusions.

Figure 3 has 192 displayed species pairs (25 genes/64 lipids) and a 139-by-64
precursor. Strict SD comparison retains 108 displayed pairs, 36 lipids and 15
displayed genes; precursor has 5,004 pairs. Native ULM/nominal p accompany
unrestricted pairing/null-bootstrap tests and distinct within-sex upper-tail
tests. Conditional tests ask whether alignment exceeds that expected from shared
sex structure. Centered networks remain descriptive. QC-CV, Fill and conditional
mapping sensitivities preserve frozen selection. See phase 05c for full details.

Four sex-zone contexts retain 434 displayed/58,656 full eligible pairs.
Conservation uses native ULM, exact within-sex SD pairing and independent
empirical-marginal null bootstrap. Whole-mouse indices are shared across features
and zones. Ordinary paired bootstrap describes stability. WD-SD delta/sign
changes and all 49 joint leave-one-SD/one-WD-out comparisons are descriptive.
Preserve predefined Mtmr11/Cyp examples. Dependent edge p-values do not treat
lipid pairs as independent mice. See phase 05b for complete assumptions.

## CLI, defaults and guards

```bash
Rscript scripts/workflows/12_wdsd_network_comparison.R --validate-only
Rscript scripts/workflows/12_wdsd_network_comparison.R --output=results/wdsd_network_comparison_reproduction
Rscript scripts/workflows/12_wdsd_network_comparison.R --report-only
python3 scripts/reports/build_wp2_report_pdf.py --input results/wdsd_network_comparison/README.md
```

Defaults: seed 20261007, pooled pairing draws 19,999, independent-marginal
bootstrap 4,999, paired stability 999; zonal pairing always all 7!=5,040.
Options: `--permutations`, `--bootstraps`, `--paired-bootstraps`, `--seed`,
`--output`, `--validate-only`, `--report-only`. Shared counts/seed route to both
phases; `--permutations` affects only pooled Monte Carlo pairings. Unknown and
duplicate options, and simultaneous validation/report modes, are rejected.

Validation checks both actual input/reference phases without saved WP2 analysis
or writes. Full runs reject any nonempty root; choose a fresh destination.
`--report-only` explicitly refreshes combined Markdown from completed root and
`followup/` results without recalculation. PDF generation is a separate Python
step using reportlab/pypdf and configured fonts, not executed by the R analysis.

## Outputs and reporting order

Default root is `results/wdsd_network_comparison/`. Root retains zonal baseline,
diet, edge/LOCO, QC and sensitivity CSV/RDS, network SIFs and plots. `followup/`
retains pooled inferences, mapping/annotation ledgers, complete contexts,
shared null/stability draws and plots. Core statistical helpers retain
`wdsd_*` names in `scripts/core/`.

Main README and `WP2_WDSD_network_comparison.pdf` present:

1. Pooled Figure 3 conservation, coverage, unconditional/conditional nulls,
   centered sensitivity and predefined Cyp examples.
2. Four sex-specific zonation networks and diet comparisons, including Mtmr11.
3. Lipid annotation overlap/resolution.
4. Shared cohort/QC and normalization-layer checks.
5. Interpretation/limitations, artifact links and reproduction.

PDF vector appendix order is pooled, zonation/diet, overlap/QC. Component
`followup/README.md` remains available; `README_zonation_original_2026-10-07.md`
archives the preceding root narrative. Historical Meeting_Summaries deliveries
are preserved. Combined prose is regenerated from machine-readable results;
links and PDF text/tables/vector figures/layout are checked before delivery.

## Verified consolidation and limitations

Actual local validation, source loading, fresh intake, CLI/default routing,
pooled-first order, report refresh and output guards passed independent checks.
A complete fresh default-count run finished both analysis phases and combined
report. All 12 primary inference/coefficient/sensitivity tables reproduced within
1e-10, including empirical p-values. Sixteen primary resampling draw/index/LOCO
banks across pooled and zonation RDSs were exactly identical to original results.
Primary repository numeric files remain unchanged. Evidence is under
`tests/outputs/wp1_wp2_migration/`; repeat actual comparisons with
`Rscript tests/test_work_package_workflows.R --wp2-output=/path/to/fresh/output`.

Focused tests cover local input hashes, intake, matching, reference replay,
pooled/zonal statistics and orchestration. Original phase plans record independent
OLS, null tail/index and scientific figure QA. Original historical hashes remain
historical records; fresh runs record local input/code paths and package manifest.
Run every modified R script immediately and update plans after changes.

Inherited limits: n=7 per sex/diet, unknown nucleus contributions/processing and
batch/litter provenance, dependent edges, conditional empirical nulls, and bulk
rather than zonal lipid measurement. Pooled conservation is strong but much is
explained by shared sex structure; displayed conditional and zonal animal-level
tests are inconclusive. Descriptive diet changes alone do not establish rewiring
or causal regulation. Reorganization does not strengthen the inferential claims.
