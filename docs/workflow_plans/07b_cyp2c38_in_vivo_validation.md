# Workflow 07b — Cyp2c38 in vivo lipidomic validation (WP1)

## Current compatibility entry (18 November 2026)

This wrapper calls only `wp1_run()` in `scripts/core/wp1_workflow.R`.
The complete public work package is [Workflow11](11_cyp2c38_in_vivo_validation.md).
Accepted phase options are `--permutations` (1999), `--bootstraps` (999),
`--seed` (20261006), `--workers` (3), `--output` and `--validate-only`.
Default output is `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/`; a full phase rejects a
nonempty directory. Validation reconstructs intake/reference coefficients without
saved WP1 output or writes. No report-only or ULM/calibration phase runs here.
The exact original whole-cell RNA, cleaned lipids, Workflow02 group/PCA export
and Workflow04 full-network SIF prerequisites are listed in Workflow11; those
are public supplied artifacts. The original execution record below is retained
as historical evidence, not a new result caused by cleanup or CSV distribution.

**Entry point:** `scripts/workflows/07b_cyp2c38_in_vivo_validation.R`.
**Helpers:** `scripts/core/cyp2c38_reference.R`, `scripts/core/cyp2c38_models.R`.
Plot and report helpers: `scripts/core/cyp2c38_plots.R`, `scripts/core/cyp2c38_report.R`.
**Scope:** October 5, 2026 meeting WP1; heterozygous in vivo Cyp2c38 KO lipidomics. This is a separate extension of Workflow 07's primary-hepatocyte siRNA validation.

**Reporting follow-up:** Workflow 07c reads these completed outputs and adds
full separate gonadal networks, native ULM scores, nominal ULM p-values and
explicit null-bootstrap p-values. Use its `ulm_followup/README.md` for the
current requested concordance reporting. This plan and the original outputs
retain the historical fixed-network analysis.

## Fixed specification, chosen before fitting KO outcomes

1. Read the integrated acquisition-mode export `data/metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv`. Preserve every annotation and a unique row ID. Exclude 16 explicit internal standards, any isotope-labelled feature, and the exact classes already excluded in Workflow 07 (BA, DGGA, HBMP, LPE-N, NATau, PEtOH, PMeOH, SL, CAR, SPB, VAE). Analytical QC columns are not animals.
2. Use log2 supplied positive intensities without another normalization or imputation. Normalization provenance, batch/litter and gap-filling details are unavailable. Retain all curated endogenous features as primary; QC CV <=30% is a predetermined sensitivity, with QC CV calculated across nine analytical QC injections.
3. Four group means: F_WT (8), F_KO (8), M_WT (7), M_KO (9). Equal-weight pooled effect is half the female plus male KO-minus-WT effects. Also estimate female, male and exploratory male-minus-female interaction. No averaging of t-statistics.
4. Joint limma model with sex-group precision weights (`arrayWeights(var.group=sex)`), intensity-trend and robust empirical Bayes. This addresses the broad sex-dependent residual scale observed before outcome fitting. Save weights and mean–variance diagnostics. Unweighted joint, separate-sex fits, QC-restricted and supplied-scale fits are prespecified sensitivities; weighted and unweighted group-mean logFC must agree.
5. Use the fixed 828 Cyp2c38 edges in `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv`. Reconstruct pooled coefficients from the original paired scaled hepatocyte RNA and log2 lipidomics (43 FCG mice), stopping if discrepancy exceeds 1e-12. Recompute within-gonad coefficients using the same feature set (22 testes, 21 ovaries). No MOFA retraining, no KO RNA, no response/FDR-driven feature selection.
6. Match by exact molecular annotation after last-pipe/class-parenthesis formatting only. Preserve chain order, slash/underscore and structural resolution. Reject duplicate joins and exclude unresolved CL reorder candidates.
7. Main concordance family: pooled effects versus pooled reference, female effects versus ovary reference, male effects versus testes reference. Signed concordance is `-cor(reference_r, KO_t)`; per-lipid consistency is `-reference_r * KO_t`. Positive values support concordance. Sex-specific effects versus pooled reference are secondary descriptions. Preserve the existing signed ULM convention: negative score supports concordance.
8. Calibrate main concordance using genotype permutations within sex, refitting weights and moderation; two-sided empirical p-values with +1 correction, BH across the three main tests. Assess combined uncertainty by independently resampling KO animals within four groups and paired FCG animals within four original genotype/gonad strata. Fixed reference selection is conditional and not bootstrapped. Save every aggregate resampling draw and invalid-draw counts. Leave-one-KO-animal-out checks quantify influence. Lipid-edge regression/ULM p-values are descriptive because species are dependent.

## Dependencies and execution

Run from the repository root. Source `scripts/support_functions.R` first. A project-local `.local/R-library` is supported; required packages include limma, statmod, ggplot2, decoupleR, digest and jsonlite. Local installation must never replace the user's system library. The existing Workflow 07 and its outputs remain independent.

The entry point supports seed/permutation/bootstrap/worker/output CLI options and writes only `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/` by default. `--validate-only` checks inputs and archive replay without writing results. Read `_meta_cross_workflow_overview.md` and the existing 03, 04 and 07 plans before changing this extension.

Original FCG chromosome/gonad groups are read from the standard Workflow 02
`results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/PCA_sample_coord_RNA_HEPA.csv`, independently checked against
`data/Suppl_Table_1_OverviewData.xlsx`: XXO=10, XYO=11, XXT=11, XYT=11. The paired
reference bootstrap retains these four group counts.

## Outputs and checks

Full differential tables and fits; complete feature/sample/QC/identity ledgers; all fixed-reference coefficients and mouse manifest; concordance/ULM/consistency tables; permutation, combined-bootstrap and leave-one-out summaries; QC, differential and predictive-validation PDFs; input SHA256/package/session provenance; result report.

Differential tables include marginal moderated 95% effect intervals (`CI.L`,
`CI.R`) and `effect_units`. Limma's `logFC` field represents a supplied-intensity
mean difference in the linear-scale sensitivity; all log2 models use log2 fold
changes. BH FDR is computed independently within each contrast. The primary
t-statistic matrix is exported explicitly.

`fixed_reference_and_matching_ledger.csv` retains two unmatched CL reorder
candidates (Alignment IDs 8059 and 8073); the second has two reference labels in
the same unordered equivalence class. Neither enters primary matching.

## Historical production verification, 2026-10-06

The full entry point and `tests/test_cyp2c38_in_vivo.R` passed. Defaults are
1,999 within-sex KO label permutations and 999 combined whole-animal bootstrap
draws, seed 20261006, three workers. All draws yielded finite aggregate scores.
The original 828 coefficients replayed within 5.56e-16; equal-sex pooled and
interaction identities and weighted/unweighted effect agreement passed.

There are zero discoveries at FDR <=0.05 in all four primary contrasts and
all prespecified sensitivity models. The primary pooled/female/male Pearson
associations are -0.06270, -0.06369 and -0.16654; their animal-level permutation
p-values are 0.848, 0.722 and 0.6215 (three-test BH FDR 0.848 throughout).
Combined bootstrap intervals span zero, and KO leave-one-mouse-out score
ranges cross zero in all three comparisons. Thus this dataset does not provide
robust in vivo validation; it does not establish absence of a biological effect.
Use `results/current/03_gonadal_sex_biased_regulatory_circuitry/cyp2c38_in_vivo/README.md` and machine-readable tables
for full interpretation. The historical siRNA result remains distinct.

Required invariants: finite log2 values, unique sample/feature keys, 8/8/7/9 animals, rank-four model, exact contrast signs, pooled and interaction identities, weighted/unweighted logFC agreement, reference numerical replay, one-to-one fixed-universe joins and resampling at whole-animal level. Test helper scripts immediately after edits and run the complete entry point immediately after creation or modification. Update this plan after verification.

## Scientific limitations

Heterozygous perturbation is not a complete deletion. Reference correlations and KO effects cannot identify direct substrates or flux. Pooled FCG associations can reflect between-gonad covariance; sex-matched reference tests estimate a different relationship and may be weaker. Permutation inference requires genotype exchangeability within sex; unrecorded family/batch structure cannot be controlled. Within-gonad correlations have 21–22 mice and bootstrap intervals are conditional on the archived feature selection. Supplied intensities appear scaled, but normalization/imputation details are unverified.
