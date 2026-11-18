# Workflow 05b - WDSD baseline conservation and Western diet perturbation (WP2)

**Compatibility phase documentation — 9 October 2026:** the active complete
work package is [12](12_wdsd_network_comparison.md), entered through
`scripts/workflows/12_wdsd_network_comparison.R`. This former entry is a thin wrapper for its
zonation/diet phase; orchestration now lives in `scripts/core/wp2_workflow.R`.
Inputs are repository-local, with byte-preserving import provenance in
`data/wp1_wp2_input_migration_manifest.csv`. The recorded production results
below remain unchanged. The combined WP2 report now presents pooled results
first, then zonation/diet, as authorized on 9 October. Earlier separate-report
instructions below document the original deliveries.

**Entry point:** `scripts/workflows/05b_wdsd_network_comparison.R`.
**Authorized scope:** Work package 2 of the 5 October 2026 meeting plan: four
original Figure 5 versus WDSD standard-diet network comparisons, followed by
four within-WDSD Western versus standard-diet comparisons. CF459 cross-cohort
robustness is WP3. Optional composition/Xenium and TF-gene replication are
outside the defined WP2 network analysis.

Read the meeting plan, Workflow 05 plan and cross-workflow overview before
editing. Reuse existing results as a frozen reference; preserve Workflow 05,
WP1, and unrelated working-tree changes. Write new results only to
`results/wdsd_network_comparison/` (or a user-specified fresh directory).

## Inputs and pairing

- WDSD CF453 lipidomics:
  `data/metabolomics/wdsd/lipidomics_cf453.csv`.
- Updated WDSD pseudobulks:
  `data/RNA/wdsd/wdsd_hep_pseudobulks.rds` and
  `wdsd_zonation_pseudobulks.rds`.
- Frozen original four context networks:
  `results/zonation/networks/cormat_{testes_portal,testes_central,Ovaries_portal,ova_central}_{SIF,full_SIF}.csv`.
- Original source matrices:
  `data/RNA/scaled_RNA_zonation.csv`, `meta_data_df_zonation.csv`,
  `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` and FCG group metadata.

Require 28 distinct CF453 mice, seven in each female/male x SD/WD group,
56 zone profiles with exactly two different zones per mouse, and exact
RNA-lipid alignment by `cf453_id`. Each mouse contributes one bulk-liver lipid
profile shared by its two RNA zones. Zones are not independent animals.
Use Seurat v5 named layer accessors; preserve supplied objects unchanged.
Record library totals/detected genes and marker consistency, but do not infer
nucleus counts or completed nucleus-level QC from aggregated pseudobulks.

## Frozen reference and preprocessing

Reconstruct all original full gene-lipid coefficients independently and compare
with the archives. Validate source branches using RNA source identity and
original lipid targets, since a gene that is also a regulon source can be
labelled TF in attributes. Exclude TF-gene branches. Original source coefficients
use scaled RNA and log2 lipid intensity; raw lipid intensities do not replay.
Original reference mice: 22 testes-bearing or 21 ovary-bearing in each zone.
Original context-specific gene and edge selection remains fixed.

Primary comparisons use the original displayed standard-network gene-lipid
edges (absolute original correlation >0.65), with all eligible full-network
edges as a complementary analysis. Never impose the original correlation
cutoff again on WDSD observations. Eligibility requires exact gene identity,
unique molecular lipid match and finite, varying features in both SD and WD.
Keep excluded/missing/invariant edges in a complete audit ledger.

Import WDSD lipid intensities and analytical controls separately. Give every
annotation row a unique feature ID; Alignment ID alone is not a key. Exclude
explicit standards, isotope labels and the exact repository-excluded classes
used by WP1. Preserve chain order, slash/underscore, ether/oxygenation and
structural resolution in formatting-only molecular matching. Primary abundances
are log2 supplied positive intensities, with no additional normalization or
imputation. QC CV<=30% and an alternative supplied RNA data layer are
prespecified sensitivities. After input inspection, before association analysis,
the supplied Fill %>=0.5 restriction was added as a separate QC sensitivity.
Do not remove outcome-driven outliers.

## Stage A - Standard-diet conservation

Map testes to recorded male and ovaries to recorded female WDSD mice within
Periportal/Pericentral contexts. Compute original versus SD coefficient-vector
Pearson association and native intercept ULM slope t/p as descriptive summaries.
Nominal edge-level p-values do not make thousands of dependent edges biological
replicates. Use the same measured edge universe in both stages.

The primary animal-level check exhaustively permutes all 7!=5,040 RNA-lipid
mouse pairings within each sex's SD group, sharing each pairing across both
zones and all features. Its null is absence of RNA-lipid pairing information,
stronger than simply zero alignment with the reference network. Report the
exact two-sided probability, including the observed identity permutation.

An independent empirical-marginal null bootstrap independently resamples RNA
mouse and lipid mouse indices (4,999 draws by default), sharing indices across
zones/features; it tests an explicit no-association null conditional on fitted
empirical margins. Use two-sided absolute ULM tails with the +1 correction.
An ordinary paired-mouse bootstrap (999 draws) describes conservation stability
and supplies no null p-value. Undefined constant-feature draws invalidate the
affected fixed-universe summary, rather than silently changing its edge set.
Report valid/invalid counts and conditional inference limitations at n=7.

## Stage B - Western diet comparison

For the same predefined pairs compute r_SD, r_WD and delta_r=r_WD-r_SD,
absolute change and sign reversal. Export a coefficient-vector scatter and
ranked descriptive edge changes in all four sex-zone contexts. Check all
7x7=49 joint leave-one-SD/one-WD mouse-out combinations, dropping the same
mouse across its two zones and all features. Record delta direction and
sign-reversal stability, ranges and per-omission ID ledgers. Large delta_r
alone does not establish reproducible rewiring. No formal differential-correlation
p-value is assigned: a naive diet-label permutation would also test changes in
the full joint distributions, including abundance, rather than specifically
equality of correlations.

Mtmr11 is a predefined illustration where measurable. Label any additional
changing-edge figure selection as descriptive. Summarize changes by gene and
lipid class while retaining individual-species tables. Correlations of zonal RNA
with bulk liver lipids are not direct measurements of zonal lipid abundance,
enzyme substrates, flux or causal regulation.

## Outputs and execution checks

Helpers: `scripts/core/wdsd_intake.R`, `wdsd_reference.R`,
`wdsd_comparison_setup.R`, `wdsd_network_statistics.R` and
`wdsd_presentation.R`. The entry point sources
`scripts/support_functions.R` first, records CLI settings/seeds, input/code
SHA256 and session versions, and saves intake/QC, pairing/molecular/edge
ledgers, complete coefficients, resampling draws/indices and LOCO summaries.
Create baseline/diet/selected-edge/QC scientific figure PDFs, a first Figure 6
draft and self-contained Markdown/PDF report with copies in Meeting_Summaries.

Run every new/modified R script immediately; use explicit namespaces. Tests
must independently check pairing, exact archived Pearson replay, strict joins,
shared mouse resampling, ULM equivalence, permutation identity/counts and
correct delta orientation. Update this plan after complete execution with actual
coverage, results and unresolved QC/provenance limitations. No outcomes are
assumed from the hoped-for conservation discussed in the meeting.

## Verified Execution - 7 October 2026

The validation-only and full entry point passed with seed 20261007, 4,999
independent-marginal null bootstraps, 999 paired/paired-reference stability
draws and all 5,040 exact pairings. Calculations use one deterministic process.
Repeated-draw row normalization was vectorized and verified numerically
identical to the original calculation, including extreme and constant inputs.

Intake: 28 paired mice (seven per sex/diet), 56 RNA zone profiles and 20,123
genes. Of 1,240 supplied lipid rows, 14 standards and 49 repository-class
rows are excluded, leaving 1,177; 1,129 satisfy QC CV<=30%. Median CV is
5.54%, maximum 67.28%. The supplied Fill % field has 89 retained rows below
0.5 despite positive exported values. All aggregate marker-order checks pass;
nucleus contributions, upstream normalization/imputation, batch and litter
provenance remain unavailable. The stored data layer is an affine replay of
scale.data, not independent normalization validation.

All 97,812 original full gene-lipid coefficients replay from scaled RNA and
log2 lipids with maximum error 8.88e-16. TF-gene branches are excluded while
Pitx3's measured gene-lipid branch is retained. Strict molecular overlap is
752 of 1,254 original lipids; 502 are unmatched. All original genes are
measured and varying in both diets. No nonfinite/invariant edges are excluded.

| Context | Full eligible / original | Displayed eligible / original | Baseline vector r | Native ULM | Nominal p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 18,800 / 31,350 | 86 / 155 | 0.2619 | 2.4870 | 0.01486 | 0.32262 | 0.31360 |
| Male pericentral | 22,560 / 37,620 | 249 / 387 | -0.02786 | -0.43797 | 0.66179 | 0.88750 | 0.88564 |
| Female periportal | 9,024 / 15,048 | 52 / 89 | -0.06307 | -0.44687 | 0.65690 | 0.76925 | 0.77346 |
| Female pericentral | 8,272 / 13,794 | 47 / 77 | 0.31110 | 2.19589 | 0.03330 | 0.27202 | 0.26460 |

Both stages therefore contain 58,656 full/434 displayed pairs. Full-network
baseline vector r ranges -0.0515 to +0.0157. No primary animal-level baseline
test supports positive conservation; QC-CV and fill sensitivities do not change
that conclusion. Permutation and bootstrap results agree closely on displayed
pairs, without evidence here of bootstrap-specific excessive conservatism.
This is inconclusive transfer evidence in this cohort, not proof of absence or
a verified low-power explanation. Neither bootstrap size nor power is calibrated.
All exact permutations and actual LOCO draws are defined. The largest invalid
bootstrap fraction is 121/4,999 (2.42%, male pericentral full), chiefly driven
by Krt4's repeated SD values; probabilities condition on valid draws.

Displayed-pair median absolute WD-SD changes are 0.358, 0.373, 0.646 and
0.795 in the above order. Of 434 displayed pairs, 56 retain sign reversal in
all 49 omissions. Full-universe medians are 0.405, 0.372, 0.426 and 0.430;
the apparent excess female change depends on the displayed selection. These
are descriptive within-cohort observations, not differential-correlation tests.

Mtmr11 has 752 full/seven displayed eligible pairs in male periportal. Its
displayed reference coefficients are all positive in the narrow 0.654-0.720
range; five SD coefficients stay positive. Vector r=-0.660 therefore describes
relative-order inversion, not universal biological sign reversal; full-vector
r=0.0203. Four displayed PE pairs retain negative delta in all 49 omissions.
Only PE(O-24:2_22:6) retains actual sign reversal throughout, but its Fill %=0.357
and it is removed by the fill sensitivity despite CV=3.80%. Three other stable
negative-delta PE pairs have Fill %=0.905. Flag this annotation limitation in
manuscript illustrations and follow-up prioritization.

The entry point exports 16 full/displayed SD/WD SIF tables, complete input and
edge ledgers, all null/stability draws and mouse indices, all 49 LOCO matrices
and omissions, gene/class summaries, seven scientific figure PDFs and the
first Figure 6 draft. The PDF builder is
`scripts/reports/build_wdsd_wp2_pdf.py`; the report and scientific figure
appendix are saved in this result directory and copied with Markdown to
`../Meeting_Summaries/2026-10-07_WP2_WDSD_network_conservation_and_diet.*`.
Reruns require a fresh `--output` directory; existing result directories are
never silently overwritten by the entry point.

Focused permanent tests cover intake, archival replay, matching, native
ULM/decoupleR equivalence, resampling, group strata and LOCO. Independent actual
output QA verifies all coefficients (maximum discrepancy 7.8e-15), all 16
exports, unchanged original Git inputs, earlier shared FCG SHA256 values and
all 19 primary/secondary exact-permutation identities and null p-value tails.
Scientific and report PDF QA is recorded under `results/wdsd_network_comparison/QA/`.
