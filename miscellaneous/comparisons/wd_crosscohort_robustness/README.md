# WP3 - CF459 WD comparison with original FCG SD and CF453 WD

Report date: 14 October 2026. Amended comparison of the original male FCG SD cohort, CF453 WD and independent CF459 WD.

## 1. Requested comparison: CF459 versus original FCG and CF453

The new question is whether the CF459 WD network is closer to the original paper's male SD network than to the CF453 WD network. **All three correlations use exactly the same frozen, strictly shared 610 lipid identities and original context-specific genes**: 15,250 periportal and 18,300 pericentral full pairs. Original displayed pairs remain secondary (69 and 202). No cohort gets a larger feature universe in this comparison.

Ages clarified by the user on 14 October 2026: CF459 male WD and original male FCG SD mice were 3 months old; CF453 male WD mice were 10 months old. These are ages at sampling, not established diet durations.

The original male cohort contains **22 testes-bearing FCG mice: 11 XXT and 11 XYT**. This is the cohort used to construct the original male Figure 5 networks; it is not the five additional wild-type males. CF453 has seven WD males and CF459 has four. The original FCG SD coefficients are recalculated from paired original RNA and log2-lipid matrices and replay the archived network before matching.

Define **delta r = correlation(CF459, original FCG SD) - correlation(CF459, CF453 WD)**. Positive delta means CF459 is closer to the original FCG SD coefficient vector; negative delta means it is closer to CF453 WD. This compares coefficient-vector alignment, not equality of gene-lipid correlations or absolute molecular abundance.

| Comparison | Pairs | Original FCG / CF459 r | CF453 / CF459 r | Original FCG / CF453 r | Delta r |
| --- | --- | --- | --- | --- | --- |
| Male periportal full | 15250 | 0.052 | -0.063 | 0.016 | 0.115 |
| Male periportal displayed | 69 | 0.220 | 0.032 | 0.102 | 0.189 |
| Male pericentral full | 18300 | 0.080 | -0.019 | 0.113 | 0.098 |
| Male pericentral displayed | 202 | 0.098 | -0.076 | 0.424 | 0.174 |

**Male periportal: CF459 is closer to original FCG SD in the observed full network (delta r = 0.115). Male pericentral: CF459 is closer to original FCG SD in the observed full network (delta r = 0.098).**

Absolute full-network alignment with original FCG remains weak despite being higher than with CF453. This is an observed, cohort-specific comparison. Original feature selection, different sample sizes and unmeasured cohort/processing differences can affect alignment. Even a positive delta does not identify an age effect or show that diet is unimportant. Age, diet, genetic background, diet duration and batch are not independently varied here; there is no within-CF459 SD arm.

[Matched-anchor coefficient figure](plots/FCG_age_diet_comparison.pdf) shows both full-network anchors beside one another, with CF459 coefficients on the same y-axis and exactly the same pairs in each zone.

### 1.1 ULM and whole-mouse evidence for the original-FCG anchor

The native ULM below regresses CF459 coefficients on the original FCG SD coefficients with an intercept. Its nominal edge-level p-value is descriptive because edges share mice, genes and lipids. The exact null enumerates the same 24 CF459 RNA-lipid pairings as the original WP3 comparison, now conditional on the fixed original-FCG anchor. The approximate null bootstrap reuses the same complete CF459 RNA and lipid index banks; no smaller p-value is selected across anchors.

| Comparison | ULM score | Nominal ULM p | Exact CF459 p (24) | Null bootstrap p (defined only) |
| --- | --- | --- | --- | --- |
| Male periportal full | 6.369 | 1.95e-10 | 0.0417 | 0.0594 |
| Male periportal displayed | 1.849 | 0.0688 | 0.3750 | 0.1882 |
| Male pericentral full | 10.802 | 4.06e-27 | 0.0417 | 0.0014 |
| Male pericentral displayed | 1.393 | 0.1651 | 0.5833 | 0.5502 |

A numerical-degeneracy audit identifies **233 of 4999 CF459 null index configurations** whose exact integer source-identity contingency makes every feature covariance mathematically zero. Tiny floating-point residuals cannot supply a defined coefficient-vector correlation. The new original-FCG-anchor calculation therefore records these vectors as undefined, without using an arbitrary coefficient cutoff. Its conditional denominator can differ from the original WP3 tables retained below. The [degeneracy ledger](fcg_comparison/null_degeneracy_ledger.csv) records every affected specification and historical/raw value; observed correlations and delta are unaffected by this null-only correction.

For transparency, the same mathematical-zero correction also gives the following audit of the original CF453-anchor null bootstrap. Source files and the historical figures are preserved; the corrected values below supersede the historical bootstrap p-values for interpretation. This correction does not alter the observed WD-versus-WD correlations or exact pairing tests.

| Comparison | Historical bootstrap p | Corrected bootstrap p | Original defined draws | Corrected defined draws |
| --- | --- | --- | --- | --- |
| Male periportal full | 0.0678 | 0.0686 | 4819 | 4766 |
| Male periportal displayed | 0.8361 | 0.8368 | 4819 | 4766 |
| Male pericentral full | 0.7656 | 0.7680 | 2153 | 2146 |
| Male pericentral displayed | 0.5439 | 0.5456 | 3084 | 3067 |

Complete corrected draws, completion bounds and predefined-gene checks are saved in the [historical-null correction audit](fcg_comparison/historical_null_correction_summary.csv). The next table gives the corrected original-FCG-anchor denominator and the separate paired-delta stability counts.

| Comparison | Null defined / requested | Null completion range | Delta defined / requested | Positive delta among defined |
| --- | --- | --- | --- | --- |
| Male periportal full | 4766 / 4999 | 0.0566 to 0.1032 | 972 / 999 | 92.3% |
| Male periportal displayed | 4766 / 4999 | 0.1794 to 0.2260 | 980 / 999 | 75.0% |
| Male pericentral full | 2146 / 4999 | 6.00e-04 to 0.5712 | 391 / 999 | 90.3% |
| Male pericentral displayed | 3067 / 4999 | 0.3376 to 0.7240 | 590 / 999 | 89.2% |

The closeness comparison uses **999 three-cohort paired whole-mouse bootstrap draws**. Each draw resamples original FCG RNA/lipids together within the 11 XXT and 11 XYT strata; CF453 and FCG are independently sampled, and the **same CF459 paired draw** enters both anchor correlations. Both zones also share indices. The positive-delta fraction describes stability conditional on defined draws; **it is not a null p-value, posterior probability or probability of an age effect**.

The original WP3 whole-mouse index banks are reused; the additional independent FCG stratified paired bank uses recorded seed 20365743.
The full frozen vector is retained in every draw. If a required gene or lipid is constant, or the complete coefficient vector is constant, the affected comparison and delta are undefined. No draw-specific edge deletion is used. Null completion ranges assign every undefined null draw to either tail extreme; these are mathematical missing-draw bounds, not confidence intervals. The n=4 CF459 limitation therefore remains visible in the new comparison.

| Zone | Omissions | Defined / planned | Minimum delta r | Maximum delta r | Positive delta among defined |
| --- | --- | --- | --- | --- | --- |
| Male periportal | joint_omission | 28 / 28 | 0.014 | 0.137 | 100.0% |
| Male periportal | CF453_only | 7 / 7 | 0.095 | 0.133 | 100.0% |
| Male periportal | CF459_only | 4 / 4 | 0.036 | 0.112 | 100.0% |
| Male periportal | FCG_only | 22 / 22 | 0.075 | 0.145 | 100.0% |
| Male pericentral | joint_omission | 14 / 28 | 0.009 | 0.111 | 100.0% |
| Male pericentral | CF453_only | 7 / 7 | 0.078 | 0.126 | 100.0% |
| Male pericentral | CF459_only | 2 / 4 | 0.023 | 0.097 | 100.0% |
| Male pericentral | FCG_only | 22 / 22 | 0.058 | 0.115 | 100.0% |

The omission audit preserves the original 28 joint WD omissions, seven CF453-only and four CF459-only omissions, and adds all 22 original-FCG single-mouse omissions. These 61 overlapping configurations are descriptive influence checks, not independent experiments. [Closeness stability figure](plots/FCG_closeness_stability.pdf) displays full-network deltas and undefined-draw counts. The [extension inference](fcg_comparison/inference.csv), [paired-draw ledger](fcg_comparison/paired_draws.csv) and [omission summary](fcg_comparison/omission_summary.csv) retain the complete full/displayed and predefined-gene results.

### 1.2 Prespecified Mtmr11 does not mirror the full-network pattern

The full-network comparison does not imply that every gene is closer to the original FCG network. The predefined periportal Mtmr11 example has a negative observed delta in both full and displayed universes; its paired delta direction is unstable. This is retained alongside the global network result rather than selecting only concordant genes.

| Example | Pairs | Original FCG / CF459 r | CF453 / CF459 r | Delta r | Positive delta among defined |
| --- | --- | --- | --- | --- | --- |
| Mtmr11 full | 610 | 0.017 | 0.063 | -0.046 | 51.1% |
| Mtmr11 displayed | 6 | -0.576 | -0.151 | -0.425 | 46.2% |

| Comparison | ULM score | Nominal ULM p | Exact CF459 p (24) | Null bootstrap p (defined only) |
| --- | --- | --- | --- | --- |
| Mtmr11 full | 0.415 | 0.6785 | 0.8750 | 0.8427 |
| Mtmr11 displayed | -1.409 | 0.2316 | 0.1667 | 0.1284 |

## 2. Original WP3 scope and observed WD network agreement

This secondary robustness analysis compares gene-lipid correlation coefficients in **seven CF453 male WD mice and four CF459 male WD mice**, separately in periportal and pericentral hepatocytes. The primary universe is every measurable common pair from the fixed original FCG full networks; the original displayed Figure 5 pairs are a secondary subset. Selection is frozen before examining either new cohort. No female network or within-CF459 diet contrast is included.

The observed full-network coefficient-vector correlations are **r = -0.063** in periportal and **r = -0.019** in pericentral hepatocytes. These quantify alignment of the two edge-coefficient vectors, not correlation of independent biological replicates. Positive values denote similarly ordered coefficients; negative values denote opposing ordering.

**Neither primary network shows substantial positive observed alignment. This comparison does not demonstrate robust positive cross-cohort conservation in the measured male WD features.** The near-zero coefficients, four-mouse resolution and unknown cohort differences limit interpretation; they do not establish absence of shared regulation.

| Comparison | Pairs | Vector r | ULM score | Same coefficient sign |
| --- | --- | --- | --- | --- |
| Male periportal full | 15250 | -0.063 | -7.807 | 47.9% |
| Male periportal displayed | 69 | 0.032 | 0.259 | 44.9% |
| Male pericentral full | 18300 | -0.019 | -2.513 | 48.7% |
| Male pericentral displayed | 202 | -0.076 | -1.073 | 46.5% |

The native ULM score is the intercept-OLS slope t-statistic for CF459 coefficients regressed on CF453 coefficients. Its sign follows coefficient-vector r. Thousands of overlapping gene-lipid edges share mice, genes and lipids; therefore its nominal p-value is descriptive and does not turn edge count into biological replication. No FDR, signed concordance score or confidence interval is reported.

[Network coefficient figure](plots/network_conservation.pdf) shows full networks first and original displayed subsets second. Axes use CF453 WD coefficients on x and CF459 WD coefficients on y, each within [-1, 1]. The dashed identity line is a visual reference; coefficient equality is not tested.

## 3. Common-feature coverage and exact pairing

| Zone | Original full pairs | Shared full pairs | Genes | Lipids | CF453 / CF459 mice |
| --- | --- | --- | --- | --- | --- |
| Male periportal | 31350 | 15250 | 25 | 610 | 7 / 4 |
| Male pericentral | 37620 | 18300 | 30 | 610 | 7 / 4 |

| Zone | Original displayed pairs | Shared displayed pairs | Genes | Lipids |
| --- | --- | --- | --- | --- |
| Male periportal | 155 | 69 | 18 | 62 |
| Male pericentral | 387 | 202 | 15 | 146 |

Strict annotation matching retains **610 shared eligible FCG lipids** across the two cohorts. Unique source feature IDs remain separate; matching uses the final pipe annotation and class formatting only. Chain order, slash/underscore, ether prefixes and oxygenation are preserved. Alignment IDs are never cross-cohort identifiers and no alias or coarser resolution is forced into the primary analysis.

The full [lipid crosswalk](lipid_crosswalk.csv) and [edge eligibility ledger](edge_eligibility_ledger.csv) document detected/excluded, missing, ambiguous, nonfinite and invariant features. No SD variance requirement is imposed on this WD-only comparison. Both zones use the same ordered mouse IDs within each cohort; each animal supplies two RNA profiles but one bulk liver lipid profile.

As a reproducibility check, CF453 WD coefficients reproduce the overlapping saved WP2 results with maximum absolute difference **5.55e-16**. WP3 reads original inputs and frozen reference selections; it does not need prior WP2 inference to calculate new CF459 coefficients.

## 4. Whole-mouse null tests and their limits

The following two-WD-cohort null tables preserve the original saved WP3 calculations. They predate the exact-zero null-degeneracy correction documented in Section 1.1 and are shown as the historical analysis, not newly corrected bootstrap estimates.
The primary exact test keeps observed CF453 coefficients fixed and enumerates all **4! = 24** complete CF459 lipid-to-RNA mouse pairings. The complementary test keeps observed CF459 coefficients fixed and enumerates all **7! = 5,040** CF453 pairings. The identity pairing and ties are included. Both tests use the two-sided absolute coefficient-vector r, equivalent to absolute ULM at a fixed edge count. These are separate conditional pairing-independence nulls, not a test of equal cohort correlations and not two independent replications.

| Comparison | Nominal ULM p | Exact CF459 p (24) | Exact CF453 p (5040) | Null bootstrap p (defined only) |
| --- | --- | --- | --- | --- |
| Male periportal full | 6.22e-15 | 0.1250 | 0.0423 | 0.0678 |
| Male periportal displayed | 0.7961 | 0.8750 | 0.7944 | 0.8361 |
| Male pericentral full | 0.0120 | 0.7917 | 0.5966 | 0.7656 |
| Male pericentral displayed | 0.2847 | 0.6667 | 0.5514 | 0.5439 |

The approximate null bootstrap independently samples complete CF459 RNA and lipid mouse vectors from their empirical margins, holding CF453 coefficients fixed (**4999** requested draws; seed **20261014**). It tests the specified RNA-lipid independence null conditional on the empirical four-mouse cohort. It is not an ordinary paired bootstrap and does not estimate a probability that the observed biological network is true.

With four mice, null resampling often has very little animal diversity. If any edge in the fixed vector becomes undefined, the entire vector draw is undefined; the code never deletes a different set of edges in each draw. The reported bootstrap p-value is conditional on defined draws. The completion range below assigns all undefined draws to either nonexceedance or exceedance; it is a mathematical missing-draw bound, **not a confidence interval or calibrated correction**.

| Comparison | Exact 24 defined | Exact 5040 defined | Bootstrap defined / requested | Completion range |
| --- | --- | --- | --- | --- |
| Male periportal full | 24 / 24 | 5040 / 5040 | 4819 / 4999 | 0.0654 to 0.1014 |
| Male periportal displayed | 24 / 24 | 5040 / 5040 | 4819 / 4999 | 0.8060 to 0.8420 |
| Male pericentral full | 24 / 24 | 5040 / 5040 | 2153 / 4999 | 0.3298 to 0.8990 |
| Male pericentral displayed | 24 / 24 | 5040 / 5040 | 3084 / 4999 | 0.3356 to 0.7186 |

A concrete source of pericentral undefined vectors is limited CF459 expression diversity: Krt4 (2 observed expression levels; 1565/4999 constant RNA resamples); Pitx3 (3 observed expression levels; 390/4999 constant RNA resamples); Shisa7 (2 observed expression levels; 1611/4999 constant RNA resamples). A resample can omit an animal carrying the distinguishing expression level and make the gene constant. Complete exact pairing permutations preserve the original expression margins and remain defined. These genes and their fixed edges are retained; no post hoc rare-gene exclusion is used to improve the result.

| Method | Draws | One unique mouse | Two or fewer unique mice | Median unique mice |
| --- | --- | --- | --- | --- |
| CF459_null_rna | 4999 | 95 | 1704 | 3.0 |
| CF459_null_lipid | 4999 | 86 | 1729 | 3.0 |
| CF453_paired | 999 | 0 | 9 | 5.0 |
| CF459_paired | 999 | 18 | 361 | 3.0 |

The four-mouse exact test is intrinsically coarse: with 24 defined configurations its smallest possible p-value is 1/24 (about 0.0417), before accounting for tied statistics. A weak or inconclusive animal-level result does not establish absence of shared regulation. Conversely, a strong nominal ULM result alone cannot overcome four biological replicates. [Null-check figures](plots/null_checks.pdf) display the full-vector nulls and defined-draw counts.

## 5. Stability and prespecified Mtmr11 example

Ordinary paired bootstrap resamples whole paired RNA/lipid mice independently within each cohort (999 requested draws). The same animal indices are used in both zones. This describes sensitivity to animal composition, without supplying an additional null-test p-value. Frozen-vector undefined draws remain visible.

| Comparison | Paired defined / requested | Positive vector r among defined | Observed direction among defined |
| --- | --- | --- | --- |
| Male periportal full | 972 / 999 | 11.6% | 88.4% |
| Male periportal displayed | 980 / 999 | 45.9% | 45.9% |
| Male pericentral full | 391 / 999 | 36.8% | 63.2% |
| Male pericentral displayed | 590 / 999 | 26.6% | 73.4% |

| Zone | Omissions | Defined / planned | Minimum r | Maximum r | Observed direction |
| --- | --- | --- | --- | --- | --- |
| Male periportal | joint_omission | 28 / 28 | -0.091 | 0.005 | 96.4% |
| Male periportal | CF453_only | 7 / 7 | -0.082 | -0.044 | 100.0% |
| Male periportal | CF459_only | 4 / 4 | -0.068 | -0.017 | 100.0% |
| Male pericentral | joint_omission | 14 / 28 | -0.052 | 0.034 | 57.1% |
| Male pericentral | CF453_only | 7 / 7 | -0.047 | 0.002 | 85.7% |
| Male pericentral | CF459_only | 2 / 4 | -0.039 | 0.021 | 50.0% |

There are **28 joint** one-CF453/one-CF459 omissions (six plus three retained mice), **seven CF453-only** omissions and **four CF459-only** omissions. These prespecified configurations are descriptive and share animals; they are not 39 independent validation experiments. [Stability figures](plots/stability_and_omissions.pdf) show the complete bank.

| Example | Pairs | Vector r | ULM score | Nominal ULM p | Null bootstrap p (defined only) |
| --- | --- | --- | --- | --- | --- |
| Male periportal full | 610 | 0.063 | 1.549 | 0.1218 | 0.4647 |
| Male periportal displayed | 6 | -0.151 | -0.306 | 0.7746 | 0.6226 |

The Mtmr11 null-bootstrap p-values in this table use the corrected CF453-anchor null bank documented in Section 1.1. Observed coefficients, native ULM and nominal ULM p are unchanged; the original saved inference and scientific figures remain historical records.

Mtmr11 periportal was chosen from the original Figure 5 context before this comparison. The [coefficient figure](plots/Mtmr11_coefficients.pdf) uses its full and displayed shared species. The three [individual-mouse examples](plots/Mtmr11_selected_pairs.pdf) were selected by the largest absolute **original FCG displayed** coefficients, with a lipid-name tie break; selection never uses CF453/CF459 agreement. They remain illustrative associations. Cohort panels have separate axes because absolute export scales are not assumed comparable.

| Lipid | Frozen FCG r | CF453 WD r | CF459 WD r |
| --- | --- | --- | --- |
| PC(32:2) | 0.720 | 0.035 | -0.964 |
| PE(O-22:1_22:6) | 0.714 | -0.484 | 0.817 |
| PS(20:4_22:6) | 0.687 | -0.341 | -0.181 |

## 6. Prespecified sensitivities and intake QC

| Sensitivity | Comparison | Pairs | Vector r | ULM score |
| --- | --- | --- | --- | --- |
| CF453_QC30 | Male periportal full | 15125 | -0.063 | -7.784 |
| CF453_QC30 | Male periportal displayed | 69 | 0.032 | 0.259 |
| CF453_QC30 | Male pericentral full | 18150 | -0.019 | -2.595 |
| CF453_QC30 | Male pericentral displayed | 198 | -0.082 | -1.151 |
| CF453_Fill50 | Male periportal full | 14875 | -0.064 | -7.865 |
| CF453_Fill50 | Male periportal displayed | 69 | 0.032 | 0.259 |
| CF453_Fill50 | Male pericentral full | 17850 | -0.018 | -2.392 |
| CF453_Fill50 | Male pericentral displayed | 192 | -0.072 | -0.995 |
| RNA_data | Male periportal full | 15250 | -0.063 | -7.807 |
| RNA_data | Male periportal displayed | 69 | 0.032 | 0.259 |
| RNA_data | Male pericentral full | 18300 | -0.019 | -2.513 |
| RNA_data | Male pericentral displayed | 202 | -0.076 | -1.073 |

CF453 QC-CV <=30% and supplied Fill >=0.5 sensitivities restrict only the CF453 annotation/QC ledger; no equivalent CF459 fields were supplied. The RNA data-layer sensitivity uses stored data expression while preserving fixed molecular membership, animal pairing and original selection. Pearson correlations are invariant to positive row scaling within a cohort, so agreement of stored data and scale.data results is expected and does not independently validate upstream normalization. Full sensitivity tables and resampling remain saved separately.

| Cohort | Male WD mice | Supplied lipid rows | Eligible lipids | Analytical QC injections |
| --- | --- | --- | --- | --- |
| CF453 | 7 | 1240 | 1177 | 7 |
| CF459 | 4 | 1463 | 1395 | Unavailable in export |

CF459 import correctly handles a UTF-8 BOM, semicolon delimiter, decimal commas, scientific notation and quoted annotation semicolons. Its 1,463 rows retain 1,395 eligible lipids after 68 repository-class exclusions, with no supplied internal-standard rows. All source annotations and 244 duplicate Alignment-ID excess rows remain auditable. No CF459 QC CV, Fill, adduct or RT/m/z evidence is invented.

Both cohorts supply named Seurat 5.3.0 RNA Assay5 counts, data and scale.data layers with 20,123 genes. Stored scale.data is primary and unchanged. CF459 scaling diagnostics use all 40 original cell-type profiles or all 16 zone profiles before hepatocyte/male subsetting. Existing data replays stored scaling to approximately 1.25 x 10^-14. Counts are finite, nonnegative integers; aggregate library sizes and marker ordering are diagnostics, not nucleus counts or proof that nucleus-level QC/gates passed.

The CF459 RNA labels f-ctrl/m-ctrl identify recorded sex; the meeting establishes WD-only exposure. Four female mice appear only in intake and scaling-roster QC, not new female networks. The supplied assay exports lack age, genetic-background, diet-duration, litter and batch fields. The age/diet descriptions and their provenance are recorded in the requested-comparison section; none establishes controlled matching or equal diet exposure. [QC and coverage figures](plots/QC_and_coverage.pdf) show the design and aggregate diagnostics.

## 7. Interpretation and reproducibility

Cross-cohort coefficient agreement supports robustness only to the extent observed in these measured, strictly matched male WD features. It does not establish equal cohort effects, eliminate age/background/batch confounding, prove causal gene regulation, or identify WD-induced change. CF459 has no SD arm. Conversely, low or unstable agreement with four mice cannot establish absence of shared biology. Report the full-network r and the displayed subset together, alongside distinct conditional null tests and animal-composition stability.

Primary results, full exclusions, bootstrap/permutation indices and omissions are retained under this cohort-specific output directory. Original FCG, siRNA and completed WP1/WP2 numerical results are preserved. Input byte hashes and remaining metadata limits are documented in data/WP3_INPUTS.md and data/wp3_input_migration_manifest.csv. The report is regenerated from saved numeric results and report_bundle.rds; no statistical rerun is needed for --report-only.

- [fcg_comparison/inference.csv](fcg_comparison/inference.csv)
- [fcg_comparison/replay.csv](fcg_comparison/replay.csv)
- [fcg_comparison/paired_summary.csv](fcg_comparison/paired_summary.csv)
- [fcg_comparison/omission_summary.csv](fcg_comparison/omission_summary.csv)
- [fcg_comparison/analysis.rds](fcg_comparison/analysis.rds)
- [crosscohort_inference.csv](crosscohort_inference.csv)
- [network_edges.csv](network_edges.csv)
- [coverage.csv](coverage.csv)
- [lipid_crosswalk.csv](lipid_crosswalk.csv)
- [edge_eligibility_ledger.csv](edge_eligibility_ledger.csv)
- [sensitivities.csv](sensitivities.csv)
- [loco_summary.csv](loco_summary.csv)
- [resampling_draw_summary.csv](resampling_draw_summary.csv)
- [predefined_pairs.csv](predefined_pairs.csv)
- [analysis.rds](analysis.rds)
- [sensitivity_analyses.rds](sensitivity_analyses.rds)

Reproduce from the repository root, using a fresh output directory for a full run:

```bash
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --output=results/wd_crosscohort_reproduction --bootstraps=4999 --paired-bootstraps=999 --seed=20261014
Rscript scripts/workflows/13_wd_crosscohort_robustness.R --report-only --output=results/wd_crosscohort_reproduction
python3 scripts/reports/build_wp3_report_pdf.py --input results/wd_crosscohort_reproduction/README.md
```

The PDF appendix begins with two new scientific vector pages for the matched original-FCG/CF453 anchors and direct-closeness stability, followed by the six original WP3 figure pages: WD network coefficients; historical exact and null checks; stability/omissions; prespecified Mtmr11 coefficients; original-reference-selected individual pairs; QC/coverage.
