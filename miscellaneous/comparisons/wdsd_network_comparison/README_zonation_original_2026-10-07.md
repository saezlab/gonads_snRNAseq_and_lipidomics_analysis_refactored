# WP2: WDSD baseline conservation and Western-diet association changes

Analysis date: 7 October 2026. This report implements work package 2 of the 5 October meeting plan.

## Main findings

2 of four displayed-edge comparisons have positive coefficient-vector association. 0 of four have positive association and exact pairing p <=0.05; 0 of four have positive association and null-bootstrap p <=0.05. These counts summarize four predefined contexts and do not define a universal conservation threshold. The largest displayed-edge point estimate is Female pericentral (vector r=0.311, ULM=2.196, nominal p=0.0333, exact pairing p=0.2720, null-bootstrap p=0.2646).

Displayed-pair median absolute SD-to-WD changes are 0.646 (female periportal) and 0.795 (female pericentral), compared with 0.358 and 0.373 in the corresponding male contexts. 56 of 434 displayed pairs remain sign reversals across all 49 joint mouse omissions. The corresponding full-universe median changes are 0.405, 0.372, 0.426, 0.430 in male PP, male PC, female PP, female PC order. The apparent female excess among displayed pairs therefore depends on the fixed feature selection. These are descriptive follow-up candidates, not evidence of sex-specific causal rewiring.

The analysis assesses whether the original Figure 5 gene-lipid coefficient patterns transfer to standard-diet CF453 mice, then describes SD-to-WD changes within CF453. Nominal ULM results are reported beside mouse-level null tests. Diet changes are descriptive and assessed for omission stability; they are not established differential correlations or causal rewiring.

## Cohort, processing and QC

The updated mouse-ID crosswalk pairs all 28 CF453 mice between modalities, with seven female SD, seven female WD, seven male SD and seven male WD mice. There are 28 whole-hepatocyte profiles and 56 paired periportal/pericentral profiles. Each mouse supplies one whole-liver lipid measurement shared by its two RNA zones; 56 profiles do not create 56 independent animals.

The meeting described approximately 10-month-old mice after nine months of dietary exposure, compared with an approximately three-month-old original cohort, and a different genetic background. These descriptions are meeting context; exact ages, strain/cross and exposure records are not present in the supplied analysis metadata.

The lipid export contains 1240 annotation rows and 7 analytical QC injections. Exclusions comprise 14 internal standards, 0 additional isotope-labeled features and 49 repository-excluded-class rows, leaving 1177 primary features. 1129 eligible features have QC CV <=30%; 48 exceed 30%. Median QC CV is 5.5% and maximum is 67.3%. The QC-CV cutoff is a sensitivity, not a primary outcome-driven filter.

Primary lipid values are log2 supplied positive intensities. No extra normalization or imputation is applied. Primary RNA values are the supplied `scale.data` layer. Row scaling the supplied `data` layer over its complete original profile roster reproduces `scale.data` numerically; this checks scaling, rather than validating the upstream normalization. PCA and aggregate marker distributions are diagnostics, not grounds for excluding animals because they weaken conservation.

The supplied `Fill %` field ranges from 0.095 to 1.000; 91 supplied rows and 89 retained rows have values below 0.5 despite positive exported intensities. Their underlying detection, missingness and possible imputation provenance is unavailable. The Fill % >=0.5 restriction was added after intake inspection and before association analysis as an additional QC sensitivity.

Expected Cyp2e1/Cyp2f2 marker ordering occurs in 56 of 56 aggregate zone profiles. 56 satisfy the conventional threshold when applied to the stored aggregate data layer. This aggregate check cannot reconstruct the original nucleus-level gate.

Per-mouse and per-zone nucleus counts are unavailable for all 28 mice. Supplied count-library totals and detected genes are reported, but are not nucleus contributions. The recorded command history contains only ScaleData; upstream nuclear QC, doublets, aggregation, missing-value treatment and normalization provenance remain unresolved.

QC source tables: [sample pairing](metadata.csv), [input QC metrics](input_qc_summary.csv), [zone marker summary](zonation_marker_summary.csv), [input limitations](input_limitations.txt). [QC figures](plots/QC_diagnostics.pdf) show lipid/RNA PCAs, analytical variability and six marker distributions.

## Frozen Figure 5 reference and molecular coverage

| Context | FCG mice | Genes | Full reference edges | Displayed reference edges | Eligible full edges | Eligible displayed edges |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 22 | 25 | 31350 | 155 | 18800 | 86 |
| Male pericentral | 22 | 30 | 37620 | 387 | 22560 | 249 |
| Female periportal | 21 | 12 | 15048 | 89 | 9024 | 52 |
| Female pericentral | 21 | 11 | 13794 | 77 | 8272 | 47 |

752 of 1,254 original molecular annotations have strict eligible matches (502 lack a strict match), giving 58656 eligible full pairs and 434 eligible displayed pairs. All originally selected context genes are measurable and varying in both diets. Here feature loss is driven by exact species matching, rather than absent genes or invariant observed measurements. Unmatched annotations are retained in the crosswalk; no lower-resolution species/chain collapse is used to force coverage.

Original Pearson coefficients are independently reconstructed from 22 testes-bearing or 21 ovary-bearing mice. All 97812 archived gene-lipid edges replay using log2 lipid intensities, with maximum absolute error 8.88e-16. Raw intensities do not replay (maximum context discrepancies 0.423, 0.539, 0.523, 0.440).

The original full exports freeze context-specific genes and 1,254 lipidomics features; primary membership is taken directly from original displayed SIF files (documented absolute original r >0.65). The 1378 TF-to-gene edges are excluded. Pitx3 is also a regulon source and carries an ATT TF label, but its measured RNA-to-lipid branch is correctly retained.

Matching preserves gene symbols and each original lipid's chain order, slash/underscore separator, ether prefix, oxygenation and structural resolution. Only the final pipe-delimited annotation and class formatting are harmonized. Formatting agreement does not independently establish molecular identity. Alignment ID is not treated as a unique key; every annotation row has a feature ID. Missing, excluded, ambiguous or invariant measurements remain in the complete audit ledgers.

The original displayed-edge and full-edge universes are frozen before WDSD correlations are examined. Eligibility requires the same measured, varying RNA/lipid features in both SD and WD. No WDSD r cutoff is imposed. Original selection of strong associations makes displayed-edge conservation conditional on the original selection; the complementary full universe reduces this particular selection restriction, while retaining the originally selected genes.

Coverage and identity audits: [reference replay](reference_replay_summary.csv), [network coverage](network_coverage.csv), [lipid crosswalk](lipid_identity_crosswalk.csv), [edge eligibility](edge_eligibility_ledger.csv).

## Stage A: standard-diet conservation

Each context compares original reference coefficients with the matching sex and zone in seven standard-diet mice. Positive coefficient-vector Pearson r and positive native ULM slope t indicate conserved direction. Negative values indicate inverse alignment. ULM p-values are nominal, because edges share genes, lipids and seven animals; the molecular edge count is not the biological sample size.

### Primary: original displayed pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 86 | 0.262 | 2.487 | 0.0149 | 0.3226 | 0.3136 |
| Male pericentral | 249 | -0.028 | -0.438 | 0.6618 | 0.8875 | 0.8856 |
| Female periportal | 52 | -0.063 | -0.447 | 0.6569 | 0.7692 | 0.7735 |
| Female pericentral | 47 | 0.311 | 2.196 | 0.0333 | 0.2720 | 0.2646 |

### Complementary: all eligible original full-network pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 18800 | 0.016 | 2.158 | 0.0309 | 0.6891 | 0.6738 |
| Male pericentral | 22560 | -0.010 | -1.450 | 0.1470 | 0.8970 | 0.8930 |
| Female periportal | 9024 | -0.042 | -3.996 | 6.51e-05 | 0.3683 | 0.3165 |
| Female pericentral | 8272 | -0.051 | -4.687 | 2.82e-06 | 0.1534 | 0.1688 |

The primary exact pairing test exhausts all 7!=5,040 RNA-lipid pairings within each sex's SD group, including the observed identity pairing. One pairing is shared across all genes, lipids and both zones. Its null is lack of mouse-pairing information, stronger than simply no alignment with this specific reference. The approximate empirical-marginal null-bootstrap sensitivity independently resamples complete RNA-mouse and lipid-mouse vectors and tests their independence conditional on the observed empirical margins. Both p-values use two-sided absolute coefficient-vector association; at fixed edge count its ordering is equivalent to absolute ULM. The bootstrap uses `(1 + exceedances)/(1 + valid draws)`; the exact permutation uses its exhaustive tail proportion.

These animal-level tests retain dependence among molecular features. They remain conditional, small-cohort tests and do not establish regulation, substrate specificity or a mechanistic pathway. A small two-sided p-value with negative alignment would indicate inverse alignment, rather than positive conservation. Neither null allows arbitrary RNA-lipid associations while asking only whether they align with the reference. Exact pairings retain all seven observed mice per margin; empirical bootstrap samples contain about 4.62 distinct mice per modality on average. Its tails can be coarse, unstable or undefined at this group size. A higher or lower bootstrap p-value does not establish conservatism; type-I error and power have not been calibrated for this WDSD bootstrap procedure.

Exact inference assumes exchangeability of RNA-lipid pairing under the independence null within each sex's SD group. Batch and litter covariates are unavailable, so pairing tests do not independently control those possible shared sources of variation. Four predefined network contexts and secondary focus-gene summaries are exploratory; nominal edge-based p-values and unadjusted mouse-level p-values do not establish multiplicity-adjusted discoveries.

An ordinary paired-mouse bootstrap resamples matched RNA/lipid mice together and describes observed conservation stability. An additional paired bootstrap resamples original FCG mice within sex-chromosome/gonad groups and WDSD paired mice independently, sharing draws across zones. Neither supplies a null p-value. Original selected edges stay fixed in these draws, so uncertainty in their original selection is not propagated. The fractions below are descriptive stability proportions, not probabilities that the biological interpretation is true.

| Context | Pairings valid / total | Null bootstrap valid / total | Paired stability valid / total | Paired same-direction fraction | Reference-paired same-direction fraction |
| --- | --- | --- | --- | --- | --- |
| Male periportal | 5040 / 5040 | 4999 / 4999 | 999 / 999 | 0.907 | 0.908 |
| Male pericentral | 5040 / 5040 | 4896 / 4999 | 980 / 999 | 0.607 | 0.602 |
| Female periportal | 5040 / 5040 | 4996 / 4999 | 999 / 999 | 0.590 | 0.581 |
| Female pericentral | 5040 / 5040 | 4999 / 4999 | 999 / 999 | 0.989 | 0.988 |

Undefined constant-feature draws invalidate a fixed-universe summary instead of silently changing its edges. Valid and invalid counts for every universe and focus gene are saved in [baseline inference](baseline_inference.csv). [Baseline figures](plots/baseline_conservation.pdf) retain original displayed and full eligible results on separate pages.

## Stage B: Western versus standard diet

The two diet groups contain different mice. Every predefined eligible pair has `delta_r = r_WD - r_SD`; absolute change and sign reversal are exported alongside both coefficients. All 49 combinations dropping one SD and one WD mouse are evaluated, with the same omission applied across both zones. Direction stability and reversal stability require all 49 combinations to be defined, rather than counting a fraction of one from a smaller valid subset.

### Primary displayed-pair association

| Context | Edges | SD-WD vector r | Descriptive ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Male periportal | 86 | 0.094 | 0.865 | 0.3895 |
| Male pericentral | 249 | 0.098 | 1.547 | 0.1231 |
| Female periportal | 52 | -0.119 | -0.847 | 0.4012 |
| Female pericentral | 47 | -0.049 | -0.329 | 0.7439 |

These native ULM statistics describe SD-WD coefficient-vector association. Their nominal p-values test the fitted edge-vector slope under independent residual assumptions; they are not differential-correlation or rewiring tests. No diet-bootstrap p-value is assigned.

### Primary displayed-pair changes

| Context | Edges | Median absolute change | Sign reversal % | Direction stable in all 49 % | Reversal in all 49 % |
| --- | --- | --- | --- | --- | --- |
| Male periportal | 86 | 0.358 | 53.5 | 29.1 | 14.0 |
| Male pericentral | 249 | 0.373 | 49.4 | 32.9 | 7.2 |
| Female periportal | 52 | 0.646 | 55.8 | 59.6 | 21.2 |
| Female pericentral | 47 | 0.795 | 63.8 | 63.8 | 31.9 |

### Full eligible-pair association

| Context | Edges | SD-WD vector r | Descriptive ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Male periportal | 18800 | 0.031 | 4.247 | 2.18e-05 |
| Male pericentral | 22560 | 0.114 | 17.261 | 2.47e-66 |
| Female periportal | 9024 | 0.010 | 0.938 | 0.3482 |
| Female pericentral | 8272 | 0.107 | 9.805 | 1.42e-22 |

### Full eligible-pair changes

| Context | Edges | Median absolute change | Sign reversal % | Direction stable in all 49 % | Reversal in all 49 % |
| --- | --- | --- | --- | --- | --- |
| Male periportal | 18800 | 0.405 | 48.4 | 32.9 | 11.8 |
| Male pericentral | 22560 | 0.372 | 46.2 | 29.0 | 9.3 |
| Female periportal | 9024 | 0.426 | 49.9 | 37.6 | 13.2 |
| Female pericentral | 8272 | 0.430 | 46.3 | 37.8 | 11.8 |

A large absolute change is a descriptive observation at n=7 per group. Stability under omissions makes that observation less dependent on one mouse; it does not make it statistically significant or prove rewiring. No differential-correlation p-value is assigned. A simple diet-label shuffle would also test diet differences in the joint expression/lipid distributions, including abundance, rather than specifically equality of correlations.

The strongest illustrative pairs below are chosen descriptively: complete 49-omission delta-direction stability first, then minimum omitted absolute delta and observed absolute delta. This figure selection is recorded; it is not used to select the analyzed universe or to infer significance.

| Context | Gene | Lipid | r SD | r WD | Change in r | Omission direction fraction |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | Mup20 | LPC(18:3/0:0) | 0.869 | -0.893 | -1.761 | 1.000 |
| Male periportal | Traf5 | CL(17:1_18:2_18:1_18:2) | 0.788 | -0.859 | -1.647 | 1.000 |
| Male pericentral | Nat8f6 | PG(18:1_22:6) | -0.883 | 0.831 | 1.713 | 1.000 |
| Male pericentral | Mup9 | LPC(15:0) | 0.627 | -0.923 | -1.551 | 1.000 |
| Female periportal | Nol4l | PE(16:1_16:1) | 0.873 | -0.811 | -1.685 | 1.000 |
| Female periportal | Setbp1 | HexCer(18:1;O2/22:0) | 0.766 | -0.862 | -1.627 | 1.000 |
| Female pericentral | Mug1 | LPI(20:4) | 0.852 | -0.859 | -1.711 | 1.000 |
| Female pericentral | Rnf170 | PG(18:1_19:1) | -0.875 | 0.872 | 1.747 | 1.000 |

[Diet coefficient figures](plots/diet_network_comparison.pdf), [selected mouse measurements](plots/changing_pairs.pdf), [complete coefficient/change table](network_edges.csv) and [LOCO summaries](loco_aggregate_summary.csv) retain the species-level evidence.

## Predefined Mtmr11 focus

Mtmr11 is present only in the frozen male periportal reference: 1,254 full original pairs and nine originally displayed pairs. 752 full pairs and 7 displayed pairs are eligible in WDSD.


### Mtmr11 baseline: original displayed pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 7 | -0.660 | -1.962 | 0.1070 | 0.1921 | 0.1978 |

### Mtmr11 baseline: full eligible pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 752 | 0.020 | 0.557 | 0.5776 | 0.9204 | 0.9323 |

The 7 original displayed Mtmr11 coefficients occupy a narrow range (0.654 to 0.720); 7 are positive originally and 5 remain positive in WDSD SD. Its coefficient-vector statistic concerns relative ordering of these associations. A negative vector association does not imply that all seven gene-lipid associations reverse sign, nor opposite biological regulation. The ordinary paired bootstrap preserves the observed vector direction in 99.3% of valid draws; including paired original-reference resampling gives 65.0% (both stability descriptions, not null p-values).

### Mtmr11 diet coefficient-vector association

| Universe | Edges | SD-WD vector r | Descriptive ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Full | 752 | 0.193 | 5.373 | 1.03e-07 |
| Displayed | 7 | 0.360 | 0.862 | 0.4280 |

### Mtmr11 diet changes and omission stability

| Universe | Edges | Median absolute change | Sign reversals | Direction stable in all 49 | Reversal in all 49 |
| --- | --- | --- | --- | --- | --- |
| Displayed | 7 | 0.586 | 5 / 7 | 4 / 7 | 1 / 7 |
| Full | 752 | 0.475 | 320 / 752 | 337 / 752 | 88 / 752 |

4 of 7 displayed Mtmr11 pairs have negative SD-to-WD delta direction in all 49 omissions. 1 of 7 remain actual sign reversals in all 49 omissions: PE(O-24:2_22:6). Persistent direction changes and persistent sign reversals are different observations; neither is a formal rewiring test.

### Displayed Mtmr11 species and analytical quality

| Lipid | r SD | r WD | Change in r | Supplied Fill | QC CV % | Reversal all 49 |
| --- | --- | --- | --- | --- | --- | --- |
| PE(17:0_20:4) | 0.862 | -0.237 | -1.099 | 0.905 | 2.9 | No |
| PE(O-20:1_22:6) | 0.057 | -0.529 | -0.586 | 0.905 | 3.2 | No |
| PE(O-22:1_22:6) | 0.200 | -0.484 | -0.684 | 0.905 | 5.0 | No |
| PE(O-24:2_22:6) | 0.367 | -0.284 | -0.650 | 0.357 | 3.8 | Yes |
| PS(20:4_22:6) | -0.233 | -0.341 | -0.108 | 0.905 | 5.7 | No |
| PC(32:2) | -0.123 | 0.035 | 0.158 | 0.905 | 9.5 | No |
| PI(38:5) | 0.706 | 0.326 | -0.380 | 0.905 | 2.6 | No |

QC qualification: PE(O-24:2_22:6) is the displayed Mtmr11 sign-reversal illustration with low supplied Fill (0.357) despite low analytical QC CV (3.8%). It is excluded by the Fill >=0.5 sensitivity, so its apparent persistent reversal has an additional input-detection/provenance limitation. 3 higher-Fill displayed PE pairs retain negative delta direction across all 49 omissions: PE(17:0_20:4); PE(O-20:1_22:6); PE(O-22:1_22:6). Their persistent delta direction does not mean all are persistent sign reversals. The low-Fill pair remains in the frozen primary analysis, with its quality information reported explicitly.

Mtmr11 was specified before examining WDSD results. Its coefficient plots compare full and displayed species without requiring their new correlations to exceed the original cutoff. Up to six species are illustrated using actual scaled RNA and log2 whole-liver lipid values from seven mice per diet, prioritizing original displayed membership and then descriptive omission stability/change. The source identities and selection rule are in [the figure ledger](plots/selected_edges_for_figures.csv).

[Mtmr11 coefficient comparison](plots/Mtmr11_coefficients.pdf) and [individual mouse measurements](plots/Mtmr11_selected_pairs.pdf) support visual inspection; they do not measure enzyme activity, lipid flux or zonal lipid abundance.

## Prespecified sensitivities

The QC-CV <=30% feature restriction was prespecified, and the supplied Fill % >=0.5 restriction was added after input inspection before association analysis. Both rerun native ULM, exact pairings and the independent-marginal null bootstrap on the explicitly restricted original pairs. Neither is independent biological validation. Different denominators can themselves change coefficient-vector statistics; outcomes do not choose the preferred filter. The supplied RNA data layer is an affine replay of the scaled layer, documented in [RNA layer sensitivity](RNA_data_layer_sensitivity.csv), rather than independent normalization evidence.


### QC CV <=30%: original displayed pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 86 | 0.262 | 2.487 | 0.0149 | 0.3226 | 0.3136 |
| Male pericentral | 245 | -0.030 | -0.467 | 0.6406 | 0.8847 | 0.8830 |
| Female periportal | 50 | 0.000 | 0.003 | 0.9978 | 0.9980 | 0.9982 |
| Female pericentral | 45 | 0.356 | 2.494 | 0.0165 | 0.1897 | 0.1732 |

### Supplied Fill % >=0.5: original displayed pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 82 | 0.253 | 2.341 | 0.0217 | 0.3516 | 0.3294 |
| Male pericentral | 236 | -0.032 | -0.486 | 0.6277 | 0.8776 | 0.8734 |
| Female periportal | 52 | -0.063 | -0.447 | 0.6569 | 0.7692 | 0.7735 |
| Female pericentral | 47 | 0.311 | 2.196 | 0.0333 | 0.2720 | 0.2646 |

The [complete sensitivity table](sensitivities.csv) records retained universes and statistics. No confidence interval, reversed concordance score or multiplicity-adjusted discovery claim is needed for this exploratory report.

## Interpretation limits and manuscript use

Baseline conservation and within-cohort diet changes answer different questions. Where baseline evidence is inconclusive, a difference from the original network cannot be attributed specifically to Western diet; age, background, processing and sampling differences remain possible. The within-CF453 SD/WD comparison addresses diet-associated coefficient changes more directly, but does not establish their causality or reproducibility beyond this cohort.

The Figure 6 draft uses the frozen original displayed pairs in its baseline row and its SD-versus-WD row. Full-universe results, animal-level uncertainty and QC remain available beside the draft. Draft labels identify sample counts, feature coverage, ULM summaries and the explicit null tests, rather than treating edge-level nominal significance as biological replication.

Recorded input limitations:

- The supplied files are processed pseudobulks; per-mouse and per-zone nucleus counts, doublet/QC summaries and composition metadata are absent.
- Only ScaleData.RNA is recorded in each command history. The upstream normalization, missing-value treatment, imputation and nuclear zonation gate cannot be established from these objects.
- Row scaling the stored data layer reproduces scale.data numerically; it supplies a transparent scaling sensitivity but no independent normalization validation.
- Lipid intensity positivity and QC CV describe the supplied export, without establishing raw-data detection, batch correction or normalization provenance.
- No animal or lipid is removed on the basis of a network association, apparent biological outlier or favorable result. The QC-CV <=30% set is a sensitivity only.
- Strict canonical labels preserve chain order, separators, ether prefixes and oxygenation. Formatting agreement is not independent chemical identification; ambiguous annotations must not be forced into matches.
- Each mouse contributes two RNA zones but one whole-liver lipid profile. Twenty-eight mice remain the biological sample count.
- Some supplied Fill % values are below 0.5 despite positive exported intensities; upstream missingness and imputation provenance are unavailable.
- Cross-cohort comparisons involve older WDSD wild-type animals and younger original FCG animals; age, background and genotype differences cannot be separated here.

## Artifacts and reproduction

- [Figure 6 draft](plots/Figure6_draft.pdf): original displayed baseline and diet comparisons across four contexts.
- [Plot manifest](plot_manifest.csv): scientific figure descriptions and file locations.
- [Baseline inference](baseline_inference.csv), [diet summaries](diet_comparison_summary.csv), [complete species edges](network_edges.csv), [LOCO summaries](loco_aggregate_summary.csv) and [sensitivity results](sensitivities.csv).
- [Mouse metadata](metadata.csv), [lipid identities](lipid_identity_crosswalk.csv), [edge eligibility](edge_eligibility_ledger.csv), [network coverage](network_coverage.csv) and [reference replay](reference_replay_summary.csv).

Recorded seed: 20261007; independent-marginal null draws: 4999; paired stability draws: 999; exact pairings per sex: 5040. Shared mouse indices and per-omission IDs are retained.

Run from the repository root, choosing a fresh output directory:

```sh
Rscript scripts/workflows/05b_wdsd_network_comparison.R --validate-only
Rscript scripts/workflows/05b_wdsd_network_comparison.R --bootstraps=4999 --paired-bootstraps=999 --seed=20261007 --output=results/wdsd_network_comparison_reproduction
```

The entry point records settings, seeds, shared mouse indices, input/code checksums and package/session versions. Original FCG results and supplied objects remain the fixed inputs.

Run audit: [settings and scope](run_manifest.json), [input/code SHA256](input_and_code_SHA256.csv), [session versions](sessionInfo.txt).
