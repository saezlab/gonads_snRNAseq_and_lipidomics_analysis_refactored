# WP2: pooled and zonation network conservation, followed by diet-associated changes

Report reorganized: 09 October 2026. The numerical analyses use the recorded 7 October 2026 run settings; this report combines both completed WP2 analyses without changing their results.

## Overview

The original displayed Figure 3 coefficient pattern has pooled FCG-to-SD vector r=0.973 and native ULM=43.628 (nominal p=1.51e-69). The overall mouse-pairing p is 0.0104 and overall null-bootstrap p is 0.0104. Those overall tests allow shared sex differences to explain agreement. When sex margins are retained, the upper-tail pairing p is 0.4256 and upper-tail bootstrap p is 0.4480. These conditional checks do not establish positive alignment beyond that expected from the shared sex structure.

2 of four displayed-edge comparisons have positive coefficient-vector association. 0 of four have positive association and exact pairing p <=0.05; 0 of four have positive association and null-bootstrap p <=0.05. These counts summarize four predefined contexts and do not define a universal conservation threshold. The largest displayed-edge point estimate is Female pericentral (vector r=0.311, ULM=2.196, nominal p=0.0333, exact pairing p=0.2720, null-bootstrap p=0.2646).

Displayed-pair median absolute SD-to-WD changes are 0.646 (female periportal) and 0.795 (female pericentral), compared with 0.358 and 0.373 in the corresponding male contexts. 56 of 434 displayed pairs remain sign reversals across all 49 joint mouse omissions. The corresponding full-universe median changes are 0.405, 0.372, 0.426, 0.430 in male PP, male PC, female PP, female PC order. The apparent female excess among displayed pairs therefore depends on the fixed feature selection. These are descriptive follow-up candidates, not evidence of sex-specific causal rewiring.

The pooled whole-hepatocyte Figure 3 analysis is presented first, followed by the four sex-specific zonation networks from Figure 5 and their standard-to-Western-diet changes. These comparisons answer different questions: the pooled pattern includes between-sex covariance, while each zonation comparison uses one sex and seven SD mice. The lipid annotation audit and shared quality controls follow both network analyses. Native ULM scores and nominal p-values are retained beside explicitly defined mouse-level null tests; diet changes remain descriptive.

## 1. Pooled whole-hepatocyte network comparison: Figure 3

### Frozen network, original selection and cohort coverage

The current revised Figure 3 A states top three correlated genes per lipid, and panel B labels that species display FULL SPECIES-LEVEL NETWORK. Its Workflow 04 integrated export has 267 edges: 192 species-level gene-lipid edges involving 25 displayed genes and 64 lipids, plus 75 TF-gene edges. This analysis uses the 192 species pairs and excludes TF edges. The figure's use of FULL refers to molecular species resolution; it is distinct from the permissive 33%-tail archive containing 769 genes, 828 lipids and 636,732 pairs. That separate archive is inventoried, rather than substituted for Figure 3.

Original selection uses MOFA Factor 1's top 2.5% in each tail, a lipid-metabolism gene filter, and the three largest absolute gene correlations per selected lipid. It does not impose a 0.65 correlation cutoff: the minimum absolute displayed coefficient is approximately 0.5847. The unthresholded precursor includes 139 selected genes by the same 64 lipids, or 8,896 pairs; it is a complementary universe preserving the original MOFA selection. Neither uses new SD correlations to choose features.

All 43 original FCG mouse pairings are retained: 21 ovaries and 22 testes, with XX/XY and group metadata preserved. Scaled whole-hepatocyte RNA and log2 original lipid intensities replay the displayed and precursor archives to maximum absolute error 9.99e-16. New SD uses 14 distinct CF453 animals, seven recorded females and seven males. There is one whole-hepatocyte RNA profile and one whole-liver lipid profile per mouse. No zone-specific RNA profiles enter this analysis.

Matching is one-to-one through the strict final molecular annotations; missing or invariant features remain recorded in the coverage ledger. The same eligible, originally selected pairs enter every observed comparison and mouse draw. A new correlation is kept even if it is weak or reverses sign. [Reference replay](followup/pooled_reference_replay.csv), [coverage](followup/pooled_coverage.csv), [complete matched coefficients](followup/matched_pooled_edges.csv) and [mouse metadata](followup/metadata.csv) retain those identities.

For the actual Figure 3 displayed network, the strict join retains 36/64 lipids, 108/192 displayed pairs and 15/25 displayed genes. The conditional compatibility sensitivity retains 39/64 lipids, 117/192 pairs and 16/25 genes. These Figure 3 denominators are separate from the 752/1,254 whole-lipid-survey overlap.

### Pooled coefficient conservation and animal-level null tests

| Universe | Pairs | Vector r | ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Figure 3 displayed | 108 | 0.973 | 43.628 | 1.51e-69 |
| Full sparse precursor | 5004 | 0.855 | 116.671 | <2.23e-308 |

| Universe | Overall pairing p | Overall bootstrap p | Within-sex pairing upper p | Within-sex bootstrap upper p |
| --- | --- | --- | --- | --- |
| Figure 3 displayed | 0.0104 | 0.0104 | 0.4256 | 0.4480 |
| Full sparse precursor | 5.00e-05 | 2.00e-04 | 0.0633 | 0.0960 |

Positive vector r / positive ULM indicate positive alignment between the original and new correlation patterns. Native ULM is the intercept-model slope t statistic. Its nominal p-value treats edge residuals as independent even though pairs share genes, lipids and mice. Edge count therefore does not represent the biological sample size. These are exploratory unadjusted summaries, not independent confirmation of every edge. A returned numerical p-value of zero is reported as <2.23e-308, indicating floating-point underflow rather than an exact zero probability.

The unrestricted null permutes complete lipid mouse vectors relative to complete RNA vectors across all 14 SD animals. The independent-marginal null bootstrap resamples RNA mice and lipid mice separately across that pooled roster. Both test broad absence of pairing information; either may reject because females and males differ similarly in both modalities. They do not isolate regulation or correspondence within sex.

The conditional null makes those rearrangements or independent resamples within the two seven-mouse sex groups. Within-sex permutations retain the observed sex means and empirical margins exactly. The independent within-sex bootstrap preserves sex membership and resamples from each sex's observed empirical distribution, while individual draw means vary. Both conditional reference-alignment distributions can consequently be centered above zero. The reported conditional p-value is the prespecified one-sided upper tail: does the positive observed alignment exceed what the shared sex structure can already generate? Lower and equal-tail diagnostics are saved without being relabeled as positive-conservation evidence.

A [descriptive sex-mean diagnostic](followup/QA/sex_mean_only_concordance.csv) keeps only the between-sex covariance in SD, normalized by observed pooled feature variation. Its alignment with the original FCG coefficient vector is 0.9745 for displayed pairs and 0.8610 for the precursor, close to observed pooled alignment 0.9733 and 0.8551. This shows why shared sex structure can generate high pooled correspondence. It supplies no new p-value and does not prove that within-sex associations are absent.

The main run uses 19999 random pairing draws per null and 4999 null-bootstrap draws. Each draw is shared across every gene and lipid. Identity pairing is saved separately as draw zero; Monte Carlo p-values use `(1 + exceedances)/(1 + valid draws)`. All tails and valid/invalid counts are in [pooled inference](followup/pooled_inference.csv). Undefined features invalidate their fixed summary vector instead of shrinking its edge set within a draw. The bootstrap's type-I error and power at seven animals per sex have not been calibrated.

Paired stability bootstraps resample matched SD animals within sex. Additional stability draws independently resample original FCG paired mice within their chromosome/gonad groups. They describe sensitivity to mouse sampling and reference uncertainty, while keeping historical feature selection fixed. They do not supply a null p-value and do not propagate uncertainty in the original feature-selection process.

[Coefficient figures](followup/plots/pooled_conservation.pdf) separate displayed and precursor universes. [Null figures](followup/plots/pooled_null_checks.pdf) show why unrestricted and sex-conditioned checks can differ.

### Removing the between-sex means

| Universe | Pairs | Vector r | ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Full sparse precursor | 5004 | 0.082 | 5.796 | 7.19e-09 |
| Figure 3 displayed | 108 | 0.269 | 2.878 | 0.0048 |

This sensitivity subtracts each cohort's female or male mean from every RNA and lipid feature before recalculating Pearson correlations on unchanged pairs. It describes a different residual network, rather than replaying Figure 3. Reported residual ULM p-values remain nominal under dependent edges; no residual-network null p-value is assigned. Pooled conservation can coexist with weak residual conservation because between-sex effects and within-sex mouse variation answer different questions.

### Four predefined Cyp genes

| Gene | Pairs | Vector r | ULM | Nominal ULM p | Within-sex pairing upper p | Within-sex bootstrap upper p |
| --- | --- | --- | --- | --- | --- | --- |
| Cyp2b9 | 33 | 0.978 | 26.310 | 9.10e-23 | 0.4825 | 0.4844 |
| Cyp2b13 | 22 | 0.971 | 18.010 | 7.91e-14 | 0.9855 | 0.5762 |
| Cyp2c38 | 9 | 0.979 | 12.689 | 4.37e-06 | 0.1964 | 0.3950 |
| Cyp2d9 | 23 | 0.984 | 25.486 | 2.80e-17 | 0.0682 | 0.3620 |

The Cyp2b9, Cyp2b13, Cyp2c38 and Cyp2d9 summaries were specified before examining new results and retain their original displayed species membership. Full precursor summaries for those genes are also saved. These secondary summaries are exploratory and their p-values are unadjusted across genes and universes. A gene's pooled coefficient pattern does not identify enzyme activity, substrate specificity or causal regulation.

[Cyp coefficient figures](followup/plots/Cyp_pooled_conservation.pdf) show all available originally displayed pairs. [Individual-mouse examples](followup/plots/selected_pooled_pairs.pdf) show one pair per predefined gene: the matched originally displayed pair with the largest absolute original FCG coefficient, with original lipid-name ordering resolving ties. The new correlation does not enter selection. Seven females and seven males are colored separately, with all mouse IDs shown; pooled and within-sex fits help expose mean-driven associations. The [figure selection ledger](followup/plots/figure_pair_selection_ledger.csv) records these choices.

## 2. Sex-specific zonation networks and diet comparison: Figure 5

### Frozen Figure 5 reference and molecular coverage

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

### Standard-diet conservation within each sex and zone

Each context compares original reference coefficients with the matching sex and zone in seven standard-diet mice. Positive coefficient-vector Pearson r and positive native ULM slope t indicate conserved direction. Negative values indicate inverse alignment. ULM p-values are nominal, because edges share genes, lipids and seven animals; the molecular edge count is not the biological sample size.

#### Primary: original displayed pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 86 | 0.262 | 2.487 | 0.0149 | 0.3226 | 0.3136 |
| Male pericentral | 249 | -0.028 | -0.438 | 0.6618 | 0.8875 | 0.8856 |
| Female periportal | 52 | -0.063 | -0.447 | 0.6569 | 0.7692 | 0.7735 |
| Female pericentral | 47 | 0.311 | 2.196 | 0.0333 | 0.2720 | 0.2646 |

#### Complementary: all eligible original full-network pairs

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

### Western versus standard diet within CF453

The two diet groups contain different mice. Every predefined eligible pair has `delta_r = r_WD - r_SD`; absolute change and sign reversal are exported alongside both coefficients. All 49 combinations dropping one SD and one WD mouse are evaluated, with the same omission applied across both zones. Direction stability and reversal stability require all 49 combinations to be defined, rather than counting a fraction of one from a smaller valid subset.

#### Primary displayed-pair association

| Context | Edges | SD-WD vector r | Descriptive ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Male periportal | 86 | 0.094 | 0.865 | 0.3895 |
| Male pericentral | 249 | 0.098 | 1.547 | 0.1231 |
| Female periportal | 52 | -0.119 | -0.847 | 0.4012 |
| Female pericentral | 47 | -0.049 | -0.329 | 0.7439 |

These native ULM statistics describe SD-WD coefficient-vector association. Their nominal p-values test the fitted edge-vector slope under independent residual assumptions; they are not differential-correlation or rewiring tests. No diet-bootstrap p-value is assigned.

#### Primary displayed-pair changes

| Context | Edges | Median absolute change | Sign reversal % | Direction stable in all 49 % | Reversal in all 49 % |
| --- | --- | --- | --- | --- | --- |
| Male periportal | 86 | 0.358 | 53.5 | 29.1 | 14.0 |
| Male pericentral | 249 | 0.373 | 49.4 | 32.9 | 7.2 |
| Female periportal | 52 | 0.646 | 55.8 | 59.6 | 21.2 |
| Female pericentral | 47 | 0.795 | 63.8 | 63.8 | 31.9 |

#### Full eligible-pair association

| Context | Edges | SD-WD vector r | Descriptive ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Male periportal | 18800 | 0.031 | 4.247 | 2.18e-05 |
| Male pericentral | 22560 | 0.114 | 17.261 | 2.47e-66 |
| Female periportal | 9024 | 0.010 | 0.938 | 0.3482 |
| Female pericentral | 8272 | 0.107 | 9.805 | 1.42e-22 |

#### Full eligible-pair changes

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

### Predefined Mtmr11 focus

Mtmr11 is present only in the frozen male periportal reference: 1,254 full original pairs and nine originally displayed pairs. 752 full pairs and 7 displayed pairs are eligible in WDSD.


#### Mtmr11 baseline: original displayed pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 7 | -0.660 | -1.962 | 0.1070 | 0.1921 | 0.1978 |

#### Mtmr11 baseline: full eligible pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 752 | 0.020 | 0.557 | 0.5776 | 0.9204 | 0.9323 |

The 7 original displayed Mtmr11 coefficients occupy a narrow range (0.654 to 0.720); 7 are positive originally and 5 remain positive in WDSD SD. Its coefficient-vector statistic concerns relative ordering of these associations. A negative vector association does not imply that all seven gene-lipid associations reverse sign, nor opposite biological regulation. The ordinary paired bootstrap preserves the observed vector direction in 99.3% of valid draws; including paired original-reference resampling gives 65.0% (both stability descriptions, not null p-values).

#### Mtmr11 diet coefficient-vector association

| Universe | Edges | SD-WD vector r | Descriptive ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Full | 752 | 0.193 | 5.373 | 1.03e-07 |
| Displayed | 7 | 0.360 | 0.862 | 0.4280 |

#### Mtmr11 diet changes and omission stability

| Universe | Edges | Median absolute change | Sign reversals | Direction stable in all 49 | Reversal in all 49 |
| --- | --- | --- | --- | --- | --- |
| Displayed | 7 | 0.586 | 5 / 7 | 4 / 7 | 1 / 7 |
| Full | 752 | 0.475 | 320 / 752 | 337 / 752 | 88 / 752 |

4 of 7 displayed Mtmr11 pairs have negative SD-to-WD delta direction in all 49 omissions. 1 of 7 remain actual sign reversals in all 49 omissions: PE(O-24:2_22:6). Persistent direction changes and persistent sign reversals are different observations; neither is a formal rewiring test.

#### Displayed Mtmr11 species and analytical quality

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

## 3. Lipid annotation overlap and unresolved molecular identity

The original 752/1,254 figure is reproducible: 752 original species (60.0%) have the same final molecular annotation in the new export. This is a deliberately strict molecular-label overlap, not an estimate of how many lipids were physically detected. 959 original labels share at least one recorded pipe-delimited alias with a retained new row; these aliases often describe only sum composition. The cumulative audit reaches 963 after also allowing chain-order and positional-resolution candidates. These broader counts are candidates, not validated additional molecular identities.

The two exports use differing combinations of molecular-species labels, sum-composition aliases, sn-position notation and library identifiers. The original exported table has 1,386 rows; its cleaned network reference retains 1,254 measured rows. WDSD supplies 1,240 rows; the established exclusion of 14 internal standards and 49 rows in repository-excluded classes leaves 1,177 primary features. Multiple measured/adduct rows are never collapsed or summed to improve overlap.

The 1,254-species denominator is the pre-existing cleaned reference, rather than every row in the original 1,386-row export. The 132 raw rows absent from that reference include 111 whose classes match the repository exclusion list. The remaining 21 are 18 SM rows with ontology ASM, two Cer_EOS rows and one Cer_EODS row. This recorded inventory establishes which rows were removed before this follow-up, but does not establish every historical filtering reason. [Raw-to-cleaned ledger](followup/QA/raw_FCG_1386_to_cleaned_1254_ledger.csv) and [excluded-class inventory](followup/QA/raw_FCG_excluded_class_summary.csv) retain the full partition.

Repository exclusions remove zero of the 752 strict matches. The later QC-CV <=30% sensitivity retains 745 of those 752 survey-wide matches, excluding 7 already matched rows. All 36 strict Figure 3 lipids pass that cutoff, so its pooled network is unchanged. The 502 initially absent exact annotations therefore cannot be attributed to these analysis exclusions.

| Category | Original rows |
| --- | --- |
| Same final molecular annotation | 752 |
| Library identifier only; structure mismatch | 7 |
| No exported annotation or same-class formula candidate | 280 |
| Positional / structural-resolution candidate | 4 |
| Same class / formula only; possible isomers | 4 |
| Shared coarse annotation alias only | 202 |
| Unresolved CL chain-order candidate | 5 |

Categories are a priority partition, so each original row appears once. A row can have more than one candidate in the detailed ledger. In particular, the 202 coarse-alias-only rows are the largest potentially resolvable group after labels already assigned to the chain-order and positional categories. The 280 rows without a recorded annotation or same-class formula candidate represent a real exported-feature difference; their absence does not independently establish absence of a chromatographic peak or abundance below detection.

The largest missing-annotation contributions are DG (76), TG (73), PE (63), PC (55), Cer (45): 312 of the 502 strict-unmatched rows. These counts include candidate and candidate-free categories; they are not experimentally confirmed detection losses.

Three concrete cautions explain why a broad match is not automatically safe:

- CL chain quartets can be written in a different underscore order, but several original or new CL rows share the same unordered quartet. CL pairing/group assignment and peak identity remain unresolved; chain sorting does not provide a globally unique rescue.
- A shared alias can describe a sum such as a total carbon / double-bond composition while the final annotations resolve different acyl chains or ether/plasmalogen structure. It supports compatibility at that coarser resolution, rather than equality of molecular species.
- Exported InChIKey or SMILES fields can agree despite conflicting or differently resolved final annotations. The CL example `CL(16:0_16:1_18:1_18:2)` versus `CL(16:0_18:1_18:1_16:2)` changes reported chain constituents. `PE(O-18:1_18:2)` versus `PE(P-18:0_18:2)` may describe the same plasmalogen composition at different structural resolution because P- notation counts the vinyl bond differently. Shared library identifiers alone do not settle these unresolved structural differences.

A separate compatibility sensitivity extends the join to 854 original/new pairs: the 752 strict pairs plus 102 recorded coarse-alias candidates. Each additional candidate has the same lipid class and adduct, precursor agreement within 5 ppm, and bidirectional degree one across the complete candidate graph. Already assigned strict features are blocked. This protects against feature reuse; it does not validate molecular identity. Those pairs remain explicitly conditional matches at their reported annotation resolution.

Mass, retention time, formula and library fields are exported diagnostics. Cross-run retention times have not been calibrated and raw spectra have not been reidentified. Resolving the remaining candidate groups requires the collaborator's original peak/annotation evidence, especially acyl-chain fragments, adduct provenance and calibrated retention matching. No experimental-detection conclusion is drawn from the annotation audit alone.

The lipid notation distinguishes molecular species, positional resolution and more general composition. The matching rules follow [updated lipid shorthand nomenclature](https://pmc.ncbi.nlm.nih.gov/articles/PMC7707175/), the [LIPID MAPS structural hierarchy](https://www.lipidmaps.org/lipid_nomenclature/rules/hierarchy) and [glycerophospholipid notation rules](https://www.lipidmaps.org/lipid_nomenclature/glycerophospholipids/pl_lpl). The safe join changes formatting only; chemically broader matches stay separate.

Audits: [complete original annotation ledger](followup/lipid_overlap_audit/original_ledger.csv), [new annotation ledger](followup/lipid_overlap_audit/new_ledger.csv), [every candidate pair](followup/lipid_overlap_audit/candidate_pairs.csv), [class-level coverage](followup/lipid_overlap_audit/class_summary.csv), [strict safe map](followup/lipid_overlap_audit/safe_map.csv), [conditional compatibility map](followup/lipid_overlap_audit/conditional_compatibility_map.csv) and [structure/identifier conflict examples](followup/lipid_overlap_audit/identifier_structure_mismatch_examples.csv). [Overlap figures](followup/plots/lipid_overlap_audit.pdf) show the loss partition and class distribution.

## 4. Shared cohort, processing, QC and sensitivity analyses

### Cohort, processing and QC

The updated mouse-ID crosswalk pairs all 28 CF453 mice between modalities, with seven female SD, seven female WD, seven male SD and seven male WD mice. There are 28 whole-hepatocyte profiles and 56 paired periportal/pericentral profiles. Each mouse supplies one whole-liver lipid measurement shared by its two RNA zones; 56 profiles do not create 56 independent animals.

The meeting described approximately 10-month-old mice after nine months of dietary exposure, compared with an approximately three-month-old original cohort, and a different genetic background. These descriptions are meeting context; exact ages, strain/cross and exposure records are not present in the supplied analysis metadata.

The lipid export contains 1240 annotation rows and 7 analytical QC injections. Exclusions comprise 14 internal standards, 0 additional isotope-labeled features and 49 repository-excluded-class rows, leaving 1177 primary features. 1129 eligible features have QC CV <=30%; 48 exceed 30%. Median QC CV is 5.5% and maximum is 67.3%. The QC-CV cutoff is a sensitivity, not a primary outcome-driven filter.

Primary lipid values are log2 supplied positive intensities. No extra normalization or imputation is applied. Primary RNA values are the supplied `scale.data` layer. Row scaling the supplied `data` layer over its complete original profile roster reproduces `scale.data` numerically; this checks scaling, rather than validating the upstream normalization. PCA and aggregate marker distributions are diagnostics, not grounds for excluding animals because they weaken conservation.

The supplied `Fill %` field ranges from 0.095 to 1.000; 91 supplied rows and 89 retained rows have values below 0.5 despite positive exported intensities. Their underlying detection, missingness and possible imputation provenance is unavailable. The Fill % >=0.5 restriction was added after intake inspection and before association analysis as an additional QC sensitivity.

Expected Cyp2e1/Cyp2f2 marker ordering occurs in 56 of 56 aggregate zone profiles. 56 satisfy the conventional threshold when applied to the stored aggregate data layer. This aggregate check cannot reconstruct the original nucleus-level gate.

Per-mouse and per-zone nucleus counts are unavailable for all 28 mice. Supplied count-library totals and detected genes are reported, but are not nucleus contributions. The recorded command history contains only ScaleData; upstream nuclear QC, doublets, aggregation, missing-value treatment and normalization provenance remain unresolved.

QC source tables: [sample pairing](metadata.csv), [input QC metrics](input_qc_summary.csv), [zone marker summary](zonation_marker_summary.csv), [input limitations](input_limitations.txt). [QC figures](plots/QC_diagnostics.pdf) show lipid/RNA PCAs, analytical variability and six marker distributions.

### Pooled molecular mapping and analytical-quality sensitivities

| Sensitivity | Universe | Pairs | ULM | Nominal ULM p | Within-sex pairing upper p | Within-sex bootstrap upper p |
| --- | --- | --- | --- | --- | --- | --- |
| Strict annotation | Sparse precursor | 5004 | 116.671 | <2.23e-308 | 0.0633 | 0.0960 |
| Strict annotation | Displayed | 108 | 43.628 | 1.51e-69 | 0.4256 | 0.4480 |
| Conditional alias map | Sparse precursor | 5421 | 122.511 | <2.23e-308 | 0.0456 | 0.0866 |
| Conditional alias map | Displayed | 117 | 47.062 | 5.63e-77 | 0.3555 | 0.4316 |
| QC CV <=30% | Sparse precursor | 5004 | 116.671 | <2.23e-308 | 0.0633 | 0.0960 |
| QC CV <=30% | Displayed | 108 | 43.628 | 1.51e-69 | 0.4256 | 0.4480 |
| Fill >=50% | Sparse precursor | 4865 | 114.991 | <2.23e-308 | 0.0823 | 0.1040 |
| Fill >=50% | Displayed | 105 | 43.662 | 2.88e-68 | 0.3405 | 0.4402 |

The 854-pair conditional molecular compatibility map is a sensitivity, rather than a revised identity claim. Its original and new coefficients correspond to the explicitly recorded, potentially different molecular annotations, so any agreement is interpreted at the supported coarse annotation resolution. The QC-CV <=30% and supplied Fill % >=0.5 restrictions test dependence on the exported quality fields. They alter the denominator and do not constitute independent biological replication. [Mapping sensitivity](followup/mapping_sensitivity_inference.csv) and [QC sensitivities](followup/pooled_QC_sensitivity.csv) preserve complete results and null definitions.

Supplied positive lipid values are log2 transformed without additional normalization or imputation. Positivity does not establish that upstream missingness or imputation was absent; the supplied Fill % field records incomplete coverage for some rows. Supplied RNA scale.data is used without changing its upstream normalization. Row scaling the stored data layer replays scale.data, which validates the scaling relationship rather than independently validating normalization. Nucleus contributions, detailed upstream nuclear QC, batch and litter covariates remain unavailable.

### Zonation quality-field and RNA-layer sensitivities

The QC-CV <=30% feature restriction was prespecified, and the supplied Fill % >=0.5 restriction was added after input inspection before association analysis. Both rerun native ULM, exact pairings and the independent-marginal null bootstrap on the explicitly restricted original pairs. Neither is independent biological validation. Different denominators can themselves change coefficient-vector statistics; outcomes do not choose the preferred filter. The supplied RNA data layer is an affine replay of the scaled layer, documented in [RNA layer sensitivity](RNA_data_layer_sensitivity.csv), rather than independent normalization evidence.


#### QC CV <=30%: original displayed pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 86 | 0.262 | 2.487 | 0.0149 | 0.3226 | 0.3136 |
| Male pericentral | 245 | -0.030 | -0.467 | 0.6406 | 0.8847 | 0.8830 |
| Female periportal | 50 | 0.000 | 0.003 | 0.9978 | 0.9980 | 0.9982 |
| Female pericentral | 45 | 0.356 | 2.494 | 0.0165 | 0.1897 | 0.1732 |

#### Supplied Fill % >=0.5: original displayed pairs

| Context | Edges | Vector r | ULM | Nominal ULM p | Exact pairing p | Null-bootstrap p |
| --- | --- | --- | --- | --- | --- | --- |
| Male periportal | 82 | 0.253 | 2.341 | 0.0217 | 0.3516 | 0.3294 |
| Male pericentral | 236 | -0.032 | -0.486 | 0.6277 | 0.8776 | 0.8734 |
| Female periportal | 52 | -0.063 | -0.447 | 0.6569 | 0.7692 | 0.7735 |
| Female pericentral | 47 | 0.311 | 2.196 | 0.0333 | 0.2720 | 0.2646 |

The [complete sensitivity table](sensitivities.csv) records retained universes and statistics. No confidence interval, reversed concordance score or multiplicity-adjusted discovery claim is needed for this exploratory report.

## 5. Interpretation and remaining input limitations

This pooled test answers whether the sex-pooled, non-zonated Figure 3 pattern transfers to a new standard-diet cohort. It is complementary to WP2's within-sex zonation analysis and should not be described as rescuing a failed zonal result: the feature selection, biological profiles and between-sex contributions differ. Differences in age, background, diet duration and processing also limit direct cohort equivalence; meeting descriptions of those differences are not fully recorded in the supplied metadata.

A defensible manuscript distinction is overall pooled pattern reproduction versus additional correspondence beyond sex margins. Neither establishes direct lipid regulation or a mechanistic pathway. Non-significant conditional tests do not prove lack of a biological relationship; they indicate that this cohort and this conditional test do not establish the specified extra correspondence. Exact-annotation gaps likewise do not prove that the new experiment failed to detect those lipids.


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
- No raw peak absence or limit-of-detection inference follows from absence of an exported annotation.
- Same sum composition, formula or precursor mass cannot distinguish acyl-chain isomers; candidates are not extra validated features.
- The supplied INCHIKEY/SMILES fields can be identical despite different reported chain/ether/positional resolution; these library identifiers are not independent chemical confirmation.
- Underscores normally omit sn-position, but both datasets have multiple CL rows with the same unordered chain quartet; CL group/pair assignment and feature identity are unresolved, so none is automatically rescued.
- Retention time and 5 ppm precursor checks are descriptive diagnostics; no cross-run retention-time calibration or peak/spectrum re-identification was performed.
- Zero-chain removal and slash/underscore conversion lose positional resolution and remain candidates only.
- Pooled standard diet compares 14 wild-type WDSD mice with 43 original FCG mice; age, genetic background and cohort differences cannot be separated.
- Low supplied Fill values may indicate upstream missingness despite positive exported intensities; peak-level absence/imputation provenance unavailable.

## 6. Artifacts, local inputs and reproduction

Workflow 12 is the single WP2 entry point. It executes pooled Figure 3 / annotation analysis first, then the four sex/zone baseline and diet comparisons, and finally writes this package report. Processed inputs are stored in the repository data tree; the data manifest records source locations and integrity checks. Shared helpers retain the original calculations and frozen historical feature membership.

Repository input provenance: [input guide](../../data/WP1_WP2_INPUTS.md) and [source/destination SHA256 manifest](../../data/wp1_wp2_input_migration_manifest.csv) document the byte-preserving copies and original archive locations.

WP2 local inputs: [whole-hepatocyte Seurat pseudobulks](../../data/RNA/wdsd/wdsd_hep_pseudobulks.rds); [paired zonation Seurat pseudobulks](../../data/RNA/wdsd/wdsd_zonation_pseudobulks.rds); [CF453 SD/WD lipidomics export](../../data/metabolomics/wdsd/lipidomics_cf453.csv); [raw FCG lipid annotation export](../../data/metabolomics/fcg/lipidomics_fcg.csv); and [cleaned FCG lipid matrix](../../data/metabolomics/Lipid_data_file_JPC_cleaned.csv). The raw FCG annotation table supplements the fixed cleaned matrix rather than replacing it.

- Pooled: [inference](followup/pooled_inference.csv), [complete matched coefficients](followup/matched_pooled_edges.csv), [coverage](followup/pooled_coverage.csv), [reference replay](followup/pooled_reference_replay.csv), [mapping sensitivities](followup/mapping_sensitivity_inference.csv) and [QC sensitivities](followup/pooled_QC_sensitivity.csv).
- Zonation/diet: [baseline inference](baseline_inference.csv), [diet summaries](diet_comparison_summary.csv), [complete species coefficients](network_edges.csv), [joint omission summaries](loco_aggregate_summary.csv), [sensitivity results](sensitivities.csv) and [Figure 6 draft](plots/Figure6_draft.pdf).
- Pairing and provenance: [all 28 mice](metadata.csv), [14 pooled SD mice](followup/SD_mouse_manifest.csv), [43 pooled FCG mice](followup/FCG_mouse_manifest.csv), [zonal reference mice](reference_mouse_manifest.csv), [strict lipid crosswalk](lipid_identity_crosswalk.csv) and [all molecular candidates](followup/lipid_overlap_audit/candidate_pairs.csv).
- Figures: [pooled figure manifest](followup/plot_manifest.csv) and [zonation/diet figure manifest](plot_manifest.csv). The combined PDF preserves all scientific figure pages and presents pooled comparisons before zonation/diet comparisons, followed by annotation and shared QC diagnostics.

Recorded seed: 20261007; pooled random pairings per null: 19999; null-bootstrap draws: 4999; paired stability draws: 999; exhaustive zonation pairings per sex: 5040. Identity and shared mouse indices remain in each phase's audit files.

Run from the repository root. A full rerun requires a fresh output directory:

```sh
Rscript scripts/workflows/12_wdsd_network_comparison.R --validate-only
Rscript scripts/workflows/12_wdsd_network_comparison.R --permutations=19999 --bootstraps=4999 --paired-bootstraps=999 --seed=20261007 --output=results/wdsd_network_comparison_reproduction
```

Refresh the combined report from completed outputs without recomputing analyses:

```sh
Rscript scripts/workflows/12_wdsd_network_comparison.R --report-only
python3 scripts/reports/build_wp2_report_pdf.py
```

Each analysis phase retains its input/code SHA256, settings, session versions, shared resampling indices and null draws. The original FCG results and supplied data are fixed inputs. The preceding zonation-only README is archived once as README_zonation_original_2026-10-07.md; the pooled component report remains in followup/README.md.

Run audits: [zonal settings](run_manifest.json), [zonal input/code SHA256](input_and_code_SHA256.csv), [zonal session](sessionInfo.txt), [pooled settings](followup/run_manifest.json), [pooled input/code SHA256](followup/input_and_code_SHA256.csv) and [pooled session](followup/sessionInfo.txt).
