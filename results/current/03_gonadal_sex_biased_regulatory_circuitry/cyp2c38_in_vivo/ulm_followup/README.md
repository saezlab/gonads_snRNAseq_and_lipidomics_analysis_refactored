# Cyp2c38 in vivo lipidomics: ULM scores, bootstrap tests and regenerated gonadal networks

Completed run: 2026-10-07 01:19:31 UTC. Every result table below comes from the saved outputs.

The regenerated testicular-network comparison gives male ULM = -5.706 and nominal ULM p = 1.59e-08; its KO-only bootstrap p is 0.4828 and KO-plus-reference bootstrap p is 0.4488. The male association has the predicted direction and a small nominal p-value, but the whole-animal test does not establish robust evidence under its contrast-null model.

The regenerated ovarian-network comparison has a positive ULM score, so its observed direction does not support the predicted inverse response to reduced Cyp2c38.

## Regenerated full-network results

The native ULM score is the slope t-statistic from a linear model with an intercept: KO-minus-WT moderated lipid t-statistics are regressed on original Cyp2c38–lipid Pearson coefficients. A **negative ULM score** supports the expected perturbation direction. ULM p-values are nominal target-level p-values; bootstrap p-values resample whole mice and retain dependence among lipid species. All reported tests are two-sided.

| KO contrast / FCG reference | Matched lipids | ULM score | ULM p-value | Bootstrap p, KO only | Bootstrap p, KO + reference |
| --- | ---: | ---: | ---: | ---: | ---: |
| Pooled KO − WT / pooled FCG | 856 | -1.930 | 0.0539 | 0.8178 | 0.8200 |
| Female KO − WT / ovarian FCG | 856 | 0.331 | 0.7406 | 0.9528 | 0.9452 |
| Male KO − WT / testicular FCG | 856 | -5.706 | 1.59e-08 | 0.4828 | 0.4488 |

All **1,254 original lipidomics features** were regenerated without MOFA selection or a correlation/significance threshold. There are **856 unique strict molecular matches** to the retained KO feature roster. The ovarian and testicular networks use 21 and 22 original paired FCG mice, respectively; pooled coefficients use all 43 mice. The 32-mouse KO cohort contains female WT=8, female KO=8, male WT=7 and male KO=9. The primary model retains 1,237 curated endogenous KO features.

## Historical feature-universe sensitivity

| KO contrast / FCG reference | Matched lipids | ULM score | ULM p-value | Bootstrap p, KO only | Bootstrap p, KO + reference |
| --- | ---: | ---: | ---: | ---: | ---: |
| Pooled KO − WT / pooled FCG | 552 | -1.473 | 0.1412 | 0.8470 | 0.8446 |
| Female KO − WT / ovarian FCG | 552 | -1.497 | 0.1350 | 0.7404 | 0.7054 |
| Male KO − WT / testicular FCG | 552 | -3.961 | 8.45e-05 | 0.6272 | 0.5848 |

The historical 828-edge MOFA-derived universe retains **552 shared species**. Its pooled and within-gonad coefficients numerically reproduce the original WP1 reference. Reporting this subset alongside the complete regenerated network makes the expanded feature universe explicit rather than attributing changed ULM scores to a different KO model.

## What the bootstrap p-values test

The earlier ordinary bootstrap resampled the observed cohorts to estimate uncertainty. Its draws were not a null distribution; the earlier reported animal-level p-values came from genotype permutations. This follow-up creates an explicit null distribution for the **ULM statistic**, with each tested KO-minus-WT mean contrast constrained to zero for every lipid while the other group-mean dimensions remain unrestricted and are estimated under that constraint. This is stronger than a null of zero ULM alignment with arbitrary nonaligned KO effects still present; the bootstrap p-values must be interpreted under this stated contrast-mean null.

For the main KO-only wild bootstrap, fit the weighted four-group model, impose the tested mean-contrast restriction, and form null means. Use unrestricted four-group residuals divided by the square root of one minus their leverage (HC2). One independent Rademacher multiplier, −1 or +1 with equal probability, multiplies each whole-mouse residual lipid vector. This retains covariance among lipids and allows differing residual scales. Refit sex precision weights, robust/trend limma moderation, the KO t-statistics and ULM in every draw. This construction follows the null-imposed wild-bootstrap framework described by [Davidson and MacKinnon, Bootstrap Methods in Econometrics](https://www.russell-davidson.research.mcgill.ca/articles/bsurvey.pdf); finite-sample validity remains an assumption assessed below.

The reference network stays fixed in the KO-only column. The KO-plus-reference sensitivity additionally resamples paired original FCG RNA/lipid mice within XXO=10, XYO=11, XXT=11 and XYT=11 and recomputes each reference network. It therefore adds network-estimation uncertainty and answers a broader question. The run uses **4999 null bootstrap draws**, **1999 within-sex genotype permutations**, seed **20261007**, and 3 workers. Empirical p-values use `(1 + number of null |ULM| values at least as large as observed)/(1 + valid draws)`. Invalid-draw counts and reproducible draw seeds are retained in the full outputs.

A null-restricted-residual bootstrap is a sensitivity because those residuals can absorb a real KO signal and change the null distribution. The original within-sex genotype permutation is also retained as an independent comparison; it assumes exchangeability of genotype labels within sex. Neither bootstrap nor permutation treats molecular edges as independent experimental animals.

| Universe / contrast | Bootstrap p, restricted residuals | Genotype permutation p |
| --- | ---: | ---: |
| Full / Pooled | 0.8105 | 0.8060 |
| Full / Female | 0.9555 | 0.9455 |
| Full / Male | 0.4885 | 0.4590 |
| Historical / Pooled | 0.8425 | 0.8400 |
| Historical / Female | 0.7275 | 0.7105 |
| Historical / Male | 0.6280 | 0.6075 |

## Why nominal ULM and whole-mouse p-values differ

The fixed fitted-null bank gives a quantitative dependence diagnostic. Under the tested contrast's zero-mean-effect scenario, the whole-mouse ULM distribution is wide and nominal independent-target ULM p-values are frequently small:

| Full-network contrast | Fitted-null draws | Null ULM SD | Fraction with nominal ULM p≤0.05 |
| --- | ---: | ---: | ---: |
| Pooled | 4999 | 7.175 | 81.5% |
| Female | 4999 | 5.085 | 73.0% |
| Male | 4999 | 7.401 | 81.4% |

These fractions describe draws from this fitted contrast-null bank; they are not empirical population-wide error rates or the nested end-to-end calibration results below. They show that many nominal ULM p-values can be small even under the imposed mean-null when shared whole-mouse lipid variation is retained.

The KO-only bootstrap p-values are 0.8178/0.9528/0.4828, compared with genotype-permutation p-values 0.8060/0.9455/0.4590 (pooled/female/male). Adding paired FCG uncertainty gives 0.8200/0.9452/0.4488. Their agreement is evidence that the discrepancy with nominal ULM p-values is not caused solely by adding reference-estimation uncertainty.

## Is the bootstrap overly conservative because the KO contrast is underpowered?

Low power and conservatism are different properties. Low power is a low probability of detecting a real effect; conservatism is a rejection probability below the nominated level when the null is true. A weak individual-lipid contrast can still contribute collectively to ULM, and its lack of significance does not automatically make the bootstrap conservative. Conversely, a large nominal ULM magnitude can arise from correlated lipids sharing noise from a few mice. Whole-mouse resampling includes that dependence.

The null-calibration audit generates globally null datasets with no female or male KO effects, using the original HC2-corrected whole-mouse residual vectors and independent Rademacher signs. It preserves the actual group sizes and estimated sex-dependent lipid covariance. Each dataset constructs its own contrast-specific null, refits the complete pipeline, and receives a fresh inner wild-bootstrap test. The nominal rejection threshold is p≤0.05. This is a calibration check within the fitted sign-symmetric residual scenario, using the same resampling family as the test. The results below report counts, rejection rates and Monte Carlo standard errors rather than assuming ideal calibration.

| Universe / contrast | Null rejections / datasets | Null rejection rate | Monte Carlo SE, percentage points | Median bootstrap p | Exact binomial p vs 5% |
| --- | ---: | ---: | ---: | ---: | ---: |
| Full / Pooled | 5 / 120 | 4.2% | 1.8 | 0.535 | 0.8352 |
| Full / Female | 10 / 120 | 8.3% | 2.5 | 0.450 | 0.0942 |
| Full / Male | 9 / 120 | 7.5% | 2.4 | 0.475 | 0.2049 |
| Historical / Pooled | 4 / 120 | 3.3% | 1.6 | 0.535 | 0.5306 |
| Historical / Female | 14 / 120 | 11.7% | 2.9 | 0.527 | 0.0028 |
| Historical / Male | 9 / 120 | 7.5% | 2.4 | 0.435 | 0.2049 |

These fitted-scenario simulations do not demonstrate excessive conservatism by the exact binomial check. Some comparisons reject above the 5% target with diagnostic p<0.05. These unadjusted checks across six comparisons do not establish universal calibration, particularly with a small number of simulated datasets.

With only 120 outer null datasets per comparison, moderate conservatism cannot be ruled out regardless of the exact binomial diagnostic. Non-rejection of the 5% target is not evidence that the true rejection rate equals exactly 5%.

Each null dataset uses 199 inner bootstrap draws. These simulations condition on the fitted residual scenario and observed feature selection. They do not rule out miscalibration from unrecorded batch/litter structure, asymmetric errors or an inaccurate fitted covariance model.

Known-effect injection curves add deliberately concordant, network-shaped KO effects at controlled strengths while retaining the actual group sizes and fitted noise. Detection fractions across the entire full-network curve are:

| Full-network contrast | Datasets per setting | λ=0 | λ=0.1 | λ=0.2 | λ=0.4 | λ=0.8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Pooled | 399 | 5.5% | 27.1% | 68.9% | 100.0% | 100.0% |
| Female | 399 | 4.8% | 26.6% | 71.4% | 99.7% | 100.0% |
| Male | 399 | 5.8% | 14.0% | 36.6% | 90.5% | 100.0% |

At λ=0.8, the median absolute injected log2 effects are pooled 0.201; female 0.148; male 0.247. Every setting's physical effect size is retained in the linked complete injection table.

Injection strength is a controlled scenario parameter: shifts oppose the centered, SD-standardized full-network correlation pattern and are scaled by the original within-sex residual SD. Detection compares each simulated absolute ULM score with the observed-data null bank's 95th percentile, held fixed; it does not run a fresh inner test for every injected dataset. These are fitted-scenario sensitivity calculations, not an estimate of post-hoc power for the unknown true biological effect. Poor detection of small injected effects supports a sample-size/noise limitation; it does not by itself diagnose an overly conservative test. These diagnostics do not justify interpreting a large bootstrap p-value as proof that Cyp2c38 has no lipid effect.

## Pairing, feature identity and interpretation

Original scaled hepatocyte Cyp2c38 expression is paired with log2 original lipid intensities; reference mice retain their recorded gonadal and chromosome groups. The separate networks describe ovarian or testicular FCG associations, while pooled associations can also contain between-gonad covariance. There is no new KO RNA measurement or MOFA refit. Molecular matching only standardizes the last pipe-delimited annotation and class-parenthesis formatting: chain order, slash versus underscore, ether prefixes, oxygenation and structural resolution remain intact. Internal standards, isotope-labelled features and repository-excluded lipid classes retain the WP1 exclusions. No animal or lipid is removed because its outcome is unfavorable.

KO models use log2 supplied positive intensities, four group means, sex-specific precision weights and robust intensity-trend empirical Bayes. The pooled effect is half the female plus male KO-minus-WT effects; t-statistics are computed from that contrast rather than averaged. Upstream normalization/gap-filling, litter and acquisition-batch provenance remain unavailable. The perturbation is heterozygous in vivo Cyp2c38 reduction; the earlier primary-hepatocyte siRNA experiment is a distinct intervention. Correlation and concordance do not establish direct substrates, flux or the consequences of complete deletion.

## Saved artifacts

- [Complete ULM and p-value summary](ULM_summary.csv), [network export manifest](network_manifest.csv), [ULM concordance figures](plots/ULM_concordance.pdf), [run provenance](provenance.json).
- [Fitted-null dependence diagnostics](null_distribution_diagnostics.csv).
- [Null-calibration summary](diagnostics/calibration_null_rejection_summary.csv) and [complete aligned-effect injection curves](diagnostics/aligned_effect_injection_summary.csv).
- [Ovaries_and_Testes, full network](networks/Cyp2c38_pooled_full_SIF.csv): 1254 molecular edges, 43 paired FCG mice.
- [Ovaries_and_Testes, fixed828 network](networks/Cyp2c38_pooled_fixed828_SIF.csv): 828 molecular edges, 43 paired FCG mice.
- [Ovaries, full network](networks/Cyp2c38_ovaries_full_SIF.csv): 1254 molecular edges, 21 paired FCG mice.
- [Ovaries, fixed828 network](networks/Cyp2c38_ovaries_fixed828_SIF.csv): 828 molecular edges, 21 paired FCG mice.
- [Testes, full network](networks/Cyp2c38_testes_full_SIF.csv): 1254 molecular edges, 22 paired FCG mice.
- [Testes, fixed828 network](networks/Cyp2c38_testes_fixed828_SIF.csv): 828 molecular edges, 22 paired FCG mice.

## Reproduce

Run from the repository root, preserving the completed output with a fresh destination:

```sh
Rscript scripts/workflows/07c_cyp2c38_ulm_resampling.R \
  --bootstraps=4999 --permutations=1999 --seed=20261007 --workers=3 \
  --output=results/validation/cyp2c38_in_vivo/ulm_followup_replay
```

The saved provenance, input hashes, R session details and all resampling draws identify the numerical run. Keep the same null, contrast orientation, molecular matching and feature universes when comparing replays.
