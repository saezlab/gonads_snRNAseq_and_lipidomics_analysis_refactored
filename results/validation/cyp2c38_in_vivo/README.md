# Cyp2c38 in vivo heterozygous knockout: lipidomics and predictive validation

Completed run: 2026-10-06 14:58:42 UTC. All numerical statements below are read from the saved result tables.

**This dataset does not establish robust in vivo validation of the fixed Cyp2c38–lipid predictions.** None of the three primary animal-calibrated concordance tests passes BH FDR 0.05, and no individual primary lipid contrast produces an FDR discovery. This is inconclusive evidence, not proof that Cyp2c38 has no lipid effect.

## Primary results

The validation statistic compares the original gene–lipid Pearson coefficient with the new KO-minus-WT moderated lipid t-statistic. Negative Pearson r has the predicted perturbation direction. The signed concordance score is **−r**, so positive scores support agreement. Bootstrap intervals below apply to that score, not directly to the displayed r.

| KO contrast / original reference | Pearson r | Signed score | Combined-bootstrap 95% score interval | Animal permutation p, two-sided | BH FDR, three tests |
| --- | ---: | ---: | --- | ---: | ---: |
| Pooled KO − WT / pooled reference | -0.063 | 0.063 | [-0.294, 0.334] | 0.8480 | 0.8480 |
| Female KO − WT / female reference | -0.064 | 0.064 | [-0.258, 0.301] | 0.7220 | 0.8480 |
| Male KO − WT / male reference | -0.167 | 0.167 | [-0.209, 0.395] | 0.6215 | 0.8480 |

Each primary comparison uses 552 fixed shared lipid features. All three intervals include zero in this run. These intervals independently resample animals in both the KO and original paired reference cohorts; they remain conditional on the established feature selection.

| Primary contrast | Features tested | FDR ≤0.05 | Higher after KO | Lower after KO |
| --- | ---: | ---: | ---: | ---: |
| Female KO − WT | 1237 | 0 | 0 | 0 |
| Male − female KO effect | 1237 | 0 | 0 | 0 |
| Male KO − WT | 1237 | 0 | 0 | 0 |
| Pooled KO − WT | 1237 | 0 | 0 | 0 |

FDR is BH-adjusted across all tested lipid features separately within each differential contrast. The interaction tests whether the male KO effect differs from the female KO effect; it is exploratory and is not a fourth predictive-validation test. Nominally small feature P values are not counted as discoveries.

## Sensitivities and influence

| Differential sensitivity | Features | Contrasts available | Maximum FDR discoveries in any contrast |
| --- | ---: | ---: | ---: |
| Unweighted joint log2 | 1237 | 4 | 0 |
| QC CV ≤30%, weighted log2 | 1215 | 4 | 0 |
| Supplied-scale weighted joint | 1237 | 4 | 0 |
| Separate-sex log2 | 1237 | 2 | 0 |

Every sensitivity model also has zero lipid discoveries at FDR ≤0.05. The QC restriction leaves 545 of the 552 shared reference features. Sensitivity concordance is descriptive: the animal-calibrated three-test family is the primary weighted-log2 analysis, and favorable secondary scores are not substitutes for its outcome.

The supplied-scale sensitivity uses the original positive numerical intensities without log2. Although limma names its coefficient column `logFC`, its values in that sensitivity are **linear mean-intensity differences**, not log fold changes. These units must not be pooled with or exponentiated as log2 effects.

| Main signed score | Range omitting one KO-cohort animal | Range omitting one paired FCG animal |
| --- | --- | --- |
| Pooled KO − WT | [-0.103, 0.130] | [0.054, 0.074] |
| Female KO − WT | [-0.017, 0.143] | [-0.033, 0.195] |
| Male KO − WT | [-0.022, 0.234] | [0.133, 0.211] |

Omitting male control **Cyp-15_M-Ctrl** changes the pooled score to -0.103 and the male sex-matched score to -0.022, reversing their observed positive signs. This is an observational influence finding, not evidence that the animal is technically invalid. It remains in all primary analyses. Female reference leave-one-out scores also cross zero. Sample identities and all omission results are retained for review.

The male sex-matched ULM score is -3.961 with edge-level p = 8.45e-05. That apparent significance treats lipid targets as independent and is **descriptive**. The whole-animal permutation p for the corresponding main comparison is 0.6215; it governs the inferential conclusion.

Secondary sex-specific effects versus the pooled reference are also reported: Female KO − WT: r = 0.277, score = -0.277, two-sided animal p = 0.1345; Male KO − WT: r = -0.243, score = 0.243, two-sided animal p = 0.4600. These comparisons have no BH value in the primary three-test family. Their direction differs between sexes; reporting both prevents selecting only the more favorable pooled-reference comparison.

## Samples, matching and preprocessing

There are **32 biological mice**: female WT 8, female KO 8, male WT 7 and male KO 9. The supplied control labels are interpreted as WT as described in the meeting. Analytical QC injections are not biological replicates. From 1319 supplied measurement rows, exclude 16 explicit internal standards, 1 additional isotope-labelled feature and 65 rows in the existing Workflow 07 excluded classes. This leaves **1237 primary features**; **1215** meet QC CV ≤30% across nine analytical QC injections.

Primary models use log2 supplied positive intensities. No new normalization, gap filling or imputation is applied. The prior normalization and processing provenance is unavailable; finite positive exported values do not prove that they were never imputed upstream. Retain complete original annotations, QC CV, feature IDs and the biological sample manifest.

The fixed prediction universe contains **828 Cyp2c38 edges** from the original global full MOFA export. Of these, **552** have unique exact KO annotation matches and 276 remain unmatched. Matching only takes the last pipe-delimited annotation and converts class/space formatting to parentheses. Chain order, slash versus underscore, ether prefixes, oxygenation and sum-versus-resolved molecular structure are preserved. Two cardiolipin reorder candidates remain excluded; no ambiguous structure is forced into a match.

Original reference measurements are already scaled hepatocyte RNA paired with log2 lipidomics from **43 FCG mice: 21 ovary-bearing and 22 testes-bearing**. Pooled coefficients reproduce all archived edges with maximum discrepancy **5.55e-16**, below the required 1e−12 tolerance. Within-gonad coefficients are recomputed for the same universe. These estimate a different association and are generally weaker than pooled coefficients, which can include between-gonad covariance. There is no KO RNA dataset and no MOFA refit.

## Models and animal-level uncertainty

A joint four-group limma design estimates F_WT, F_KO, M_WT and M_KO means. Female effect is F_KO−F_WT; male effect is M_KO−M_WT; pooled effect is **0.5×(female effect+male effect)**; interaction is male effect−female effect. The pooled estimand gives the sexes equal weight despite unequal male genotype counts. It is not an average of t-statistics. Sex-group precision weights from `arrayWeights(var.group=sex)` are 1.4594 for female and 0.6852 for male samples. These weights affect precision, not the pooled estimand. Weighted/unweighted group-mean coefficients agree to 2.13e-14.

Models use robust, intensity-trend empirical Bayes moderation. Full tables contain effects, confidence intervals, moderated t, P values and BH FDR. Separate-sex sensitivity fits moderate each sex independently. Per-lipid signed consistency is `−reference_r × KO_t`; it is an agreement score, not a direct causal or enzyme-substrate estimate.

For inference, **1999 genotype permutations within sex** preserve the original counts and refit precision weights and moderation. Two-sided empirical p values use the +1 numerator/denominator correction; BH adjusts only the three primary tests. **999 combined whole-animal bootstraps** independently resample KO mice within four sex/genotype groups and paired FCG mice within XXO=10, XYO=11, XXT=11 and XYT=11. Reference RNA and lipids remain paired. Established feature selection is held fixed. Percentile 2.5%/97.5% score limits are conditional intervals.

For the three main comparisons, invalid permutation draws are 0/0/0 and invalid bootstrap draws are 0/0/0 (pooled/female/male). Seed 20261006 and 3 workers are recorded; each draw has a saved seed. Whole-animal leave-one-out checks refit KO models or recompute reference correlations as appropriate.

Permutation inference assumes genotype labels are exchangeable within sex. Unavailable litter/batch metadata cannot be controlled. Small cohorts, heterozygous protein reduction, upstream processing uncertainty and dependent lipid species limit interpretation. A null result here neither establishes biological equivalence nor rules out smaller, context-dependent or indirect effects. This in vivo heterozygous experiment is distinct from Workflow 07’s earlier siRNA knockdown in primary hepatocytes. Neither experiment alone identifies direct substrates, metabolic flux or the full impact of complete deletion.

For context, the [existing Workflow 07 siRNA scores](../global_network_validation_scores.csv) contain an on-target Cyp2c38 row with ULM score -4.00009 and p = 7.05e-05. This is also an edge-based descriptive p value. Different interventions and experimental settings mean comparing these p values does not quantify relative validation strength.

## Artifacts

- Inputs and QC: [sample manifest](sample_manifest.csv), [complete annotation/QC ledger](feature_annotation_and_QC_ledger.csv), [log2 analysis matrix](log2_analysis_matrix.csv), [reference matching ledger](fixed_reference_and_matching_ledger.csv), [paired FCG mouse manifest](FCG_paired_reference_mouse_manifest.csv), [contrast matrix](contrast_matrix.csv).
- Differential results: [complete primary table](differential_primary_all.csv), [t-statistic matrix](t_statistic_matrix.csv), [discovery summary](differential_discovery_summary.csv), [pooled](differential_pooled_KO_vs_WT.csv), [female](differential_female_KO_vs_WT.csv), [male](differential_male_KO_vs_WT.csv), [interaction](differential_interaction_M_minus_F.csv), [all sensitivities](differential_sensitivities_all.csv), [saved model fits](limma_fits.rds).
- Validation: [main/secondary concordance](concordance_primary_and_secondary.csv), [sensitivity concordance](concordance_sensitivities.csv), [per-lipid predictive consistency](per_lipid_predictive_consistency.csv), [sensitivity consistency](per_lipid_consistency_sensitivities.csv), [descriptive signed ULM](signed_ULM_descriptive_scores.csv), [descriptive lipid-class summary](lipid_class_descriptive_summary.csv).
- Uncertainty: [all permutation draws](permutation_draws.csv), [all combined-bootstrap draws](combined_bootstrap_draws.csv), [KO omission checks](leave_one_KO_mouse_out.csv), [FCG omission checks](leave_one_FCG_mouse_out.csv).
- Figures: [QC/PCA and model diagnostics](plots/QC_PCA.pdf), [differential analysis](plots/differential.pdf), [predictive validation](plots/concordance.pdf).
- Provenance: [run settings/checks/package versions](provenance.json), [input/code SHA256 manifest](input_and_code_SHA256.csv), [R session information](sessionInfo.txt).

## Reproduction

Run from the analysis repository root with the recorded packages available. Use a fresh replay destination to preserve this completed run:

```sh
Rscript scripts/workflows/07b_cyp2c38_in_vivo_validation.R \
  --permutations=1999 --bootstraps=999 --seed=20261006 --workers=3 \
  --output=results/validation/cyp2c38_in_vivo_replay
```

Run the meaningful model/sign/invariant checks separately:

```sh
Rscript tests/test_cyp2c38_in_vivo.R
```

The workflow’s `--validate-only` option checks intake and original-reference reconstruction without fitting KO effects. The saved SHA256 manifest and provenance file identify the exact inputs, code, software versions and resampling settings. Retain matching/reporting conventions when interpreting a replay.
