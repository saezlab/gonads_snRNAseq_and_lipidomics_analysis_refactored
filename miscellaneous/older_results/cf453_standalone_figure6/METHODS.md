# CF453 Figure 6 methods

## Design and inputs

The standalone analysis uses 28 paired CF453 wild-type mice, seven per sex/diet group. Each mouse supplies periportal and pericentral RNA pseudobulks and one whole-liver lipid measurement. Analytical QC injections are not biological replicates. Input files are the local CF453 lipid CSV and the updated whole-hepatocyte/zonation Seurat objects. All matching uses explicit mouse IDs; zones remain paired within mouse. No other cohort supplies model parameters, feature selection or test statistics.

The provenance audit confirms data = ln(1 + 10000 × supplied counts / profile library), and scale.data is a per-gene z-score of data. Maximum normalization replay errors were 4.441e-16 and 8.882e-16 for hepatocyte and zone objects, respectively. This establishes the stored transformation, not the original count-summing provenance.

## RNA model and testing families

Genes require positive supplied counts in at least seven of 56 zone profiles and expression SD >10^-8 across profiles, independent of lipid outcomes; 18,486 pass. Primary expression is the supplied unscaled natural-log normalized data. An eight-cell sex × diet × zone design is fit with limma, a common duplicateCorrelation mouse block (estimated correlation 0.71605), and mean-variance-trend, robust empirical Bayes. RNA effect units are differences in log-normalized expression, not log2 fold changes. Model-based plotted ranges are nominal 95% intervals from moderated standard errors, without multiplicity adjustment.

Primary sex-by-zone contrasts are (male PP−PC)−(female PP−PC) separately in SD and WD; the three-way test is interaction_WD−interaction_SD. BH adjustment pools all eligible genes and all three primary contrasts. Four simple zonation contrasts form a separate family; the additional families below remain distinct.

| family | features | contrasts | tests | q_005 |
| --- | --- | --- | --- | --- |
| primary_interactions | 18,486 | 3 | 55,458 | 1,337 |
| secondary_zonation | 18,486 | 4 | 73,944 | 13,908 |
| secondary_diet_by_zone | 18,486 | 4 | 73,944 | 2,459 |
| secondary_diet_zonation | 18,486 | 2 | 36,972 | 534 |
| secondary_sex_by_zone | 18,486 | 4 | 73,944 | 4,490 |

Gene-context membership requires both diet-specific primary interaction q≤0.05 and the corresponding within-sex/diet zonation q≤0.05. Positive PP−PC effects map to PP and negative effects to PC; multi-context memberships are retained. A context name denotes membership under this rule, not absence of expression or association in another context.

## Lipid model and quality sets

Internal standards, isotope-labelled annotations and documented repository-excluded classes are removed, leaving 1,177 broad features. Primary analysis requires analytical QC CV≤30% (1,129 features); the further Fill≥0.5 sensitivity retains 1,058. Supplied positive intensities are log2 transformed without additional normalization or imputation. A four-cell sex × diet limma model tests male−female within diet, WD−SD within sex, and their interaction. Lipid coefficients are differences of log2 abundance. The two sex contrasts share one BH family, two diet contrasts another, and the interaction a third.

| family | features | contrasts | tests | q_005 |
| --- | --- | --- | --- | --- |
| lipid_sex | 1,129 | 2 | 2,258 | 1,095 |
| lipid_diet | 1,129 | 2 | 2,258 | 1,692 |
| lipid_interaction | 1,129 | 1 | 1,129 | 267 |

Chemical annotations retain feature identity and molecular resolution; no cross-cohort overlap restriction, abundance-based outlier exclusion or forced annotation merging is used. QC and Fill sensitivities change feature eligibility and denominators, not the selected RNA programs.

## Associations, summaries and selection

Within each sex/zone, candidate genes are the union of SD and WD nominations. Genes and lipids must vary in both diets. Pearson correlations use seven matched mice per sex/diet and identical pair universes across diets. Nominal correlation probabilities use t with five degrees of freedom. One BH family covers all 2,750,244 primary pair/context tests; none pass q≤0.05. These probabilities remain conditional on the RNA selection. Broad-only lipid pairs retain coefficients and omission estimates but have no primary-family q-value.

The exploratory display rule is observed |r|≥0.65 with an active RNA nomination in that diet. A stable displayed edge additionally has all seven leave-one-mouse-out correlations defined and of the observed sign; each omission uses six mice. Omission stability does not require every omitted correlation to exceed 0.65.

999 ordinary paired-mouse bootstrap draws are generated separately within each sex/diet, using the same draw indices across both zones and every RNA/lipid feature. RNA and its matched whole-liver lipid vector travel together. Percentile ranges (2.5th–97.5th), sign fractions and valid-draw counts describe fixed-selection stability; they are not null-test p-values or calibrated post-selection intervals. Constant-feature draws are undefined and reported. WD-minus-SD delta-r draws combine independent diet-group resamples; positive fractions are descriptive stability summaries.

For each sex/zone, representative genes rank by the number of distinct stable displayed lipid partners across diets, then the smallest diet-specific interaction q, then gene symbol. At most one representative with positive stable degree is selected. Up to six partners stable in at least one diet are ranked by their maximum absolute correlation across diets, then lipid ID. Both diets and all corresponding coefficients are retained. The leading partner supplies the individual-mouse example; this is a data-selected illustration.

C assigns each lipid to the sex/zone contexts containing any qualifying displayed edge among diet-active candidates. Exactly one context denotes a unique membership, more than one denotes shared, and zero denotes unassigned. The class denominator is every eligible measured lipid in that class, including unassigned lipids. C_pair_density instead counts |r|≥0.65 pairs over all eligible union-candidate-gene × lipid pairs within each class and sex/zone, using the same pairs in both diets; it does not require an active nomination in both diets. Different contexts may have different union sizes, but each context's SD/WD comparison uses the same denominator. Stable fractions, FDR fractions and mean/median |r| are also exported.

The display-threshold audit uses 0.50, 0.60, 0.65, 0.70 and 0.80 without selecting a preferred threshold from outcomes. Seed: 20261030; requested bootstrap draws per sex/diet: 999.

## Interpretation

This is an observational analysis of sex associations in one CF453 cohort; it does not isolate gonadal hormone or sex-chromosome effects. One whole-liver lipid vector is paired with both RNA zones per mouse: RNA-zone associations do not measure spatial lipid concentrations, enzyme activity, substrate specificity or flux. Nucleus contributions, original nuclear QC/gating, count aggregation, batch/litter and assignment provenance remain unverified. Seven mice per sex/diet limit correlation precision. Representative selection uses these same data, and omission/bootstrap checks condition on the fixed RNA and partner selections. Stable signs, large correlations and bootstrap ranges do not provide independent replication or confirm causal regulation.
