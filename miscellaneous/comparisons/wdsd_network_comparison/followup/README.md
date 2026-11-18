# WP2 follow-up: lipid annotation overlap and pooled Figure 3 conservation

Analysis date: 7 October 2026. This addendum addresses the two requested follow-ups to WP2. Original results and the preceding Figure 5 / diet report remain unchanged.

## Main findings

The original 752/1,254 figure is reproducible: 752 original species (60.0%) have the same final molecular annotation in the new export. This is a deliberately strict molecular-label overlap, not an estimate of how many lipids were physically detected. 959 original labels share at least one recorded pipe-delimited alias with a retained new row; these aliases often describe only sum composition. The cumulative audit reaches 963 after also allowing chain-order and positional-resolution candidates. These broader counts are candidates, not validated additional molecular identities.

The original displayed Figure 3 coefficient pattern has pooled FCG-to-SD vector r=0.973 and native ULM=43.628 (nominal p=1.51e-69). The overall mouse-pairing p is 0.0104 and overall null-bootstrap p is 0.0104. Those overall tests allow shared sex differences to explain agreement. When sex margins are retained, the upper-tail pairing p is 0.4256 and upper-tail bootstrap p is 0.4480. These conditional checks do not establish positive alignment beyond that expected from the shared sex structure.

## 1. Why does the lipid overlap appear limited?

The two exports use differing combinations of molecular-species labels, sum-composition aliases, sn-position notation and library identifiers. The original exported table has 1,386 rows; its cleaned network reference retains 1,254 measured rows. WDSD supplies 1,240 rows; the established exclusion of 14 internal standards and 49 rows in repository-excluded classes leaves 1,177 primary features. Multiple measured/adduct rows are never collapsed or summed to improve overlap.

The 1,254-species denominator is the pre-existing cleaned reference, rather than every row in the original 1,386-row export. The 132 raw rows absent from that reference include 111 whose classes match the repository exclusion list. The remaining 21 are 18 SM rows with ontology ASM, two Cer_EOS rows and one Cer_EODS row. This recorded inventory establishes which rows were removed before this follow-up, but does not establish every historical filtering reason. [Raw-to-cleaned ledger](QA/raw_FCG_1386_to_cleaned_1254_ledger.csv) and [excluded-class inventory](QA/raw_FCG_excluded_class_summary.csv) retain the full partition.

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

Audits: [complete original annotation ledger](lipid_overlap_audit/original_ledger.csv), [new annotation ledger](lipid_overlap_audit/new_ledger.csv), [every candidate pair](lipid_overlap_audit/candidate_pairs.csv), [class-level coverage](lipid_overlap_audit/class_summary.csv), [strict safe map](lipid_overlap_audit/safe_map.csv), [conditional compatibility map](lipid_overlap_audit/conditional_compatibility_map.csv) and [structure/identifier conflict examples](lipid_overlap_audit/identifier_structure_mismatch_examples.csv). [Overlap figures](plots/lipid_overlap_audit.pdf) show the loss partition and class distribution.

## 2. Which Figure 3 network was compared?

The current revised Figure 3 A states top three correlated genes per lipid, and panel B labels that species display FULL SPECIES-LEVEL NETWORK. Its Workflow 04 integrated export has 267 edges: 192 species-level gene-lipid edges involving 25 displayed genes and 64 lipids, plus 75 TF-gene edges. This analysis uses the 192 species pairs and excludes TF edges. The figure's use of FULL refers to molecular species resolution; it is distinct from the permissive 33%-tail archive containing 769 genes, 828 lipids and 636,732 pairs. That separate archive is inventoried, rather than substituted for Figure 3.

Original selection uses MOFA Factor 1's top 2.5% in each tail, a lipid-metabolism gene filter, and the three largest absolute gene correlations per selected lipid. It does not impose a 0.65 correlation cutoff: the minimum absolute displayed coefficient is approximately 0.5847. The unthresholded precursor includes 139 selected genes by the same 64 lipids, or 8,896 pairs; it is a complementary universe preserving the original MOFA selection. Neither uses new SD correlations to choose features.

All 43 original FCG mouse pairings are retained: 21 ovaries and 22 testes, with XX/XY and group metadata preserved. Scaled whole-hepatocyte RNA and log2 original lipid intensities replay the displayed and precursor archives to maximum absolute error 9.99e-16. New SD uses 14 distinct CF453 animals, seven recorded females and seven males. There is one whole-hepatocyte RNA profile and one whole-liver lipid profile per mouse. No zone-specific RNA profiles enter this analysis.

Matching is one-to-one through the strict final molecular annotations; missing or invariant features remain recorded in the coverage ledger. The same eligible, originally selected pairs enter every observed comparison and mouse draw. A new correlation is kept even if it is weak or reverses sign. [Reference replay](pooled_reference_replay.csv), [coverage](pooled_coverage.csv), [complete matched coefficients](matched_pooled_edges.csv) and [mouse metadata](metadata.csv) retain those identities.

For the actual Figure 3 displayed network, the strict join retains 36/64 lipids, 108/192 displayed pairs and 15/25 displayed genes. The conditional compatibility sensitivity retains 39/64 lipids, 117/192 pairs and 16/25 genes. These Figure 3 denominators are separate from the 752/1,254 whole-lipid-survey overlap.

## Pooled coefficient conservation

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

A [descriptive sex-mean diagnostic](QA/sex_mean_only_concordance.csv) keeps only the between-sex covariance in SD, normalized by observed pooled feature variation. Its alignment with the original FCG coefficient vector is 0.9745 for displayed pairs and 0.8610 for the precursor, close to observed pooled alignment 0.9733 and 0.8551. This shows why shared sex structure can generate high pooled correspondence. It supplies no new p-value and does not prove that within-sex associations are absent.

The main run uses 19999 random pairing draws per null and 4999 null-bootstrap draws. Each draw is shared across every gene and lipid. Identity pairing is saved separately as draw zero; Monte Carlo p-values use `(1 + exceedances)/(1 + valid draws)`. All tails and valid/invalid counts are in [pooled inference](pooled_inference.csv). Undefined features invalidate their fixed summary vector instead of shrinking its edge set within a draw. The bootstrap's type-I error and power at seven animals per sex have not been calibrated.

Paired stability bootstraps resample matched SD animals within sex. Additional stability draws independently resample original FCG paired mice within their chromosome/gonad groups. They describe sensitivity to mouse sampling and reference uncertainty, while keeping historical feature selection fixed. They do not supply a null p-value and do not propagate uncertainty in the original feature-selection process.

[Coefficient figures](plots/pooled_conservation.pdf) separate displayed and precursor universes. [Null figures](plots/pooled_null_checks.pdf) show why unrestricted and sex-conditioned checks can differ.

## Removing the between-sex means

| Universe | Pairs | Vector r | ULM | Nominal ULM p |
| --- | --- | --- | --- | --- |
| Full sparse precursor | 5004 | 0.082 | 5.796 | 7.19e-09 |
| Figure 3 displayed | 108 | 0.269 | 2.878 | 0.0048 |

This sensitivity subtracts each cohort's female or male mean from every RNA and lipid feature before recalculating Pearson correlations on unchanged pairs. It describes a different residual network, rather than replaying Figure 3. Reported residual ULM p-values remain nominal under dependent edges; no residual-network null p-value is assigned. Pooled conservation can coexist with weak residual conservation because between-sex effects and within-sex mouse variation answer different questions.

## Four predefined Cyp genes

| Gene | Pairs | Vector r | ULM | Nominal ULM p | Within-sex pairing upper p | Within-sex bootstrap upper p |
| --- | --- | --- | --- | --- | --- | --- |
| Cyp2b9 | 33 | 0.978 | 26.310 | 9.10e-23 | 0.4825 | 0.4844 |
| Cyp2b13 | 22 | 0.971 | 18.010 | 7.91e-14 | 0.9855 | 0.5762 |
| Cyp2c38 | 9 | 0.979 | 12.689 | 4.37e-06 | 0.1964 | 0.3950 |
| Cyp2d9 | 23 | 0.984 | 25.486 | 2.80e-17 | 0.0682 | 0.3620 |

The Cyp2b9, Cyp2b13, Cyp2c38 and Cyp2d9 summaries were specified before examining new results and retain their original displayed species membership. Full precursor summaries for those genes are also saved. These secondary summaries are exploratory and their p-values are unadjusted across genes and universes. A gene's pooled coefficient pattern does not identify enzyme activity, substrate specificity or causal regulation.

[Cyp coefficient figures](plots/Cyp_pooled_conservation.pdf) show all available originally displayed pairs. [Individual-mouse examples](plots/selected_pooled_pairs.pdf) show one pair per predefined gene: the matched originally displayed pair with the largest absolute original FCG coefficient, with original lipid-name ordering resolving ties. The new correlation does not enter selection. Seven females and seven males are colored separately, with all mouse IDs shown; pooled and within-sex fits help expose mean-driven associations. The [figure selection ledger](plots/figure_pair_selection_ledger.csv) records these choices.

## Molecular mapping and analytical-quality sensitivities

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

The 854-pair conditional molecular compatibility map is a sensitivity, rather than a revised identity claim. Its original and new coefficients correspond to the explicitly recorded, potentially different molecular annotations, so any agreement is interpreted at the supported coarse annotation resolution. The QC-CV <=30% and supplied Fill % >=0.5 restrictions test dependence on the exported quality fields. They alter the denominator and do not constitute independent biological replication. [Mapping sensitivity](mapping_sensitivity_inference.csv) and [QC sensitivities](pooled_QC_sensitivity.csv) preserve complete results and null definitions.

Supplied positive lipid values are log2 transformed without additional normalization or imputation. Positivity does not establish that upstream missingness or imputation was absent; the supplied Fill % field records incomplete coverage for some rows. Supplied RNA scale.data is used without changing its upstream normalization. Row scaling the stored data layer replays scale.data, which validates the scaling relationship rather than independently validating normalization. Nucleus contributions, detailed upstream nuclear QC, batch and litter covariates remain unavailable.

## Interpretation and reproduction

This pooled test answers whether the sex-pooled, non-zonated Figure 3 pattern transfers to a new standard-diet cohort. It is complementary to WP2's within-sex zonation analysis and should not be described as rescuing a failed zonal result: the feature selection, biological profiles and between-sex contributions differ. Differences in age, background, diet duration and processing also limit direct cohort equivalence; meeting descriptions of those differences are not fully recorded in the supplied metadata.

A defensible manuscript distinction is overall pooled pattern reproduction versus additional correspondence beyond sex margins. Neither establishes direct lipid regulation or a mechanistic pathway. Non-significant conditional tests do not prove lack of a biological relationship; they indicate that this cohort and this conditional test do not establish the specified extra correspondence. Exact-annotation gaps likewise do not prove that the new experiment failed to detect those lipids.

Recorded input limitations:

- The supplied files are processed pseudobulks; per-mouse and per-zone nucleus counts, doublet/QC summaries and composition metadata are absent.
- Only ScaleData.RNA is recorded in each command history. The upstream normalization, missing-value treatment, imputation and nuclear zonation gate cannot be established from these objects.
- Row scaling the stored data layer reproduces scale.data numerically; it supplies a transparent scaling sensitivity but no independent normalization validation.
- Lipid intensity positivity and QC CV describe the supplied export, without establishing raw-data detection, batch correction or normalization provenance.
- No animal or lipid is removed on the basis of a network association, apparent biological outlier or favorable result. The QC-CV <=30% set is a sensitivity only.
- Strict canonical labels preserve chain order, separators, ether prefixes and oxygenation. Formatting agreement is not independent chemical identification; ambiguous annotations must not be forced into matches.
- The complete WDSD intake contains 28 mice with two RNA zones and one whole-liver lipid profile each. This pooled follow-up uses 14 standard-diet mice, whole-hepatocyte RNA and no zone profiles.
- No raw peak absence or limit-of-detection inference follows from absence of an exported annotation.
- Same sum composition, formula or precursor mass cannot distinguish acyl-chain isomers; candidates are not extra validated features.
- The supplied INCHIKEY/SMILES fields can be identical despite different reported chain/ether/positional resolution; these library identifiers are not independent chemical confirmation.
- Underscores normally omit sn-position, but both datasets have multiple CL rows with the same unordered chain quartet; CL group/pair assignment and feature identity are unresolved, so none is automatically rescued.
- Retention time and 5 ppm precursor checks are descriptive diagnostics; no cross-run retention-time calibration or peak/spectrum re-identification was performed.
- Zero-chain removal and slash/underscore conversion lose positional resolution and remain candidates only.
- Pooled standard diet compares 14 wild-type WDSD mice with 43 original FCG mice; age, genetic background and cohort differences cannot be separated.
- Low supplied Fill values may indicate upstream missingness despite positive exported intensities; peak-level absence/imputation provenance unavailable.

Run from the repository root, using a fresh output directory:

```sh
Rscript scripts/workflows/05c_wdsd_lipid_overlap_and_pooled_network.R --validate-only
Rscript scripts/workflows/05c_wdsd_lipid_overlap_and_pooled_network.R --permutations=19999 --bootstraps=4999 --paired-bootstraps=999 --seed=20261007 --output=results/wdsd_network_comparison_followup_reproduction
```

The workflow retains all shared resampling indices and null draws, original feature membership, molecular candidate and exclusion ledgers, input/code SHA256, settings and package/session versions. [Plot manifest](plot_manifest.csv) lists the figure pages. The original analyses and preceding WP2 report are preserved.
