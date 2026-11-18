# Independent scientific interpretation review

Reviewed on 7 October 2026 against production `paired_network_inputs.rds`, the complete coefficient and inference CSVs, the frozen Workflow 04 species export, and the draft follow-up README. This review does not alter observations, sample inclusion or network membership.

## Numerical and selection checks

The script [independent_interpretation_numeric.R](independent_interpretation_numeric.R) independently recomputes every one of the 5,004 strict precursor pairs using `stats::cor`, including pooled SD/FCG and sex-centered SD/FCG coefficients. The maximum discrepancy from the production coefficients is 1.32e-15. Ordered RNA/lipid mouse identities agree with metadata. The 108 retained displayed pairs are the exact originally displayed pairs, with no new SD correlation threshold.

The original Figure 3 species display contains 192 edges, 25 genes and 64 lipids. It selects three original absolute correlations per lipid from the frozen 139-gene by 64-lipid sparse precursor. The minimum original displayed absolute correlation is 0.5846677; five original edges are below 0.65. Applying an additional 0.65 filter would change the requested historical network. Strict matching retains 36/64 lipids, 108/192 displayed pairs and 15/25 displayed genes. The conditional compatibility sensitivity retains 39/64 lipids, 117/192 pairs and 16/25 genes. The broader 752/1,254 annotation-overlap survey has a different denominator; the 769-gene by 828-lipid permissive archive is also a different reference universe.

## Why the pooled result and sex-conditioned result differ

The pooled displayed-vector association is 0.9732661, and the sparse-precursor association is 0.8551496. Sex-conditioned permutation-null medians are 0.9727288 and 0.8419613, respectively. Thus the sex-conditioned null does not center near zero: strong pooled concordance is compatible with the between-sex structure even when individual RNA–lipid pairing information within sex is removed.

A deterministic diagnostic independently verifies this point. For gene expression x and lipid abundance y, retaining the actual pooled variance denominators, the expected coefficient under uniformly random within-sex pairing is

```
r_sex_mean_only = sum_s n_s (mean(x_s) - mean(x)) (mean(y_s) - mean(y)) /
                  sqrt(sum_i (x_i - mean(x))^2 * sum_i (y_i - mean(y))^2)
```

This expectation follows because the pairing permutation preserves each sex's means and total sums of squares, while the within-sex cross-product has expectation zero. Correlating those expected SD coefficients with the frozen FCG coefficients gives 0.9745132 for the displayed pairs and 0.8610483 for the precursor. These values describe `cor(reference, E[SD coefficient])`; they are not the expected value of the nonlinear vector-association statistic and supply no additional p-value. Exact values are saved in [sex_mean_only_concordance.csv](sex_mean_only_concordance.csv).

It is defensible to describe strong reproduction of the sex-pooled Figure 3 pattern in the new standard-diet cohort. The reported conditional tests do not establish additional alignment beyond shared sex structure. Those statements answer different biological questions and can both be true. Removing the sex means produces a different residual network; its modest overall vector associations do not disprove sex-mediated biology or a gene's regulatory role. Pooled associations likewise do not establish substrates, enzymatic activity, flux or causality.

The historical MOFA selection emphasizes Factor 1 and therefore sex-associated features. Keeping that selection fixed avoids selecting on the new SD outcomes, but it does not remove shared sex structure from the replication question. All four predefined Cyp summaries and both edge universes are secondary exploratory comparisons; unadjusted edge-level ULM probabilities should not be interpreted as independent biological replication of thousands of pairs. There are 14 SD mice, not 108 or 5,004 experimental units.

## Annotation interpretation

The supplied original FCG annotation export contains 1,386 rows; the existing cleaned reference contains 1,254. All cleaned labels are uniquely present in the supplied export. Of the 132 exported rows absent from that reference, 111 have classes appearing in the repository's exclusion list: CAR 49, SL 24, DGGA 10, PMeOH 8, LPE-N 7, PEtOH 5, SPB 3, BA 2, HBMP 1, NATau 1 and VAE 1. The other 21 consist of 18 SM rows with recorded ontology ASM and three Cer rows with ontology Cer_EOS (two) or Cer_EODS (one). None of these 132 has an explicit internal-standard annotation. These are recorded annotation categories; the historical reason for each exclusion was not reconstructed. The [complete denominator ledger](raw_FCG_1386_to_cleaned_1254_ledger.csv), [class/ontology summary](raw_FCG_excluded_class_summary.csv) and [reproducible inventory script](raw_FCG_denominator_audit.R) preserve this distinction. The strict map remains unchanged.

Exact annotation overlap remains 752 against both all 1,240 supplied new rows and the 1,177 retained rows. Intake exclusions therefore do not explain the 502 originally unmatched annotations. Restricting the complete exact map to QC CV<=30% retains 745, because seven exact-matched lipids exceed the threshold; all 36 strict Figure 3 lipids pass that QC restriction. This additional QC loss does not explain why those 502 labels were unmatched in the primary analysis.

The literal recorded-alias inventory covers 959 original labels; cumulative coverage including positional candidates reaches 963. These are annotation candidates, not additional proven molecular identities. The separate 854-row sensitivity includes 102 bidirectionally unique coarse-alias/same-adduct/5 ppm candidates, blocked from reusing any strict-matched feature. It compares potentially different molecular annotations under a compatibility assumption. It must not be presented as a corrected count of identical molecular species or as class-summed abundance analysis.

LIPID MAPS distinguishes unknown versus known sn positions, sum composition versus constituent chains, and CL DG-moiety resolution. The CL chain-sorting caution is additionally grounded in these actual exports: multiple distinct original/new rows share one unordered chain quartet. That observed collision prevents a unique automatic molecular-feature rescue. Nomenclature alone does not prove that each recorded underscore order encodes a resolved CL pairing; those pairings remain uncertain.

The O/P distinction needs careful wording. An O- annotation can have less structural resolution than a P- annotation, and the P- vinyl bond is not counted in the chain DBE. Consequently, `PE(O-18:1_18:2)` and `PE(P-18:0_18:2)` are not necessarily different chemical entities. Equality of exported InChIKey/SMILES cannot independently settle those resolution differences. The CL example `CL(16:0_16:1_18:1_18:2)` versus `CL(16:0_18:1_18:1_16:2)` does contain different reported chain constituents despite identical supplied identifiers. Prefer “conflicting or differently resolved annotations” to asserting that every identifier-matched label pair is chemically contradictory.

Absence of an exported annotation cannot establish absence of a raw peak or an abundance below detection. Cross-run retention times and raw MS/MS assignments have not been revalidated. A targeted annotation reconciliation with the collaborator should begin with the candidate ledgers and preserve unresolved one-to-many mappings.

## Wording clarifications communicated and incorporated in the report

- Within-sex permutations preserve the observed group means exactly. Within-sex independent replacement bootstraps preserve sex membership and draw from each sex's empirical margins, but individual draw means vary. They should not both be described as preserving observed means exactly.
- State both annotation-survey coverage (752/1,254) and Figure 3 coverage (36/64 lipids, 108/192 pairs, 15/25 displayed genes).
- Qualify the complete WDSD intake's 28 mice separately from the 14 SD animals used in this whole-hepatocyte follow-up; no zonal RNA profiles enter the pooled comparison.
- Use structural-resolution language for O/P pairs, without automatically equating or declaring them different molecular identities.

## Verified primary sources

- [Liebisch et al., updated shorthand notation, Journal of Lipid Research 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7707175/) — source for the structural-resolution hierarchy and underscore/slash distinction. The same publication is available through [the journal](https://www.jlr.org/article/S0022-2275%2820%2960017-7/fulltext).
- [LIPID MAPS phospholipid and lysophospholipid rules](https://www.lipidmaps.org/lipid_nomenclature/glycerophospholipids/pl_lpl) — directly states unknown/known sn positions, CL DG-moiety examples, lyso zero-chain notation and O/P bond conventions.
- [LIPID MAPS structural hierarchy](https://www.lipidmaps.org/lipid_nomenclature/rules/hierarchy) — explains sum-composition, molecular-species and sn-position levels.
- [Fiehn Lab, UC Davis MS adduct calculator](https://www.fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator) — supports the ion-adduct shifts used for the neutral-mass diagnostic. Those mass checks are diagnostic rather than peak re-identification.

The independent check files preserve the original inputs. No manuscript wording should imply that the conditional non-significant tests prove absence of a biological relationship, or that the large pooled ULM statistics prove direct regulation.
