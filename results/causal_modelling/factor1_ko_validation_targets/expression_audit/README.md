# Expression evidence for the Factor 1 TF-to-lipid KO proposal

Reproduce from the repository root:

```sh
Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R --validate-only
Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R --output=NEW_DIR
```

`gene_expression.csv` keeps baseline count evidence distinct from RNA loadings and inferred activity. Full celltype/sex and control-library tables, source hashes and session details accompany it. No new fit or perturbation prediction is calculated.

- FCG samples are 43 unique mice (21 ovarian/22 testicular), excluding the separate five WT additions. Counts are finite nonnegative integers.
- FCG counts are per-mouse/per-celltype pseudobulks; nonzero across mice is not single-nucleus prevalence and does not exclude ambient or contaminating RNA.
- CPM uses total supplied counts per library. FCG and cultured bulk RNA are distinct assays/states; their CPM ratio is not a measured culture-induced fold change.
- Four historical neg2 RNA libraries establish detection in a prior culture system, not biological independence, protein activity, functional isoform, future expression or knockdown efficiency.
- Factor 1 RNA loadings, footprint TF activities and MOON inferred scores are separate fields; none is baseline abundance.
- No RNA expression filter was applied during MOON fitting; frozen RNA TF-target sign coherence does not establish mediator expression or cell autonomy.
- Pld2 and Pparg/Cebpb/Vdr context genes outside the selected graph do not acquire model-supported paths merely by inclusion in this expression audit.
- The upcoming experimental system is unspecified. Historical culture sex/timing should not be imposed on the new experiment.
