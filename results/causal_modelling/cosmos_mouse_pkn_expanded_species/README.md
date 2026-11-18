# Expanded measured species reaction PKN

Separate exploratory expansion of the frozen mouse `pkn_lipid_rules`, built on 6 October 2026. See `docs/cosmos_lipid_species_expansion.md` for methods, interpretation and reproduction.

Chain-compatible PKN: 331,032 edges and 419 candidate measured product features. Glycan-conditional PKN: 331,472 edges and 431 candidates, with separate glucose/galactose identities. The 3,230 new hypotheses preserve complete directional mouse Rhea chemistry. They are not independent experimental evidence.

No MOON scoring or input-value assignment was performed. The product ledger is not an input vector. Positive participation does not enforce concentrations, flux or joint substrate availability. Unfiltered reachability is not the selected workflow's RNA-filtered preparation or final solution.

```bash
python3 scripts/pkn/build_expanded_lipid_rules.py --output /tmp/expanded_species_replay
python3 scripts/pkn/validate_expanded_lipid_rules.py \
  --output /tmp/expanded_species_replay --report /tmp/expanded_species_validation.json
Rscript scripts/pkn/validate_expanded_lipid_rules.R /tmp/expanded_species_replay
```

Choose an empty replay directory.
