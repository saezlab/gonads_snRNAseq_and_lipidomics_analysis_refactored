# Factor 1 TF to lipid KO validation prioritization

Purpose: derive ten ranked experimental hypotheses from the frozen conservative
TF→lipid fits, combining exact path signs, reaction specificity, source-level
biological auditing and actual expression evidence. This is a read-only proposal
workflow, not a new MOON fit or an intervention simulation.

## Inputs

- Four/six-layer selected network and evidence under
  `results/causal_modelling/factor1_tf_to_lipid_expanded_species/`.
- Full-precision frozen Factor 1 inputs and the one-to-one lipid mapping ledger.
- Actual FCG counts `data/RNA/fcg_pseudobulks.rds` and historical culture counts
  in `data/RNA/legacy/merged_counts_matrix_250508_VH00211_433_AACVFM3HV_with_gene_symbols_and_groups.csv`.
- Reviewed ranking/exclusions `support/cosmos_mouse/factor1_ko_validation_review.json`
  and primary literature `support/cosmos_mouse/ko_validation_literature_sources.json`.

## Entry points

`scripts/moon/ancillary/audit_factor1_ko_expression.R` checks true integer
counts and separates abundance evidence from Factor 1 loadings/activity scores.
It exports 38 genes, 43 FCG hepatocyte mice and four historical culture libraries.
`--validate-only` writes nothing; `--output=NEW_DIR` refuses nonempty destinations.
The full run, validation mode, help and overwrite checks passed on 13 October 2026.

`scripts/moon/ancillary/propose_factor1_ko_validation.py` reads both networks,
removes the declared unsupported edges/catalyst instances in an audit copy,
enumerates signed paths and verifies ten unique reviewed target sets. Local
direction is −1 times the path-sign product. It outputs ranked sets, per-lipid
directions, exact paths, exclusions and source manifests. A full run and
validation-only run passed; no source graph/score/PKN changes occur.

## Result and limitations

Output root: `results/causal_modelling/factor1_ko_validation_targets/`.
Report: `docs/cosmos_factor1_ko_validation_targets.md`. Ten designs use eight
biological targets; four additional constituent singles are required to
interpret every pair/triplet. Directions are conditional local effects, not
concentration/flux predictions or synergy estimates.

Reject Ppara→Plpp1, Cebpa→Pld1 and Gata6→Mapk1 as source-supported experimental
chains because of cited paralog/direction mismatches. Reject Enpp2 diacyl-PC→PA
assignment. Exclude the retracted RXRA/PLD1 citation; the remaining VDR/RXRA
keratinocyte evidence requires a hepatocyte context gate. Prioritize cPLA2
LPC18 readouts rather than all PUFA-retaining generated LPC rules. PLPP1 pool
accessibility and PLC physiological PIP2 versus measured PI remain conditional.

The proposal does not recompute scores after edge exclusions. Historical MOON
figures and fitted graphs remain computational artifacts with these newly
identified biological-source limitations. The primary experiment choices are
therefore chemical anchors and conditional mediation tests, not ten validated
upstream TF pathways. Any repair of the PKN and refitting is separate work.
