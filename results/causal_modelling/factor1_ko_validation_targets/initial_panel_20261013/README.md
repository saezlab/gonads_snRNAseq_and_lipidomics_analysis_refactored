# KO validation proposal for Factor 1 TF to lipid hypotheses

Read [the ranked report](../../../docs/cosmos_factor1_ko_validation_targets.md). Ten reviewed sets prioritize direct biochemical anchors and conditional combinations. Rank 1 Pld1, rank 2 Pld1/Pld2 and rank 3 Pla2g4a are the first pilot choices. Ranks 9–10 RXRA require a receptor-context and PLD1-mediator gate. This is a source-audited proposal, not a new MOON fit or KO simulation.

`prioritization/` contains ranked_target_sets.csv, per-readout expected_lipid_directions.csv, exact audited_paths.csv/path_audit.json, unsupported-edge exclusions, review snapshot, validation summary and source hashes. `expression_audit/` contains actual-count expression/context evidence from 43 FCG hepatocyte pseudobulks and four historical culture libraries.

Known source problems: Ppara→Plpp1 cites PPARG; Cebpa→Pld1 cites CEBPB; Gata6→Mapk1 reverses the cited kinase/TF relationship; Enpp2 diacyl-PC→PA is unsupported. The retracted RXRA/PLD1 citation is excluded; the remaining VDR/RXRA evidence is contextual. Original fitted graphs/scores remain frozen and are not recalculated after these exclusions.

The ranking deliberately repeats targets to interpret combinations. Four additional constituent singles bring the ten designs to 14 distinct target-containing conditions; estimating a full three-way interaction requires an additionalpair. No additivity, synergy or steady-state lipid direction is guaranteed.
