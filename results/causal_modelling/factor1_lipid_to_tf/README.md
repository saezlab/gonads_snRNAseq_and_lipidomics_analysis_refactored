# Factor 1 lipid-to-TF analysis

See the [analysis report](../../../docs/cosmos_factor1_lipid_to_tf.md) for
interpretation and the [active workflow guide](../../../scripts/moon/README.md)
for reproduction commands.

The selected exploratory network is in
[`refinement/selected_all_tfs/`](refinement/selected_all_tfs/): **455 signed edges,
233 nodes, LPC(18:0/0:0) and 64 observed TF activities**. It uses the separately
labelled acid/base-normalized precursor PKN, compressed scoring at pruning
distance 12/scoring layers 12, and strict dual cutoffs 1.65/1.10.

The original PKN has no directly supported measured signalling lipids or
lipid-only routes to TF observations. The positive sensitivity connects
LPC-associated remodelling to unmeasured prostanoids and Ptger3. It requires an
unmeasured arachidonoyl-CoA donor and regenerates LPC; it does not establish
direct LPC receptor activation or net flux. Esr1 and Dbp are supplied downstream
observations, not independently predicted controls.

`analysis_summary.json` records the **coarse** search (72 configurations,
16 fits, 752 threshold evaluations, no compact candidates). The subsequent
`refinement/refinement_summary.json` records 1,431 finer evaluations and
26 compact all-TF candidates. These are successive stages, not contradictory
counts. The 1.65/1.10 choice is a representative, not a biological optimum.

- `coverage/`: exact mapping, direct regulation and reaction-shortcut audit.
- `runs/`: frozen fitted states, compression maps, scores and fitting ledgers.
- `refinement/`: finer extraction, selection and complete graph/evidence exports.
- `score_audit/`: independent final-graph rescoring; all 16 retained tables match
  exactly, with lipid sign exclusions and expected-TF status preserved.
- `normalized_rules_topology/`: restoring the measured-species rule layer gives
  identical prepared graphs and compression for every tested configuration.
- `refinement/selected_all_tfs/reaction_evidence_audit/`: all 51 reaction nodes,
  31 Rhea accessions and the two prostanoid–Ptger3 records audited in full.

Dependency manifests and byte-preserved source snapshots retain execution
provenance. Active ancillary scripts live in `scripts/moon/ancillary/`; snapshots
inside result folders are audit records, not maintained execution entry points.
