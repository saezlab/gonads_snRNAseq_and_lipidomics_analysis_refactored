# Actual MOON preparation and compression equivalence

Run `Rscript scripts/moon/ancillary/audit_normalized_rules_pruning.R --analysis=ANALYSIS_DIRECTORY --topology=TOPOLOGY_AUDIT_DIRECTORY --output=NEW_EMPTY_DIRECTORY` after the maintained `scripts/moon/ancillary/audit_normalized_rules_topology.py` topology audit. Occupied output directories are refused. The scripts stored beside these results are byte-preserved historical snapshots, not the maintained entry points.

This audit calls the same shared moon_prepare and moon_safe_compress functions as the lipid-to-TF runner, under its frozen software/input contract. It performs no MOON scoring and does not overwrite any existing fit.

For both TF observation panels (all 353 mapped TFs and the 49 with absolute activity above two), and pruning steps 6, 8 and 12, the normalized curated and rule-expanded graphs produce identical canonical signed pruned edges and identical upstream/downstream numeric vectors. The compressed scoring edge sets and complete compression maps are also identical. The pruned states match every available saved normalized-curated prepared state.

The 149 restored measured-species inputs have no route to a protein or observed TF under the lipid-only precursor policy. The shared alternating forward/backward pruning removes those disconnected additions before scoring. Restoring them therefore cannot change these fits: the identical scoring inputs, protected nodes, compression and shared RNA/regulon evidence imply the same deterministic MOON/coherence computation for each fixed scoring-depth/threshold configuration. This is equivalence of the computation supplied to scoring, not a claim that a second fit was performed.

pruning_equivalence.csv records exact equality checks, SHA256 hashes, edge/node/input counts and comparisons against saved states. Edge hashes use canonical sorted source/sign/target tab-delimited lines without a trailing newline. Numeric-vector hashes sort by node and use 17 significant digits. pruning_equivalence_ledgers.csv preserves the actual alternating-pruning iterations for both starting graphs. Dependency hashes and sessionInfo record the implementation and environment.

The full rule-expanded normalized graph remains a derived diagnostic copy, and the 155 accepted measured identifiers are unchanged. No original PKN or fit was modified.
