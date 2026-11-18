# Preserving measured-species rules during chemical normalization

The maintained scripts are `scripts/moon/ancillary/audit_normalized_rules_topology.py` and `scripts/moon/ancillary/audit_normalized_rules_pruning.R`. The two script copies in this result directory are byte-preserved snapshots whose hashes remain in the original dependency manifests; use the maintained versions for a replay.

From the project root, choose new output directories:

```sh
python3 scripts/moon/ancillary/audit_normalized_rules_topology.py --analysis results/causal_modelling/factor1_lipid_to_tf --output /tmp/normalized_rules_topology_replay
Rscript scripts/moon/ancillary/audit_normalized_rules_pruning.R --analysis=results/causal_modelling/factor1_lipid_to_tf --topology=/tmp/normalized_rules_topology_replay --output=/tmp/normalized_rules_pruning_replay
```

Both maintained scripts refuse occupied output directories and hash their dependencies. No original PKN, fit, or accepted measurement mapping is changed. Migration replays reproduced all ten topology result files and all fourteen pruning/equivalence CSVs byte-for-byte; only the reproduction documentation, script paths and dependency manifests differ.

The original rule-expanded graph is projected through exactly one frozen `chemical_node_aliases.csv` lookup per endpoint. Projecting the original curated graph reproduces the prebuilt pH7.3 curated graph exactly. No accepted input ID is changed. The lipid-only precursor policy is then applied.

The normalized rules restore all 155 accepted input nodes, versus six in the curated graph, and add 898 precursor edges. They do not change any reachable lipid–TF pair, shortest distance, upstream candidate, or signed edge on an input-to-terminal-TF path corridor. Only LPC(16:0/0:0) and LPC(18:0/0:0) reach observed TFs. The other 149 rule-mapped inputs remain without a path to a protein under this policy. Thus omission of the measured-species rules affects input-node presence but does not explain the limited upstream signaling candidates in this particular topology. This audit is not a MOON refit and does not claim equality of scores from independently rerun pipelines.

`cleanup_audit.csv` and `conflicting_signed_pairs.csv` report deduplication, selfloops and opposing signs. `comparison.csv` gives the exact pair/distance/corridor comparison for all 353 mapped TFs and the 49 TFs with absolute activity above two, made terminal separately. `input_reachability.csv` preserves all 155 inputs under both policies. `endpoint_reachability.csv` records each reachable pair; `topology_summaries.csv` also reports distance limits 6, 8, 12 and unbounded. These are topological distances, not iterative MOON pruning or scoring depth.

`pkn_lipid_rules_ph73_precursor.csv` is a derived diagnostic copy of the normalized, filtered rule-expanded graph, not a replacement for the frozen PKN. `added_rule_precursor_edges.csv` isolates its additional edges relative to the normalized curated precursor graph.
