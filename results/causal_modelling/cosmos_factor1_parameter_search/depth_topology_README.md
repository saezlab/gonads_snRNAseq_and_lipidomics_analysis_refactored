# Factor 1 depth topology audit

**Depth 5 already covers every directed TF–lipid pair reachable in the full
terminal-observation PKN. Larger depths add graph context and cycles, without
adding reachable input pairs or shortening their routes.** This is a topology
result, not a comparison of MOON scores or biological plausibility.

| Depth | Edges | Nodes | TFs | Lipids | Reachable TF–lipid pairs |
|---:|---:|---:|---:|---:|---:|
| 4 | 18,423 | 2,226 | 45 | 72 | 3,240 |
| 5 | 35,036 | 5,114 | 49 | 74 | 3,626 |
| 6 | 44,127 | 6,673 | 49 | 74 | 3,626 |
| 7 | 46,161 | 7,284 | 49 | 74 | 3,626 |
| 10 | 47,262 | 7,484 | 49 | 74 | 3,626 |
| 14 | 47,285 | 7,497 | 49 | 74 | 3,626 |

The complete 4–14 range is in `depth_topology_audit.csv`. Depth 6 exactly
reproduces the previously saved pruned PKN by signed edge identity.

Across depths 5–14, the all-pairs shortest-path distribution is unchanged:
107 pairs need three edges, 1,223 need four, 1,644 need five, 650 need six, and
two need seven. The seven-edge pairs are Pitx3 to LPC(14:0/0:0) and
LPC(17:0/0:0). Edge counts include explicit enzyme-reaction nodes and are not
counts of biochemical reactions. Every retained TF reaches at least one lipid
within five edges; every retained lipid has at least one TF within five edges.

Increasing depth from 4 to 5 restores Dlx3, Hoxd13, Pitx3, Rorc,
LPC(14:0/0:0), and LPC(17:0/0:0), adding 386 reachable pairs. Increasing
depth from 5 to 14 adds 12,249 edges; 12,244 of those lie on a directed cycle
at their depth of introduction. It adds no new reachable pairs, no shorter
pairwise paths, and no additional edges on the union of endpoint-minimal
routes. Cyclic edges may still represent real feedback or alternative routes;
topology alone cannot dismiss them or establish mechanistic support.

Depths 5, 6, and 7 therefore provide a focused sensitivity range: 5 reaches all
accessible inputs, 6 reproduces the pilot, and 7 accommodates every pair's
shortest length numerically. Depths 8–14 predominantly test expansion of cyclic
context. No claim is made that these depths optimize inference quality.

## Method and limits

Inputs are unchanged: 49 PKN-present TF activities with `abs(activity) > 2` and
155 mapped nonzero lipid weights. All outgoing edges from those 155 observed
lipids are removed first. The existing runner's alternating forward/backward
depth-limited induced-subgraph pruning is then repeated to a fixed point.
There is no compression, RNA coherence filtering, MOON scoring, sign-consistency
selection, or expression-presence filtering in this audit. Edge signs remain
stored, but reachability and cycle membership ignore signs.

This pruning depth is **not** a hard maximum TF–lipid path length: each retained
node is near some TF and some lipid, so an individual pair can require more
than the depth parameter. The CSV explicitly reports such pairs. Likewise,
these shortest paths are not MOON inference levels.

The unrestricted terminal-observation graph can connect only 74 of 155 mapped
lipids; the remaining 81 cannot be recovered by increasing this pruning depth.
This limit depends on the selected TF seeds, PKN directionality, and the
terminal convention. It does not imply those lipids lack biological regulation.

`endpoint_minimal_route_edges` is the union of shortest routes from any TF to
each lipid and from each TF to any lipid; it is not an enumeration of all
shortest paths for every pair. Cycle edges are edges whose endpoints lie in
the same nontrivial strongly connected component.

## Reproduce and inspect

From the repository root:

```sh
Rscript scripts/moon/audit_factor1_depth_topology.R
```

Dependencies are the existing `cosmosR` and `igraph` R packages. Output files:

- `depth_topology_audit.csv`: counts, coverage, shortest-length summaries,
  cycle metrics, and incremental changes.
- `depth_shortest_path_distribution.csv`: complete pairwise length histograms;
  `Inf` denotes unreachable pairs among all 49 × 155 candidate combinations.
- `depth_tf_lipid_shortest_paths.csv`: individual pair distances at each depth.
- `depth_input_coverage.csv`: retained inputs and nearest opposite-input distance.
- `depth_topology_source_manifest.csv`, `depth_topology_session_info.txt`,
  `depth_topology_run.log`: provenance, software, and execution record.
