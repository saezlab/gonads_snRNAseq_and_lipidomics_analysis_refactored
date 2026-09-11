# MOON compatibility audit for the mouse lipid PKN

Audit date: 2026-09-08. Scope: source inspection and synthetic numerical/topology tests only; no biological MOON scores were computed for this audit.

Implemented compatibility layer: `scripts/moon/moon_safe_helpers.R`. Its direct
`Rscript` execution runs synthetic contract tests; all passed immediately after
creation. `moon_safe_compress()` expands protection to every parent in a signed
child-signature group containing a protected node, then invokes the installed
compressor. It exports a complete original-to-compressed map and verifies that
mapping reproduces the compressed graph. `moon_safe_decompress_scores()` uses
that complete map, avoiding the package's no-merge decompression failure.
`moon_extract_paths()` applies the requested two thresholds only to inferred
scores, exempts observed nonzero values from the secondary threshold, never
promotes raw observed weights to primary nodes, makes observations sinks and
unconditionally enforces signed upstream-to-observation reachability. Tests
include a real unprotected merge alongside a protected duplicate group,
permuted edge order, zero compression, low-weight observations, high-weight
observations, incorrect upstream signs, disconnected coherent components and
dead-end cycles. The main runner and biological outputs are outside this
sub-audit's ownership.

The installed `cosmosR` is version **1.21.2**, GitHub commit **b2ed195738479c4ad204cace37fd11ae2fd02bb8**. Installed `decoupleR` is **2.17.0**. Function bodies of `moon`, compression/decompression, coherence filtering, both reducers, expression filtering, and neighborhood helpers were compared programmatically with `/tmp/fcgl-pkn-review/saezlab-cosmosR-b2ed195/R/` and were identical.

References inspected:

- [Current COSMOS_basic example](https://github.com/saezlab/COSMOS_basic/blob/main/scripts/net_compr_MOON.R), fetched from its raw GitHub URL on the audit date.
- [MOON and reduction source at the installed commit](https://github.com/saezlab/cosmosR/blob/b2ed195738479c4ad204cace37fd11ae2fd02bb8/R/decoupleRnival.R).
- [Compression source](https://github.com/saezlab/cosmosR/blob/b2ed195738479c4ad204cace37fd11ae2fd02bb8/R/compress_same_children.R).
- [Expression-filter source](https://github.com/saezlab/cosmosR/blob/b2ed195738479c4ad204cace37fd11ae2fd02bb8/R/filter_pkn_expressed_genes.R).
- The checkout's `AGENTS.md`, data-to-PKN mapping guide and latent-factor guide, and installed COSMOS MOON skill references.

## What transfers from the tutorial

Its overall order is appropriate: prepare TF activities and metabolite observations, clean/prune the PKN, optionally compress, score backward from measured metabolites, remove RNA-incoherent TF-target edges, rescore to convergence, decompress, and extract a consistent subnetwork. The numbers and identifier transformations cannot be copied uncritically.

The example uses `n_steps = 6`, ULM, upstream TF `abs(score) > 2`, and double-threshold extraction at primary 1.7 and secondary 1. It also adds cytosolic and mitochondrial copies of each metabolite. Its own iteration loop has a maximum of nine actual scores; when stopped by the cap, the final filtered graph can differ from the graph used for the last scores. Always rescore the final graph when reaching a cap.

## PKN contract and expression filtering

The requested RDS is a data frame with 254,589 rows and columns in **`source, interaction, target`** order. It has 33,953 `Enzyme__` nodes, 9,673 `Gene__` nodes, 9,116 `Metab__` nodes and 13 `UniProt__` nodes. These are graph identifiers, not the historical identifiers expected by several package conveniences.

- Pass only `source`, `target`, `interaction` to `meta_network_cleanup()`. Require finite numeric signs in `{-1,+1}` and nonempty endpoints. The cleaner averages all other columns for duplicated source-target pairs, so annotations must be kept outside the three-column working table. Conflicting signs are discarded after averaging.
- **Before compression explicitly reorder as `source,target,interaction`.** `compress_same_children()` reads/replaces the first two columns positionally. Passing the RDS's original order renames duplicate parent sources but leaves incoming target identifiers unchanged, producing broken connections. A synthetic `A→L, B→L, U→A` test reproduced this failure.
- Do not call `filter_pkn_expressed_genes_fast(names(raw_gene_symbol_RNA), pkn)` unchanged. It recognizes direct membership, names containing `Metab`/`orphanReac`, and historical `Gene123__...` patterns. It does not interpret `Gene__mouseSymbol` or `Enzyme__hash`.
- Prefixed RNA names rescue directly matched `Gene__` nodes, but **do not rescue any reaction/enzyme nodes**. Either use a documented type-aware context filter or omit this filter for the first compatibility pilot. Preserve reaction nodes while auditing the presence of their associated mouse catalysts separately. Unmapped `UniProt__` nodes need an explicit treatment.
- `prepare_metab_inputs()` only adds `Metab__` and compartment suffixes. It performs no chemical mapping. Build inputs directly from the measured-feature crosswalk and actual PKN node IDs; do not add prefixes twice or invent additional compartments. Compartment copies, if used, repeat a single whole-tissue measurement and are not independent observations.
- TF input names, TF-regulon sources/targets and RNA consistency-vector names must all use exactly the same `Gene__` domain. Preserve gene symbols in separate annotation columns.

## Compression safety: an additional counterexample

The skill recommends placing protected parents last, then checking that no protected node was renamed. **The check is essential; ordering alone is insufficient.**

Synthetic example:

```r
net <- data.frame(source = c("A", "B", "P"), target = "L", interaction = 1)
z <- cosmosR::compress_same_children(net, sig_input = c(P = 1), metab_input = c(L = 1))
z$duplicated_signatures
# A, B, AND protected P all map to "parent_of_L1".
```

Here two unprotected parents already duplicate the child signature. The implementation identifies that signature and then replaces every parent sharing it, including `P`, even though `P` occurs last. This is established by execution against the installed package.

Protect upstream and downstream observations plus every TF-regulon source/target needed by the coherence loop. Require all protected nodes that existed before compression to retain their IDs and require none in `names(duplicated_signatures)`. On failure, use the uncompressed PKN. A separately documented safe implementation could prohibit all groups containing a protected node, but modifying the package is unnecessary for a pilot.

Additional decompression failure: when no nodes were merged (`duplicated_signatures = NULL`), `decompress_moon_result()` silently drops all nonleaf original nodes. Synthetic `A→B→L` with scores for A, B and L decompressed to L alone. Therefore skip the helper entirely if compression merged no nodes. If it did merge nodes, verify the decompressed ID set explicitly and use `source_original` as the new `source`; the function leaves virtual IDs in its existing `source` column.

## Exact scoring semantics and small-sample behavior

`moon()` renames the sign column to `mor` and removes outgoing edges of every downstream input before any scoring. Thus observations become scoring sinks. The original PKN may contain outgoing lipid regulatory edges, but these outgoing edges do not explain other nodes during that run. Upstream lipid-to-protein regulation can still occur through unmeasured lipid nodes. An analysis of feedback from the *measured* lipids would require a separately defined direction or observation role.

The first layer uses all supplied downstream values as the decoupleR matrix. Each subsequent layer uses only the immediately previous inferred layer; a source is removed from further scoring after it is scored. `n_layers` is an upper bound, not a guarantee. A first scoring layer is attempted even if `n_layers` is zero. The `while` guard checks for more than one edge and more than one target overlap, not adequate statistical support or nonzero variance.

- `minsize = 1` is hardcoded at every layer. It cannot be changed by an argument to `moon()`.
- `upstream_input` does **not** anchor the numerical propagation, force TFs into results, or weight their inference by the magnitude of the TF activity. After the propagation it removes overlapping inferred nodes whose signs disagree with supplied upstream activities. Reachability pruning and final seed anchoring provide the explicit TF connection.
- Nonzero downstream values are appended as level 0 at the default `downstream_cutoff = 0`. That cutoff controls their inclusion in the returned table; it does not exclude them from initial scoring.
- ULM returns a t-like score, not a p-value. `moon()` discards decoupleR's p-values. `n_perm` has no effect for ULM. `wmean` hardcodes two permutations internally, regardless of `n_perm`; the reported weighted mean itself is the nonnormalized score. `norm_wmean` does use `n_perm`.
- `return_levels` is ignored: levels are always returned.

In installed decoupleR, ULM computes `r = cor(mor_matrix, measurement_matrix)` and `score = r * sqrt((N - 2) / ((1-r+1e-20)*(1+r+1e-20)))`. The coefficient vector is zero for observations not connected to that source. Consequently the result is relative to the whole layer's measurement background, not simply the signed sum of its directly connected observations.

Confirmed numerical examples:

| Synthetic input | Result |
|---|---|
| Seven positive observations `L1..L7 = 1..7`, `P1→L1` and `P7→L7`, both positive | `P1 = -1.732051`, `P7 = +1.732051` |
| Two observations 1 and 2, one positively connected parent for each | Both parent ULM scores are 0; p-value calculation warns of NaNs |
| One observation 2, positive parent→observation | Parent ULM score is NA |
| Three constant observations 1,1,1 | Parent ULM scores are NA; zero-standard-deviation warning |
| One observation 2 with `statistic="wmean"` | Parent weighted-mean score is 2 |

Require finite scores at **every realized layer** before interpreting or extracting. NA/Inf output is not a successful ULM result. A later layer receiving NA can fail the decoupleR input checks. Report sparse-layer failure explicitly; a documented `wmean` sensitivity can test basic signed propagation but has a different statistic and should not be labeled the tutorial ULM result.

The score baseline changes when the observation set changes, so comparing selected Factor 1 lipids against all mapped Factor 1 weights is substantive sensitivity analysis, not merely a computational convenience. The prior build's “seven Factor 1 lipids” referred to the **64 lipid subset in the saved correlation network**. It does not establish that only seven of all lipid-view Factor 1 weights can enter the new PKN. Resolve the full model-weight coverage before committing to that limitation.

## TF-target coherence

`filter_incohrent_TF_target()` expects scores in `source,score,level` order, regulons with `source,target,mor`, and named RNA values in the same identifier domain. It renames its second merged score column positionally; keep the canonical score-table column order rather than inserting annotations before score.

It removes an edge when `sign(inferred_TF_score * RNA_target * regulon_mor) < 0`. Missing RNA is excluded by an inner join. Zero RNA does not cause removal. It does not inspect the PKN edge sign itself, perform a pre-score pass, translate IDs, use protein evidence, or remove weak RNA support automatically.

Use the RNA **feature weights** as a separately labeled consistency vector. Do not substitute TF mRNA weights for TF activities. If RNA weights are thresholded, declare that threshold and whether below-threshold targets are omitted/zeroed; those operations disable the current sign check rather than removing edges.

This is a consistency screen using the same RNA factor signature that generated TF footprints; it is not independent validation. Retain unfiltered and filtered scoring runs, removed edge IDs and per-iteration logs if this comparison is practical. After a bounded loop, score the exact last graph again if no fixed point was reached. If decompression is used, project the final coherence exclusions onto the uncompressed PKN before extraction; the tutorial's extraction against the initial uncompressed graph can otherwise reintroduce rejected edges.

## Extraction and raw-weight thresholds

Level 0 retains **raw lipid factor weights**, whereas inferred layers have ULM scores. A single numeric threshold applied to both domains is not interpretable without a stated scaling choice. In particular, copying the example's secondary threshold of 1 can remove every lipid observation when feature weights are small.

Preferred pilot choice: preserve the original weights; apply candidate ULM cutoffs to inferred nodes, retain measured level-0 nodes explicitly, then prune the consistent graph to TF-to-observation reachability. Alternatively use the package reducer at cutoff 0 as a structural baseline and implement the separate observed/inferred selection step transparently. Do not replace true weights by artificial large scores merely to pass a reducer without clearly marking that temporary extraction-only representation. Positive global rescaling leaves ULM inferred scores invariant, but changes which raw level-0 inputs pass common thresholds; centering can change the biological signs and is not an innocuous repair.

Neither reducer guarantees the claimed final graph structure in all cases. Execution reproduced both failures:

- Single-threshold reducer keeps `U→X→Y→X` alongside `U→L`, although the X/Y cycle has no route to level-0 L. Recursive leaf pruning cannot remove a cycle.
- Double-threshold reducer keeps fully disconnected coherent `X→Y` alongside `U→L`. It only performs its connectivity step after finding an inconsistent edge, so no connectivity filtering runs when all edges initially pass.

On the **final sign-consistent graph**, independently calculate nodes reachable forward from retained intended TFs and nodes reachable backward from retained measured observations. Keep their intersection and repeat the edge/node consistency checks. Every exported node must lie on a directed TF-to-observation walk and every retained edge must satisfy `sign(score[source]*score[target]) == interaction`. For a more focused narrative export, prefer explicitly enumerated bounded TF-to-lipid routes or a monotonic score-level subgraph; reachability alone can still include detours through cycles that are on valid walks.

Keep measured lipids as terminal observations in the explanatory extraction if following the scorer's sink convention. Otherwise outgoing edges restored from the original PKN can create apparent mechanisms that were not part of MOON's scoring graph.

`get_moon_scoring_network()` requires one matching score row for the requested node and internal IDs matching the graph actually scored. It explains eligible downstream topology by score levels, not a unique causal path. `translate_res()` is HMDB-specific, positional, and converts `Gene...` labels to `Enzyme...`; avoid it for the new ChEBI/MeasuredLipid/Gene domains. Use the build's node annotation table to add display labels while retaining immutable internal IDs.

## Suggested bounded pilot and sensitivities

These are analysis choices, not package defaults or significance claims:

1. Verify sign-oriented Factor 1 weights and derive TF activities from the matching hepatocyte RNA view. Audit all lipid-view weights against exact and measured-rule PKN nodes; report the 64-feature subset separately. No superclass match should be silently treated as an exact molecular species.
2. Start with separate forward and reverse pruning radii of 6 and `n_layers = 6`, following the example. Compare radius/depth 4 and 8 if the coverage or inferred TF set changes materially. Reaction representation introduces explicit enzyme nodes, so graph distance is not a count of elementary biochemical reactions.
3. Start with footprint TF `abs(score)>2` as an exploratory threshold, retaining its target counts and any footprint p/FDR fields separately. Compare at least one lower/higher threshold if the inferred network relies on marginal TFs. TF score magnitude is not a calibrated causal effect.
4. Attempt compression only with complete protection checks and matching-column order; use an uncompressed run if checks fail. Preserve mapping metadata and verify any decompression is lossless for scored nodes.
5. Run exact package ULM first and inspect finite scores, realized levels, observation coverage and inferred/signed-compatible upstream TF overlap. Record direct target counts and the number of distinct measured lipids supporting candidate routes, since an inferred score may depend on only one lipid or many copies of one measurement.
6. For extraction, use level-specific selection: retain all nonzero observed weights and compare inferred absolute-score cutoffs 0, 0.5, 1 and 1.7 (or a prespecified subset). Require final TF-to-observation reachability. Report edge/node overlap and whether named hypotheses recur across thresholds.
7. If sparsity makes ULM unstable or undefined, report this result and add a `wmean` sensitivity only as a different scoring model. If only a handful of independent observed lipids remain, leave-one-lipid-out checks are more informative than increasing permutation counts for ULM.
8. Preserve scores, actual scoring PKN, uncompressed coherent PKN, inputs, dropped IDs, all thresholds, logs, provenance and display annotations. Separate class-template/species-rule edges and direct curated regulation in interpretation.

Even a successful result is a signed prior-knowledge hypothesis aligned with a latent covariance axis. It does not establish that TFs cause the lipid changes, that every reaction's substrates/cofactors are available, or that the graph implies feasible metabolic flux. These limits are especially relevant for rule-instantiated lipid species and compartment-expanded measurements.
