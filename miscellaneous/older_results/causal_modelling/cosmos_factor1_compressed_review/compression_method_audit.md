# Compression method audit

Audit date: 2026-09-09. Installed cosmosR: **1.21.2**, commit **b2ed195738479c4ad204cace37fd11ae2fd02bb8**. This audit reads the saved radius-6 analysis and installed function bodies. It does not modify scoring helpers or rerun the MOON parameter search.

## What “uncompressed” actually meant

The previous uncompressed candidates were **scored without same-children compression**. In `scripts/moon/search_factor1_parameter_space.R`, the uncompressed branch calls `moon_safe_compress(prepared$net, node_ids(prepared$net))`. The wrapper protects every original node, finds no eligible duplicate group, and therefore **never calls `cosmosR::compress_same_children()`**. Its returned scoring graph is exactly the original pruned graph; this identity was tested on the saved RDS objects.

This differs from a compressed run whose scores are subsequently decompressed into original biological IDs. Such an export can legitimately be named `pkn_final_uncompressed.rds`, or have individual original nodes in GraphML, while its **scores were calculated on merged groups**. The run manifest and `pkn_final_scoring.rds`, rather than the export name, establish which graph was scored.

## Why the concern is justified

The authors’ [COSMOS_basic guidance](https://github.com/saezlab/COSMOS_basic#network-compression) explicitly identifies compression as protection against repeatedly propagating the same downstream information: identical-child parents acquire redundant inferred scores and can act as repeated observations in the next upstream fit. The [example pipeline](https://github.com/saezlab/COSMOS_basic/blob/main/scripts/net_compr_MOON.R) therefore compresses before scoring and decompresses afterward. It presents an uncompressed run as a sensitivity analysis.

[decoupleR’s ULM documentation](https://saezlab.github.io/decoupleR/reference/run_ulm.html) specifies that activity is the t-statistic from regressing molecular readouts on signed regulator weights, including zero weights for unrelated features. Thus compression changes the effective observations and background for upstream fits. **A higher uncompressed score is not automatically stronger evidence.** It is also incorrect to call every higher uncompressed value an artifact: changing the regression background, propagation levels, and subsequent RNA-coherence deletions can increase or decrease particular scores. This latter distinction follows from the implemented fits and filtering, not from an independent biological validation.

## Actual compression and protection counts

The frozen radius-6 topology has 49 TF inputs and 74 lipid inputs. Counts below are before the first coherence iteration.

| Representation | Scoring edges | Scoring nodes | Original nodes merged | Merged groups | Explicitly protected IDs | Protection after whole-group expansion |
|---|---:|---:|---:|---:|---:|---:|
| Previous uncompressed mode | 44,127 | 6,673 | 0 | 0 | 6,673 | 6,673 |
| Previous compressed mode | 36,569 | 4,331 | 2,720 | 378 | 2,591 | 2,777 |
| Inputs-only compression diagnostic | 33,191 | 3,882 | 3,393 | 602 | 123 | 123 |

The diagnostic last row applies the same safe one-pass compressor but protects only actual upstream/downstream input IDs. **It is not a scored candidate or a recommendation to bypass RNA coherence.** Compared with the previous compressed mode, it merges 673 additional original nodes; these include 20 regulon sources and 484 regulon targets (the categories can overlap).

The previous compressed mode protects **every regulon source and target present in the PKN**, as well as input nodes. This preserves identifiers expected by the unchanged RNA coherence helper, which joins TF score IDs, regulon IDs, RNA target IDs and PKN endpoint IDs. It is a conservative compatibility choice, but it also preserves substantial repeated signal. Restricting protection without changing how coherence is applied would silently miss regulon checks after original nodes are renamed.

## Safety checks and remaining redundancy

The installed [compression source](https://github.com/saezlab/cosmosR/blob/b2ed195738479c4ad204cace37fd11ae2fd02bb8/R/compress_same_children.R) uses positional columns and row-ordered child signatures. Its exclusion of input names while identifying duplicates does not exclude protected names from the subsequent whole-group replacement. The local wrapper correctly canonicalizes `source,target,interaction`, sorts signed children, excludes an entire signature group if it contains a protected parent, and verifies the complete original-to-group projection. The saved graph passes these invariants. The wrapper’s synthetic contract tests also pass. **No new bug was found in its stated one-pass, protected-ID behavior.**

One-pass compression does not remove all possible redundancy. Merging child nodes can make previously different parents equivalent; later RNA-edge deletions can also change equivalence classes. The current compressed initial graph still has **257 exact signed-child duplicate groups involving 857 parents**: 243 groups contain a protected ID and 14 do not. At the final depth-5 coherent scoring graph, there are **260 groups involving 871 parents**, with 246 groups touching protected IDs. The inputs-only one-pass diagnostic still has 60 newly formed duplicate groups involving 245 parents. These counts describe structural equivalence, not independent estimates of score inflation.

## Defensible compressed alternative

Prioritize solutions scored with compression and assess lower extraction thresholds there; absolute cutoffs are not directly comparable across representations. The previously emphasized uncompressed networks should be treated as representation sensitivities.

A justified stronger-compression sensitivity can protect only measured input identities while preserving coherence as follows: compress the current original-ID PKN; score it; decompress inferred scores; apply the existing RNA-coherence rule to original-ID edges; then recompress after any deletions. This lets merged regulon nodes retain their original RNA evidence without wholesale protection. Record convergence using original edge deletions, save each mapping and scored graph, and explicitly record any input that becomes disconnected. Use only original edges whose mapped signed counterparts survive the final scoring graph when extracting a network; compression cleanup can eliminate contradictory incoming edges after endpoints merge.

This changes the compression/coherence implementation deliberately and therefore belongs in a separately labelled sensitivity run. Retaining one-pass compression per iteration preserves the package’s compression policy; recursive compression to a fixed point would be an additional methodological change. Neither approach validates a mechanism by itself.

## Reproduction and outputs

Run from the project root:

```sh
Rscript results/cosmos_factor1_compressed_review/compression_method_counts.R
Rscript scripts/moon/moon_safe_helpers.R
```

The first command passed and generated `compression_method_counts.csv`, `compression_final_graph_redundancy.csv`, and `additional_nodes_merged_with_inputs_only_protection.csv` in this directory. The second passed the protected-group, canonical-order, complete-decompression, threshold, sign, sink, and connectivity contracts. No helper or scoring script was changed during this audit.
