# COSMOS/MOON reproducibility and Methods fact audit

> **Current reproduction, verified 2026-09-15:** `Rscript scripts/moon/run_factor1_tf_to_lipid.R` refits the selected compressed model from the frozen mapped inputs and PKN using `scripts/moon/config/factor1_tf_to_lipid.json`. Its completed output in `results/causal_modelling/factor1_tf_to_lipid/` matches frozen C03 edge/node identities, compression mapping, score IDs/levels and numerical scores within 1e-10, and reproduces expected-TF exclusions, evidence exports and plots. See `reference_comparison.json` and the run/dependency/output manifests. Use `--validate-only` for setup checks without output and a new `--output` directory for replay. It does not rebuild the PKN or mapped inputs. The audit below records the earlier frozen-state review; its historical source hashes and path literals are preserved. Original scripts now reside in `scripts/moon/archive/2026-09/`, while current commands are documented in `scripts/moon/README.md`.

Audit scope: read the relocated saved artifacts, relevant analysis source, package metadata and archived manifests. No MOFA model or MOON state was refitted, and no historical analysis script or result was changed in this audit. Paths below are relative to the repository root unless specified otherwise. The current artifact root is `results/causal_modelling/`.

## Selected analysis identity

The selected baseline is **`cosmos_factor1_compressed_review/threshold_search/C03`**, not `sensitivity_threshold_search/C03`. Candidate identifiers are local to each search directory.

| Item | Verified value | Relocated evidence |
|---|---|---|
| Underlying scoring state | `radius6_depth4_compressed` | `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth4_compressed/run_state.rds` |
| Topology pruning radius | 6 steps in each directed pruning operation | Same state and `run_manifest.json` |
| Maximum MOON layers | 4; realized maximum 4 | Same state and `coherence_ledger.csv` |
| Compression | Same signed children; safe input/regulon protection enabled | `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/topology_radius_6/compressed/compression.rds` |
| Primary / secondary cutoffs | Strict absolute inferred-score cutoffs `>1.85` / `>0.55` | `threshold_search/C03/extraction_audit.json` |
| Selected network | 153 original-ID nodes, 388 signed directed edges | `threshold_search/C03/network_ATT.csv`, `network_SIF.csv`, `solution.rds` |
| Upstream TFs | Ahr, Hnf4a, Ppara, Tbx2 | `threshold_search/selected_endpoints_and_margins.csv` |
| Downstream observations | DG(20:0_18:1), DG(24:0_18:2), DG(24:1_18:2) | Same endpoint table |

The full local path prefix for rows abbreviated with `threshold_search/` is `miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/`.

## Factor 1 input construction

1. The source was the saved named column `Factor1` in `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/main_model_RNA_lipid_weight.RData`, verified against native factor 1 in `results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/RNA/mofa_res_12factor.hdf5`. The model was initialized with 12 factors and retained 10; it contains 43 samples. Extraction used names, not factor column position. The saved hepatocyte and lipid loadings matched native HDF5 loadings with maximum absolute differences of 4.885e-15 and 5.107e-15, respectively.
2. The input set contains 20,123 hepatocyte RNA feature loadings, 1,254 lipid feature loadings and 367 inferred hepatocyte TF footprint activities. TF activities are inferred from target-gene loading patterns, not the TF's own RNA loading.
3. TF footprints were reproduced with `decoupleR::run_ulm(..., minsize = 20)` from the hepatocyte RNA loading vector. The regulon has 37,262 signed rows and combines the repository's `support/TF_regulons.RData` CollecTRI network with `support/Gonadal_Hormone_Genes.csv`, which replaces the original Ar and Esr1 regulons. The hormone table's `_spec` suffixes and uppercase AR/ESR1 source names were normalized to Ar/Esr1. Recomputed scores reproduce the saved activities to maximum absolute error approximately 2.26e-12. Raw footprint p-values and BH-adjusted p-values are stored, but **TF input selection used `abs(activity) > 2`, not a p-value cutoff**.
4. Factor 1 was jointly oriented toward testes-positive sample coordinates. The original coordinate already met this convention: means were -0.7422907155 for ovaries and +0.7085502284 for testes, so the multiplier was +1. The same convention applies to RNA loadings, lipid loadings and TF activities. No differential-expression t-statistics or direct testes-minus-ovaries effects replaced these values. Factor 1 explains 8.906% of hepatocyte RNA variance and 29.958% of lipid variance; its saved association with gonadal type is R²=0.969836. These three quantities must not be conflated.
5. Lipids were mapped with exact reported identities. Explicit measured-species rule nodes were preferred; otherwise a single unambiguous exact-resolution ChEBI identity was accepted. This gave **155 distinct measured features/nodes: 149 measured-species rule nodes and six ChEBI nodes**. The other 1,099 observations were excluded: 1,098 lacked an accepted exact node and one had multiple exact ChEBI alternatives. No accepted feature required compartment copying, and no two accepted features mapped to the same node. Accepted compartment labels are `unknown`.
6. No lipid-class proxy, compatible but more specific structure, protonation substitution, feature averaging, z-scoring or rescaling of the mapped lipid weights was used. The previous correlation network's 64 selected lipid features were annotated only and did not restrict the new input set. All original weights are preserved.
7. Of 367 TF activities, 353 have PKN nodes. Fifty-one activities pass `abs(activity)>2`; 49 also have PKN nodes. These 49 form the initial upstream vector and all remain after baseline pruning. Of 155 mapped lipid nodes, 74 survive baseline directed topology pruning. No RNA expression-presence filter was imposed, because a nonzero factor loading does not establish expression detection.

Evidence: `miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/{factor1_inputs.rds,factor1_inputs_raw.rds,README.md,factor1_model_audit.json,factor1_hepatocyte_tf_activities.csv,factor1_hepatocyte_tf_activities_mapped.csv,factor1_lipid_node_mapping.csv,hepatocyte_tf_inference_regulons.csv,input_source_manifest.csv}`. Source implementations are `scripts/moon/archive/2026-09/prepare_factor1_inputs.R`, `scripts/moon/archive/2026-09/audit_factor1_model.py` and `scripts/moon/archive/2026-09/map_factor1_inputs.R`.

Factor stability evidence is descriptive agreement across already fitted factor-number choices, not bootstrap, donor-resampling or repeated-seed validation. The TF activity inference and RNA coherence checks use the same RNA loadings; the latter is internal consistency rather than independent experimental validation.

## Baseline PKN processing and MOON scoring

The selected PKN is `results/causal_modelling/cosmos_mouse_pkn/pkn_lipid_rules.rds`; the companion PKN audit supplies its resource provenance. The scoring runner retained the signed `source`, `interaction`, `target` representation, cleaned it, verified signs of -1 or +1, and enforced unique signed source-target relationships and no self-loops. Outgoing edges of mapped measured lipid nodes were removed so observations were terminal.

Directed forward reachability from the upstream TF vector and reverse reachability from the downstream lipid vector were each limited to radius 6. Induced edges were retained, inputs were restricted to present nodes, and this operation was repeated until edge count stabilized. The pruning ledger is:

| Iteration | Edges | Nodes | Upstream TFs | Downstream lipid nodes |
|---|---:|---:|---:|---:|
| 1 | 44,392 | 6,791 | 49 | 74 |
| 2 | 44,127 | 6,673 | 49 | 74 |
| 3 | 44,127 | 6,673 | 49 | 74 |

The baseline compressor canonicalizes each parent's signed child set. It protects measured input identifiers and all TF-regulon source/target identifiers present in the graph; if a protected parent belongs to a same-child signature group, the entire group is protected. The safe wrapper calls `cosmosR::compress_same_children()` for the eligible groups and verifies the complete original-to-compressed graph projection and unchanged protected identifiers. Initially 2,720 original nodes were merged into 378 virtual groups, reducing 6,673 original nodes/44,127 edges to 4,331 scoring nodes/36,569 edges. There are 2,591 directly protected identifiers, expanded to 2,777 by the signature-group rule. This **is actual network compression**.

`cosmosR::moon()` used `statistic="ulm"` and `n_layers=4`, with the unscaled 74 lipid weights as downstream measurements and the 49 TF activities supplied upstream. The ULM score is calculated against the preceding layer's complete node background, including zero coefficients for non-targets. A node is scored at its first eligible layer and is not rescored when more distant layers are reached. The upstream vector applies a final sign-agreement filter to overlapping scored TF rows; it is not a hard upstream-seeded optimization constraint. Default downstream cutoff 0 preserves the nonzero observations at level 0. ULM does not use the `n_perm` argument.

After each MOON call, the package helper `filter_incohrent_TF_target()` compared inferred TF signs, original signed regulons and hepatocyte RNA loadings. It removes a matching TF-target pair if `sign(TF_score × RNA_target_loading × regulon_sign) < 0`. The script allowed at most ten iterations and declared convergence when the exact ordered signed edge set no longer changed. Baseline C03's scoring state converged after three iterations:

| Iteration | Compressed edges scored | Edges removed | Returned compressed score rows | Returned upstream TFs |
|---|---:|---:|---:|---:|
| 1 | 36,569 | 7,628 | 1,889 | 19 |
| 2 | 28,941 | 1,195 | 1,711 | 13 |
| 3 | 27,746 | 0 | 1,680 | 13 |

The final compressed graph has 27,746 edges and projects to 35,304 original-ID edges. Its 1,680 returned compressed score rows decompress to 2,207 original-ID score rows. The baseline does **not** include the separate initial RNA-coherence or repeated inputs-only recompression sensitivities. Its original compression mapping remains fixed during the coherence loop.

The historical implementation source for the baseline sequence is `scripts/moon/archive/2026-09/search_factor1_parameter_space.R`; the current selected refit uses `scripts/moon/run_factor1_tf_to_lipid.R` and `scripts/moon/lib/workflow.R`, retaining the same compression/decompression and extraction implementation in `scripts/moon/moon_safe_helpers.R`. Saved evidence is in `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/topology_radius_6/` and `.../runs/radius6_depth4_compressed/`.

## Threshold search, final extraction and export

The later **compressed-only** search reused the converged compressed scoring states at depths 4, 5, 6, 7, 8, 10 and 12. Depths 8, 10 and 12 were exact aliases of depth 7 for scores, PKN and inputs, yielding four distinct scored states. The extraction-only coarse grid used primary cutoffs 0.5–3.0 by 0.1 and secondary cutoffs 0–2.0 by 0.1, restricting secondary≤primary. Local refinement used 0.025 increments within ±0.1 of threshold/coverage Pareto candidates and near-cap 500–999-edge anchors. There were 2,376 distinct scored-state/threshold combinations, including empty results; 253 had 1–499 edges and 249 additionally retained at least two upstream TFs and two downstream observations.

C03 is a selected interpretable representative, not a calibrated statistical optimum. Its 388 signed edges are identical at primary 1.85 or 1.875 combined with secondary 0.55 or 0.575. Secondary 0.55 leaves Hnf4a 0.0516453 above the cutoff, whereas secondary 0.60 leaves it only 0.0016453 above the cutoff. These margins describe a fixed-score extraction plateau, not score uncertainty.

The final extractor is the repository's **guarded `moon_extract_paths()` dual-threshold implementation**, not an unmodified call to the package's dual-threshold reducer. It requires:

- inferred nodes above the strict secondary absolute-score cutoff;
- at least one endpoint of every candidate edge above the strict primary absolute-score cutoff;
- upstream input nodes with inferred signs agreeing with supplied footprint signs;
- edge-wise sign agreement, `sign(score_source × score_target) == interaction`;
- every retained node reachable from an intended upstream TF and able to reach a measured level-0 observation in the selected directed graph;
- level-0 weights exempt from inferred-node score cutoffs, ineligible as primary nodes, and terminal.

For C03, threshold selection generated 1,630 candidate edges; 1,029 passed score-sign checks; directed endpoint-connectivity pruning retained 388 edges/153 nodes. No manual edge trimming was used to satisfy the strict `<500` target. Materialized candidates were re-extracted with the uncached helper and checked against their frozen full-precision scoring states. Archived `selected_validation.csv` reports true for exact size/grid agreement, final-PKN edge membership, exact state score identity, strict secondary/primary conditions, score-sign agreement and seed-to-observation connectivity.

Provenance export preserves internal identifiers and joins exact source/sign/target triples to local PKN evidence. C03 has 393 evidence records for its 388 edges, three reaction records and six edges with measured-species rule hypotheses. `network.graphml` includes distinct `factor1_weight`, `moon_score` and upstream activity attributes; raw observed loadings and inferred MOON scores must not be described as a single comparable measurement scale.

For route illustrations the exporter selects one shortest path per reachable TF–observation pair, breaks shortest-path ties by maximum minimum absolute inferred score and then deterministic identifiers, and excludes level-0 weights from the score minimum. C03 has 12 reachable TF–lipid pairs. Its chosen paths contain **7–9 edges**, all longer than six; neither a pruning radius of six nor four MOON layers is a hard route-length cap. The separate catalyst-only audit retains all 12 pairs; the ≤6-edge and catalyst-only ≤6-edge route audits retain none. Thus the manuscript should not describe this network as connecting TFs to lipids within four or six reactions/edges.

Relevant files:

| Function | Script | Main relocated evidence |
|---|---|---|
| Compressed threshold search and SIF/ATT/GraphML | `scripts/moon/archive/2026-09/search_compressed_threshold_review.R` | `.../cosmos_factor1_compressed_review/threshold_search/{search_manifest.json,grid_complete.csv,selected_candidates.csv,state_validation.csv,selected_validation.csv}` |
| Selected network and extraction audit | Same script | `.../threshold_search/C03/{solution.rds,network_SIF.csv,network_ATT.csv,network.graphml,extraction_audit.json}` |
| Local cutoff robustness and matched-policy comparison | `scripts/moon/archive/2026-09/audit_compressed_candidate_neighborhood.R` | `.../threshold_search/C03_local_threshold_neighborhood.csv`; `.../sensitivity_threshold_search/matched_threshold_sensitivity.csv` |
| Provenance and route audit | `scripts/moon/archive/2026-09/export_factor1_network.py` | `.../threshold_search/C03/network_evidence_network_summary.json`; `candidate_paths_network*.csv`; evidence CSV/JSONL files |
| Route illustrations | `scripts/moon/archive/2026-09/plot_factor1_network.R` called by exporter | `.../threshold_search/C03/figures/candidate_paths_network_network.pdf` |
| Complete scores and sign-filter explanations | `scripts/moon/archive/2026-09/export_compressed_score_retention.R` | `.../cosmos_factor1_compressed_review/score_retention/radius6_depth4_compressed_scores_all.rds` and `.csv` |

In the table, `.../` means `results/causal_modelling/`.

## Focus TF diagnostic facts relevant to interpretation

The complete baseline depth-4 score archive has 2,238 original-ID rows before upstream sign filtering: 2,207 are returned and 31 are omitted for sign conflict. Esr1, Nfic and Dbp survive input selection, directed pruning and compression, and are present in the final scoring graph. Their positive raw scores (+0.1317963, +0.7185897 and +0.3451557) conflict with their negative supplied footprints (-7.0887486, -4.3492313 and -7.6608848), so they are omitted **after scoring**. These TFs must not be described as pre-score pruning losses or recovered mechanisms. Full diagnostics are under `miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/tf_audit/`.

The published-style baseline can mention sensitivity checks without substituting them into its Methods. At matching depth 4 and cutoffs 1.85/0.55, broader compression with inputs-only protection and repeated recompression yields the exact same 388 edges; initial RNA coherence yields 368 edges with the same four TFs/three DGs and signed-edge Jaccard 0.913924. These are methodological checks rather than independent biological validation.

## Software and provenance status

The archived baseline session records R 4.6.1, cosmosR 1.21.2, decoupleR 2.17.0 and igraph 2.3.3 on aarch64 macOS Tahoe 26.6.2. Current installed cosmosR package metadata gives commit `b2ed195738479c4ad204cace37fd11ae2fd02bb8`, also recorded in the previous analysis report. The pinned COSMOS_basic template is commit `7108c6b849c712b1282440dc122c9d783f563da8`, `scripts/net_compr_MOON.R`. Python input audit records Python 3.12.13, numpy 2.5.2 and h5py 3.16.0. These are source/audit environments, not claims that the entire repository has a locked executable environment.

Read-only MD5 verification after relocation:

| Artifact | Current MD5 | Comparison |
|---|---|---|
| Relocated `pkn_lipid_rules.rds` | `be8da3886f19f8b121a389e442113edd` | Matches archived baseline manifest |
| Relocated `factor1_inputs.rds` | `679f0f336499ec59106eb2107bb49ec4` | Matches archived baseline manifest |
| Relocated depth-4 compressed `run_state.rds` | `3f035b52f990919fc51c0c450edd7a2f` | Matches archived compressed-search manifest |
| `scripts/moon/moon_safe_helpers.R` | `86be2ca41edc65841f7674209c8b36fd` | Matches both archived manifests |
| Current `search_factor1_parameter_space.R` | `93b4571b8a92414497aeee05730e98ea` | Differs from baseline manifest's `06126fb2a0ea772ebf86d828fdcba3f5` |
| Current `search_compressed_threshold_review.R` | `e906362cd4c39f1be09039cb1749be81` | Differs from search manifest's `508aa9adaa79e473988f4d88ca458b9d` |

All five original input-source files named in `inputs/input_source_manifest.csv` still exist at their original non-COSMOS paths and match their archived hashes: the weight RData, TF regulon RData, hormone-target CSV, MOFA HDF5 and original correlation RData. Hash equality establishes unchanged bytes, not correctness of the underlying scientific assumptions.

## Historical relocation and source-replay caveats

The relocated outputs are intact and readable, but historical commands are **not an immediately runnable end-to-end replay** in the new directory layout.

- Core preparation, mapping, scoring, threshold-search and audit scripts still hardcode old `results/cosmos_*` locations. The old directories no longer exist. Executing them unchanged can fail on reads or create a new old-style output tree. Existing report/README links also point to old locations.
- Frozen `run_state.rds` retains a relative `run_dir` such as `miscellaneous/older_results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth4_compressed` and an absolute `topology_reference` containing the old location; both resolve to nonexistent paths. CSV selection tables and manifest hash keys retain the same historical path strings. These should be treated as archived provenance and resolved to the relocated files by an adapter, not silently assumed current.
- The final exploratory search script hashes differ from the source hashes saved with the original run. The numerical PKN/input/state objects and the safe helper still match, but claiming byte-identical source replay requires the historical script snapshots or a reviewed change record in addition to path migration.
- `scripts/moon/archive/2026-09/export_factor1_network.py` offers explicit `--results`, `--pkn` and `--input-mapping` arguments. All three must be supplied for the relocated baseline. Its default mapping fallback points to the old pilot folder, and currently calls `.exists()` on `None` if no mapping candidate exists. Supplying the correct mapping avoids that failure. The exporter writes derived evidence/GraphML/plots; it was not rerun as part of this read-only artifact audit.

The following export command was **documented but not executed in the earlier audit**; its script location is updated to the archive. It is a historical record, not a current reproduction instruction. Current exports are produced by the main workflow in its new output directory:

```sh
python3 scripts/moon/archive/2026-09/export_factor1_network.py \
  --results miscellaneous/older_results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03 \
  --pkn results/causal_modelling/cosmos_mouse_pkn \
  --input-mapping miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_lipid_node_mapping.csv \
  --stem network
```

The immutable relocated state can be audited without refitting by loading its RDS path explicitly, resolving the saved scoring PKN and compression map from the current artifact root, applying `moon_extract_paths()` to the frozen scores at 1.85/0.55 and comparing the result to saved C03. The subsequently completed frozen-state audit is now archived as `scripts/moon/archive/2026-09/audit_manuscript_causal_model.R`, with successful outputs under `results/causal_modelling/manuscript_audit/`. Its historical self-path is preserved, so it is not an active execution entry point; use the verified current workflow described above.

Current reproducibility statement: “Archived inputs, compression maps, converged scoring states, selected networks, provenance and validation reports are available. The current selected-model workflow refits compressed MOON from frozen mapped inputs and the assembled PKN and reproduces the archived result, including score and graph comparisons. Historical preparation, parameter-search and sensitivity scripts remain preserved separately; this refit does not rebuild the input bundle or PKN, or rerun every original exploration.”
