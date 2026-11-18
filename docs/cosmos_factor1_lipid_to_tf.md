# Factor 1 lipid-to-TF mechanistic hypotheses

Analysis completed on 15 September 2026 using the existing compressed MOON
pipeline, frozen Factor 1 measurements and mouse prior-knowledge resources.

The analysis produces a **455-edge, 233-node exploratory network connecting
LPC(18:0/0:0) to 64 TF observations through a prostaglandin–Ptger3 branch**.
This requires the existing acid/base-normalized PKN sensitivity. The original
PKN, under the frozen accepted-input mapping, has no directly supported measured
signalling lipid and no lipid-only precursor route to the observed TFs.
Consequently, the positive result is a
conditional lipid-remodelling hypothesis, not evidence that LPC directly
activates a receptor or that its abundance causes the TF signature.

The selected network is in
[`refinement/selected_all_tfs/`](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/):
[GraphML](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/network.graphml),
[signed edges](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/network_SIF.csv),
[annotated nodes](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/network_evidence_nodes.csv),
[illustrative routes](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/network_candidate_paths.svg)
and [selection record](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/selection.json).
The original TF-to-lipid manuscript result is unchanged.

**Coverage correction, 15 September 2026:** the subsequent
[regulatory measurement audit](cosmos_lipid_regulatory_measurement_coverage.md)
identified specific fatty-acid candidates omitted by the abbreviation-based
mapper. FA(20:4) and FA(22:6) were measured and are compatible with AA and DHA,
respectively, but those specific isomers were not unambiguously assigned or used
as observations. All coverage statements below refer to the frozen mapping;
they are not exhaustive estimates of available biological knowledge. This
correction does not change the fitted scores, topology or selected network.

## Inputs and biological direction

The analysis uses the same jointly oriented, testes-positive Factor 1 as the
forward analysis. Upstream measurements are the 155 accepted species-level
lipid loadings: six ChEBI-mapped lipids and 149 explicitly represented
measured-species reaction hypotheses. Their magnitudes are retained without
rescaling; they are abundance associations rather than ligand potencies. The
current MOON implementation applies upstream sign compatibility to inferred
scores; it does not model lipid dose responses.

The main observation panel contains all 353 finite, nonzero PKN-present
hepatocyte TF footprint activities. A separate sensitivity retains the 49 with
absolute activity greater than two. TF observations are terminal: their outgoing
edges are removed before pruning, their activities are supplied at level zero,
and they cannot seed the primary inferred-score threshold. Hepatocyte RNA
loadings remain separate evidence for internal TF–target coherence.

These are view-specific feature weights and derived TF footprints, not sample
factor coordinates. Factor 1 accounts for 8.906% of hepatocyte RNA variation and
29.958% of lipid variation in the previously audited model. Earlier stability
checks were descriptive comparisons across fitted factor-number choices, not
resampling validation. The data are cross-sectional; no time ordering or
hepatocyte subcellular ligand exposure is measured.

## Why a different PKN policy was necessary

The original signed PKN has 254,589 edges, but only 77 chemical→protein
regulatory edges from 49 chemicals; the lipid subset comprises 38 edges from
23 lipid identities. **None of the 155 accepted lipid inputs has a direct
regulatory edge.** Among all 1,254 lipid features considered during mapping,
the frozen crosswalk also supplies no exact current curated regulatory match,
including explicit acid/base equivalents. Most of the PKN's scale therefore does not translate
into measured signalling-lipid coverage.

Unrestricted reaction traversal appears to connect many measured lipids to
Nr3c1. The shortest examples pass through H⁺ and an unrelated steroid reaction.
Removing common currency metabolites reduces this problem but still permits
longer cofactor shortcuts. These graphs were retained as diagnostics and were
not used to nominate signalling lipids. The precursor policy permits only
lipid→reaction→lipid steps, curated lipid→protein regulation and protein→protein
edges. Catalyst identities remain reaction annotations; no protein→reaction
traversal or reversed PKN edges are introduced.

| PKN policy | Outcome |
| --- | --- |
| Original PKN, direct chemical regulation only | No accepted lipid input with an outgoing regulatory edge |
| Original PKN, lipid-only precursor traversal | No measured-lipid→observed-TF path, including unbounded reachability |
| Existing acid/base-normalized curated PKN, lipid-only traversal | Two potentially connected measured lipids: LPC(16:0/0:0) and LPC(18:0/0:0) |
| Normalized PKN with measured-species rules restored | Same connected inputs, pruned signed graphs, input vectors and compression maps at every tested pruning/panel configuration |

The normalized sensitivity uses the frozen, explicit ChEBI alias/acid–base
crosswalk; it does not assign lipid-class proxies or infer new chemical
identities. The prebuilt curated variant contains only six of the accepted
input nodes. A separate rule-preserving audit restores all 155 nodes and adds
898 lipid-precursor edges, but the additional 149 inputs still have no path to a
protein. They disappear before scoring. Exact pruning and compression
equivalence was verified for all six panel/distance combinations and all 16
saved prepared states. Thus omission of the rule layer does not explain the
positive or negative results here.

Coverage and shortcut examples are in
[`coverage/`](../results/causal_modelling/factor1_lipid_to_tf/coverage/);
the rule-preserving comparison is in
[`normalized_rules_topology/`](../results/causal_modelling/factor1_lipid_to_tf/normalized_rules_topology/).
Class-level LPC–Trpv2 and more-specific DG–Trpc2 mapping opportunities remain
excluded from exact input assignment; those receptors also lack onward
protein→TF connectivity in the current graph. Cached STITCH/BRENDA alternatives
were audited separately and were not incorporated.

## Fitting, extraction and selection

The initial search crossed three graph policies, two TF observation panels,
pruning distances 6/8/12 and scoring layers 4/6/8/12: 72 configurations in total.
Sixteen had sufficient topology and observations to fit; all converged. The
remaining configurations lacked mapped upstream nodes, directed topology or
enough observations for ULM. The normalized all-TF graph at pruning distances
6/8 reached only Rxra, which is insufficient for ULM scoring.

The 16 fits used safe identical-signed-child compression with input/regulon
protection and a fixed compression map during RNA-coherence filtering. The
coarse grid evaluated 752 valid dual-threshold pairs; none met the compact
network criterion. A targeted refinement of the seven fits with nonempty
coarse networks evaluated a further 1,431 pairs at 0.025 spacing. All extractions
passed validation, and 26 met the criterion of fewer than 500 edges with at
least one upstream lipid and two observed TFs. All 26 came from the all-TF
panel; no strong-TF sensitivity met the size criterion. This was a targeted
refinement, not an exhaustive continuous search or statistical model selection.

| Selected setting or result | Value |
| --- | --- |
| Prior | Acid/base-normalized curated PKN; lipid-only precursor policy |
| Observation panel | All 353 PKN-present TF activities; 208 retained after pruning |
| Pruning distance / scoring layers | 12 / 12 |
| Initial pruned graph / compressed graph | 9,896 / 8,223 signed edges |
| RNA-coherence convergence | Three passes; 8,060 final compressed scoring edges |
| Strict primary / secondary cutoffs | Absolute inferred score > 1.65 / > 1.10 |
| Selected original-ID graph | 455 signed edges; 233 nodes |
| Node composition | 162 gene products, 51 reaction instances, 20 metabolites |
| Measurements connected | One lipid; 64 observed TF activities |
| Actual chemical→protein transitions | PGF₂α→Ptger3 and PGD₂→Ptger3 |

Selection favours a rounded cutoff with a compact, interpretable output. It is
not a biological optimum. At secondary 1.10, primary 1.625 gives 462 edges,
1.65 gives 455, 1.675 gives 425, and 1.70 eliminates the network. At primary
1.65, secondary 1.15 gives 442 edges with the same 64 TFs. The upstream lipid's
absolute inferred score is only 1.177648, creating a strict upper bound on its
secondary cutoff. Thresholds therefore cannot be increased arbitrarily while
retaining an explanation. The smallest nonempty strong-TF network in the
targeted refinement still has 1,974 edges.

Pruning distance and scoring layers are separate settings. The selected LPC
is scored at level 11. Shortest exported routes to individual selected TFs
span 11–17 graph edges; neither setting is a hard maximum on the final displayed
walk length. The exported paths are illustrative shortest routes, not every
possible path or an independently fitted minimal causal model.

## What the selected network actually proposes

Its central branch is:

```text
measured LPC(18:0/0:0)
  → acyl remodelling through PC(18:0/20:4)
  → arachidonate (not assigned a measurement) → PGH₂ → PGF₂α / PGD₂
  → Ptger3 → protein signalling → observed TF activities
```

The first reactions are RHEA:37480, annotated to Lpcat3, and RHEA:40520,
annotated to Pla2g10/Pla2g4a/Pla2g6. They consume an unmeasured
arachidonoyl-CoA donor to acylate LPC and subsequently release arachidonate while
regenerating LPC. Adding their full equations cancels LPC and PC:

```text
arachidonoyl-CoA + H₂O → arachidonate + CoA + H⁺
```

Thus the arachidonate chain does not originate in the stearoyl chain of measured
LPC(18:0). This is a potential acyl-remodelling/availability mechanism conditional
on donor pools and enzyme context. The graph does not establish that LPC is
rate-limiting, the direction of a steady-state abundance effect, or net flux.
Prostaglandin synthesis also requires omitted oxygen/reducing equivalents.
The complete reaction and regulatory audit is documented in
[the selected-network evidence note](cosmos_factor1_lipid_tf_selected_evidence.md).

The frozen Guide to Pharmacology annotations identify PGF₂α and PGD₂ as positive
Ptger3 agonists, with pKi values 7.12 and 6.55 respectively. These are affinity
annotations, not functional EC50 values. Their attached primary study measured
binding to mouse prostanoid receptors expressed in CHO cells; it does not
validate a mouse hepatocyte response. The inferred activating signs follow the
database's agonist annotations, and the physiological relevance of these
non-exclusive receptor interactions remains conditional on concentrations and
cellular context. [Kiriyama et al., 1997](https://pubmed.ncbi.nlm.nih.gov/9313928/).

The selected LPC weight is −0.37735 and its inferred MOON score is −1.17765.
The scores for arachidonate, PGH₂, PGF₂α/PGD₂ and Ptger3 are approximately
−4.021, −2.658, −1.735 and −2.250. These indicate signed alignment with the
testes-positive Factor 1 signature. They are not measurements of the mediator
concentrations, probabilities of causality or independent tests of receptor
activation. Although 64 TF observations are connected, all routes cross the
same Ptger3 entry point: these are not 64 independent signalling discoveries.

## Expected TFs and the alternative DHA–Rxra route

Esr1 and Dbp are present in the selected network at level zero, with supplied
activities −7.08875 and −7.66088. Illustrative paths end through
Ptger3→Egfr→Esr1 and Ptger3→Egfr→Jak1→Dbp. Their presence is compatible
with their measured footprints, but is **not independent recovery or prediction
of their activities**, because they were used as observations. Nfic has no
retained topology in the all-TF policy and is absent from the selected network.
Removing other TF terminal boundaries in the strong-TF sensitivity permits Nfic
to be an observation, but that sensitivity supplies no compact solution.

The normalized graph also contains an LPC(16:0)→DHA→Rxra topology. In the
all-TF fit, however, LPC(16:0)'s inferred score (−0.11976) opposes its positive
measured weight (+0.08679), so it is correctly excluded by the upstream sign
filter. Rxra's supplied footprint is weakly negative (−0.74957). In the
strong-TF sensitivity, this weak observation is omitted and Rxra instead becomes
a positively inferred intermediate (~+0.60), allowing positive LPC(16:0).
That sign reversal reflects observation-set sensitivity, not convergent support.
DHA and Rxra are absent from the selected 455-edge network. The corresponding
reaction cycle also requires unmeasured DHA-CoA and regenerates LPC; see the
[DHA–Rxra evidence audit](cosmos_factor1_lipid_tf_evidence.md).

The score audit independently recomputed every retained score in all 16 runs on
the actual final compressed scoring graph, with **maximum difference zero**.
LPC(18:0) has two preceding-layer compressed children supporting its score;
LPC(16:0) has one. Compression prevents duplicated graph neighborhoods from
inflating support, but it does not create independent biological evidence or
remove correlation among overlapping TF footprints. Per-input loss stages and
score support are retained in
[`score_audit/`](../results/causal_modelling/factor1_lipid_to_tf/score_audit/).

## Reproduce and audit

The main entry point is
[`scripts/moon/run_factor1_lipid_to_tf.R`](../scripts/moon/run_factor1_lipid_to_tf.R),
with [configuration](../scripts/moon/config/factor1_lipid_to_tf.json) and
[workflow plan](workflow_plans/moon_factor1_lipid_to_tf.md). From the repository
root, use a fresh output directory:

```sh
Rscript scripts/moon/run_factor1_lipid_to_tf.R --validate-only --include-ph73
Rscript scripts/moon/run_factor1_lipid_to_tf.R --include-ph73 \
  --output=results/causal_modelling/factor1_lipid_to_tf_replay
Rscript scripts/moon/ancillary/refine_lipid_tf_thresholds.R \
  --results=results/causal_modelling/factor1_lipid_to_tf_replay
Rscript scripts/moon/ancillary/audit_lipid_tf_scores.R \
  --analysis=results/causal_modelling/factor1_lipid_to_tf_replay
```

The refinement command reads saved fits, performs the finer search, extracts
the documented 1.65/1.10 representative and calls the provenance exporter.
It does not refit scores or trim the graph manually. Omit `--include-ph73` to
reproduce only the original-PKN feasibility results; no subsequent refinement
is warranted when there are no nonempty fits. Additional normalized-rule and
reaction-evidence audit commands are provided in their output READMEs.

The software contract checks cosmosR 1.21.2 at commit
`b2ed195738479c4ad204cace37fd11ae2fd02bb8`, decoupleR 2.17.0, igraph 2.3.3 and
jsonlite 2.0.0; R 4.6.1 was used. Python ancillary tools use the standard library.
No live external resources are required to replay the frozen analysis.

| Audit question | Files under `results/causal_modelling/factor1_lipid_to_tf/` |
| --- | --- |
| What was measured and how was it oriented? | `input_decision.json`, `run_manifest.json`, frozen input bundle referenced there |
| Why are most lipids unavailable as signalling inputs? | `coverage/summary.json`, `measured_mapping_coverage.csv`, `example_paths.csv` |
| Which configurations failed or were empty? | `run_summary.csv`, `threshold_grid.csv`, `analysis_summary.json` |
| Can compression, scores and coherence be reproduced? | `runs/<run_id>/run_state.rds`, `compression.rds`, final scoring/original-ID PKNs and ledgers |
| What changed in finer extraction? | `refinement/threshold_grid.csv`, `refinement_summary.json`, selected `threshold_neighborhood.csv` |
| Are selected edges and measurement roles valid? | Selected `extraction_audit.json`, `network_signaling_audit.json`, SIF/ATT and GraphML |
| Where did each edge and reaction come from? | Selected `network_evidence.csv`, `network_evidence.jsonl.gz`, `network_candidate_path_steps.csv` and evidence note |
| Why were LPC16, Nfic or Rxra absent? | `score_audit/lipid_input_audit.csv`, `focus_tf_audit.csv`, `two_lpc_score_summary.csv` |
| Does restoring normalized measured-species rules change the fit? | `normalized_rules_topology/pruning_equivalence.csv` and its exact signed-edge/compression hashes |

The coarse `analysis_summary.json` correctly reports zero eligible candidates;
the later `refinement/refinement_summary.json` reports the 26 finer-grid
candidates. All solution edges were checked against the final fitted PKN for
sign agreement, observed-TF terminal status and upstream-to-observation
connectivity. The exporter additionally checks exact signed prior evidence,
retaining original endpoints and explicit normalization provenance. The 45
previously regenerated forward-analysis artifacts were checksum-verified as
unchanged during this work.
