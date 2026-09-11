# Esr1, Nfic and Dbp retention audit

All three TFs were scored and then removed by the installed `cosmosR::moon()` upstream-sign filter. They were **not pruned out before scoring**. This holds for every audited compressed radius-6 run at depths 4, 5, 6, 7, 8, 10 and 12, for the uncompressed radius-6/radius-7 depths 4–6 and radius-12 depth-6 controls (42 TF/run comparisons), and for eight additional compressed sensitivities (24 more TF/run comparisons): **66 focus TF/run comparisons in total**.

The diagnostic recomputed the saved final scoring PKNs with `upstream_input = NULL` solely to expose scores hidden by the final sign filter. These unconstrained scores were not used to generate a replacement solution. All saved returned scores were reproduced within 1e-10 from the exact RDS inputs. The three focus TFs' raw scores were already sign-incompatible before the iterative RNA-coherence procedure.

| TF | Factor 1 footprint activity | Initial compressed MOON score | Final compressed MOON score | First scoring level | Active children in that layer |
|---|---:|---:|---:|---:|---|
| Esr1 | -7.08874857626001 | +0.131796253954624 | +0.131796253954624 | 3 | Lpl |
| Nfic | -4.34923126849243 | +0.165143185835450 | +0.718589712937770 | 4 | Sp1, Nr3c1, Src |
| Dbp | -7.66088476896073 | +0.369642797826770 | +0.345155683867094 | 4 | Esr1 |

These compressed values are identical across all seven audited scoring depths. Factor 1 is oriented testes-positive; the negative footprint activities therefore correspond to the opposite direction of this latent axis. The positive lipid-driven MOON scores disagree with these supplied signs. A generic role of ESR1 in hepatocyte physiology does not establish that this particular Factor 1 lipid signature supports the same inferred direction through the measured-species PKN.

## Exact stage checks

All three nodes are in the raw PKN, pass `abs(TF_activity) > 2`, remain in the radius-6 pruned PKN and upstream vector, are protected from compression with unchanged identifiers, and remain in the final scoring graph. No outgoing edges from these three TFs were removed by the post-score coherence loop. Their initial sign conflict prevents them from being among the returned TF rows used by that loop; other TFs' filtering can nevertheless change their layer background.

| TF | Raw outgoing edges | Pruned outgoing edges | Final compressed outgoing edges | Outgoing coherence removals | Incoming coherence removals at depth 6 |
|---|---:|---:|---:|---:|---:|
| Esr1 | 536 | 302 | 300 | 0 | 48 |
| Nfic | 68 | 38 | 38 | 0 | 6 |
| Dbp | 34 | 15 | 15 | 0 | 3 |

Esr1's 302 to 300 outgoing-edge change is caused by compressed target identities collapsing duplicate projected edges, not by removal of the Esr1 node. The stage audit preserves the distinction between original-node edges and compressed scoring edges.

## Why the inferred signs disagree

MOON's ULM score is a regression/correlation statistic over the complete preceding layer, with zero coefficients for non-targets. It is not the uncentered signed mean of a TF's target scores.

- **Esr1:** only Lpl contributes as an immediate child at its first eligible layer. Lpl's score is -0.108966299655943, compared with a background mean of -0.236806209456588 across 50 preceding-layer nodes. An activating target that is less negative than the layer background produces a weak positive ULM score. Removing that single child leaves no informative coefficient variation, so this score has no multi-target support.
- **Nfic:** three activating children occur in the 231-node preceding layer: Sp1 (+0.534331447984020), Nr3c1 (+0.131796253954624), and Src (-0.055203197214175). Each is above the layer mean (-0.235269038951820), giving a positive score. Removing any one of those three layer features leaves the reconstructed score positive.
- **Dbp:** only Esr1 contributes at the first eligible layer. Esr1's +0.131796253954624 score is above the same -0.235269038951820 layer mean. Because `moon()` applies upstream-sign filtering after all layers have been computed, Esr1's raw score contributes to Dbp before both rows are omitted from the returned results.

All five active edges have transcription provenance in the local OmniPath GRN evidence. Esr1→Lpl is absent from the original TF inference regulon, which contains a custom receptor footprint. Provenance presence is not proof that the measured mouse lipid species determine the TF's mechanism in this tissue.

## Which tuning is reasonable?

Increasing `n_layers` beyond 4 cannot rescue these three nodes in the present topology. MOON scores each node at its first eligible layer and removes it as a candidate source for subsequent layers; it does not revisit the node when longer paths become available. Topology radius 7 and the broad radius-12 uncompressed control retain the same positive sign conflicts. Lowering extraction thresholds also cannot reinstate a row that the upstream sign filter has removed.

A defensible extra sensitivity is a **uniform initial RNA-coherence filter for all measured TF activity inputs**, applied before compression, with the existing Factor 1 activity and RNA signs. Restrict deletion to signed triples confirmed as transcription edges in the PKN evidence and matching the input regulon's sign, and remove a triple only when `sign(TF_activity * RNA_target_weight * edge_sign) < 0`. This is an explicit consistency prior and uses no new independent evidence. It was tested as a comparison with the original analysis, not tuned separately for Esr1, Nfic or Dbp.

All active scoring edges listed above already agree with the corresponding TF-activity × target-RNA sign. A pre-score filter would therefore act on other edges and layer backgrounds, rather than directly deleting the observed conflicting path. Among the pruned original outgoing edges, the signed input-regulon overlap and potential conflicts are Esr1 3/302 (1 conflict), Nfic 36/38 (18 conflicts), and Dbp 14/15 (4 conflicts); confirm transcription provenance before removal. The sparse Esr1 overlap follows the existing custom footprint and is not justification for silently changing its input mapping.

Do not drop the upstream sign boundary, flip a TF activity, selectively remove Lpl or another inconvenient target, or lower the input-activity cutoff solely to force these controls into a network. The first two would hide the disagreement; the latter two do not fix the limited informative downstream support.

## Additional compressed sensitivities completed

The following independent checks retain same-signed-children compression and the accepted input mapping. The topology checks were scored to a maximum allowed depth of 8. The RNA-coherence and inputs-only compression policies were each tested at depths 4 and 8 to compare processing policies at matching depths; all settings reach these TFs at levels 3 and 4.

| Sensitivity | Esr1 raw score | Nfic raw score | Dbp raw score | Focus TFs returned after sign filtering |
|---|---:|---:|---:|---:|
| Topology radius 4 | +0.098596951 | +0.628108591 | +0.334129168 | 0/3 |
| Topology radius 5 | +0.105523814 | +0.686601741 | +0.347830517 | 0/3 |
| Topology radius 7 | +0.131796254 | +0.718003517 | +0.344692799 | 0/3 |
| Topology radius 12 | +0.089731404 | +0.914216833 | +0.369042842 | 0/3 |
| Radius 6, initial RNA-coherence filter, depth 4 | +0.132701939 | +0.729162439 | +0.352345858 | 0/3 |
| Radius 6, initial RNA-coherence filter, depth 8 | +0.132701939 | +0.729162439 | +0.352345858 | 0/3 |
| Radius 6, protect measured inputs only during compression, depth 4 | +0.131796254 | +0.699337886 | +0.334276675 | 0/3 |
| Radius 6, protect measured inputs only during compression, depth 8 | +0.131796254 | +0.699337886 | +0.334276675 | 0/3 |

Across the final frozen graphs of the original and added compressed checks, Esr1 remains between +0.0897 and +0.1327, Nfic between +0.6281 and +0.9142, and Dbp between +0.3341 and +0.3690. No setting resolves their sign conflict, and every added setting retains the same informative-child counts: Esr1 one, Nfic three, Dbp one.

The initial RNA-coherence sensitivity evaluated **all 49 initial measured TF inputs**, before any MOON sign filtering. It removed 1,750 edges from 48 TFs. Independent reconstruction confirmed exactly the expected deletion set and that every deleted signed triple has both transcription provenance and an agreeing signed pair in the original TF regulon. It was applied once before topology pruning/compression; the usual bounded inferred-score coherence loop then ran to convergence. The focus TFs' retained active child edges remain compatible with their original RNA signs.

The inputs-only compression sensitivity protects all 49 measured TFs and 74 measured lipids, permits other regulon nodes to merge, decompresses scores for original-identifier coherence filtering, and recompresses after every coherence iteration. At depth 8 it converged in five iterations, with 3,415 original nodes merged into 606 groups and 24,054 final scoring edges corresponding to 34,228 original-node edges. At depth 4 it converged in three iterations, with 3,391 nodes merged into 600 groups and 25,039 scoring edges corresponding to 35,275 original-node edges. Independent validation checked every saved iteration's original-to-compressed projection, the final saved score/checkpoint identity, and the final coherence fixed point. More complete compression and rescoring therefore leave the focus TFs' exclusion reason unchanged, including in the matched-depth comparison.

All eight new runs passed independent checks of final graph projection, protected input identity, final coherence convergence, score-copy agreement, positive focus-TF raw scores, and post-score omission. The raw diagnostic CSV scores were sufficient to reconstruct each focus score to within 1e-10 without refitting a MOON run. The matched-depth additions leave the score ranges above unchanged.

## Files and reproduction

- `focus_tf_stage_audit.csv`: exact exclusion stage and graph/input checks per TF/run.
- `all_upstream_tf_sign_audit.csv`: corresponding sign audit for all supplied TFs.
- `focus_tf_score_background.csv`: reconstructed ULM scores and target/background summaries.
- `focus_tf_layer_influence.csv`: per-feature covariance and leave-one-layer-feature-out score diagnostics; these are descriptive recalculations, not complete graph rescoring.
- `focus_tf_outgoing_edges_by_stage.csv`: raw, pruned and final original-node outgoing edges, RNA signs and input-regulon agreement.
- `exact_transcription_evidence.csv`: OmniPath GRN evidence for matching signed transcription triples.
- `audit_full_precision.rds` and `*_diagnostic_scores.rds`: full-precision numeric results; CSVs are display/interchange copies.
- `installed_moon_source.txt`, `installed_coherence_source.txt`, `session_info.txt`: implementation and environment record.
- `additional_sensitivity_tf_audit.csv`, `additional_sensitivity_validation.csv`, `prescore_filter_validation.csv`: independent validation of the eight additional compressed sensitivities.

Run from the project root:

```sh
Rscript scripts/moon/audit_factor1_focus_tf_retention.R
Rscript scripts/moon/validate_factor1_compressed_tf_sensitivities.R
```
