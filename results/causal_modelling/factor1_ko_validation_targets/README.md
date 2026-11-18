# KO validation proposal for Factor 1 TF to lipid hypotheses

Read [the ranked report](../../../docs/cosmos_factor1_ko_validation_targets.md).
The revised ten-set allocation alternates five enzyme-only and five
TF-containing designs. Every pair tests a possible alternate producer of the
same lipid output; no current triplet has sufficient backup rationale.
Six singles/four pairs use nine genes. Add Pld2, Plcb1 and Plcg1 singles for
13 distinct target-containing conditions before controls and replication.

**Current exports: `alternating_prioritization/`.** Ranked sets, prespecified
lipid directions, exact paths, source-corrected hypothesis edges, exclusions,
review snapshot, validation summary and source hashes are included.
`expression_audit/` is reused: true counts from 43 FCG hepatocyte pseudobulks
and four historical culture libraries. Initial `prioritization/` outputs and
`initial_panel_20261013/` source snapshots are preserved but superseded.

Known limitations: Ppara→Plpp1 cites PPARG; Cebpa→Pld1 cites CEBPB;
Gata6→Mapk1 reverses the cited kinase/TF relationship; Enpp2 diacyl-PC→PA is
unsupported. Corrected Cebpb/Pparg prefixes are separately labelled unscored
biochemical follow-ups, not selected/sign-coherent Factor 1 explanations.
RXRA's remaining source requires VDR/context availability. TF/PLD2 designs
require PLD1 suppression and retained PLD2 before compensation interpretation.
No adaptive hepatic compensation, stronger pair effect or net steady-state
direction is guaranteed. Original graphs/PKNs/scores remain frozen; no new fit.


The updated shareable PDF is
[`output/pdf/cosmos_factor1_ko_validation_targets.pdf`](../../../output/pdf/cosmos_factor1_ko_validation_targets.pdf).
It includes ten full-route annotations with KO targets highlighted, reaction
chemistry and separate model/unscored/unsupported upstream status. Full path
records and exact reaction participants are in
`route_annotation_revision_20261013/full_tf_route_audit.json` and
`full_tf_to_lipid_routes.csv`; prior report/PDF snapshots are retained there.
`pdf_export_audit.json` records text/link and page-render verification.
