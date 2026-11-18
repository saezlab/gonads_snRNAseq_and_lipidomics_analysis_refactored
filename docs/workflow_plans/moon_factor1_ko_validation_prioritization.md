# Factor 1 TF to lipid KO validation prioritization

Purpose: propose ten experimental designs from source-screened frozen TF→lipid
networks and explicitly separate literature-derived unscored TF extensions.
The revised allocation alternates five enzyme-only and five TF-containing
sets. Every pair addresses a plausible alternate producer of the same lipid;
this is not a MOON refit, flux model or quantitative intervention simulation.

## Inputs

- Frozen four/six-layer selected ATT/SIF files under
  `results/causal_modelling/factor1_tf_to_lipid_expanded_species/`.
- Full-precision Factor 1 inputs and one-to-one measured-lipid mapping ledger.
- Existing expression audit from true FCG counts `data/RNA/fcg_pseudobulks.rds`
  and culture counts in
  `miscellaneous/older_inputs/legacy_rna/merged_counts_matrix_250508_VH00211_433_AACVFM3HV_with_gene_symbols_and_groups.csv`.
- Explicit ranking, roles, exclusions, regulatory prefixes and compensation
  gates in `support/cosmos_mouse/factor1_ko_validation_review.json`.
- Verified sources in `support/cosmos_mouse/ko_validation_literature_sources.json`.

## Entry points and invariants

`ancillary/audit_factor1_ko_expression.R` under `scripts/moon/` separates
abundance evidence from Factor 1 loadings/activity scores. Its unchanged audit
contains 38 genes, 43 FCG hepatocyte mice and four historical culture libraries.
The full run, validation mode, help and overwrite guard previously passed.

`ancillary/propose_factor1_ko_validation.py` under `scripts/moon/` enumerates
screened frozen paths and separately prefixes Cebpb→Pld1 and Pparg→Plpp1.
These prefixes never mutate the PKN or graph: their score layers/full-path
MOON score are null; downstream frozen-path evidence is labelled separately.
Cebpa/Ppara scores cannot be transferred to corrected paralogs. The script
checks original score layers/signs, strict category alternation, unique sets,
shared requested outputs for every pair member, no serial co-targets, null
corrected-TF scoring provenance and unchanged original graph hashes. Local
KO direction is −1 times the directed path-sign product, regardless of the
baseline score sign. Validate immediately after edits with `--validate-only`;
full exports require a fresh `--output` directory and refuse overwrites.

## Current results and experimental controls

The revised full export and validation-only run passed on 13 October 2026:
27 set/readout rows and 52 ranked paths (44 frozen, eight source-corrected).
Negative checks rejected wrong category alternation, a serial Rxra/Pld1 pair
and nonempty output reuse. All eight frozen graph file hashes and 14 current
source-manifest hashes were verified; no expression script changes were needed.

Report: `docs/cosmos_factor1_ko_validation_targets.md`.
Output root: `results/causal_modelling/factor1_ko_validation_targets/`.
Current tables are in `alternating_prioritization/`; initial sources are in
`initial_panel_20261013/` and initial exports remain in `prioritization/`.
Existing `expression_audit/` is reused. Ranked sets, per-readout directions,
exact paths, exclusions, source-corrected hypothesis ledger, review snapshot,
validation summary and source hashes make the revised panel auditable.

The ranking is Pld1; Rxra; Pld1/Pld2; Cebpb; Pla2g4a; Cebpb/Pld2; Plpp1;
Pparg; Plcb1/Plcg1; Rxra/Pld2. Six singles/four pairs use nine genes. Add Pld2,
Plcb1 and Plcg1 singles for 13 distinct target-containing conditions.
Compensation means potentially retained functional capacity; it does not
assert adaptive hepatic upregulation. No triplet is justified by current
backup evidence. All Pld2 sets use six-layer path support.

## Biological limitations

Reject Ppara→Plpp1, Cebpa→Pld1 and Gata6→Mapk1 as source-supported chains
because of cited paralog/direction mismatches; reject Enpp2 diacyl-PC→PA.
RXRA evidence excludes a retracted paper and retains a VDR/vitamin-D gate.
The source-corrected Cebpb/Pparg paths are unscored follow-up hypotheses,
not selected or sign-coherent Factor 1 explanations. Pparg is absent from
both frozen PKNs. Cebpb gene loss affects multiple isoforms; transfer from
adipocytes or endothelial cells requires an actual mediator response.
Both TF/Pld2 pairs require TF-induced PLD1 reduction and retained PLD2
function before compensation interpretation. PLPP1 pool accessibility and
PLC physiological PIP2 versus measured PI remain unresolved.

Ranking balances experimental categories rather than assigning equal
confidence. TF→enzyme→lipid claims require the mediator response as well as
lipids. Existing graphs/scores were fitted before source screening and remain
frozen. PKN correction and rescoring are separate analysis.

## Full TF-to-lipid route annotations and PDF

The 13 October route annotation revision adds one entry per ranked set with
fitted upstream TF context, intermediaries, targeted steps and measured-panel
endpoint. `route_annotation_revision_20261013/full_tf_route_audit.json` and
`full_tf_to_lipid_routes.csv` retain full signed node IDs, source/input status,
KO positions and complete chemical records; 222 CSV records include both
four/six-layer paths, local subpaths and unscored extensions. Prior report/PDF
and renderer/QA snapshots are preserved in this revision's `before/` folder.
The rank/config, original networks/scores and existing local-direction exports
are unchanged.

The report distinguishes input Nr2f1/Nr4a2→RXRA context from the intermediary
RXRA footprint. NR2F1's source sign is questioned; NR4A2/RXRA association and
VDR/RXRA-mediated PLD1 regulation require partner-context compatibility.
Hoxa10→receptor→Gnai1 inhibits PLD2 and is six-layer model context only;
Hoxa10→Ednra evidence is endometrial, not hepatocyte-specific. PLCβ1's fitted
TF prefix requires six layers, while PLCγ1 has four-layer upstream routes.
No supplied-TF prefix survives declared source exclusions for Pla2g4a or
Plpp1. Enzyme-only experiments cannot establish unperturbed upstream links.
Cebpb/Pparg regulatory prefixes remain unscored and sign-discordant.

Primary endpoint chemistry is explicit. LPC18 local alternatives include
glycerol transacylation, whose coproduct is MAG rather than free fatty acid.
The prioritized PC(18:0_20:4) hydrolysis is therefore distinguished from all
other LPC18 rules. PLPP1's PA precursors differ from the PLD PA3 panel and
cannot be chained into one measured linear conversion.

`scripts/reports/render_cosmos_ko_validation_pdf.py` renders the Markdown to
`output/pdf/cosmos_factor1_ko_validation_targets.pdf`. The current report has
nine pages, with one landscape ranking page and two landscape route pages.
The renderer was run immediately after edits. Poppler rendering and full
page-by-page visual inspection are required; PDF text/link checks are stored
in `pdf_export_audit.json`. No R script edits or new MOON fit are involved.

## 16 November 2026 repository cleanup

Historical causal result folders now live under
`miscellaneous/older_results/causal_modelling/`. Active configurations and
consumer defaults explicitly use their frozen inputs/reference states there.
The per-file old/new path and SHA256 ledger is
`docs/provenance/causal_cleanup_20261116.csv`; the active/frozen/archive dependency
map is `docs/CAUSAL_ANALYSIS_MAP.md`. No numerical settings, graph policy, mapping,
score or fitted result changed. Saved-artifact extraction/reference replay passed
without fitting or scoring. Full-fit version gates remain unchanged; the installed
decoupleR 2.9.7 differs from the pinned 2.17.0 and blocks full-fit validation.
