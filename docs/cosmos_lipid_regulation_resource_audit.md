# Audit of the promised and implemented lipid-to-protein resources

The cache-audit sections below were completed 15 September 2026 by reading the frozen downloads, evidence tables, source manifests, builder and assembled graphs. No cache source was refreshed, no input mapping was changed, and no biological analysis or MOON fit was rerun. The subsequent directionality-audited overlay section documents a separately frozen raw-STITCH retrieval. This note concerns resource admission and provenance; measured-feature coverage is a separate audit.

## Correction to the initial scope description

The [feasibility document](cosmos_lipid_pkn_feasibility.md#proposed-network-layers) listed **preferred resources for a proposed network**, including BRENDA, MRCLinksDB, CellPhoneDB and an additional STITCH tier. That table was not a record of the completed build. The actual conservative analysis PKN contains **GtoPdb-derived chemical-to-protein regulation only**. Cached BRENDA and STITCH were retained in the separate exploratory PKN; MRCLinksDB was inspected and rejected; CellPhoneDB was not integrated. No fresh BRENDA rebuild was implemented.

Thus a description claiming that the analysis PKN incorporated all four originally named regulation resources would be incorrect. The implemented restrictions are recorded in [regulation_build_notes.md](../results/causal_modelling/cosmos_mouse_pkn/regulation_build_notes.md) and [regulation_manifest.json](../results/causal_modelling/cosmos_mouse_pkn/regulation_manifest.json), and enforced in code. They limit the resource that was actually tested; they do not establish that those databases have no potentially useful lipid knowledge.

## Downloaded records and their disposition

Counts below distinguish downloaded records, retained **evidence rows**, and unique graph edges. An evidence row is not necessarily an independent experiment or a lipid-specific interaction.

| Resource | What was frozen | Retained regulation evidence | Admission to the analysis PKN |
| --- | --- | --- | --- |
| GtoPdb | Four CSVs from release 2026.2, published 15 June 2026, frozen 8 September; 24,599 interaction records, including 1,226 mouse-target records and 125 mouse-target metabolite records | 83 metabolite evidence rows with direct GtoPdb ChEBI cross-references; 20 rows from 19 records with secondary PubChem chemical mappings | Only the 83 direct-mapping rows were eligible for the curated tier; the 20 secondary-mapping rows were exploratory |
| BRENDA | 786 rows inside the downloaded OmniPath-metabo mouse cache; no separate fresh BRENDA download | 772 fallback rows; 731 distinct chemical/UniProt/sign triples before assembly conflict handling | Exploratory only; 14 cached targets failed validation against the frozen mouse UniProt export |
| STITCH | 40,428 rows inside the same mouse cache; no separate fresh STITCH download | 40,360 fallback rows; 25,495 distinct chemical/UniProt/sign triples before assembly conflict handling | Exploratory only; 68 cached targets failed mouse UniProt validation |
| MRCLinksDB | 641 rows inside the mouse cache: 429 labelled ligand–receptor and 212 labelled transport | None | All 641 rejected because the upstream adapter's positive sign was unconditional, not an independently checked activation effect |
| CellPhoneDB | No records in the frozen mouse cache, no separate download in the source manifests, and no adapter invocation in the local regulation builder | None | Not implemented in this build; the recorded rationale was a human-only nonpeptidic adapter with assumed positive effects, which was not silently projected to mouse |

The GtoPdb supporting CSVs contain 13,856 ligand-mapping rows, 13,856 ligand rows and 3,383 target/family rows. PubChem synonyms were downloaded for 18 missing metabolite identifiers. These support chemical mapping, not new biological interactions. Of the 125 mouse-target metabolite records, 23 were excluded: nine lacked a ChEBI mapping, four had unresolved/multisubunit targets, eight had chemical-connectivity discrepancies and two had stereochemical/protonation discrepancies. The 102 accepted source records become 103 evidence rows because one secondary mapping has two chemical alternatives. A further 179 nonmetabolite pharmacological rows were saved separately and are not read by the PKN integrator.

The full raw mouse cache contains 235,302 records from eight resource labels, with no CellPhoneDB label. Its component timestamps are 21 April 2026; retrieval was 8 September 2026. These timestamps are **not BRENDA or STITCH release identifiers**. The service export does not supply the underlying release ledger, original chemical-mapping alternatives or per-edge literature citations for those fallback rows. [Cached resource counts](../results/causal_modelling/cosmos_pkn_audit/cached_resource_audit.csv); [cached response and source hashes](../results/causal_modelling/cosmos_pkn_audit/manifest.json); [regulation source manifest](../results/causal_modelling/cosmos_mouse_pkn/regulation_manifest.json).

## What the assembled graphs actually contain

The integrator places only `evidence_tier == "curated_mouse_current"` in `lipid_curated`; all other regulation tiers enter `exploratory`. The `lipid_rules` variant inherits the curated regulation layer and adds metabolic reaction hypotheses, not extra lipid-to-protein regulatory evidence. Opposing signs for the same assembled node pair are removed from each exported graph, while raw evidence is preserved. [Admission logic, lines 215–226](../scripts/pkn/assemble_mouse_pkn.py#L215); [variant and conflict export, lines 242–264](../scripts/pkn/assemble_mouse_pkn.py#L242).

| Exported graph | Distinct GtoPdb-supported regulatory edges | Distinct BRENDA-supported regulatory edges | Distinct STITCH-supported regulatory edges | Total distinct chemical-to-protein edges |
| --- | --- | --- | --- | --- |
| `pkn_lipid_curated.csv` | 77 | 0 | 0 | 77 |
| `pkn_lipid_rules.csv` | 77 | 0 | 0 | 77 |
| `pkn_exploratory.csv` | 97 | 707 | 23,552 | 24,342 |

Resource columns in the exploratory row overlap and therefore should not be summed. Counts were obtained by matching `edge_evidence.jsonl.gz` to the exact signed triples present in each exported CSV, preserving resource identity.

The 83 curated GtoPdb evidence rows initially represent 81 unique signed triples. Four rows concern two conflicting pairs, ATP–Kcnj8 and GTP–Kcnj8, with both agonist and antagonist signs. Removing those pairs leaves 79 supporting evidence rows and 77 unique graph edges. The raw count of 83 must therefore not be reported as the final number of regulatory graph edges.

Across all tiers, [regulation_summary.json](../results/causal_modelling/cosmos_mouse_pkn/regulation_summary.json) reports 41,235 evidence rows, 26,245 unique pre-assembly chemical/UniProt/sign triples, 4,343 chemical IDs and 3,225 protein accessions. Those are broad chemical-regulation counts, **not counts of curated lipid signaling in the analysis PKN**. The final PKN's approximately 250,000 total edges primarily describe protein networks and reaction projections; its total size does not quantify the lipid-ligand regulation layer.

## What “mouse” and “signed” mean here

For GtoPdb, the builder selects records whose target species is explicitly Mouse, then requires a single target accession present in the frozen mouse UniProt export. It does not perform human-to-mouse orthology projection for this layer. However, **a mouse-target database annotation is not proof that every cited paper tested mouse hepatocytes**. The separate [DHA–RXR audit](cosmos_factor1_lipid_tf_evidence.md) illustrates that distinction. Original references and assay fields are preserved; a database action label and a paper's directly measured outcome must be interpreted separately.

For cached BRENDA and STITCH, the downloaded cache is labelled organism 10090 and retained target accessions are validated as mouse. The builder itself fills `organism` and `source_organism` with 10090. Given missing per-edge provenance, these fields do not establish native-mouse experimental evidence independently of the cache's species assignment. All 772 BRENDA and 40,360 STITCH fallback rows have empty literature-reference and assay-description fields.

GtoPdb positive/negative signs come from explicit pharmacological action labels, with an explicit activator/inhibitor type used only when action is blank. Binding, neutral, mixed, biphasic and unknown effects are not assigned signs. STITCH rows must contain an activation/inhibition action consistent with their sign; the original actor orientation cannot be reconstructed from this cache. BRENDA signs are inherited from the signed cached allosteric records. None of these entries is a direct lipid-to-transcript expression edge. [Sign parsing, lines 293–313](../scripts/pkn/build_lipid_regulation.py#L293); [GtoPdb admission, lines 364–444](../scripts/pkn/build_lipid_regulation.py#L364); [cache admission, lines 446–488](../scripts/pkn/build_lipid_regulation.py#L446).

## Precise implementation pointers

- [build_lipid_regulation.py, lines 31–36](../scripts/pkn/build_lipid_regulation.py#L31): source paths and the four GtoPdb filenames. Several frozen paths predate the move into `results/causal_modelling`; this audit reads the relocated artifacts rather than executing that builder.
- [Lines 78–90](../scripts/pkn/build_lipid_regulation.py#L78): `refresh()` downloads GtoPdb and optional PubChem mappings. There is no fresh BRENDA, MRCLinksDB or CellPhoneDB retrieval in this function.
- [Lines 451–458](../scripts/pkn/build_lipid_regulation.py#L451): the cache loop considers BRENDA, STITCH and MRCLinksDB, then rejects MRCLinksDB before other admission tests.
- [Lines 476–488](../scripts/pkn/build_lipid_regulation.py#L476): explicit cached fallback tiers and missing-provenance caveats.
- [Lines 505–509](../scripts/pkn/build_lipid_regulation.py#L505): separate evidence, sign conflicts, pharmacological supplement, exclusions and mapping-review exports.
- [Lines 541–550](../scripts/pkn/build_lipid_regulation.py#L541): recorded species, sign, fallback and CellPhoneDB limitations. The CellPhoneDB explanation is a documented upstream-adapter assessment; no local CellPhoneDB records were independently imported or evaluated in this audit.
- [regulation_exclusions.csv](../results/causal_modelling/cosmos_mouse_pkn/regulation_exclusions.csv): every recorded exclusion, including all 641 MRCLinksDB rows and the 82 failed cached mouse-target validations.
- [build_manifest.json](../results/causal_modelling/cosmos_mouse_pkn/build_manifest.json): frozen raw-file, build-code and output hashes. It contains GtoPdb, PubChem, ChEBI and the shared cached response, but no fresh BRENDA or CellPhoneDB download.

A corrected description is: **“The conservative mouse PKN integrated current GtoPdb metabolite–protein activity annotations with signaling, transcription and lipid-reaction resources. Cached BRENDA/STITCH regulation was archived in an exploratory variant; MRCLinksDB associations and the human CellPhoneDB nonpeptidic adapter were not admitted to the conservative analysis.”**

## Independent check of current upstream adapters

Read-only source inspection on 15 September 2026 verified OmniPath-metabo commit `ffe066a3fe3a87f390a05cfc5b84bdc3a62761e9` and pypath commit `33f37fbaab59993d24f5c32bb4b3e7587085bcf4`. The observations below describe those upstream implementations and possible improvements, **not additional resources imported into the PKN used here**. No adapter was executed to rebuild biological evidence, and no PKN or MOON result changed.

- **MRCLinksDB:** the adapter assigns `mor=1` unconditionally after classifying targets as ligand–receptor or transport. Its separate transporter adapter also assigns `+1`. Consequently, these exported positive signs cannot be treated as independently checked activation effects. The relationships may still support future action-specific curation; excluding the assumed signs does not establish biological absence. [Ligand–receptor/transport adapter](https://github.com/saezlab/omnipath-metabo/blob/ffe066a3fe3a87f390a05cfc5b84bdc3a62761e9/omnipath_metabo/datasets/cosmos/resources/mrclinksdb.py#L111-L125), [dedicated transporter adapter](https://github.com/saezlab/omnipath-metabo/blob/ffe066a3fe3a87f390a05cfc5b84bdc3a62761e9/omnipath_metabo/datasets/cosmos/resources/mrclinksdb.py#L194-L204).

- **CellPhoneDB:** the nonpeptidic adapter returns no records for organisms other than human. It parses labels such as `Adenosine_byNT5E_and_SLC29A1` into a chemical ligand and producing machinery; the machinery remains metadata, not enzyme-to-metabolite edges. Receiving proteins/complexes are resolved separately, and the resulting nonpeptidic edges receive an unconditional `+1`. Thus neither biosynthetic machinery labels nor the emitted sign should be mistaken for independently established chemical identity and functional activation. A mouse extension would require explicit orthology, ligand identity, action and complex checks. [Parsing semantics](https://github.com/saezlab/omnipath-metabo/blob/ffe066a3fe3a87f390a05cfc5b84bdc3a62761e9/omnipath_metabo/datasets/cosmos/resources/cellphonedb.py#L26-L41), [species restriction](https://github.com/saezlab/omnipath-metabo/blob/ffe066a3fe3a87f390a05cfc5b84bdc3a62761e9/omnipath_metabo/datasets/cosmos/resources/cellphonedb.py#L200-L230), [emitted sign](https://github.com/saezlab/omnipath-metabo/blob/ffe066a3fe3a87f390a05cfc5b84bdc3a62761e9/omnipath_metabo/datasets/cosmos/resources/cellphonedb.py#L267-L285).

- **BRENDA:** the adapter genuinely distinguishes activators (`+1`) and inhibitors (`−1`). Its underlying pypath parser retrieves the bulk textfile after license acceptance and retains PubMed references and reaction constants, but the OmniPath-metabo wrapper does not copy those fields into its interaction records. Fresh BRENDA integration with retained provenance is therefore a reasonable improvement, not something shown to be impossible because SOAP credentials are unavailable. The current official bulk download is release 2026.1 under CC BY 4.0; the large data files were not downloaded in this audit. [Signed adapter](https://github.com/saezlab/omnipath-metabo/blob/ffe066a3fe3a87f390a05cfc5b84bdc3a62761e9/omnipath_metabo/datasets/cosmos/resources/brenda.py#L48-L83), [bulk retrieval code](https://github.com/saezlab/pypath/blob/33f37fbaab59993d24f5c32bb4b3e7587085bcf4/pypath/inputs/brenda/_main.py#L108-L116), [underlying evidence fields](https://github.com/saezlab/pypath/blob/33f37fbaab59993d24f5c32bb4b3e7587085bcf4/pypath/inputs/brenda/_main.py#L293-L306), [official download page](https://www.brenda-enzymes.org/download.php).

- **STITCH:** the current wrapper filters directed activation/inhibition records by score, but then normalizes every retained chemical–protein pair to chemical→protein. In pypath, `directed` means that original partner A is acting; partner A may be a protein. Therefore that wrapper alone does not guarantee chemical-as-actor orientation and can reverse a directed protein→chemical record. A fresh signed import should explicitly require the original acting source to be chemical and preserve its action, direction and evidence scores. This code-level concern does not quantify affected rows in the frozen cache, whose original actor information is absent. [Filtering and orientation normalization](https://github.com/saezlab/omnipath-metabo/blob/ffe066a3fe3a87f390a05cfc5b84bdc3a62761e9/omnipath_metabo/datasets/cosmos/resources/stitch.py#L268-L313), [underlying actor flag and entity order](https://github.com/saezlab/pypath/blob/33f37fbaab59993d24f5c32bb4b3e7587085bcf4/pypath/inputs/new_stitch/_actions.py#L38-L53).

These checks reinforce the distinction between **absence from the tested conservative graph** and absence from available biological knowledge. They identify concrete integration and provenance work still needed before claiming comprehensive signed lipid-to-protein coverage.

## Directionality-audited STITCH v5 overlay

A separate exploratory overlay now replaces every cache-derived STITCH record
from the original exploratory build with raw records imported directly from the
mouse STITCH v5 action and detailed-link files. It is saved in
[`results/causal_modelling/cosmos_mouse_pkn_stitch_directional/`](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/).
The original 40,428/40,360 cache counts and the original 303,925-edge
exploratory graph above remain historical facts about the frozen base build;
they were not edited in place.

A raw record is admitted only if original partner A is a `CIDm`/`CIDs`
chemical, partner B is a mouse `ENSMUSP` protein, `a_is_acting=t`, `action` is
exactly `activation` or `inhibition`, the action score is at least 700, and its
detailed-link text-mining score is below 700. The sign comes solely from
`action` (+1 for activation and −1 for inhibition); `mode` is retained only as
provenance. Chemical mapping requires one current ChEBI target after explicit
alternate-ID canonicalisation through the frozen UniChem source mapping;
protein mapping requires one exact `UniProt_AC` alias in the frozen mouse
UniProt export. The current STRING v12 alias snapshot is used for this
identifier step because a matched v5 alias archive is no longer publicly
available; that limitation is retained per edge and does not alter the raw
action orientation.

The overlay contains 7,481 retained raw evidence rows, representing 4,336
unique projected signed triples; 4,243 survive graph-wide opposing-sign
quarantine. It replaces all 40,360 cache STITCH rows and exports a
284,602-edge signed exploratory graph. The complete evidence table preserves
raw endpoints, acting flag, signed action, score, detailed-link channels,
mapping basis and hashes. Inspect
[`stitch_directional_evidence.csv`](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/stitch_directional_evidence.csv),
[`stitch_directional_exclusions.csv`](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/stitch_directional_exclusions.csv),
[`stitch_directional_mapping_review.csv`](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/stitch_directional_mapping_review.csv),
[`conflicts_exploratory_stitch_directional.csv`](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/conflicts_exploratory_stitch_directional.csv),
and the [source manifest](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/stitch_directional_manifest.json).
[`stitch_directional_validation.json`](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/stitch_directional_validation.json)
verifies the source-direction predicates, removal of cached STITCH rows,
graph projection, conflict policy and immutability of the base PKN; the
[R validation](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/r_validation.json)
confirms that `cosmosR::meta_network_cleanup()` leaves the overlay unchanged.

This verifies source-record directionality. It does not establish a measured
lipid's molecular identity, hepatocyte context, direct experimental causality
or a current resource release: the mouse STITCH v5 files were last modified on
16 August 2016. The overlay is not included in `pkn_lipid_curated` or
`pkn_lipid_rules`. It is used only in the separately documented relaxed
lipid-to-TF sensitivity after a conditional measurement-to-chemical mapping
audit; it is not used in the manuscript TF-to-lipid or original reaction-network
lipid-to-TF fits. That use remains exploratory and does not convert the overlay
into a current or conservative resource.
