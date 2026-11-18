# Measured-lipid coverage of the lipid-to-protein regulatory prior

Audit completed on 15 September 2026. This examines the frozen measurements, chemical crosswalk, downloaded source records and assembled PKNs. The audit itself adds annotation audits, not a biological fit: it did not change measurement assignments, PKN edges, MOON scores or selected networks. A later, separately documented relaxed lipid-to-TF sensitivity uses the directionality-audited STITCH overlay with its own conditional mapping ledger.

The limited coverage has **three distinct causes: incomplete implementation of the initially proposed resource collection, incomplete chemical candidate generation, and genuine differences between the structural resolution of measurements and regulatory evidence**. The long LPC precursor routes therefore should not be interpreted as evidence that the databases contain no shorter regulatory hypotheses for our measured lipid panel. Earlier coverage conclusions were too broad when they treated the frozen mapping as an exhaustive chemical search.

## What actually entered the scored PKN

The quoted BRENDA/MRCLinksDB/CellPhoneDB/STITCH description appears under **Proposed network layers** in the [initial feasibility document](cosmos_lipid_pkn_feasibility.md#proposed-network-layers). It describes the intended collection, not the completed regulatory layer. The scored `pkn_lipid_rules` contains 254,589 total edges, but only **77 chemical-to-protein regulatory edges from 49 chemicals**, including **38 lipid-to-protein edges from 23 lipid identities**. All 77 come from GtoPdb. Most of the network's size comes from protein interactions, transcriptional regulation and reaction projections.

| Resource | Frozen material actually available | Status in the scored PKN |
| --- | --- | --- |
| GtoPdb | Current release 2026.2; explicit mouse-target activity annotations | 77 chemical→protein edges, including 38 lipid→protein edges |
| BRENDA | 786 cached records; 772 pass the local export gate | No edges in the scored PKN; 707 surviving edges supported by BRENDA in the separate exploratory graph |
| STITCH | Historical cache: 40,428 records; 40,360 pass the original local export gate | No edges in the scored PKN. The original exploratory graph had 23,552 cache-STITCH-supported edges, but its source orientation cannot be reconstructed. A separate raw directionality-audited overlay exists and has not been used for this measurement-coverage screen. |
| MRCLinksDB | 641 cached records: 429 ligand–receptor and 212 transport | Excluded because the adapter assigned positive signs unconditionally |
| CellPhoneDB | No records imported into the frozen mouse collection | Not implemented; its human nonpeptidic adapter was not projected to mouse |

BRENDA/STITCH edge counts include chemicals beyond lipids. They are not independent experiments; resource-supported edge sets can overlap. The cache's April 2026 timestamp is not a verified BRENDA/STITCH release identifier. Cached rows lack original chemical-mapping alternatives and per-edge references, and STITCH's original acting partner cannot be recovered from this export. These are real provenance limitations, but retaining the resources outside the scored graph also made its regulation coverage narrower. The raw directionality-audited STITCH overlay is a separate, unscored sensitivity resource; this cache-based coverage audit does not estimate its measured-lipid coverage. [Full source-admission audit and upstream implementation checks](cosmos_lipid_regulation_resource_audit.md).

## Denominators and what a match means

The cleaned input file and Factor 1 lipid-weight table contain **1,254 lipid features**. This is the denominator for this audit; the manuscript's earlier 1,386-lipid headline is not the size of this particular input table. The existing mapper accepted **155 inputs**: six ChEBI-mapped lipids and 149 measured-species reaction-hypothesis nodes. A historical gene–lipid correlation matrix selected **64 lipids**; that subset is not the complete set of available Factor 1 weights. The two candidates reaching TFs under the prior normalized precursor policy were **LPC(16:0/0:0) and LPC(18:0/0:0)**.

We keep three mapping relationships separate:

- **Same reported identity:** matching identifiers, explicit alternate identifiers or explicit conjugate acid/base relationships. This does not resolve any structural ambiguity already present in the measured label.
- **Compatible specific identity:** the regulatory chemical could fit the measurement but has more structural detail, such as double-bond positions, stereochemistry or chain arrangement. This is an identity hypothesis, not an accepted measurement assignment.
- **Broad lipid class:** the resource describes a class or ancestor of the feature. It does not establish that each measured member regulates that protein, nor justify copying one edge to every species.

“Direct regulatory match” here means that the matched chemical is the source of a chemical→protein resource record, without traversing reaction nodes. It does **not** guarantee direct physical binding or a single biochemical step.

## Coverage of the inputs actually used

There are **0/155 accepted inputs with an exact outgoing chemical→protein regulatory edge** in the scored PKN. Neither of the two LPC candidates has a species-specific regulatory match in the inspected BRENDA, STITCH, MRCLinksDB or GtoPdb records. CellPhoneDB is unassessed in this collection.

Both LPC candidates have a **class-level LPC→Trpv2** opportunity in GtoPdb. That one generic source term matches 53 measured LPC features, but represents one resource edge, not evidence for 53 independently tested species. Under the previously audited graph, Trpv2 lacks the onward protein→TF connectivity needed for this analysis. Likewise, the compatible DG→Trpc2 candidates below do not establish a usable TF route under that graph. These statements concern the frozen topology and should not be generalized to absence of biological signaling.

BRENDA's existing class/ancestor mappings reach 163/1,254 features (13.0%), but only **four generic chemical terms and seven distinct resource edges**: phosphatidylserine→Dgkq/Smpd3/Smpd5, phosphatidic acid→Acot12, phosphatidylglycerol→Smpd5 and cardiolipin→Smpd3/Smpd5. They cover 22 accepted inputs and two of the historical 64 features. They cannot supply species-specific upstream measurements without a separate class-level model and an explicit aggregation assumption.

## Regulatory identity opportunities beyond the original mapper

The original SwissLipids search in [map_lipids.py](../scripts/pkn/map_lipids.py) searched abbreviation fields. Some important structured records have **empty abbreviations** despite carrying names and ChEBI identifiers. The audit found, for example, empty abbreviations for arachidonate (`SLM:000000296`), oleate (`SLM:000000418`) and DHA anion (`SLM:000001084`). The frozen candidate list therefore missed useful alternatives.

The companion [resolution audit](../results/causal_modelling/lipid_regulatory_resolution_audit/) follows explicit ChEBI subclass relationships beneath the measured composition terms, allowing only explicit alternate IDs and conjugate acid/base equivalence. Simple fatty-acid candidates must additionally have a compatible elemental formula; three oxidized/nitrated FA(18:2) candidates fail that check. This search finds **14 features with specific regulatory-chemical candidates, including 11 with candidates absent from the frozen crosswalk**. It only screens the 128 features already carrying reported-resolution ChEBI identifiers, so it is not an exhaustive remapping of all 1,254 measurements. All 14 are outside the historical 64-lipid correlation subset; their Factor 1 weights nevertheless exist.

| Measured feature | Regulatory chemical and protein | Evidence and implication |
| --- | --- | --- |
| FA(20:4) | Arachidonic acid/arachidonate → Trpv4 | One positive GtoPdb edge already survives in the scored PKN. Its chemical is more specific than the measured 20:4 term. The primary experiment describes indirect activation through CYP-derived EETs, not AA binding directly to TRPV4. |
| FA(22:6) | DHA → Rxra | One positive GtoPdb edge already survives in the scored PKN. DHA is a specific isomer compatible with the measured 22:6 term; the measurement does not establish that identity. |
| FA(20:4) | Arachidonic acid → Acsl4 | One cached BRENDA edge survives in the exploratory graph; identity and original-evidence limitations remain. |
| FA(18:1) | Oleic acid → Lpin1 | Cached BRENDA has both positive and negative signs; neither survives conflict removal. Finding the mapping does not yield a usable signed edge. |
| FA(16:0) | Hexadecanoate/palmitic acid → numerous proteins | The frozen crosswalk already has an acid/base match to the historical STITCH cache, but the input was excluded for multiple ChEBI IDs. Those IDs are nested specificity levels, not unrelated structures; that exclusion merits review. The cache gives 325 records collapsing to 113 signed chemical–protein pairs across 109 proteins before graph conflict handling, but their source orientation is not auditable. |
| DG(18:0_20:4), DG(38:4) | Specific stearoyl–arachidonoyl glycerol (SAG) → Trpc2 | Existing compatible GtoPdb candidates; unresolved chain positions/double-bond identities prevent treating them as exact species assignments. |
| LPS(18:0) | Specific stearoyl lysophosphatidylserine → Gpr34 | Existing compatible candidate, but GtoPdb's secondary chemical-mapping tier is exploratory, not in the scored PKN. |

For AA–TRPV4, see [Watanabe et al., 2003](https://pubmed.ncbi.nlm.nih.gov/12879072/). The [DHA–RXR evidence audit](cosmos_factor1_lipid_tf_evidence.md) separates the GtoPdb mouse-target annotation from original assays, including mouse-brain evidence and human RXRα structural evidence; none establishes a causal effect in this hepatocyte dataset.

Recovering a chemical candidate also does not establish a strong Factor 1 signal.
The saved oriented input weights are approximately +0.000148 for FA(20:4),
−0.07655 for FA(22:6), and +0.004190 for FA(16:0). These are unscaled loadings,
not significance statistics or receptor activities; this coverage audit did not
calculate a new MOON score for these candidate assignments.

Across the frozen and additional candidate searches, GtoPdb has **four measured features with compatible identities for primary-tier regulatory chemicals**: FA(20:4), FA(22:6) and the two DG features. This is 0.32% of the panel and remains conditional coverage, not four accepted regulatory measurements. LPS(18:0) adds one secondary-tier opportunity. The historical cache screen identified compatible specific STITCH candidates for **13 features** (12 fatty acids and PC(36:2), 1.04% of the panel), but its source orientation cannot be reconstructed and it is not directionality-audited coverage. FA(16:0) also has the earlier acid/base match and must not be counted twice. The raw source overlay reports its own identifier and graph audits in [`cosmos_mouse_pkn_stitch_directional`](../results/causal_modelling/cosmos_mouse_pkn_stitch_directional/). It is now paired with a separate long-form relaxed mapping ledger for the exploratory lipid-to-TF sensitivity; that conditional use does not modify this coverage audit, the conservative PKN, or the original lipid-to-TF fit. These small percentages quantify the inspected frozen evidence plus this targeted mapping audit, not the maximum coverage achievable with a fresh full mapping effort.

MRCLinksDB adds potentially relevant associations for fatty acids, LPS and DGs, but its numerical +1 values cannot be used as verified activity signs. CellPhoneDB contributes no assessable coverage estimate until its records, chemical identities, effects and mouse projection have actually been examined. Counts across resources and mapping categories overlap and must not be summed.

The [consolidated per-measurement table](../results/causal_modelling/lipid_regulatory_coverage_overview/all_measured_lipid_regulatory_coverage.csv)
deduplicates these overlaps. Across all inspected tiers, **16/1,254 features
(1.28%) have exact-at-reported-resolution or compatible specific-chemical
opportunities**, while 216 features have broad-class annotations. These are
annotation opportunities, including uncertain cached evidence and conflicting
records, not 16 confirmed signaling lipids. Only one of the 155 accepted inputs
has a compatible specific-chemical opportunity (DG(18:0_20:4)); 28 have class-level
annotations. The two LPC candidates have class annotations only. The table
retains every feature, its original exclusion, Factor 1 weight, resource tier,
mapping category and newly discovered candidate status.

## Consequences for the long routes and a defensible next step

The fitted LPC routes remain valid descriptions of what the supplied graph and accepted input mapping allow. They do not establish that LPC is the only plausible signaling entry point. In particular, **FA(20:4) and FA(22:6) were measured, while their specific AA/DHA mediator identities were not unambiguously assigned or used as observations**. Describing AA/DHA simply as absent from the lipid measurements was too categorical. The acyl-CoA donor pools required by the fitted reaction routes are still unmeasured.

The evidence supports improving candidate generation and reviewing the resource tiers before expanding the network-depth search further. A defensible follow-up would preserve the present fit as a reference, recover specific fatty-acid candidates from ontology/name fields, and compare explicit identity assumptions in a separately labelled sensitivity. Analytical annotation or standards would be preferable to silently equating FA(22:6) with DHA or FA(20:4) with arachidonic acid. BRENDA can be rebuilt from current bulk data with references retained. The raw-source STITCH overlay now verifies its original acting partner, but remains a historical v5 exploratory sensitivity and still needs a dedicated measurement-identity audit. MRCLinksDB and CellPhoneDB need interaction-specific effect and species checks rather than automatic positive signs. Generic class evidence could support a separate class-level question, but should not be duplicated over every measured species.

This audit establishes implementation and mapping limitations. It does not yet establish which newly identified candidates will have useful directed TF routes or coherent MOON scores, and it does not validate the selected LPC network biologically.

## Reproduce and inspect

From the project root, use new empty output directories:

```sh
python3 scripts/moon/ancillary/audit_direct_lipid_regulatory_mapping.py \
  --output results/causal_modelling/lipid_regulatory_mapping_audit_replay
python3 scripts/moon/ancillary/audit_regulatory_resolution_candidates.py \
  --output results/causal_modelling/lipid_regulatory_resolution_audit_replay
python3 scripts/moon/ancillary/summarize_lipid_regulatory_coverage.py \
  --frozen-audit results/causal_modelling/lipid_regulatory_mapping_audit_replay \
  --resolution-audit results/causal_modelling/lipid_regulatory_resolution_audit_replay \
  --output results/causal_modelling/lipid_regulatory_coverage_overview_replay
```

All three scripts use Python's standard library, refuse occupied destinations and save dependency SHA256 hashes. They were executed against the frozen resources; a separate agent checked the ontology interpretation, formula exclusions, compartment-preserving graph lookup and the AA/DHA records. Network/model inputs were not edited.

- [Frozen mapping audit](../results/causal_modelling/lipid_regulatory_mapping_audit/): `measurement_inventory.csv` lists all 1,254 inputs and exclusions; `feature_edge_matches.csv` links features to individual resource records; `coverage_counts.csv` separates the four panels and mapping relationships; `feature_reason_funnel.csv` also includes unmatched features. Palmitate ontology and SwissLipids records are exported separately.
- [Additional specificity audit](../results/causal_modelling/lipid_regulatory_resolution_audit/): `specific_regulatory_candidates.csv` gives the chemical identity hypotheses, exact ChEBI subclass statements, elemental-formula checks, evidence IDs, original-mapper omissions and graph-edge availability. `resource_summary.csv` counts unique candidate features within each evidence tier.
- [Consolidated coverage overview](../results/causal_modelling/lipid_regulatory_coverage_overview/): one row per measured lipid, per-resource counts and unique feature unions across all four input panels; this joins the two detailed audits without changing them.
- [Resource implementation audit](cosmos_lipid_regulation_resource_audit.md): source downloads, exclusions, cache limits and precise build-code pointers.

The source inputs are `results/causal_modelling/cosmos_mouse_pkn/mapping_*.csv`, `regulation_edges.csv`, the explicit ChEBI alias/acid–base tables and both PKN CSVs; the raw mouse cache is `results/causal_modelling/cosmos_pkn_audit/raw/om-mouse-pkn.json.gz`. The resolution screen additionally reads frozen `support/cosmos_mouse/raw/regulation_chebi.obo.gz` and `uniprot_mouse.tsv.gz`. No download is needed to replay these coverage counts.
