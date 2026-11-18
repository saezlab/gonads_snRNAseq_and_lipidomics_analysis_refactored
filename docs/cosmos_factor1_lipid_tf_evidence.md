# Factor 1 lipid-to-TF analysis: DHA–RXR and acyl-remodeling evidence audit

Evidence reviewed on 15 September 2026. This document audits existing prior knowledge relevant to a separately labelled sensitivity analysis using acid/base-normalized chemical identities. It does not add PKN edges, change the exact mapped LPC measurements, report inferred scores or establish a causal effect in hepatocytes.

The DHA–RXR connection has experimental support, but the frozen database record does not provide a mouse hepatocyte assay or a quantitative potency. The accompanying LPC-to-DHA route represents transfer and release of an existing DHA acyl chain. It cannot establish synthesis of DHA from the measured LPC species.

**Measurement-resolution clarification:** the dataset contains FA(22:6), which
is compatible with DHA but does not establish its specific double-bond/isomer
identity. The initial mapper missed this conditional candidate because the
structured SwissLipids record has an empty abbreviation. DHA was therefore not
assigned an observation in the fitted networks; it should not be described as
unconditionally absent from the measurement panel. The DHA-CoA donor remains
unmeasured. See the [mapping and resource audit](cosmos_lipid_regulatory_measurement_coverage.md).

## Frozen PKN provenance

The existing regulation edge is `reg_9244aa3ac6805d503e77` in [regulation_edges.csv](../results/causal_modelling/cosmos_mouse_pkn/regulation_edges.csv). Its chemical source is DHA acid, `CHEBI:28125`; its protein target is mouse Rxra, UniProt `P28700`; and its sign is positive. The source class `curated_mouse_current` reflects the **target-species annotation of the GtoPdb export**. It should not be interpreted as a statement that each attached publication experimentally tested mouse RXRα in liver.

The frozen [GtoPdb interaction export](../support/cosmos_mouse/raw/regulation_gtopdb_interactions.csv.gz), version 2026.2, published 15 June 2026, specifies:

| Field | Frozen value |
| --- | --- |
| Ligand | Docosahexaenoic acid; ligand ID 1051 |
| Target | Retinoid X receptor-α; target ID 610 |
| Target species and protein | Mouse; Rxra; P28700 |
| Type and action | Agonist; Agonist |
| Endogenous | True |
| Selectivity | Not Determined |
| Concentration, potency and assay description | No value supplied |
| Attached PubMed references | 10529898; 11690641 |

The [official GtoPdb ligand page](https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=1051&tab=biology) also lists mouse RXRα agonism without a potency value. The live website required login when directly opened during this audit; the frozen export remains the reproducible evidence source, with the indexed official page providing corroboration.

**Record identifier caveat:** `GtoPdb:interaction_record:20691` is the PKN builder's one-based enumeration of the parsed frozen interaction export. [build_lipid_regulation.py](../scripts/pkn/build_lipid_regulation.py) constructs it with `enumerate(interactions, start=1)`. It is not an official persistent GtoPdb interaction accession and should not be used to construct a live interaction URL. The source-file checksum, export version, ligand/target identifiers and row fields together identify the underlying record. The positive sign is assigned from the explicit agonist action, not from binding alone.

## What the attached references actually establish

**PMID 10529898 — Desvergne and Wahli, 1999.** *Peroxisome proliferator-activated receptors: nuclear control of metabolism* is a review in *Endocrine Reviews*. It is not an original DHA–RXR activation experiment. Its presence in the database reference field must not be presented as an independent mouse RXRα assay. [Publisher record](https://academic.oup.com/edrv/article/20/5/649/2530847); [PubMed](https://pubmed.ncbi.nlm.nih.gov/10529898/).

**PMID 11690641 — Lampen, Meyer and Nau, 2001.** *Phytanic acid and docosahexaenoic acid increase the metabolism of all-trans-retinoic acid and CYP26 gene expression in intestinal cells* examined retinoid responses in **human Caco-2 intestinal cells**. DHA alone did not alter CYP26 expression or retinoic-acid metabolism. Combining DHA with all-trans-retinoic acid or the synthetic RAR ligand Am580 increased CYP26 induction relative to the RAR ligand alone. This supports context-dependent retinoid signaling; it is not direct evidence for an isolated mouse RXRα response or hepatic regulation. The paper also describes phytanic-acid PPAR transactivation experiments, which should not be misattributed to DHA–RXR. Exact DHA treatment concentrations were not verified from the accessible primary text and are deliberately not supplied here. [Primary paper](https://www.sciencedirect.com/science/article/pii/S0167478101003050); [PubMed](https://pubmed.ncbi.nlm.nih.gov/11690641/); DOI: [10.1016/S0167-4781(01)00305-0](https://doi.org/10.1016/S0167-4781(01)00305-0).

## Additional primary evidence for the existing DHA–RXR hypothesis

These papers provide interpretive context for an existing edge. They were **not added to the frozen PKN**, and should not be described as its original attached provenance.

**Mata de Urquiza et al., 2000, PMID 11118147.** *Docosahexaenoic acid, a ligand for the retinoid X receptor in mouse brain* identified DHA by purification and mass spectrometry of an RXR-activating factor from **adult mouse brain**. Activation was assessed in cell-based assays. This is the foundational mouse-brain study; its PMID is absent from the edge's two attached references. The accessible abstract does not resolve the exact isoform-specific assay conditions or concentration series, and the publisher's full text was inaccessible during this audit. No exact assay dose or mouse RXRα potency is inferred from the abstract. [Primary-study record](https://pubmed.ncbi.nlm.nih.gov/11118147/); DOI: [10.1126/science.290.5499.2140](https://doi.org/10.1126/science.290.5499.2140).

**Egea, Mitschler and Moras, 2002.** *Molecular recognition of agonist ligands by RXRs* directly resolved DHA bound to recombinant **human RXRα ligand-binding domain**, expressed in *Escherichia coli*. The protein was incubated with ninefold ligand excess and a GRIP1 coactivator peptide; the crystal complex adopted an active conformation. Table 2 reports an approximate DHA dissociation constant of **50–100 µM**, not nanomolar. This is a value listed in the structural paper, not a mouse hepatocyte EC50 or a potency supplied by the frozen GtoPdb record. The experiment supports physical binding and an agonist-compatible conformation under its preparation conditions. It does not establish receptor occupancy, heterodimer-dependent transcription or an effective concentration in this liver dataset. [Primary paper, including Table 2 and Methods](https://academic.oup.com/mend/article/16/5/987/2741678); DOI: [10.1210/mend.16.5.0823](https://doi.org/10.1210/mend.16.5.0823).

## Reaction equations and net cancellation

The following equations, chemical identities and mouse catalyst annotations were checked in the frozen [rhea_reactions.jsonl](../results/causal_modelling/cosmos_mouse_pkn/rhea_reactions.jsonl). The reaction accessions link to Rhea; the frozen records are the evidence inspected for this audit.

| Reaction | Directed equation | Mouse catalysts annotated in the frozen record |
| --- | --- | --- |
| [RHEA:37476](https://www.rhea-db.org/rhea/37476) | LPC(16:0) + DHA-CoA → PC(16:0/22:6) + CoA | Lpcat1, Lpcat4 |
| [RHEA:41232](https://www.rhea-db.org/rhea/41232) | PC(16:0/22:6) + H₂O → LPC(16:0) + DHA⁻ + H⁺ | Pla2g2a, Pla2g2d, Pla2g2e |

Here, the Rhea LPC is specifically 1-hexadecanoyl-sn-glycero-3-phosphocholine (`CHEBI:72998`). The PC is specifically 1-hexadecanoyl-2-docosahexaenoyl-sn-glycero-3-phosphocholine (`CHEBI:74963`), with the all-cis DHA double-bond positions specified by Rhea. DHA-CoA is `CHEBI:74298`, CoA is `CHEBI:57287`, and released docosahexaenoate is `CHEBI:77016`. These are precise reaction identities, not evidence that the measurements resolve every positional detail. Catalyst alternatives do not imply a required multi-enzyme complex.

Summing the two equations cancels the LPC and PC intermediates:

```text
DHA-CoA + H₂O → DHA⁻ + CoA + H⁺
```

The DHA chain therefore originates in a pre-existing **DHA-CoA pool**, which is unmeasured in this proposed route. LPC(16:0) is regenerated. This is an acyl-transfer/release cycle with net thioester hydrolysis, not DHA synthesis from the LPC palmitoyl chain. Substrate availability could conditionally influence turnover, but the graph alone does not establish rate limitation, the sign of a steady-state abundance response or net flux.

## Chemical normalization and interpretation limits

The frozen [ChEBI acid/base crosswalk](../results/causal_modelling/cosmos_mouse_pkn/regulation_chebi_acid_base_relations.csv) explicitly relates DHA acid (`CHEBI:28125`, charge 0) to docosahexaenoate (`CHEBI:77016`, charge −1) as conjugate acid/base forms, using ChEBI version 255. The regulatory and reaction resources use these different identifiers. A separately labelled sensitivity analysis may use this recorded relationship to normalize identifiers while retaining the existing exact mapped LPC inputs. It is an identity-normalization assumption, not a newly established causal interaction or a measured change in protonation state. Such normalization must not imply equal membrane permeability, nuclear access or assay exposure of the two forms.

The defensible hypothesis is **LPC-associated acyl remodeling linked to release of an RXR-active fatty acid, conditional on availability of a DHA-containing precursor pool and the relevant enzymatic and receptor context**. It remains uncertain whether this route operates in these hepatocytes or explains their Factor 1 association. Factor weights do not provide absolute ligand concentrations, receptor occupancy, nuclear availability, heterodimer context or reaction flux. Changes in bulk lipid abundance do not directly measure the free intracellular DHA available to RXR.

Accordingly, any network using this route should distinguish a measured LPC entry node, unmeasured metabolic intermediates and a literature-supported receptor interaction. It should not claim that LPC is a direct RXR agonist, that LPC is converted into a newly synthesized DHA chain, or that DHA caused the observed hepatocyte transcriptional changes.
