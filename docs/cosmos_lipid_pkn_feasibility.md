# Lipid-aware COSMOS PKN: feasibility and source audit

Assessed 2026-09-08. This is a completed feasibility audit and a candidate protein
backbone, **not a completed lipid PKN or a COSMOS analysis**.

## Recommendation

Build a mouse-specific PKN using current OmniPath signaling/CollecTRI, an
explicit Mouse-GEM metabolic backbone, SwissLipids/Rhea lipid reactions, and
separately curated lipid-to-protein regulation. Use lipid reaction rules to
extend species coverage only in a clearly identified sensitivity network.
Do not use the live prebuilt COSMOS network unchanged.

This is feasible at the resource and software level. Its ability to explain
the particular measured lipids and selected Cyp genes is still unestablished.
More edges do not establish molecular specificity or causal identifiability.

## What was checked

- Local README, agent guide, workflow dependency overview, saved RNA–lipid
  networks, lipid DE table, LINEX2 adjacency, and the existing meta-network.
- OmniPath-metabo main at `0d7bb56ac7c8b6e5ad0f19f54874ef84d4e7f84c`, including
  the requested AGENTS.md, actual builders, formatter, and resource processors.
- COSMOSR master snapshot `b2ed195` (short commit ID), AGENTS.md and relevant
  source. COSMOSR supports an external signed PKN; the packaged PKN is not
  mandatory. Current package guidance defaults to MOON when no method is chosen.
- Live mouse COSMOS service, fresh mouse OmniPath PPI and CollecTRI exports,
  SwissLipids download metadata and enzyme table, and exact mouse UniProt
  gene-name queries for the seven validation genes.

Downloads and SHA-256 hashes are retained in
`results/causal_modelling/cosmos_pkn_audit/manifest.json`; compressed source responses are in
`results/causal_modelling/cosmos_pkn_audit/raw/`. These freeze the observed responses, not all
underlying database releases. The service does not supply a full build ledger.

## Important implementation findings

**The development plan is behind the implementation.** The current builder
contains GEM, KEGG, iMM1415, Recon3D, OmniPath and CellPhoneDB processors in
addition to the originally proposed resources. This is not evidence that every
processor or species configuration has been validated. In particular, explicitly
set `gem={'gem': 'Mouse-GEM'}`: the GEM processor defaults to Human-GEM and says
its `organism` argument does not filter a model. Do not assume
`build(organism=10090)` alone safely selects every resource.

**The live cache has unusable protein-layer signs.** The downloaded mouse network
contains 235,302 raw records, with these resource counts:

| Resource | Records | Observed sign state |
|---|---:|---|
| Mouse-GEM metabolism | 115,434 | All +1 |
| Mouse-GEM transport | 27,282 | All +1 |
| STITCH | 40,428 | Positive and negative |
| CollecTRI | 28,414 | All 0 |
| OmniPath signaling | 22,034 | All 0 |
| BRENDA | 786 | Positive and negative |
| MRCLinksDB | 641 | All +1 |
| TCDB | 283 | All +1 |

These are records, not unique mechanistic edges or lipid-specific interactions.
The mouse cache components have April 21, 2026 modification dates. Keeping only
±1 would discard the entire signaling and GRN layers. Fresh direct OmniPath
exports have signs: 39,849 GRN records and 32,468 PPI records. The current source
parses `True`/`False` signs correctly, whereas its `directed` attribute still
compares against `'1'`. A historical boolean-parsing problem is a plausible
explanation for the cache, but its exact build history has not been verified.

The audit produces a **53,562-edge candidate protein backbone** from fresh
exports. It uses explicit direction and exactly one supported sign, excludes
unsigned pairs, and quarantines 6,194 pairs with conflicting signs across either
layer. This conservative conflict policy is an analysis choice. Raw references
and layer provenance remain in a separate evidence table; an activity-versus-
expression distinction is still needed when integrating the layers.

**Rhea is not a ready-to-enable default resource.** A `rhea_reactions()` helper
exists, but is absent from the main `PROCESSORS` registry. It requires `lipinet`
or supplied parsed nodes and returns `source,target,weight`, unlike the current
Interaction interface. Its equation parser truncates both sides to the smaller
number of participants before pairing them. That can discard reactants/products
and cannot preserve reaction stoichiometry. Do not use it unchanged to create
the lipid reaction layer.

**Identifier formatting needs an adapter.** The service returns raw ChEBI and
UniProt entities, with some ChEBI values prefixed and others numeric strings.
`format_pkn()` creates compartmental metabolite nodes, reaction/enzyme instances,
and protein-to-instance connectors. Format the assembled network once so the
connector filter sees its protein layers. Export only
`source,interaction,target` for R, retaining all provenance elsewhere. Check
reverse-node handling (`_rev` versus some R helpers' `_reverse`) and use an
explicit node table for expression filtering. HMDB display helpers must not be
used to translate ChEBI nodes.

## Lipid coverage in this repository

The current `gonade_type_ttop.csv` contains **1,254 unique feature IDs**, rather
than the 1,386 headline count in older project descriptions. This is the actual
audited result-table universe; it does not establish why other features were
excluded. Some IDs represent non-lipid metabolites.

The saved LINEX2 adjacency contains 1,129 nodes and 7,552 nonzero matrix entries
(3,776 undirected pairs). Of the 1,254 DE features, 910 have exact node-label
matches and 853 have at least one connection. The Factor 1 correlation matrix
has 139 genes and 64 lipid columns; only 22 of those lipid labels are connected
in this LINEX2 matrix. Name reconciliation could change these overlaps.

These counts measure **label overlap and adjacency**, not coverage by signed,
enzyme-supported reactions. The adjacency is symmetric and lacks enzymes,
reaction IDs, substrate/product roles and evidence. Its values cannot be
reinterpreted as activation signs. The class-conversion support tables also
lack the detail needed for a mechanistic signed network.

No exact ChEBI/SwissLipids coverage percentage is claimed: the feature-to-entity
mapping has not yet been built. The generated inventory deliberately marks
every mapping as `not_yet_mapped`.

## Proposed network layers

**Implementation status, audited 15 September 2026:** this table records the
original proposal. The completed conservative PKN uses GtoPdb for chemical→protein
regulation; cached BRENDA/STITCH were placed in a separate exploratory variant,
MRCLinksDB was excluded because its signs were assumed, and CellPhoneDB was not
integrated. See the [resource audit](cosmos_lipid_regulation_resource_audit.md) and
[measured-lipid coverage audit](cosmos_lipid_regulatory_measurement_coverage.md),
which also documents specific fatty-acid candidates missed by the initial mapper.

| Layer | Preferred resources | Admission rule |
|---|---|---|
| Protein signaling and TF regulation | Fresh OmniPath and CollecTRI mouse exports | Directed, signed evidence; preserve layer, references and orthology provenance |
| Metabolic and transport context | Explicit Mouse-GEM; curated transporter resources | Preserve reaction IDs, compartments, reversibility and gene–protein–reaction rules |
| Curated lipid reactions | SwissLipids enzyme/Rhea links, Rhea participants and UniProt catalysts | Match chemistry at its supported structural resolution and taxon |
| Lipid-specific reaction expansion | Reconstructed LINEX2/Rhea/Reactome reaction rules | Separate rule-derived edges; require transformation, enzyme evidence and specificity constraints |
| Lipid → protein regulation | BRENDA, MRCLinksDB, CellPhoneDB non-peptidic interactions; signed STITCH as an additional tier | Explicit activation/inhibition or documented ligand effect; exclude unsigned binding from signed core |
| Lipid identity/hierarchy | SwissLipids, ChEBI, LIPID MAPS; Goslin normalization | Annotation/mapping relations, not automatic signed causal edges |

SwissLipids' live download metadata reports September 8, 2026 for the inspected
enzyme table. It contains 5,714 rows, including 927 mouse rows, with UniProt,
gene, Rhea, reaction text and evidence tags. These are enzyme/reaction annotation
rows, not 927 independently established mouse lipid species interactions. The
`SwissLipids ID` in this table is an SLP protein entry, not an SLM lipid ID.
Its full lipid and lipid-to-UniProt downloads offer hierarchy and cross-mapping
fields and must be joined at their stated mapping level.

SwissLipids is a particularly useful extension because it joins structural
hierarchy to curated reactions and enzymes. A hierarchy-projected enzyme
association is not automatically a direct regulatory lipid→protein edge.
Similarly, a STITCH edge involving a lipid-binding enzyme does not by itself
identify substrate turnover or a particular lipid species' effect on activity.

## Representation and scientific safeguards

1. Preserve each original feature ID and reported identification level. Parse
   class, chains, sum composition, ether/plasmalogen state and oxygenation.
   Keep `_` (unordered chains) distinct from `/` (position-specified notation)
   unless the measurement documentation establishes otherwise. Existing naming
   helpers can replace underscores with slashes, so punctuation alone is not
   proof of experimental positional resolution. Do not resolve HexCer to a
   specific sugar isomer without evidence.
2. Store exact matches, compatible ambiguous matches, and class/ancestor matches
   separately. Never copy a generic class edge to every measured species and
   count those as curated molecular interactions. Avoid duplicating one
   measurement over alternative structures as independent observations.
3. Preserve a reaction-level master graph, including all participants, catalyst
   rules and compartments. Derive the signed COSMOS projection from it.
   An enzyme complex requires its component logic; flattening components to
   interchangeable enzymes loses that requirement.
4. Declare metabolic edge semantics. The inspected GEM adapter assigns +1 to
   substrate→enzyme and enzyme→product edges. These encode a qualitative
   reaction traversal, not experimental substrate activation or a prediction
   that greater enzyme activity raises substrate abundance. A consumption-based
   projection would instead need explicit negative effects on substrates, with
   its own assumptions and sensitivity analysis. Neither projection is a flux
   balance or kinetic model.
5. Retain true compartments in the prior. Whole-liver lipid abundance does not
   identify compartment-specific concentrations, and hepatocyte RNA does not
   make bulk lipid changes exclusively hepatocyte-derived. Pooling or expanding
   observations over compartments must remain an explicit later choice.
6. Keep the measured correlation networks outside the PKN. They provide
   hypothesis priorities and comparisons; adding them as causal edges would
   make their mechanistic explanation circular. RNA abundance is not enzyme
   activity. TF footprints and RNA coherence are possible later inputs, after
   the biological direction and mapping are agreed.

## Consequences for the study's candidate genes

An exact mouse UniProt gene-name audit, including reviewed and unreviewed hits,
finds cached incident records for Cyp2b9 (STITCH), Cyp2c38 (GEM and protein
layers), Sult3a1 (GEM), and Abcb1a (several resources). Cyp2b13, Cyp2d9 and Mtmr11
have no matching incident records in this snapshot. Only Cyp2c38 and Abcb1a
occur as exact gene symbols in the fresh protein exports examined. This is an
annotation-coverage result, not evidence of biological inactivity; accession
history, alternative mappings and reaction curation can change it.

Consequently, a broader network may suggest upstream pathways more readily
than explain every validated gene's substrate specificity. Do not transfer
human CYP reactions to expanded mouse paralogs indiscriminately. Any targeted
literature-curated additions should have explicit evidence and remain separate
from the unbiased resource-only benchmark.

## Concrete build and acceptance plan

1. Freeze all source releases/commits and parameters. Rebuild mouse resources
   explicitly; retain the April cache only as a benchmark. Check each resource's
   organism behavior rather than enabling all defaults.
2. Build and review the lipid identity crosswalk, reporting exact, ambiguous,
   class-only and unmapped coverage both globally and for the 64 Factor 1 lipids.
3. Implement a lossless Rhea/SwissLipids reaction adapter and recover enzyme-
   annotated LINEX2 reaction records. The adjacency alone is insufficient.
4. Assemble separate networks: current packaged PKN benchmark; refreshed core;
   curated lipid extension; additional rule-expanded lipid extension. Keep
   STITCH inclusion and evidence thresholds as sensitivity settings.
5. Test ID uniqueness, species, signs/conflicts, all reaction participants,
   compartment transitions, complex rules, reverse-instance identity and
   protein-to-reaction connectors. Include a one-substrate/two-product fixture
   to expose the existing Rhea parser's truncation failure.
6. Report measured-feature and candidate-gene coverage, directed connectivity,
   and paths from appropriate RNA-derived regulators to lipid observations.
   Compare edge tiers and path stability. Connectivity alone is insufficient
   if paths are dominated by ubiquitous metabolites or generic lipid classes.

The next implementation milestone is a **mapped, provenance-preserving lipid
PKN with coverage comparisons**, before selecting any MOON inputs or scoring.
No direction-dependent MOON run or final input mapping was performed here.

## Reproduction and outputs

The historical probe is preserved at `miscellaneous/quick_tests/code/audit_cosmos_pkn.py`. It retains its original path assumptions and is not a current reproduction command; use `docs/CAUSAL_ANALYSIS_MAP.md` for active entry points.
It uses the frozen compressed downloads and current local result files, verifies
download checksums and validates output sign/uniqueness invariants. It was run
successfully after implementation. No existing R workflow was modified or rerun.

Outputs in `results/causal_modelling/cosmos_pkn_audit/`: `summary.json`, resource sign audit,
lipid inventory, candidate protein backbone, quarantined sign conflicts,
protein evidence, mouse SwissLipids enzyme records, and candidate-gene audit.
The backbone is deliberately labelled a component, not a complete COSMOS PKN.

## Primary sources

- [Requested OmniPath-metabo plan](https://github.com/saezlab/omnipath-metabo/blob/main/AGENTS.md)
- [Pinned builder](https://github.com/saezlab/omnipath-metabo/blob/0d7bb56ac7c8b6e5ad0f19f54874ef84d4e7f84c/omnipath_metabo/datasets/cosmos/_build.py)
- [Pinned resource implementations](https://github.com/saezlab/omnipath-metabo/tree/0d7bb56ac7c8b6e5ad0f19f54874ef84d4e7f84c/omnipath_metabo/datasets/cosmos/resources)
- [COSMOSR agent guidance](https://github.com/saezlab/cosmosR/blob/b2ed195/AGENTS.md)
- [Live COSMOS cache status](https://metabo.omnipathdb.org/cosmos/status)
- [SwissLipids downloadable resources](https://www.swisslipids.org/api/index.php/downloadData)
- [SwissLipids hierarchy and reaction query examples](https://beta.sparql.swisslipids.org/)
- [Rhea reaction knowledgebase](https://www.rhea-db.org/)
- [LIPID MAPS reaction explorer](https://www.lipidmaps.org/resources/tools/reactions)
