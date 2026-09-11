# Mouse COSMOS PKN build

Build date: 2026-09-08. This implements the first network build proposed in
`docs/cosmos_lipid_pkn_feasibility.md`. It constructs priors and mapping
candidates; no omics scores, direction-dependent MOON analysis or final input
vectors have been selected.

## Deliverables and intended use

Outputs are in `results/cosmos_mouse_pkn/`. Each `pkn_*.csv` has exactly
`source,interaction,target`, with numeric signs ±1. Matching `.rds` files contain
the same three-column R data frames. `r_validation.json` records validation
against the installed `cosmosR::meta_network_cleanup()` implementation.

| Export | Interpretation |
|---|---|
| `pkn_core` | Fresh mouse OmniPath signaling/CollecTRI plus current Mouse-GEM, with supported flux direction and single-gene catalyst alternatives |
| `pkn_lipid_curated` | Core plus direction-annotated mouse Rhea/SwissLipids chemistry and fresh, structurally screened GtoPdb metabolite–protein activity records |
| `pkn_lipid_curated_ph73` | Optional chemical-identity sensitivity view of the curated network; ChEBI aliases and explicit acid/base relations support normalization to a Rhea pH7.3 target |
| `pkn_lipid_rules` | Curated network plus explicitly hypothetical species-level reactions instantiated from mouse-supported templates and matched measured acyl chains |
| `pkn_exploratory` | Adds species rules, uncertain-direction Rhea possibilities, enzyme-complex enabling approximations, secondary chemical mappings, unverified legacy accessions and cached BRENDA/STITCH evidence |

Use `pkn_lipid_curated` for the initial resource-based analysis and keep the
exploratory graph separate. All graphs are qualitative models, including the
curated one. The component statistics and exclusions are in
`build_summary.json`; final graph and mapping coverage is in
`coverage_summary.json`. Counts of graph edges are not counts of independent
experiments or biochemical reactions.

The completed build contains:

| Network | Signed edges | Nodes |
|---|---:|---:|
| Updated core | 231,386 | 46,114 |
| Curated lipid expansion | 251,998 | 52,157 |
| Curated expansion plus species rules | 254,589 | 52,755 |
| Broad exploratory sensitivity network | 303,925 | 70,049 |

The curated expansion has exact ChEBI matches for only **7 of 1,254 measured
features and none of the 64 Factor 1 lipids**. It has compatible structural
candidates for 181 features and hierarchy matches for 750; these do not count
as exact coverage. The separate rule extension connects 149 measured features,
including seven Factor 1 lipids. Thus the first build is useful for testing
mechanistic hypotheses, but does not solve comprehensive molecular-species
coverage.

One explicit rule hypothesis is
`PC(20:0_22:6) → Ptdss1-associated reaction → PS(20:0_22:6)`, from mouse catalyst
Q99LH2 and directional template RHEA:45089. This is chain-compatible extension
of the curated class reaction, not direct evidence that Ptdss1 processes those
measured species. The complete reaction also retains the non-lipid participants.

## Sources and versions

- Mouse-GEM snapshot `fc769ac1da421c4e4e9e1cf2a9aec3fe4f02ef6a`, fetched directly
  from the model repository. The SBML contains 12,987 reactions and 8,454
  compartment-specific metabolite records. GPR trees and flux bounds are retained.
- Mouse OmniPath/CollecTRI responses fetched September 8, 2026, using the direct
  service rather than the unsigned April COSMOS protein cache. Original
  UniProt accessions, references and dataset labels are preserved.
- Rhea release 142, September 2, 2026, RDF plus reviewed UniProt crosswalk;
  SwissLipids mouse enzyme annotations and September 8 lipid hierarchy snapshot.
- GtoPdb 2026.2 native mouse metabolite actions. Direct signed actions are
  retained; source mappings are screened against chemical structure where
  identifiers permit it. Secondary PubChem mappings remain exploratory.
- UniProt mouse accession/primary-gene export, September 8, 2026. Multiple
  accessions for one primary gene collapse to a gene-centric node. Isoform
  suffixes map to the canonical entry's gene where available; isoform-specific
  functional differences are therefore not modeled. Multi-locus protein
  entries remain protein nodes. Unverified accessions are exploratory only.
- ChEBI release 255, September 1, 2026, for explicit conjugate acid/base
  relationships, alternative IDs and lipid ancestry. The Rhea-synchronized
  ChEBI hierarchy is also preserved separately.
- Cached BRENDA/STITCH from the audited April 21 COSMOS service components are
  included only in `pkn_exploratory`. These have weaker reconstructable
  provenance than fresh resources and are not represented as newly refreshed.

`build_manifest.json` hashes source downloads, local feature inputs, build code
and outputs. The individual resource manifests contain URLs, queries and
versions. Frozen raw sources are in `support/cosmos_mouse/raw/`; the direct
OmniPath and cached fallback responses also depend on
`results/cosmos_pkn_audit/raw/`. Preserve these when moving the build.

## Reaction representation

The master reaction JSONL files preserve all substrates/products, stoichiometry,
reaction references, compartments and catalyst evidence. Mouse-GEM also retains
the complete GPR tree. Rhea preserves symbolic coefficients and generic compound
identities rather than silently discarding them.

The scoring projection represents a reaction as:

```text
mouse gene product ──(+1 enables)───────────┐
                                          v
lipid substrate ──(+1 participation)──> enzyme/reaction instance ──(+1 formation)──> lipid product
```

All non-lipid participants remain in the master and projected reactions too.
Each reaction/catalyst/direction combination has its own stable instance. The
protein-to-instance connector lets signaling connect to metabolism while
preserving reaction identity. An enzyme's transcript is not assumed to measure
the activity of any of its instances.

These positive metabolic edges describe **qualitative reaction traversal**.
They do not establish substrate activation of a protein, imply that more enzyme
raises substrate concentration, or implement steady-state flux/kinetic balance.
Regulatory activation/inhibition is a different edge mechanism in the evidence
table. Direct lipid-to-transcript regulatory evidence was not established by
this build; lipid-to-gene-product edges describe protein activity.

The conservative graphs omit:

- GPR AND complexes from their signed projection; their logic remains intact
  in the master. The exploratory graph uses positive member-to-instance
  connectors, an enabling approximation that cannot enforce an AND condition.
- Reactions without annotated catalysts, exchange/sink reactions without both
  sides, and blocked Mouse-GEM directions.
- Rhea UN annotations when direction is unspecified. The exploratory graph
  tests both orientations without claiming physiological reversibility. A
  directional annotation supersedes an UN annotation for the same reaction
  master and catalyst, avoiding redundant hypothetical reverse paths.
- Rhea transport reactions whose relative compartments cannot be assigned.
  Mouse-GEM transport retains its actual compartments.
- Generic Rhea protein/polymer participants and symbolic coefficients from the
  conservative projection; their master records and exploratory edges remain.

Rhea LR/RL means an annotated reaction orientation, not proof of net flux in
these livers. Independently sourced overlapping reactions are retained with
provenance and may share chemistry; path or enrichment analysis must account for
this redundancy rather than treating every instance as independent evidence.

## Chemical identity and lipid resolution

Metabolite IDs look like `Metab__CHEBI:12345_c`. The `unknown` suffix is an
explicit absence of compartment information and is never silently equated with
cytoplasm. Compartment-free observations are not duplicated into model
compartments during this build. Non-ChEBI model entities retain their namespace.

Gene products use `Gene__MouseSymbol`; multi-locus/unresolved entries retain
`UniProt__accession`. Instances use `Enzyme__stable_hash`.
Use `nodes.jsonl` for type, compartment, gene and catalyst information.
Legacy COSMOS helpers that parse numeric Gene prefixes, HMDB names or allowed
compartment suffixes are not compatible adapters for this naming scheme.
`moon()` accepts the signed network contract, but later filtering and input
construction must use the explicit node table. Use direct named vectors after
the biological mapping is agreed, not HMDB-specific convenience helpers.

The optional pH7.3 graph preserves original identities in the master. A merge
requires a ChEBI alt-ID relation, or a Rhea target supported by an explicit
ChEBI conjugate component (or curated Rhea mapping). The acid/base relationship
is curated; selection of the major species at pH7.3 is usually computational.
`chemical_node_aliases.csv` records both. Compartments remain separate and all
eligible raw signs are reconciled again after normalization, including evidence
that was already quarantined in the unnormalized graph.

`mapping_features.csv` retains all 1,254 feature IDs from the current lipid DE
table. `mapping_candidates.csv` distinguishes exact reported-resolution names,
compatible structural candidates, and hierarchy/ancestor annotations.
`mapping_chebi_nodes.csv` also retains broader Rhea protonation candidates with
their origin; these are not all admitted into the conservative normalized graph.

Exact reported resolution is not necessarily a fully resolved molecule.
Underscore/slash, oxygenation, ether/plasmalogen and sugar-identity distinctions
are preserved. Goslin's attempted reinterpretation of 26 two-moiety CL names as
dilysocardiolipin is rejected; these keep class annotations pending clarification
of the original identification convention.

The mapping has 877 exact SwissLipids reported-resolution matches, but only 128
features with an exact-resolution ChEBI cross-reference. For the 64 Factor 1
lipids these counts are 41 and 3. A cross-reference alone does not mean that a
reaction involving that entity exists in the PKN. See `feature_coverage.csv`
for actual network membership. Exact, compatible and hierarchy coverage counts
overlap; they are not mutually exclusive categories and are not additive.
Hierarchy matches are never copied into signed species-level edges.

## Measured-species reaction hypotheses

The separate `lipid_rules_reactions.jsonl` extends eight directional Rhea
templates after checking current LINEX2 lipid definitions at commit
`04e0099655a6b7a8cb57aaeca7b1301b3524079d`. It yields 149 instantiated reactions
covering 149 measured features (seven Factor 1 lipids), associated with 33 mouse
catalysts. The transformations are DG↔PA, PC→DG/PA/PS, DG→PC/PE and PI→DG.

Only matching, fully reported diacyl-chain compositions are paired. Sum-only
lipids, ethers, oxidation variants, sphingolipids and coupled multi-lipid
transformations are excluded from this first rule implementation. All other
participants and stoichiometries from each source reaction are preserved.
Contradictory reported sn ordering is rejected; unreported positions remain
unresolved. Identical acyl compositions do not establish that the enzyme accepts
the species, so the entire extension is a **substrate-specificity hypothesis**.

Rule-derived measured entities use `MeasuredLipid:original_feature_name` inside
the metabolite namespace and remain distinct from curated ChEBI species.
The rule pairs table, master records and manifest retain exact matches,
template participants, references and assumptions. Actual graph coverage is
reported separately as `species_rule_connected`, never counted as curated exact
coverage. These rules use feature identities, not correlation or abundance
signs, so they do not encode the observed RNA–lipid association as prior truth.

## Quality controls and limitations

The build excludes unsigned/undirected protein interactions, self edges,
opposing signs for a source–target pair and unverified mouse catalyst instances
from conservative graphs. Conflict reconciliation happens after gene identity
collapse. Every retained enzyme instance must retain its declared catalyst
connectors. Evidence tables retain omitted conflicts and their references.

The source audit quarantines MRCLinksDB's assumed +1 associations and chemical
cross-reference discrepancies in current GtoPdb data. These include an
incorrect PGE2 structure and stereochemical discrepancies for other ligands.
Affected records, structural keys and candidate corrections are preserved in
`regulation_mapping_review.csv`; no guessed replacements enter the graph.

The reaction parser's tests cover unequal side counts, repeated RDF subjects,
symbolic/numeric coefficients and direction variants. Projection tests cover
GPR AND/OR logic, bounds, uncertain direction, unverified catalyst downgrade,
gene-level sign conflicts and compartment separation. R validation requires
that the installed COSMOS cleaner leave every exported edge unchanged.

This is a biochemical prior, not a liver-specific causal result. Bulk lipid
abundance cannot be assigned exclusively to hepatocytes or liver compartments
from these resources. Mouse-native sources can still contain annotations
inferred from other species; provenance does not make them mouse-liver
experiments. No fitted correlation signs are used to create prior edges.

Candidate-gene coverage remains uneven. `focus_gene_coverage.csv` and
`focus_reaction_examples.csv` give reproducible lookups for the seven study
genes. Missing CYP paralogs or Mtmr11 cannot be rescued by indiscriminate human
orthology transfer. Candidate-specific literature curation would be a separate,
explicit evidence addition.

## Rebuild

The builders use frozen inputs by default. Source refresh options are explicit
and may change scientific results; archive a new manifest when refreshing.

```bash
python3 scripts/pkn/build_mouse_gem.py
python3 scripts/pkn/build_rhea_lipid_reactions.py
python3 scripts/pkn/build_lipid_regulation.py
# Use the pinned mapping environment or install the packages listed below.
uv run --with pygoslin==2.2.5 --with scipy==1.13.1 --with numpy==2.0.2 python scripts/pkn/map_lipids.py
python3 scripts/pkn/build_lipid_rules.py
python3 scripts/pkn/assemble_mouse_pkn.py
python3 scripts/pkn/audit_mouse_pkn.py
Rscript scripts/pkn/validate_export.R
python3 scripts/pkn/write_build_manifest.py
```

Mapping was executed with the isolated `/tmp/fcgl-pkn-goslin` environment.
Other current builders use Python's standard library. R requires `cosmosR` and
`jsonlite`. The project R analysis workflows were not modified or rerun.

## Primary resources

- [Mouse-GEM pinned model](https://github.com/SysBioChalmers/Mouse-GEM/tree/fc769ac1da421c4e4e9e1cf2a9aec3fe4f02ef6a)
- [OmniPath-metabo](https://github.com/saezlab/omnipath-metabo)
- [OmniPath interactions service](https://omnipathdb.org/interactions)
- [Rhea downloads](https://www.rhea-db.org/help/download)
- [SwissLipids downloads](https://www.swisslipids.org/api/index.php/downloadData)
- [Guide to Pharmacology downloads](https://www.guidetopharmacology.org/download.jsp)
- [ChEBI](https://www.ebi.ac.uk/chebi/)
