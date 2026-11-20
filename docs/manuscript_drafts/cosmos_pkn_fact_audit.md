# PKN facts verified for the manuscript draft

Audit date: 2026-09-09. These are drafting notes, not replacement results. All
paths below are relative to the project root. The historical manifests retain
the original `results/cosmos_*` paths; their current location is
`results/causal_modelling/cosmos_*`. SHA-256 values for `pkn_lipid_rules.csv`,
`pkn_lipid_rules.rds`, `pkn_lipid_curated.csv`, `build_summary.json`,
`edge_evidence.jsonl.gz` and `nodes.jsonl` still match the original build manifest.

## Scale and defensible terminology

| Quantity | Verified value | Provenance |
|---|---:|---|
| Updated core signed graph | 231,386 edges; 46,114 nodes | `results/causal_modelling/cosmos_mouse_pkn/build_summary.json`, core |
| Curated lipid-expanded graph | 251,998 edges; 52,157 nodes | Same, lipid_curated |
| Graph used for MOON, with measured-species hypotheses | **254,589 edges; 52,755 nodes** | Same, lipid_rules; independently counted CSV triples and unique endpoints |
| Gene-product nodes in graph used for MOON | 9,673 | Same, lipid_rules.node_types |
| Additional unresolved/multi-locus protein nodes | 13 | Same |
| Metabolite nodes (compartment-specific entities) | 9,116 | Same |
| Catalyst/reaction/direction instance nodes | 33,953 | Same |
| Mouse-GEM source model | 12,987 reaction records; 8,454 compartment-specific metabolite records | `mouse_gem_manifest.json` |
| Mouse-associated lipid Rhea master export | 2,826 reaction masters; 5,360 directional/undefined variants; 896 catalysts | `rhea_summary.json`; these are master-export counts, before conservative projection |
| Added measured-species rules | **149 reaction records, 8 directional templates, 33 mouse catalysts, 149 measured features** | `lipid_rules_summary.json` |
| Added graph structure from measured-species rules | 2,591 edges; 449 catalyst/reaction instances; 149 measured-lipid nodes | Difference of variant counts plus node-table/evidence intersection |
| Fresh metabolite-to-protein activity evidence | 83 rows before graph reconciliation; **77 unique final edges, 79 retained evidence rows** | `regulation_summary.json`; final graph/evidence intersection |
| Direct chemical-to-transcript expression edges | **0** | `regulation_summary.json` |

Resource-attributed final signed edges were independently counted by joining
`edge_evidence.jsonl.gz` to `pkn_lipid_rules.csv`, admitting only evidence with
`first_variant` in core/lipid_curated/lipid_rules:

| Resource | Unique final signed edges | Retained evidence rows |
|---|---:|---:|
| OmniPath GRN | 34,604 | 34,897 |
| OmniPath signaling | 18,743 | 19,183 |
| Mouse-GEM | 178,665 | 178,665 |
| Rhea | 20,535 | 20,655 |
| Guide to Pharmacology | 77 | 79 |
| Measured-species rules | 2,591 | 2,591 |

These resource counts are **not additive**: 626 signed triples have both GRN and
signaling attribution. Graph edges and reaction-instance nodes are also not
counts of independent experiments or distinct biochemical reactions. The
33,953 reaction-instance nodes comprise 29,820 Mouse-GEM, 3,684 Rhea and 449
rule-derived instances. Isoenzyme alternatives, permitted directions and
independent source annotations can instantiate multiple nodes for related
chemistry.

Suggested Results wording: “We assembled a mouse prior-knowledge network
integrating signed signaling and transcriptional interactions with metabolic
reactions and metabolite–protein regulation. The curated foundation contained
251,998 directed edges; an explicitly annotated extension of measured-species
reaction hypotheses increased the network to 254,589 edges connecting 52,755
nodes.” This reflects the resource's scale without presenting all edges as
experimentally validated, or the species-rule extension as curated substrate
specificity.

## Frozen versions and assembly

- Direct OmniPath mouse exports, retrieved 2026-09-08: `collectri` for the GRN;
  `omnipath,ligrecextra` for signaling. Query URLs and checksums are in
  `results/causal_modelling/cosmos_pkn_audit/manifest.json`; raw exports are in
  its `raw/` directory. These are mouse exports, but source annotations can
  involve orthology; do not call every annotation a mouse experiment.
- Mouse-GEM commit `fc769ac1da421c4e4e9e1cf2a9aec3fe4f02ef6a`; exact SBML and
  mapping snapshots are in `support/cosmos_mouse/raw/`.
- Rhea release 142, 2026-09-02, using reaction RDF, directional identifiers and
  reviewed UniProt crosswalks; SwissLipids mouse enzyme and lipid hierarchy
  snapshots retrieved 2026-09-08. See `rhea_manifest.json`.
- GtoPdb 2026.2, published 2026-06-15, retrieved 2026-09-08: native mouse
  metabolite action records with explicit signs and screened ChEBI mappings.
  Only the `curated_mouse_current` tier enters the graph used for MOON.
- ChEBI release 255, 2026-09-01, supports identity, ancestry and acid/base
  audits. The optional pH7.3-normalized graph is a separate export, not the
  `pkn_lipid_rules` graph used for this analysis.
- Mouse UniProt accession/primary-gene snapshot retrieved 2026-09-08 validates
  catalyst accessions and maps single-gene proteins to gene-product nodes.
- Measured-species rule definitions were checked against LINEX2 commit
  `04e0099655a6b7a8cb57aaeca7b1301b3524079d`, dated 2024-10-02. Rhea 142 supplies
  current directional chemistry and mouse catalysts. Rules transform DG↔PA,
  PC→DG/PA/PS, DG→PC/PE and PI→DG. Only matching reported diacyl compositions
  are paired; unreported sn positions remain unresolved. Sum-only species,
  ether/oxidized lipids, sphingolipids and coupled multi-lipid transformations
  were excluded from this first rule extension.

The assembly uses signed signaling and transcriptional edges, plus explicit
reaction instances with positive catalyst→reaction, substrate→reaction and
reaction→product connectors. Full participants, stoichiometries, compartments,
source references and GPR logic are retained in master JSONL records. The
signed projection omits AND-complex approximations, undefined Rhea directions,
unresolved transport compartments, generic/symbolic Rhea participants,
unverified catalysts, boundary reactions, blocked directions, unsigned or
undirected interactions, self loops and opposing signs for the same pair.
Complex approximations, uncertain directions and cached STITCH/BRENDA evidence
remain in the separate exploratory export and do not enter `pkn_lipid_rules`.

Positive metabolic connectors represent qualitative reaction traversal. They
do not assert substrate activation of an enzyme, steady-state flux balance or
a direct effect on metabolite concentration. Direct chemical regulatory edges
refer to protein activity; effects on transcription require downstream GRN
connections. No observed correlation or abundance sign was used to create the
species-rule prior.

## Coverage denominators and limits

The build's mapping universe is 1,254 features from
`results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gonade_type_ttop.csv`, not all 1,386 assay features. The
curated graph contains exact-resolution chemical matches for only seven of
these features; the rule extension connects 149 measured features separately.
The build's “64 Factor 1 lipids” denotes the earlier correlation-network
subset, not the later full set of Factor 1 lipid weights used for MOON. Avoid
importing its seven-of-64 rule coverage into the MOON input description.
Compatible/hierarchy mapping candidates do not establish exact PKN membership
and must not be added to exact coverage counts.

The full regulation table's 26,245 unique edges and 41,235 evidence rows include
cached exploratory data. Neither number describes the fresh curated
metabolite-to-protein layer used in this analysis. MRCLinksDB's 641 assumed
positive associations were excluded entirely.

## Script/data pointers

Assembly entry points are `scripts/pkn/build_mouse_gem.py`,
`build_rhea_lipid_reactions.py`, `build_lipid_regulation.py`, `map_lipids.py`,
`build_lipid_rules.py` and `assemble_mouse_pkn.py`. Auditing/export validation
uses `scripts/pkn/audit_mouse_pkn.py`, `test_pkn.py`, `validate_export.R` and
`write_build_manifest.py`. The existing scripts and manifests still refer to
historical output paths after the user's directory move; any claimed executable
reproduction command must address this relocation explicitly.

The main audit bundle is
`results/causal_modelling/cosmos_mouse_pkn/{build_manifest.json,build_summary.json,coverage_summary.json,edge_evidence.jsonl.gz,nodes.jsonl,mouse_gem_reactions.jsonl,rhea_reactions.jsonl,lipid_rules_reactions.jsonl,lipid_rules_pairs.csv,feature_coverage.csv}`.
Detailed build interpretation is in `docs/cosmos_mouse_pkn_build.md`.

## Primary-resource citations checked during this audit

The [OmniPath resource](https://omnipathdb.org/) documents combined signaling
and GRN knowledge and the lack of tissue specificity. It identifies Türei et
al. (2021), “Integrated intra- and intercellular signaling knowledge for
multicellular omics analysis”, *Molecular Systems Biology* 17:e9923.

The [Mouse-GEM repository](https://github.com/SysBioChalmers/Mouse-GEM) identifies
Wang et al. (2021), “Genome-scale metabolic network reconstruction of model
animals as a platform for translational research”, *PNAS* 118:e2102344118,
[doi:10.1073/pnas.2102344118](https://doi.org/10.1073/pnas.2102344118).
Use frozen local SBML counts above, not the current repository README's model
summary, which differs from the snapshot.

[Rhea direction documentation](https://www.rhea-db.org/help/reaction-side-direction)
explains distinct undefined, left-to-right, right-to-left and bidirectional
identifiers; merely listing left and right sides does not specify flux.

The [LINEX2 repository](https://github.com/LipiTUM/LINEX2_package) identifies
Rose and Koehler et al. (2023), “Lipid network and moiety analysis for revealing
enzymatic dysregulation and mechanistic alterations from lipidomics data”,
*Briefings in Bioinformatics*, bbac572,
[doi:10.1093/bib/bbac572](https://doi.org/10.1093/bib/bbac572), and identifies
Rhea and Reactome as data sources. Our species-rule implementation uses the
specific Rhea-derived templates described above; it is not a generic import
of every LINEX2 edge.

The [GtoPdb description](https://www.guidetopharmacology.org/about.jsp) supports
describing it as an expert-curated ligand–activity–target resource. Version
2026.2 is verified from the frozen native file headers and manifest, rather
than inferred from the live website.
