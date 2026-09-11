# Mouse chemical-to-protein regulation build

The regulation layer estimates signs of **protein activity**, including receptor,
channel, enzyme, and transporter modulation. It contains no direct
chemical-to-transcript expression edges. Connecting receptor activity to RNA
requires the separate signed protein and transcriptional network.

## Reproduce

From the project root, run `python3 scripts/pkn/build_lipid_regulation.py`.
The default consumes frozen local inputs and uses only the Python standard
library. `--refresh` downloads public GtoPdb and PubChem inputs;
`--refresh-chebi` downloads the current ChEBI ontology. The mouse UniProt
validation export is supplied by the parent build. Source URLs, versions,
checksums, and limitations are in `regulation_manifest.json` and
`regulation_chebi_manifest.json`.

## Evidence tiers

| Tier | Evidence rows | Meaning |
|---|---:|---|
| `curated_mouse_current` | 83 | GtoPdb 2026.2 metabolite actions curated for mouse, with direct GtoPdb ChEBI cross-references passing available structure checks |
| `curated_mouse_current_secondary_chemical_mapping` | 20 | Same mouse action evidence; chemical mapping instead follows GtoPdb PubChem CID to PubChem ChEBI synonyms |
| `cached_brenda_fallback` | 772 | Explicit activator/inhibitor records from the frozen OmniPath-metabo mouse service cache |
| `cached_stitch_fallback` | 40,360 | Explicit activation/inhibition records from that cache; causal orientation cannot be independently reconstructed from the export |

The 41,235 evidence rows represent 26,245 distinct signed chemical/protein
edges. A resource row count is not the number of independent experiments.
There are 988 chemical/protein pairs with conflicting signs across retained
evidence. `regulation_sign_conflicts.csv` preserves those records; the parent
builder applies its own graph conflict policy. The main curated assembly uses
only `curated_mouse_current`; all other tiers remain exploratory.

Of the 103 fresh metabolite records, 56 map to direct lipid descendants in the
September ChEBI ontology. This does not establish overlap with a measured lipid
species, and broad lipid classes do not supply acyl-chain specificity.

Nonmetabolite GtoPdb small molecules are retained separately in
`regulation_pharmacological_supplement.csv` (179 rows), rather than included in
the metabolic evidence table. GtoPdb classification as Metabolite is not itself
a lipid annotation or proof of endogenous exposure in liver.

## Sign and mapping safeguards

- Only native mouse GtoPdb target records are used. The accession must exist in
  the current mouse UniProt export. Human-only CellPhoneDB interactions are not
  projected to mouse in this build.
- Binding, neutral, mixed, biphasic, and unknown action labels are excluded.
  Explicit Activator/Inhibitor Type can provide a sign when Action is blank.
  Partial agonists, biased agonists, and antagonists retain contextual caveats.
- All 641 MRCLinksDB rows are excluded from signed regulation because the
  upstream adapter assigns `+1` unconditionally, including transport relations.
- ChEBI and GtoPdb InChIKeys are compared whenever both exist. The build found
  incorrect source cross-references, including PGE2 assigned to CHEBI:116390,
  which is an unrelated chlorinated aromatic chemical. Ten metabolite records
  and sixteen pharmacological records were quarantined for connectivity or
  stereochemical/protonation discrepancies. Review records preserve the original
  source mapping and any exact-InChIKey ChEBI candidates, without repairing it
  silently. Generic classes lacking an InChIKey remain explicitly unverified.
- PubChem synonym mappings can be broad or ambiguous. For example,
  lysophosphatidylserine maps to both a generic class and a stearoyl species;
  alternatives remain separate with a shared mapping group and count.
- The cached BRENDA/STITCH export does not preserve original chemical IDs,
  mapping alternatives, or per-edge literature references. No missing evidence
  was fabricated. The cache components were modified on 2026-04-21 and
  downloaded on 2026-09-08; this is not a fresh BRENDA/STITCH rebuild.

## Chemical identity support

The frozen ChEBI release 255 (2026-09-01) provides 8,794 distinct explicit
conjugate acid/base pairs and 19,425 `alt_id` cross-references. These are exported
in `regulation_chebi_acid_base_relations.csv` and
`regulation_chebi_aliases.csv`. No identity relation is a causal PKN edge.

`regulation_chebi_source_crosswalk.csv` annotates 4,343 retained chemical
identifiers with names, formulae, charge, direct lipid hierarchy membership,
aliases, conjugate components, and optional Rhea pH7.3 targets with their origin.
There are 811 sources with explicit conjugate forms and 632 direct lipid
descendants. Current ChEBI lacks 57 cached chemical IDs; they remain transparent
in this annotation table and should not be treated as validated identities.

Examples supported by explicit ChEBI relations are sphingosine 1-phosphate
CHEBI:37550 ↔ CHEBI:60119, sphingosine CHEBI:16393 ↔ CHEBI:57756, PGD2
CHEBI:15555 ↔ CHEBI:57406, and arachidonic acid CHEBI:15843 ↔ CHEBI:32395.
Their Rhea choice of the major microspecies at pH7.3 is computational even
though the underlying acid/base relation is curated. No enantiomer, tautomer,
functional-parent, compartment, or lipid-class identity merge is performed by
this regulation builder.

## Validation

The builder was executed after every change. It checks evidence-ID uniqueness,
signed actions, mouse protein membership, valid ChEBI identifier syntax,
exclusion of MRCLinksDB from regulation, an expected fresh sphingosine
1-phosphate action, and known acid/base pairs. It also asserts that the
sphingosine enantiomer relation is not admitted to the acid/base table.

Independent inspection confirmed the parent's current chemical normalization
map is idempotent and all 7,492 nonidentity pH7.3 mappings connect forms within
explicit ChEBI conjugate components. The original chemical IDs remain in the
primary graph; normalization is a separate sensitivity view.
