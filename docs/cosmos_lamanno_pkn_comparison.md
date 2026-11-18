# Comparison with the mouse brain lipid atlas prior knowledge

Audit date: 2 October 2026. This comparison concerns the deposited resources for Fusar Bassini et al., *The lipidomic architecture of the mouse brain*, Nature, DOI [10.1038/s41586-026-11050-0](https://www.nature.com/articles/s41586-026-11050-0), and our frozen mouse `pkn_lipid_rules` resource. It does not alter the PKN or rerun MOON.

The lipid reaction component can be recovered, including Supplementary Table 10 and a broader annotation-derived lipid graph. The complete signed signaling-to-lipid graph used by CORNETO, together with its assembly code, was not located in the cited code and data deposits. The recovered resources identify useful extensions to our measured-species reaction hypotheses. They do not show that our underlying Rhea chemistry is absent or obsolete: most listed enzymes are already present, while our species-instantiation rules are deliberately narrower.

## Recovered resources and their scope

The supplementary-table guide explicitly identifies Table 10 as used for both enzyme-lipid correlation and CORNETO. The reaction-matrix methods describe LINEX2 followed by manual curation; the CORNETO methods additionally describe a Rhea draft refined with paired LC-MS acyl-chain information. These descriptions are complementary, but the available table does not expose the complete assembly or signed projection. [Supplementary methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-11050-0/MediaObjects/41586_2026_11050_MOESM1_ESM.docx), [supplementary tables](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-11050-0/MediaObjects/41586_2026_11050_MOESM3_ESM.zip).

| Recovered component | Contents established by this audit | Interpretation |
| --- | --- | --- |
| `LaManno_SuppTable10.csv` | 350 unique directed reagent-product rows; 113 lipid labels; 15 classes; 30 printed enzyme symbols; 23 class transformations | Published lipid reaction component. The 113 connected labels are a subset of the 172 annotated lipids described in the study. |
| `csv/corereactions_wenzymes.csv` in the data deposit | Byte-identical to the supplied and published Table 10 | Confirms the deposited identity of the table. |
| `csv/reactions_linex_now.csv` and `mixed/graph_from_annotations.pickle` | 4,562 directed lipid pairs; 675 graph nodes, of which 611 have an edge; 3,620 edges carry an annotation containing `RHEA:` | Broader deposited lipid draft; correspondence to the final CORNETO PKN has not been established. |
| Complete CORNETO PKN | Not located in the cited deposits | Its total size, exact signs, enzyme-lipid projection and complete lipid mapping cannot be independently compared. |

The 350 table rows expand to 522 row-enzyme annotations because some rows name alternative enzymes. These are not 522 independently established reactions. The broader graph has pooled enzyme, UniProt, Rhea and Reactome annotations, and the ordered lipid pairs exactly match `reactions_linex_now.csv`. All 675 nodes are typed as lipids; it is not a signaling graph. Its sum-composition labels have inferred representative molecular strings, which must not be interpreted as resolved measurements.

The search covered the published supplements, the cited [lipidbrainatlas repository](https://github.com/lamanno-epfl/lipidbrainatlas), its [archived code](https://zenodo.org/records/20702555), the [other-data bundle](https://zenodo.org/records/15379565), and EUCLID. The cited lipidbrainatlas snapshot is commit `3fd1f0db7763b781689f39e3cf835c5635560f31`; its [pregnancy notebook](https://github.com/lamanno-epfl/lipidbrainatlas/blob/3fd1f0db7763b781689f39e3cf835c5635560f31/009b-PregnancyAnalysis.ipynb) loads the core reaction table but does not provide the final CORNETO graph assembly. This is a retrieval limitation, not proof that the full graph does not exist elsewhere.

## Comparison with our reaction layer

Our complete selected PKN contains 254,589 signed projected edges and 52,755 nodes. Its underlying Rhea layer uses release 142, dated 2 September 2026: 5,360 reaction variants, 2,826 reaction masters and 896 catalyst accessions. These whole-network counts include signaling and generic metabolism and cannot be compared directly with Table 10's lipid-pair count.

| Comparable aspect | Table 10 | Our measured-species rule layer |
| --- | --- | --- |
| Instantiated lipid-pair hypotheses | 350 | 149 |
| Lipid classes represented | 15 | 6: DG, PA, PC, PE, PI, PS |
| Directed class transformations | 23 | 8 |
| Main scope | Lysolipid acylation/hydrolysis, headgroup changes and sphingolipid conversions | Chain-preserving diacyl headgroup conversions supported by directional mouse reaction templates |
| Structural resolution | Predominantly summed compositions for multi-chain lipids | Exact reported two-ester-chain multisets; unresolved sn positions remain uncertain |
| Reaction provenance | Reagent-product-enzyme table; no Rhea identifiers or complete stoichiometry in Table 10 | Rhea templates, catalysts, participants and stoichiometry retained in master records |

Of the 30 printed Table 10 enzyme symbols, **29 are already active in our PKN**, 27 have active Rhea reaction instances, and four have measured-species rule instances: `Plcb1`, `Pld1`, `Pld2` and `Smpd1`. `Glb1` and `Ugt8a` have active Mouse-GEM instances. The sole unmatched printed symbol is `Psd`; `Pisd` and its PS-to-PE reaction are present, suggesting a nomenclature issue that should be checked before import.

Table 10 contains 143 lyso-to-diacyl acyl-addition rows, 147 diacyl-to-lyso hydrolysis rows, and 24 sphingolipid-conversion rows. These families are outside our current species-instantiation scope, although their generic chemistry is already represented. Examples in our frozen Rhea records include:

| Reaction family | Example existing mouse-supported chemistry |
| --- | --- |
| PC hydrolysis | `Pla2g2a`, RHEA:38780, PC to LPC plus fatty acid |
| Lysolipid acylation | `Mboat1`, RHEA:12938, LPC plus acyl-CoA to PC; `Agpat1`, RHEA:19710, LPA to PA |
| PE and PS interconversion | `Ptdss2`, RHEA:27607, PE to PS; `Pisd`, RHEA:20829, PS to PE |
| Sphingolipid conversions | `Smpd1`, RHEA:19254, SM to Cer; `Sgms1`, RHEA:18766, Cer plus PC to SM plus DG |
| Glycosphingolipid conversion | `B4galt5`, RHEA:31496, glucosylceramide to lactosylceramide |

The Pemt-mediated PE-to-PC route also exists in our chemistry as three methylation steps; Table 10 collapses it to one pair. Gene coverage alone does not establish that every enzyme-pair assignment in the table agrees with our records. The enzyme inventory and transformation table identify where further reconciliation is needed.

## Mapping to our 1,254 lipid measurements

Mappings below compare reported nomenclature, not chemical identity or validated enzyme specificity. Conditional projections preserve class, oxygenation and ether state, but sum chain lengths and double bonds. They do not resolve isomers, sn positions or HexCer glycan identity.

| Mapping measure | Table 10 | Broader deposited lipid draft |
| --- | ---: | ---: |
| Measurements with identical normalized reported shorthand | 46 | 125 |
| Additional measurements requiring sum-composition projection | 91 | 267 |
| Total possible measurement matches | **137** | **392** |
| Distinct matched lipid labels | 90 | 268 |
| Maximum measurements sharing one matched label | 3 | 3 |
| Possible matches appearing on the product side | Not used as a validated product endpoint count | 377 |

The different normalization definitions are recorded in the script: Table 10 uses our parsed canonical names, while the broader draft also admits identical original labels after parentheses and limited class-synonym normalization. Neither uses the draft's inferred representative molecular string as an identity match.

Our previous direction-aware audit accepted **128 active measured product endpoints**, comprising 123 species-rule products and five exact chemical mappings. Only 32 of those overlap the 137 possible Table 10 matches; 105 table matches lie outside our accepted endpoint set, and 96 accepted endpoints lie outside Table 10. This is descriptive overlap between sets with different resolution and evidence requirements, not a fair coverage ratio. The union of 233 is not a new validated endpoint count.

The historical 64-lipid correlation-network subset has eight possible matches to Table 10 and 16 to the broader draft. That subset flag must not be confused with availability of Factor 1 weights: all 1,254 mapped measurements have Factor 1 loadings.

## Why summed labels cannot inherit every reaction

For example, a Table 10 reaction from LPC 16:0 to PC 36:1 is potentially consistent with PC(16:0_20:1), but cannot simply be assigned to PC(18:0_18:1): the latter does not contain the retained 16:0 chain. Matching the total carbon and double-bond count conceals this distinction.

Across the 288 table reaction-to-measurement assignments where a resolved diacyl measurement could inherit a lyso/diacyl transformation by sum matching, **172 lack the required retained lyso chain**, affecting 51 measured features. This does not invalidate the original summed-composition reaction; it invalidates that particular naive assignment to a resolved measurement. It also explains why nominal lipid-label coverage can substantially exceed defensible species-reaction coverage.

Other annotations merit review rather than automatic import: the printed `Psd` symbol; the LPC-to-LPC rows assigned to `Lpcat3/4`; and the PG-to-DG row assigned to `Plcb1`. These are audit flags, not demonstrated biochemical errors. Table 10 also contains no direct lipid-to-protein regulatory signs, so it cannot independently validate the relaxed STITCH signaling hypotheses or the selected LPS candidate.

## Modeling implications

Their CORNETO model combines OmniPath/SIGNOR signaling, TF-enzyme links and curated lipid reactions. Inputs are receptor activity proxies derived partly from lipid correlations and local receptor expression; outputs are pregnancy lipid changes. Expression weights, graph-dropout alternatives and permutations support model selection. This differs from our MOFA Factor 1 measurements, TF footprints and compressed MOON scoring. [Supplementary methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-11050-0/MediaObjects/41586_2026_11050_MOESM1_ESM.docx).

The table does not reveal how substrate consumption, product formation and enzyme regulation were signed in their final model. Our positive metabolic participation edges represent qualitative reaction traversal, not steady-state flux or the claim that enzyme activation increases substrate concentration. Comparing modeling algorithms would require the same audited graph and compatible input semantics.

The supported next step is to **extend our measured-species reaction rules**, using these resources to prioritize acyl-remodeling and sphingolipid families. Each proposed assignment should preserve retained chains, donors, co-products, mouse catalyst evidence, direction and unresolved structural ambiguity. PE-to-PC should preserve its multi-step chemistry. Such additions should remain a separate exploratory overlay until audited. Requesting the authors' full signed PKN and assembly/mapping code would resolve the remaining comparison gap.

## Reproduction and audit files

The new standalone audit is `scripts/pkn/compare_lamanno_lipid_network.py`; source snapshots and their manifest are under `support/cosmos_mouse/comparison/lamanno_2026/`. Outputs are under `results/causal_modelling/lamanno_2026_pkn_comparison/`:

- `summary.json`: counts, scope and mapping limitations.
- `table10_measured_mapping_candidates.csv` and `draft_measured_mapping_candidates.csv`: feature-level candidate assignments and resolution flags.
- `table10_chain_compatibility.csv`: retained-chain checks for individual lyso/diacyl assignments.
- `table10_enzyme_coverage.csv` and `reaction_family_comparison.csv`: overlap with our active PKN and species rules.
- `table10_vs_accepted_product_endpoints.csv`: descriptive comparison with the frozen product-endpoint audit.
- `dependency_hashes.csv`: SHA256 hashes of the script, comparison sources and frozen local inputs.

From the repository root:

```bash
python3 scripts/pkn/compare_lamanno_lipid_network.py --self-test
python3 scripts/pkn/compare_lamanno_lipid_network.py --output /tmp/lamanno_pkn_comparison_replay
```

The replay destination must be empty. Python 3.9 or later and its standard library are sufficient. The source manifest checks all deposited comparison-file hashes before computation. Table 10's SHA256 is `40da8d72675d3b25a9e491c381b405def77d09d016184fc6b01517f641fb9033`.

The original compressed graph member and an inert export reader are preserved with the source snapshot; exports can be regenerated in a temporary copy of that directory using `read_graph_annotations_inert.py`. The reader replaces pickle-referenced external classes with data-only stand-ins. The comparison itself reads only the exported CSV files. Source discovery and modeling interpretation are documented here; the standalone script reproduces the quantitative local comparison.
