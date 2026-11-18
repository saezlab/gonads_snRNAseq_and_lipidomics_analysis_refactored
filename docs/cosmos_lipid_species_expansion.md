# Expanded measured lipid reaction hypotheses

Built and audited on 6 October 2026. This separate exploratory extension of the frozen mouse `pkn_lipid_rules` follows the [mouse brain lipid atlas comparison](cosmos_lamanno_pkn_comparison.md). It increases candidate measured reaction-product coverage from **128 to 419 features**, or **431** with explicitly conditional glycan identities. This demonstrates that the previous deficit was partly caused by narrow species-instantiation rules. It does not establish enzyme specificity or agreement with Factor 1 activity signs.

The manuscript PKN, input vectors, selected networks and MOON scores remain unchanged. No new measurement vector or MOON fit was generated.

## Construction and structural gates

The builder uses the frozen mouse-supported Rhea 142 export and a reviewed catalog of 33 directional templates. It does not import the brain study's lipid-pair table as signed regulatory evidence. Original orientation, catalysts, references, participants and coefficients are retained. Every substitution declares its side/index and source chemical. Unsupported directions and transport are excluded.

| Added family | Instantiated hypotheses |
| --- | ---: |
| Acyl-CoA-dependent lysolipid acylation | 668 |
| Diacyl lipid hydrolysis | 1,215 |
| Sphingomyelin hydrolysis | 44 |
| Coupled sphingomyelin synthesis and annotated reverse | 1,054 |
| Glycosphingolipid transformations with conditional sugar identity | 54 |
| PE and PS interconversion | 42 |
| Three separate Pemt methylation steps | 153 |
| Total | 3,230 |

These are generated alternatives, not independent experiments. Templates, alternative catalysts and coupled lipid pairs can support overlapping hypotheses.

- **Lysolipids:** the diacyl feature must contain the retained lyso chain. The other chain becomes an explicitly unmeasured fatty-acid product or acyl-CoA donor. Reported lyso sn positions, including `LPC(chain/0:0)`, must agree with template specificity; unreported positions remain conditional.
- **Sphingolipids:** sphingobase and N-acyl identities are conserved separately. Fixed sphingosine templates require `18:1;O2`; unresolved 4E/stereochemistry remains hypothetical. Generic two-R templates allow variable base chains within the audited unmodified O2 scope. Unsupported O3 bases and oxygenated N-acyl chemistry are excluded.
- **Sphingomyelin synthesis:** Cer/SM substitutions require a joint chain-compatible measured PC/DG pair; donor and co-product remain explicit.
- **Pemt:** PE-to-PC retains three reactions and two typed, unmeasured chain-specific intermediates. SAM, SAH and proton participants are preserved.
- **Glycans:** HexCer is not silently equated with both glucose and galactose identities. Alternatives have separate `CandidateMeasuredLipid:` nodes carrying the source measurement ID and glycan condition. No measurement vector is duplicated across these nodes.

Sum expansion, ethers/plasmalogens, general TG/CL reactions and unsupported modifications remain outside scope. The only measured Hex2Cer is summed, so the two audited HexCer-to-Hex2Cer templates produce zero exact-chain instances. No directional mouse Rhea Hex2Cer-to-HexCer template was admitted.

## Mapping and topology

| Property | Baseline | Chain-compatible | Including glycan alternatives |
| --- | ---: | ---: | ---: |
| Signed projected edges | 254,589 | 331,032 | 331,472 |
| Active graph nodes | 52,755 | 66,621 | 66,731 |
| Candidate measured product features | 128 | 419 | 431 |
| Reachable within four directed edges | 72 | 339 | 349 |
| Reachable within six directed edges | 125 | 389 | 401 |
| Reachable within eight directed edges | 125 | 419 | 431 |

Reachability uses the same frozen 49 TF seed IDs on the unfiltered graph. It ignores activity signs, RNA expression/coherence, joint substrate availability and scoring. Generic cofactors can contribute permissive paths. These are topology diagnostics, not the selected-workflow pruning count or validated signaling explanations.

The chain-compatible products span 15 classes, including 47 LPC, 17 LPE, 27 Cer and 17 SM features. Glycan alternatives add ten HexCer and two Cer product features. The historical 64-lipid correlation subset has 24/28 candidate products; that flag does not describe Factor 1 weight availability, which covers all 1,254 features.

Positive metabolic edges describe catalyst/substrate participation and product formation, not concentrations or flux. Synthetic acyl pools and methylated intermediates have no measurements or asserted exact ChEBI identities. New species records use `chemically_balanced=null`: formal source chemistry and matched chains are retained, but independent atom-balance testing was not performed.

## Why the brain study appeared more connected

**Dataset-specific curation matters.** Their Rhea draft was manually refined using paired LC-MS information. Their table includes reaction families omitted by our earlier species rules. Our new 128-to-419 gain directly tests this explanation: low coverage was not solely a frontier of biological knowledge.

**Their lipid units and panel differ.** Their study analyzes 172 annotated lipids, with many multi-chain labels at sum resolution; Table 10 connects 113 labels. Our 1,254-feature panel contains many resolved chains, plus 180 TGs and 89 cardiolipins, classes absent from Table 10. Summed labels pool molecular possibilities that our retained-chain checks distinguish. [Study](https://www.nature.com/articles/s41586-026-11050-0), [published supplements](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-11050-0/MediaObjects/41586_2026_11050_MOESM3_ESM.zip).

**Upstream constraints differ.** Their receptor activity proxies partly derive from correlations with lipid changes and local receptor expression; ours use RNA-derived TF footprint signs. Receptor-to-enzyme routes can use post-translational signaling without a TF transcriptional boundary. This can make reconstruction more permissive; its exact contribution has not been quantified. [Supplementary methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-11050-0/MediaObjects/41586_2026_11050_MOESM1_ESM.docx).

**Connected labels and selected targets are different stages.** Our original sequence was 1,254 weights → 155 mapped observations → 128 active product endpoints → 74 observations after selected topology preparation → three DGs in the 388-edge solution. The final number additionally reflects score/sign agreement, RNA coherence, dual thresholds and the requested compact size. It cannot be compared directly with their 113 reaction-table labels. Esr1, Nfic and Dbp survived topology/compression but were omitted when lipid-driven scores disagreed with RNA-derived activity signs.

**Algorithms and reporting scope differ.** CORNETO can optimize a parsimonious fitting subnetwork; MOON computes recursive ULM activity scores before extracting sign-consistent edges. A selectable path may still yield weak MOON scores when informative children are sparse or backgrounds unfavorable. The atlas reports per-lipizone models over 46 lipizones, whereas our primary result is one Factor 1 signature. [CORNETO framework](https://www.nature.com/articles/s42256-025-01069-9).

Compression did not cause initial mapping loss: the original 74 retained observations were protected, and stronger compression reproduced the selected network. Keep compression in a subsequent comparison. The complete signed atlas PKN remains unavailable, so its final connected-lipid fraction and substrate/product sign semantics cannot be compared fairly.

## Reproduction and checks

From the repository root, using an empty replay destination:

```bash
python3 scripts/pkn/expanded_lipid_matching.py --self-test
python3 scripts/pkn/build_expanded_lipid_rules.py --output /tmp/expanded_species_replay
python3 scripts/pkn/validate_expanded_lipid_rules.py \
  --output /tmp/expanded_species_replay --report /tmp/expanded_species_validation.json
Rscript scripts/pkn/validate_expanded_lipid_rules.R /tmp/expanded_species_replay
```

Python uses its standard library; R requires `cosmosR` and `jsonlite`. The reviewed catalog is `support/cosmos_mouse/expanded_species_templates.json`; input/script hashes are in `build_manifest.json`.

Outputs under `results/causal_modelling/cosmos_mouse_pkn_expanded_species/` include:

- `pkn_expanded_chain.csv/.rds` and `pkn_expanded_glycan_conditional.csv/.rds`.
- `expanded_species_reactions.jsonl`: full master reactions, substitutions and assumptions.
- `extension_edge_evidence.jsonl.gz` and `nodes.jsonl`: mechanisms and conditional chemical identities.
- `measured_product_candidate_inventory.csv`: proposed feature-to-product candidates, not an input vector.
- `template_instantiation_audit.csv`: includes zero-match templates.
- `focus_tf_unfiltered_reachability.csv`: topology-only Esr1/Nfic/Dbp diagnostics.
- `independent_validation.json` and `r_validation.json`.

All baseline edges and baseline node records are retained. Independent validation passed 118 check categories: template participants, coefficients, catalysts, references/directions, sn/base/oxygen gates, glycan separation, synthetic acyl formulas, Pemt intermediates and exact graph/evidence projection. Both graphs passed the installed COSMOS cleaner unchanged. Fresh construction reproduced all plain-text outputs and uncompressed evidence exactly.

Use the chain-compatible resource first for a subsequent compressed MOON analysis, after auditing new input assignments. Glycan alternatives need a separate sensitivity so conditional interpretations are not counted as independent measurements.

The compressed TF→lipid follow-up was completed on 6 October 2026; see
[`cosmos_factor1_expanded_moon_rerun.md`](cosmos_factor1_expanded_moon_rerun.md).
It preserves frozen mappings, adds 268 one-to-one observations and compares
four/six-layer fits. A generic-water bridge made the initial stronger-cutoff
Cer topology unsuitable. The preferred cofactor-filtered, consecutive-layer
sensitivity has 138 edges, eight TFs and 13 measured features at 0.75/0.75.
This follow-up does not alter the resource builder's original inventories or
full chemical masters.
