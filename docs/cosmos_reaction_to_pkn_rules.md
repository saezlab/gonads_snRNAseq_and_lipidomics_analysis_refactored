# Verified reaction-to-COSMOS projection rules

The completed follow-up is documented in
[`cosmos_factor1_expanded_moon_rerun.md`](cosmos_factor1_expanded_moon_rerun.md).
An initial Ar/Arnt→Cer topology did use water as an artificial bridge. The
preferred conservative sensitivity blocks generic cofactor substrate traversal,
retains full masters, and displays consecutive scoring-layer paths: 138 edges,
eight TFs and 13 measured features at 0.75/0.75. The rules below describe the
complete resource projection; the subsequent working-graph policy is explicit
and separate.

Implementation audit performed on 6 October 2026 against the current project scripts, the frozen expanded reaction records, the exported edge evidence and installed `cosmosR` 1.21.2. This is a code/data description, not a claim that the resulting lipid concentration mechanisms are established.

## Two distinct operations

**Species instantiation** creates a candidate chemical reaction from a generic directional Rhea template. **Graph projection** subsequently represents that complete reaction as ordinary signed `source, interaction, target` edges. These operations should not be conflated: a defensible structural substitution does not make the resulting graph a kinetic or stoichiometric model.

The expanded species builder checks compatible measured chains and modifies only declared participant slots. Other participants, coefficients, catalysts and provenance remain in the complete reaction record. Lyso acylation/hydrolysis conserves the retained chain and names the missing acyl donor/free-fatty-acid product as an unmeasured hypothesis. Sphingolipid rules conserve the base and N-acyl identities separately; sphingomyelin synthesis additionally requires a chain-compatible PC/DG pair. Pemt remains three methylation reactions with two typed unmeasured intermediates. Glycan interpretations are separate conditional identities. No additional measured value is assigned by the builder.

Code: `scripts/pkn/build_expanded_lipid_rules.py:147` and `scripts/pkn/expanded_lipid_matching.py`. Full records: `results/causal_modelling/cosmos_mouse_pkn_expanded_species/expanded_species_reactions.jsonl`.

## Projection into the signed graph

For a reaction `A + B -> C + D`, catalysed by mouse enzyme E, the code generates a reaction-specific node R and these edges:

| Edge | Sign | Meaning in this construction |
| --- | ---: | --- |
| E -> R | +1 | Catalyst enables this reaction hypothesis |
| A -> R | +1 | Substrate participates in the reaction |
| B -> R | +1 | Other substrate/cofactor participates in the reaction |
| R -> C | +1 | Product formation |
| R -> D | +1 | Other product formation |

R is named `Enzyme__<hash>`, but denotes an **enzyme–reaction–direction instance**, not a second physical protein. A single gene can therefore connect to many distinct reaction instances. The expanded hash is determined by resource, generated reaction ID and catalyst identity; the baseline hash also contains traversal orientation. Alternative catalyst genes have separate instances, so the graph does not imply that all alternative enzymes must be present together. Mouse UniProt accessions are mapped to primary mouse symbols where possible; unresolved mouse accessions retain their identity.

All metabolic projection signs are +1. Negative signs elsewhere in the PKN come from signed transcriptional/signaling/lipid–protein regulation; they are not assigned by converting a substrate into a depleted metabolite. In particular, the code does **not** create `enzyme -> substrate` edges with -1. Nor does it create a direct signed enzyme–lipid concentration relationship independent of the reaction node.

Code: `scripts/pkn/assemble_mouse_pkn.py:140` and `scripts/pkn/build_expanded_lipid_rules.py:256`.

## Actual measured LPC example

Expanded record `EXPANDED_LIPID_RULE:03c3d71b154cae067166`, supported by directional `RHEA:12938`, contains:

`LPC(16:0) + hypothetical ACoA(18:0) -> PC(16:0_18:0) + CoA`

Each participant coefficient is 1. The added chain is 18:0; its acyl-CoA availability is unmeasured. The reported lipid notation does not establish its sn position, double-bond positions, stereochemistry or the enzyme's exact species specificity.

The mouse catalyst alternatives in this record project to `Lpcat3`, `Mboat1` and `Mboat2`. For the Lpcat3 alternative the **exact exported triplets** are:

```text
Gene__Lpcat3                                      +1  Enzyme__4816473df5b47346
Metab__MeasuredLipid:LPC(16:0)_unknown             +1  Enzyme__4816473df5b47346
Metab__HypotheticalAcyl:ACoA(18:0)_unknown         +1  Enzyme__4816473df5b47346
Enzyme__4816473df5b47346                          +1  Metab__MeasuredLipid:PC(16:0_18:0)_unknown
Enzyme__4816473df5b47346                          +1  Metab__CHEBI:57287_unknown
```

The two other catalysts have their own reaction-instance IDs and the same chemical participants. This supports alternatives to investigate; it does not demonstrate that each enzyme acts on these exact measured lipid species in hepatocytes.

As a second example, record `EXPANDED_LIPID_RULE:05c00a47285c73e76cb2`, supported by `RHEA:19254`, preserves `SM(18:1;O2/26:0) + H2O -> Cer(18:1;O2/26:0) + phosphocholine + H(+)`. Its mouse enzyme alternatives are Enpp7 and Smpd1–5. The graph retains water, phosphocholine and proton rather than collapsing this into an isolated SM -> Cer pair.

## Direction, compartments and provenance

- Rhea LR and RL records use the explicit RDF substrate/product sides, so RL is already oriented correctly before projection. The measured-species extension accepts only an annotated LR or RL traversal, with no invented reverse. `reversible=false` here means only that this accession permits one annotated traversal; it is not proof of thermodynamic irreversibility.
- The baseline admits an explicitly BI-annotated reaction in both directions. UN records are explicitly marked direction-uncertain and placed only in the exploratory sensitivity tier; they are absent from the selected baseline plus chain-compatible extension. Mouse-GEM directions follow its allowed flux bounds.
- Mouse-GEM OR/AND gene–protein rules are expanded into alternative catalyst groups. A multi-gene complex is retained in master/node evidence but restricted to the exploratory tier, because an ordinary signed graph does not enforce its AND requirement. The expanded Rhea rules use single-gene alternative catalysts, not inferred complexes.
- Chemical IDs include their compartment: `Metab__<chemical>_<compartment>`. `unknown` remains distinct from cytosol/other compartments; no synthetic bridges are added. Rhea transport records with unresolved compartments are excluded. The measured Rhea-derived extension uses unknown localization and consequently cannot establish organelle-specific mechanisms.
- Full participants, stoichiometric coefficients, source versions, template IDs, references and substitutions remain in master records and edge evidence. The MOON input itself has only three columns: each coefficient is reduced to the +1 participation/formation sign, not used as a numerical flux weight. Self loops are excluded; duplicate triplets are consolidated; opposing signs for the same endpoint pair are excluded by baseline construction/cleanup.

Direction source: `scripts/pkn/build_rhea_lipid_reactions.py:217`; baseline tiering: `scripts/pkn/assemble_mouse_pkn.py:140`–`213`; compartment IDs: `scripts/pkn/assemble_mouse_pkn.py:105` and expanded builder `:277`.

## What MOON can and cannot infer from this representation

The complete chemical record is a hyperedge: the reaction requires its substrates together. Its ordinary graph projection allows traversal through any individual incoming edge. Installed `cosmosR::moon()` does not inspect the complete record or enforce joint substrate availability, mass balance, reversibility constraints, reaction rates or compartment-specific flux. An enzyme can acquire a score from a product observation even when other substrates are unmeasured. An extracted path can omit a co-substrate or co-product even though the retained master record contains it.

MOON uses the exported signs as `mor` in recursive `decoupleR::run_ulm()` calls. A high reaction-instance score means its downstream score pattern agrees with the positive prior; it is not a measured catalytic rate. More enzyme activity need not mean higher steady-state product abundance: synthesis, consumption, transport and substrate availability jointly determine concentration. Accordingly, a TF-to-reaction-to-lipid route is an exploratory pathway hypothesis with prior/sign agreement, not an identified concentration or flux mechanism.

The installed function also removes edges **whose source is a supplied downstream observation** before scoring. Observed lipids therefore behave as terminal readouts within this fit; a measured lipid cannot also propagate its value through outgoing substrate-participation edges in that same scoring call.

Verification: installed `cosmosR` version 1.21.2, `body(cosmosR::moon)`; local current source `/Users/dugourd/Documents/cosmosR/R/decoupleRnival.R:35`.

## Risks to audit in the expanded rerun

1. **Redundant alternatives:** the chain-compatible extension contains 3,176 generated reaction hypotheses and 13,412 catalyst-specific reaction instances. Before topology pruning these have 2,088 distinct signed output sets; 13,195 instances lie in 1,871 duplicate-output groups. These are alternatives, not independent biochemical observations. Safe same-children compression should be retained with protected input/regulon identifiers. It removes exact downstream duplicates; it does not remove every correlation caused by partial/shared outputs or turn generated reactions into independent evidence.
2. **Cofactor shortcuts:** in the full chain-compatible graph, `H2O_unknown` has 151 incoming and 10,842 outgoing edges, `H(+)_unknown` has 11,381 incoming and 301 outgoing edges, and `CoA_unknown` has 2,740 incoming and 211 outgoing edges. A permissive path can bridge otherwise unrelated reactions through these shared chemicals. These counts do not establish that a selected solution uses such a shortcut. Audit actual selected routes and scoring support; report those that depend on such bridges. A separate sensitivity that removes non-informative currency-metabolite traversal would need a declared rationale and should retain the full masters.
3. **Conditional chemistry:** chain compatibility does not establish enzyme substrate specificity or positional identity, and synthetic acyl pools have no measured values. Avoid interpreting every alternative catalyst/path as equally supported in the liver. Glycan alternatives should remain a separate sensitivity and should not be counted as independent lipid measurements.
4. **Coverage versus fit:** an increased number of structurally mapped endpoints changes the downstream evidence and ULM background. Improved raw reachability need not improve upstream TF sign agreement or compact-network quality. Report mapping, topology, scoring and final selection as separate stages.

The duplication and cofactor counts above were computed from the frozen `pkn_expanded_chain.csv` plus `extension_edge_evidence.jsonl.gz`, before expression filtering, boundary handling, topology pruning or compression. They are diagnostic risks, not claims of artifacts in the eventual selected network.
