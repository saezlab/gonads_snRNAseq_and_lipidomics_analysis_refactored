# Selected Factor 1 lipid-to-TF network: prostanoid and reaction evidence

Audited on 15 September 2026 against the existing [selected all-TF solution](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/selection.json). This is an evidence audit, not a refit or PKN expansion. It distinguishes annotated chemistry and receptor binding from hypotheses about this liver dataset.

The selected 455-edge, 233-node network contains **two chemical-to-protein regulatory transitions: PGF2α→Ptger3 and PGD2→Ptger3**. Both prostanoids are unmeasured intermediates in the selected hypothesis, rather than the measured LPC entry species. The previously audited DHA–RXR route is not the receptor route selected here. [Selected regulatory transitions](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/network_regulatory_transitions.csv).

## What supports the receptor transitions

Both transitions inherit positive signs from the frozen GtoPdb 2026.2 **Full agonist** action annotation for mouse EP3/Ptger3, UniProt P30557, target ID 342. The frozen export supplies binding affinities, not functional-response EC50 values:

| Ligand | GtoPdb ligand ID | Export row enumeration | pKi | Original Ki in the frozen export | Original PKN edge |
| --- | --- | --- | --- | --- | --- |
| PGF2α | 1884 | 8293 | 7.12 | 75 nM | `reg_159e4e6ac413b95f846b` |
| PGD2 | 1881 | 8278 | 6.55 | 280 nM | `reg_2b78f8e6bc97e51127e9` |

These original Ki values were read directly, not back-calculated from rounded pKi. The export's assay-description fields are empty. Row numbers are the builder's enumeration of this frozen export, **not official persistent GtoPdb interaction IDs**. [Frozen export](../support/cosmos_mouse/raw/regulation_gtopdb_interactions.csv.gz); [recovered complete records](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/reaction_evidence_audit/selected_regulatory_evidence.csv); official [PGF2α](https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=1884&tab=biology) and [PGD2](https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=1881&tab=biology) records.

Both records cite **Kiriyama et al., 1997, PMID 9313928**. This primary study expressed mouse prostanoid receptors in CHO cells and measured competitive radioligand displacement. It supports receptor-binding specificity in a heterologous cell system, rather than an experiment in hepatocytes or sex-biased liver physiology. The accessible primary abstract describes binding measurements; it does not establish the quantitative functional efficacy of these two ligands at EP3. The positive activity sign therefore remains explicitly attributed to database curation, while Ki is interpreted strictly as binding affinity. Full agonism, receptor coupling, effective exposure and downstream transcription in this dataset are not independently demonstrated by that binding result. [Primary study](https://pubmed.ncbi.nlm.nih.gov/9313928/); DOI: [10.1038/sj.bjp.0701367](https://doi.org/10.1038/sj.bjp.0701367).

The selected regulatory sources use the conjugate-base identifiers `CHEBI:57404` and `CHEBI:57406`, whereas the original regulatory records use acid identifiers `CHEBI:15553` and `CHEBI:15555`, respectively. These were connected by the existing frozen chemical-normalization crosswalk. This is a stated identity-normalization sensitivity, not new evidence for biological signaling or measured protonation. [Selected edge provenance](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/network_regulatory_transitions.csv).

## Complete chemistry behind the principal route

The following are abbreviated versions of the full frozen Rhea equations. Catalyst lists include **all annotated alternatives for each listed reaction in the assembled mouse resource**, not just the first catalyst shown on a representative path. An enzyme name identifies a possible catalyst; its inclusion does not measure that enzyme's activity.

| Rhea reaction | Abbreviated complete equation | Frozen mouse catalyst annotations |
| --- | --- | --- |
| [37480](https://www.rhea-db.org/rhea/37480) | LPC(18:0) + arachidonoyl-CoA → PC(18:0/20:4) + CoA | Lpcat3 |
| [40520](https://www.rhea-db.org/rhea/40520) | PC(18:0/20:4) + H₂O → LPC(18:0) + arachidonate + H⁺ | Pla2g4a, Pla2g6, Pla2g10 |
| [23729](https://www.rhea-db.org/rhea/23729) | Arachidonate + AH₂ + 2 O₂ → PGH2 + A + H₂O | Ptgs1, Ptgs2 |
| [42597](https://www.rhea-db.org/rhea/42597) | Arachidonate + 2 O₂ → PGG2 | Ptgs1, Ptgs2 |
| [42601](https://www.rhea-db.org/rhea/42601) | PGG2 + AH₂ → PGH2 + A + H₂O | Ptgs1, Ptgs2 |
| [10601](https://www.rhea-db.org/rhea/10601) | PGH2 → PGD2 | Hpgds |
| [45314](https://www.rhea-db.org/rhea/45314) | PGH2 + NADPH + H⁺ → PGF2α + NADP⁺ | Akr1b1, Akr1b7 |
| [12894](https://www.rhea-db.org/rhea/12894) | PGH2 → PGE2 | Ptges, Ptges2 |
| [24510](https://www.rhea-db.org/rhea/24510) | PGE2 + NADPH + H⁺ → PGF2α + NADP⁺ | Cbr1 |

Here, AH₂/A denotes the reduced/oxidized donor pair specified generically by Rhea. The overall PGH2 reaction and its two component steps describe overlapping chemistry; their coexistence is not independent experimental support. Similarly, the direct PGH2→PGF2α reduction and the route through PGE2 are alternatives. NADPH and the oxygen/reducing requirements remain necessary biochemical conditions even though the lipid-only traversal graph does not carry all nonlipid participants.

The initial two reactions cancel LPC and PC:

```text
Arachidonoyl-CoA + H₂O → arachidonate + CoA + H⁺
```

Thus the arachidonate chain comes from an **unmeasured arachidonoyl-CoA pool**, while LPC(18:0) is regenerated. The route links LPC-associated acyl remodeling to possible mobilization of an existing polyunsaturated acyl chain. It is not synthesis of arachidonate from the LPC stearoyl chain, and the LPC abundance weight alone does not establish net arachidonate or prostanoid production. The annotated reactions are chemically plausible conditional steps; the network does not solve their joint stoichiometry, flux or substrate limitations.

## Complete audit and reproduction

The selected network contains **51 reaction-instance nodes corresponding to 31 distinct Rhea accessions**. The [complete reaction table](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/reaction_evidence_audit/selected_reaction_evidence.csv) covers every one, including additional phospholipid, glycerolipid and N-acyl routes outside the representative prostanoid path. It preserves complete equations, all selected instance IDs, all frozen mouse catalyst genes and UniProt annotations, original participant stoichiometries, incident selected edges and annotation provenance. Participant presence in the graph is recorded separately from measurement status; neither presence nor omission measures biochemical availability.

Reproduce the tables without changing the fit:

```bash
python3 scripts/moon/ancillary/audit_lipid_tf_selected_reactions.py \
  --results=results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs \
  --output=NEW_EMPTY_DIRECTORY
```

The [audit script](../scripts/moon/ancillary/audit_lipid_tf_selected_reactions.py) refuses occupied output directories. [Dependency hashes](../results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs/reaction_evidence_audit/dependency_hashes.csv) identify the exact selected network, frozen PKN records, export and script inspected.

The defensible biological statement is a **conditional hypothesis linking LPC-associated remodeling, arachidonate-derived prostanoids and an EP3-connected signaling network consistent with supplied TF footprints**. It does not identify LPC as a direct EP3 ligand, establish exposure of hepatocytes to the inferred prostanoids, or prove that these reactions caused the observed TF activities. Those TF activities are supplied downstream observations; fitting a route to them is not independent validation of their regulation.
