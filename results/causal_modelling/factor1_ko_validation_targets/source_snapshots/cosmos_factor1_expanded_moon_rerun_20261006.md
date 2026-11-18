# Factor 1 TF→lipid MOON with expanded reaction hypotheses

Executed and independently audited on 6 October 2026. This is a separate
exploratory sensitivity of the manuscript TF→lipid analysis. The frozen
manuscript result and its replay assertions are preserved.

## Outcome

The preferred conservative result contains **138 edges, 94 nodes, eight TF
inputs and 13 measured lipid features**, at strict absolute inferred-score
cutoffs **0.75/0.75**, pruning six steps and scoring four layers. Six scoring
layers gives 177 edges, 117 nodes and nine TFs, with the **same 13 lipid
features and the same cutoffs**. Four layers is the more compact representative
within this tested grid. Neither choice is a calibrated statistical optimum.

The eight TFs are Ahr, Gata6, Klf10, Nr2f1, Nr4a2, Ppara, Rora and Tbx2.
The readouts comprise three DG, seven LPC and three PA features. The six-layer
sensitivity adds Hoxa10. All 138 primary-network edges follow decreasing,
consecutive MOON scoring levels, all 60 displayed reaction instances retain
their catalyst connectors, and the network contains no generic cofactor nodes
or bridges. The exact signed network, GraphML, 36 representative TF/readout
routes, complete reaction masters and participant audit are available in:

`results/causal_modelling/factor1_tf_to_lipid_expanded_species/currency_v2_expanded_layers4_lower_thresholds/expanded_network_audit/`.

These are low-to-moderate numerical consistency hypotheses. Some route
bottlenecks are close to the cutoff: Ppara's inferred score is −0.7566, while
Mapk1 on the Gata6→LPC(20:5) example has −0.7554. Cutoffs are not calibrated
P values, and increased mapping coverage does not demonstrate causal control.

## Coverage and controlled comparisons

The raw chain-compatible PKN has 331,032 edges and 66,621 nodes, including
3,176 new chain-compatible reaction hypotheses. It is an additive extension of
the frozen 254,589-edge measured-species PKN. The separate glycan-conditional
tier was not scored here. Its glucose/galactose alternatives are not treated
as extra observations.

There are 1,254 Factor 1 lipid features. The old frozen input bundle mapped
155 of them; 128 were active reaction-product endpoints in the old PKN. The
expanded resource has 419 candidate measured product features. The accepted
input ledger preserves all 155 old node identities and values and adds 268
previously unmapped features once each, giving **423 unique observations**.
It includes 418 chosen active products and five preserved non-product
mappings. These numbers answer different mapping questions and should not be
interchanged. In particular, LPC(20:0/0:0) keeps its frozen ChEBI identity,
although its newly generated measured-node alias is the product candidate.

Five features have both a frozen ChEBI product alias and a generated measured
product alias. Each receives only one loading. Values are read at full
precision from the frozen RDS, using the original feature ID and common
Factor 1 orientation (+1; testes-positive coordinate). CSV exports are checked
within 1e-13, and original mapped R doubles are restored verbatim.

| Graph and observations | Post-preparation lipid inputs | Scoring layers | Extraction | Result |
|---|---:|---:|---|---|
| Frozen manuscript graph and inputs | 74 | 4 | Historical topology, 1.85/0.55 | 388 edges, four TFs, three DG features |
| Expanded graph, frozen155 observations | 140 | 4 | Same topology/cutoffs | Empty |
| Expanded graph, expanded423 observations | 371 | 4 | Same topology/cutoffs | 617 edges, five TFs, eight Cer features |
| Expanded graph, expanded423 observations | 371 | 4 | Higher-cutoff topology, 3/1 | 101 edges, Ar/Arnt and three Cer features; unsuitable water bridge |
| Currency-filtered v2 graph, expanded observations | 369 | 4 | Consecutive score support, 1.85/0.55 | Empty |
| Currency-filtered v2 graph, expanded observations | 369 | 4 | Consecutive score support, 0.75/0.75 | **138 edges, eight TFs, 13 lipid features** |
| Currency-filtered v2 graph, expanded observations | 369 | 6 | Consecutive score support, 0.75/0.75 | 177 edges, nine TFs, same 13 lipid features |

The final comparison changes the graph, observation background and extraction
policy; it is not an otherwise identical numerical replacement of the original
network. The same frozen TF activities, RNA coherence signature, regulon,
joint factor orientation and safe compression policy are used throughout.

## Why the apparently stronger Cer network was rejected

The initial 101-edge result passed signed-edge and input-to-observation
connectivity checks, but every one of its six representative TF/Cer paths
passed through water. For example:

`Ar → Cyp1b1 → RHEA:48717 reaction → H2O → sphingomyelin-hydrolysis instance → Cer`.

The Cyp reaction transforms a hydroperoxy fatty acid to an oxo fatty acid and
releases water. That water node then connects to sphingomyelin hydrolysis.
The acyl oxylipid does not become the ceramide. A shared water co-product is
therefore providing a permissive graph connection between unrelated chemistries.
The selected water node has 26 incoming and 23 outgoing edges. Blocking its
outgoing substrate connectors removes every TF→Cer pair in that representative.

Those paths also have score levels `[4,4,3,2,1,0]`: the TF→Cyp edge is a
same-level edge, not the child relationship used to calculate that TF's score.
The initial stronger-cutoff networks are retained as diagnostics, with exact
provenance, but should not be presented as convincing Cer mechanisms.

The conservative sensitivity removes **37,666 exact positive
substrate-participation connectors** from explicitly declared generic cofactor
nodes to reaction instances. It blocks their use as reaction-to-reaction
bridges while retaining complete master reactions, their coefficients and
chemical participants. The panel covers water/protons, CoA, nucleotide and
redox cofactors, phosphate/diphosphate and SAM/SAH, including Mouse-GEM aliases
H+/Pi/PPi. Fatty-acid/acyl pools, measured lipids and second messengers such as
cAMP/cGMP are excluded from this removal policy. Four legitimate signed
chemical→protein regulatory edges from ATP/GDP are preserved.

This leaves a 293,366-edge working PKN. It is a conservative topological
sensitivity, not a claim that cofactors are biologically irrelevant. It cannot
evaluate specific cofactor-mediated regulation without appropriate evidence.
The full resource remains unchanged. Currency policy v1 and the alias-complete
v2 give exactly identical prepared graphs, compression maps and scores here;
the additional alias edges were already absent after preparation.

## Preferred network and reaction interpretation

Examples with direct recursive score support are:

- `Ppara → Plpp1 → PA-dephosphorylation instance → DG(16:0_22:6)`;
  scoring levels 3→2→1→0.
- `Nr2f1 → Rxra → Pld1 → PC-headgroup-hydrolysis instance → PA(18:1_18:2)`;
  levels 4→3→2→1→0.
- `Gata6 → Mapk1 → Pla2g4a → PC-acyl-hydrolysis instance → LPC(20:5)`;
  levels 4→3→2→1→0.

The primary network uses 54 distinct master reaction hypotheses represented by
60 catalyst-specific instances: 42 newly expanded PC-hydrolysis hypotheses
supported by directional RHEA:15802 and 12 baseline measured-species masters.
The new hypothesis layer contributes LPC production routes. These routes assume
compatible chains and enzyme species specificity; they do not establish exact
sn positions, double-bond locations, stereochemistry or substrate availability.

The 13 selected readouts are:

| Class | Measured features |
|---|---|
| DG | DG(16:0_22:6), DG(18:0_20:4), DG(18:0_22:6) |
| LPC | LPC(18:0/0:0) via frozen CHEBI:73858, LPC(18:0), LPC(19:1), LPC(19:1/0:0), LPC(20:4), LPC(20:5), LPC(22:4) |
| PA | PA(16:0_18:2), PA(16:1_18:2), PA(18:1_18:2) |

LPC(18:0) and LPC(18:0/0:0), and the two LPC(19:1) annotations, may overlap
chemically. They are separate original assay features with their own loadings;
no loading was copied to multiple nodes. Nevertheless, count **13 measured
features**, rather than 13 unequivocally resolved chemical species or
independent biochemical observations. Correlated readouts and overlapping
reaction hypotheses can still affect score interpretation after compression.

## Expected TFs

Esr1, Nfic and Dbp survive preparation and retain their identities through
compression and the final coherence-filtered PKN. In both final scoring depths,
their pre-upstream-sign-filter scores are positive, whereas their supplied TF
activities are negative:

| TF | Supplied activity | Inferred score | Final classification |
|---|---:|---:|---|
| Esr1 | −7.0887 | +0.1312 | Upstream sign conflict |
| Nfic | −4.3492 | +0.8996 | Upstream sign conflict |
| Dbp | −7.6609 | +0.2862 | Upstream sign conflict |

They are removed by MOON's final upstream sign filter, not pre-pruned. Nfic/Dbp
had weak negative scores in the unfiltered expanded graph, but that agreement
did not survive the cofactor policy. These outcomes do not refute their liver
biology: the tested prior, abundance-to-reaction projection and single factor
signature may inadequately represent their mechanisms. No signs or inputs were
changed to force expected TFs into the solution.

## Reaction-to-PKN rules

The detailed code/data explanation is in
[`cosmos_reaction_to_pkn_rules.md`](cosmos_reaction_to_pkn_rules.md). Two operations
are distinct:

1. **Chemical instantiation:** substitute structurally compatible measured
   chains into declared slots of a complete directional Rhea reaction, retaining
   the remaining substrates/products, coefficients, catalysts and evidence.
2. **Graph projection:** for `A + B --E--> C + D`, generate an enzyme–reaction–
   direction instance R and positive edges `E→R`, `A→R`, `B→R`, `R→C`, `R→D`.

The +1 signs mean enabling, participation or product formation in this
construction. They do not demonstrate substrate activation of the enzyme or
positive concentration control. No enzyme→substrate −1 depletion edges are
invented. Alternative catalysts have separate instances; the graph does not
require all alternative enzymes together. Unknown compartments remain distinct,
and new reverse directions are not invented.

Full reactions are chemical hyperedges. Ordinary MOON triplets do not enforce
joint co-substrate availability, stoichiometric balance, thermodynamics or flux.
The final displayed graph omits 236 of 356 master participant/catalyst connectors,
largely because co-substrates/co-products are not scored readouts or because of
terminalization/pruning. All selected catalyst links are present, but that is
not full reaction feasibility. Complete masters and omission bookkeeping remain
available for audit.

## Methods, software and reproducibility

Use the separate ancillary runner, not the manuscript entry point with modified
reference assertions. Installed versions are cosmosR 1.21.2 at
`b2ed195738479c4ad204cace37fd11ae2fd02bb8`, decoupleR 2.17.0, igraph 2.3.3 and
jsonlite 2.0.0; individual sessions and dependency hashes are saved. Python
preparation/export use the standard library.

Observed downstream lipids are made terminal before fitting, matching installed
`moon()` behavior. Adding observed lipids can therefore remove precursor
traversal as well as add evidence. After terminalization, directed forward and
reverse pruning are repeated until stable. Compression groups identical signed
children while protecting input and regulon identifiers. A complete original→
compressed map is verified. ULM scoring and RNA TF-target coherence are repeated
until no edges change; the final four-layer fit converges in four passes. RNA
is used for coherence, with no additional expression filter introduced here.

The final 43,212-edge prepared topology compresses to 34,578 edges. Its final
coherence-filtered scoring graph has 25,165 edges. The complete fixed mapping,
original-ID graph, scores, pruning/coherence ledgers and removed edges remain
together. Threshold extraction uses this same final graph, then retains only
edges with `source_level = target_level + 1`, strict score/sign agreement and
TF→terminal-observation connectivity.

An initial 75 valid higher-cutoff combinations was followed by a union grid of
15 primary and 14 secondary values,165 valid pairs per final depth. Empty results
and invalid pairs are retained. Eligible representatives require 0<edges<500,
at least two TFs and two lipid readouts; ranking prioritizes secondary then
primary cutoff, then readout/TF coverage. Both final depths have12 eligible
combinations. The final 0.75/0.75 cutoff is a within-grid choice, not a claim that
higher cutoffs can never yield any single-TF/single-readout route.

To reproduce in a new destination, from the repository root:

```bash
python3 scripts/moon/ancillary/prepare_expanded_factor1_mapping.py --validate-only
python3 scripts/moon/ancillary/prepare_expanded_factor1_mapping.py --output=/tmp/moon_expanded_replay/inputs

Rscript scripts/moon/ancillary/run_expanded_factor1_tf_to_lipid.R \
  --mode=expanded --mapping=/tmp/moon_expanded_replay/inputs/input_mapping.csv \
  --currency-policy=no_generic_substrate \
  --currency-edges=support/cosmos_mouse/moon_currency_substrate_edges_v2.csv \
  --layer-consistent --pruning-steps=6 --layers=4 \
  --output=/tmp/moon_expanded_replay/fit4

Rscript scripts/moon/ancillary/reextract_expanded_factor1_thresholds.R \
  --fit=/tmp/moon_expanded_replay/fit4 --layer-consistent \
  --output=/tmp/moon_expanded_replay/thresholds4

python3 scripts/moon/ancillary/export_expanded_factor1_network.py \
  --results=/tmp/moon_expanded_replay/thresholds4 \
  --pruned-pkn=/tmp/moon_expanded_replay/fit4/pkn_final_uncompressed.csv

Rscript scripts/moon/ancillary/plot_expanded_factor1_routes.R \
  /tmp/moon_expanded_replay/thresholds4/expanded_network_audit \
  /tmp/moon_expanded_replay/figures
```

All new nonempty fit/extraction/export directories are refused. For six layers,
repeat the fit with `--layers=6` and new destinations. For the full-graph control,
use `--mode=frozen_control --currency-policy=full`; omit layer-constrained
extraction for the original topological comparator. The original-cutoff
comparison is saved separately under each fit's `fixed_thresholds/`.

The currency list is reproduced with
`scripts/moon/ancillary/build_currency_traversal_policy.py --alias-complete`
and explicit new output filenames (see `--help`). Its v2 edge SHA256 is
`fae23006954fd2570dc20f30598ca2b6120e52ee5c9041a382a43d99e217e00e`.

`scripts/moon/ancillary/audit_expanded_factor1_rerun.R` independently verifies
the saved seven fits and seven threshold extractions: exact full-precision input
assignments, explicit currency-policy edge subtraction, signed graph/compression
projection, convergence, score IDs/levels, extraction signs/connectivity,
consecutive-layer policy and selection. `independent_validation.json` passed.
GraphML round-trips preserve the exact signed primary/sensitivity networks.

The result-root `selected_network.json` identifies the preferred and sensitivity
outputs; `inputs/mapping_audit.json` records mapping decisions. Every fit retains
its input RDS, script snapshot, configuration, software and dependency hashes.
Every extraction retains its fit reference, grid and source snapshot. Evidence
exports include exact signed provenance, masters, participant completeness,
cofactor audits and path tables. The baseline PKN/input/helper hashes still match
the frozen configuration.
