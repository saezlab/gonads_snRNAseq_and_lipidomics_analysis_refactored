# Ranked knockout and knockdown targets for TF to lipid validation

The revised panel **alternates five enzyme-only sets with five TF-containing
sets**. Every pair tests a plausible alternative route that could preserve
the same lipid output after a single perturbation. It contains six singles and
four pairs; no current triplet has a sufficiently specific compensation
rationale to merit inclusion.

This is a balanced experimental allocation, not a global ordering of evidence
strength: the TF hypothesis at rank 2 is more conditional than the enzyme pair
at rank 3. Ranking within each category combines chemical specificity, primary
evidence, expression and interpretability. The assumed system is primary mouse
hepatocytes, consistent with previous validations; the forthcoming experimental
system has not been confirmed.

An important limit concerns the TF arm. RXRA is an intermediary in the frozen
selected graph. **CEBPB and PPARG are source-corrected, unscored follow-up
hypotheses**, rather than selected MOON explanations of Factor 1. Their
regulatory actions motivate perturbation tests, but their Factor 1 activity
signs do not support the corresponding activating branches as current
sign-coherent explanations. Five TF-containing sets test three unique TFs
and two enzyme modules, not five independent fitted TF mechanisms.

## Source audit and evidence boundaries

The frozen conservative four-layer network contains 138 edges and 94 nodes;
the six-layer sensitivity contains 177 edges and 117 nodes. Computational
sign/layer consistency does not guarantee biological source accuracy:

- **Ppara→Plpp1:** the source supports **PPARG**, whereas PPARA and PPARD did
  not induce LPP1 in that endothelial-cell experiment. The Ppara edge is
  excluded; PPARG→PLPP1 is a separately labelled candidate extension.
  [Huang et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4052932/).
- **Cebpa→Pld1:** the source supports **CEBPB**, explicitly not CEBPA, in mouse
  3T3-L1 adipocyte differentiation. The Cebpa edge is excluded; CEBPB→PLD1 is
  a separate candidate extension.
  [Gao et al., 2010](https://pubmed.ncbi.nlm.nih.gov/19950202/).
- **Gata6→Mapk1:** the cited mechanism is ERK phosphorylation of GATA6 and
  GATA6 activation of NOX1, rather than GATA6 transcriptional activation of MAPK1.
  [Adachi et al., 2008](https://pubmed.ncbi.nlm.nih.gov/18454176/).
- **Enpp2 as diacyl-PC→PA catalyst:** verified autotaxin chemistry supports
  lysophosphatidylcholine→lysophosphatidic acid, not the selected diacyl-PC→PA
  assignment. These catalyst instances are excluded.
  [Hausmann et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3064516/).
- **Rxra→Pld1:** reliable evidence involves vitamin-D/VDR–RXRα promoter
  regulation in human keratinocytes. Another supporting reference was retracted.
  Hepatocyte transfer remains conditional.
  [Unretracted study](https://pubmed.ncbi.nlm.nih.gov/17433303/),
  [retraction notice](https://pubmed.ncbi.nlm.nih.gov/39888011/).

Original graphs, PKNs and scores remain frozen. Path calculations remove the
declared unsupported edges in an audit copy. For Cebpb and Pparg, a separately
labelled literature-derived regulatory edge is prefixed to existing
enzyme-to-lipid paths. These added edges receive **no MOON score or score
layer**; no Cebpa/Ppara score is transferred to a corrected paralog. This
revision neither refits the PKN nor establishes that these extended paths
explain the observed Factor 1 pattern.

## Prespecified lipid readouts

| Panel | Exact measured features | Interpretation |
|---|---|---|
| **PA3** | PA(16:0_18:2), PA(16:1_18:2), PA(18:1_18:2) | Candidate PC→PA products |
| **DG2** | DG(16:0_22:6), DG(18:0_20:4) | Candidate PLPP1-associated PA→DG products |
| **DG3** | DG2 plus DG(18:0_22:6) | Candidate PLC-associated DG readouts |
| **LPC18** | LPC(18:0), LPC(18:0/0:0) | Possibly overlapping chemical annotations; not independent replication |

LPC18 is the primary cPLA2 panel because hydrolysis of an sn-2 arachidonoyl
18:0/20:4 PC leaves an 18:0 LPC. Other generated LPC rules retaining the PUFA
remain exploratory specificity hypotheses.
[Primary specificity study](https://pubmed.ncbi.nlm.nih.gov/7496992/).

Arrows below describe **conditional local production/control effects** after
functional target loss. They do not predict effect sizes, synergy or guaranteed
changes in bulk steady-state abundance. Directions follow regulatory action
and chemistry rather than reversal of a Factor 1 loading or MOON score.

## Revised alternating ranking

| Rank | KO/KD set | Category | Primary expectation | Rationale and necessary qualification |
|---:|---|---|---|---|
| **1** | **Pld1** | Enzyme single | **PA3 ↓** | Strongest direct PC→PA anchor with hepatic perturbation evidence; exact acyl compositions remain hypotheses. |
| **2** | **Rxra** | TF single | **PA3 ↓ if PLD1 decreases** | Tests the contextual RXRA→PLD1 branch. Establish receptor context and PLD1 response; Vdr is weakly expressed in historical cultures. |
| **3** | **Pld1 + Pld2** | Enzyme pair | **PA3 ↓** | Block parallel PLD isoforms: PLD2 could retain PC→PA production after Pld1 loss. Compare both singles; compensation is anticipated rather than established in these hepatocytes. |
| **4** | **Cebpb** | TF single | **PA3 ↓ if PLD1 decreases** | Tests source-corrected CEBPB→PLD1. Requires transfer from the adipocyte context and net PLD1 suppression. Not a selected, sign-coherent MOON TF route. |
| **5** | **Pla2g4a** | Enzyme single | **LPC18 ↓** | Tests PC acyl hydrolysis using the most credible retained-chain readout. |
| **6** | **Cebpb + Pld2** | TF + enzyme pair | **PA3 ↓, conditional on the CEBPB→PLD1 response** | If Cebpb loss lowers PLD1, remaining PLD2 could preserve PA output. Pld2 loss blocks this proposed backup. |
| **7** | **Plpp1** | Enzyme single | **DG2 ↓** | Direct anchor for PA→DG; requires enzyme access to the relevant PA pool. |
| **8** | **Pparg** | TF single | **DG2 ↓ if PLPP1 decreases and accesses that pool** | Tests source-corrected PPARG→PLPP1. Evidence is endothelial; Pparg is absent from the frozen PKN and is not a selected MOON prediction. |
| **9** | **Plcb1 + Plcg1** | Enzyme pair | **DG3 ↓, conditional on signalling-pool contribution** | Block alternative PLC inputs capable of DAG production. Activity in the same accessible pool and exact-chain contribution remain unproven. |
| **10** | **Rxra + Pld2** | TF + enzyme pair | **PA3 ↓, conditional on the RXRA→PLD1 response** | If Rxra loss lowers PLD1, PLD2 could maintain PA. The pair tests this bypass and requires both constituent singles. |

Pld2 is supported by the six-layer sensitivity rather than the primary
four-layer graph. All Pld2-containing designs use six-layer enzyme-to-lipid
paths. Corrected TF prefixes remain unscored in either setting.

## Compensation rationale for every pair

**Pld1/Pld2** is the strongest design. Either enzyme can raise PA in primary
mouse hepatocytes, and overlapping PA species profiles have been demonstrated
in other mammalian cells. This supplies a specific alternative producer of the
same output. It does not prove adaptive PLD2 upregulation after Pld1 loss or
equal access to every measured PA pool.
[Primary-hepatocyte study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4319019/),
[molecular-species comparison](https://pubmed.ncbi.nlm.nih.gov/11736663/).
Hepatocyte-specific Pld1 loss also has reduced PA/LPA in published mouse liver
lipidomics; the exact three chains here remain our hypotheses.
[Hepatic perturbation study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10149370/).

A mouse forebrain study provides a functional precedent: single PLD KOs retained
substantial PLD-dependent transphosphatidylation, whereas double KO nearly
abolished it without expression upregulation of the remaining isoform.
Residual capacity can therefore mask a single perturbation without adaptive
RNA/protein induction. This is brain evidence, not proof of the same response
in our hepatocytes.
[Santa-Marinha et al., 2020](https://www.med.uminho.pt/en/Research/Selected%20Publication/Santa-Marinha-et-al-2020.pdf).

**Cebpb/Pld2 and Rxra/Pld2** apply the same logic upstream: reduced
TF-dependent PLD1 production might leave PLD2-mediated PA production intact.
To interpret these pairs as preventing compensation, establish that the TF
single lowers functional PLD1 and PLD2 remains available. If the TF already
suppresses both PLDs, PLD2 is not the proposed retained backup. Measure the
mediator and both enzyme proteins/activities alongside lipids. A pair's PA
decrease alone could be driven entirely by Pld2 loss and would not validate
the TF link.

**Plcb1/Plcg1** tests alternative receptor-linked signalling inputs converging
on DAG. This is a lower-confidence backup hypothesis: physiological PLC
signalling principally uses phosphoinositides such as PIP2, whereas the PKN
instantiates measured unphosphorylated PI. Both routes must operate in the
experimental context; other PLC isoforms remain. Do not require measured
bulk PI accumulation as a validation criterion.
[PLCβ1 substrate study](https://pubs.acs.org/doi/10.1021/bi961606w),
[purified PLCγ1 study](https://pubmed.ncbi.nlm.nih.gov/9915811/).

The earlier Pld1/Pla2g4a combination joined different outputs, while
Rxra/Pld1 and Mapk1/Pla2g4a were serial mediation designs. They are replaced
in this compensation-focused ranking. A stronger pair effect or synergy is
not required: the question is whether a remaining route preserves production
after a single intervention.

**No triplet is currently recommended.** Independent RXRA and CEBPB promoter
effects do not show that these TFs substitute for each other in hepatocytes.
Adding both TFs to Pld2 would assume an unverified regulatory backup. A later
triplet should follow evidence that a third route retains the same product
after the corresponding pair, with single and constituent-pair comparisons.

## Interpretation of the TF arm

**RXRA, ranks 2 and 10.** Its cited mechanism requires vitamin-D/VDR–RXRα
regulation in HaCaT keratinocytes. Transfer requires the appropriate
partner/ligand context and a demonstrated PLD1 response. Low historical Vdr
expression is a feasibility limitation. RXRA controls broad nuclear receptor
physiology, so a general lipid phenotype does not establish this specific route.
The retracted cervical-cancer source is excluded.
[Promoter study](https://pubmed.ncbi.nlm.nih.gov/17433303/),
[hepatic RXRA study](https://pmc.ncbi.nlm.nih.gov/articles/PMC85811/).

**CEBPB, ranks 4 and 6.** The source supports PLD1 transcription in early
mouse adipocyte differentiation, with PLD1 rather than PLD2 induction in that
context. This motivates the alternative-PLD2 design but does not establish
hepatocyte selectivity. Whole-gene loss affects CEBPB isoforms with different
regulatory actions, making net PLD1 response a necessary measurement.
[Primary promoter study](https://pubmed.ncbi.nlm.nih.gov/19950202/).

**PPARG, rank 8.** The corrected source supports PLPP1 transcription in human
endothelial cells; mouse hepatocyte transfer must be demonstrated. Even after
PLPP1 suppression, DG2 depends on accessible PA. PLPP enzymes have
extracellular/luminal catalytic topology, and hepatic PLPP1 literature more
directly supports LPA regulation. Measure LPA as a context endpoint if available.
[Regulatory study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4052932/),
[PLPP topology](https://pubmed.ncbi.nlm.nih.gov/18215144/),
[hepatocyte PLPP1 study](https://pubmed.ncbi.nlm.nih.gov/28273928/).

Factor 1 evidence warrants restraint. Cebpb's footprint is negative (−2.22),
whereas PA3 measurements and inferred PLD1 activity are positive. Pparg's
footprint is weakly positive (+1.20), whereas DG2 and inferred PLPP1 activity
are negative. Their added activating routes are therefore **not sign-coherent
explanations of the frozen Factor 1 measurements**. Rxra's weak footprint
(−0.75) differs from its positive intermediate MOON score; it was not an
upstream seed. These observations motivate biochemical follow-up and
falsification tests rather than validation of an already fitted causal
explanation. Refitting a corrected PKN would be separate analysis.

## Expression, controls and experimental sequence

Actual-count checks cover 43 FCG hepatocyte pseudobulks and four historical
negative-control cultured-hepatocyte RNA libraries. Mean culture CPM values
are Pld1 38.5, Pld2 46.6, Rxra 125.0, Cebpb 31.3, Pla2g4a 8.1, Plpp1 14.8,
Pparg 25.2, Plcb1 15.9 and Plcg1 93.1; all nine are detected. Vdr is low
(0.46 CPM). Detection is feasibility evidence, not protein/activity evidence
or proof of donor independence. Confirm engagement in the new preparations.

The panel has **nine unique genes, six singles and four pairs**. Add constituent
singles **Pld2, Plcb1 and Plcg1**: the ten designs then require **13 distinct
target-containing conditions**, before non-targeting controls, rescues and
replication. Every pair requires its two single arms and matched perturbation
burden. Measure the proposed backup after the primary single; an RNA increase
is not required for residual catalytic capacity to matter.

A balanced first wave can test Pld1 and Pld2 chemistry in parallel with Rxra
and Cebpb singles, with PLD1/PLD2 expression, protein or activity as early
readouts. Advance TF/Pld2 pairs when the proposed mediator response is observed.
Pla2g4a tests the separate LPC branch; Plpp1/Pparg and the PLC pair require
their compartment/signalling checks. Thus TF hypotheses receive an early test,
while combinations remain grounded in observed biology.

A TF→enzyme→lipid claim requires TF engagement, the predicted mediator response
and a compatible lipid response. Downstream enzyme rescue under TF loss adds
mediation evidence, with expression/activity controls for the rescue itself.
Monitor hepatocyte identity and cell health to distinguish a local branch
from general state changes. KO and partial KD differ functionally; transcript
loss alone does not establish loss of enzyme activity.

Prespecify panels and contrasts, distinguish them from exploratory lipidome-wide
results and control multiplicity. Bulk abundance integrates production,
consumption, transport, compartmentalization and feedback. A local production
decrease need not lower the total pool. Disagreement tests whether the branch
dominates that condition rather than automatically refuting its chemistry.

## Provenance and reproduction

Current exports are in
`results/causal_modelling/factor1_ko_validation_targets/alternating_prioritization/`:

- `ranked_target_sets.csv`: ranks, categories, compensation rationales and gates.
- `expected_lipid_directions.csv`: prespecified set/readout expectations.
- `audited_paths.csv` and `path_audit.json`: signed paths, separating frozen
  support from source-corrected unscored prefixes.
- `source_corrected_hypothesis_edges.csv`: literature-derived candidate edges.
- `excluded_path_edges.csv`, `review_snapshot.json`, `validation_summary.json`
  and source manifests: exclusions, choices and reproducibility checks.

Existing `expression_audit/` is reused. Initial sources and outputs are preserved
in `initial_panel_20261013/` and `prioritization/`; they are superseded by this
ranking. Historical MOON graphs and scores are unchanged.

Choices and references reside in
`support/cosmos_mouse/factor1_ko_validation_review.json` and
`support/cosmos_mouse/ko_validation_literature_sources.json`. From the repository
root, using a fresh destination:

```bash
python3 scripts/moon/ancillary/propose_factor1_ko_validation.py --validate-only
python3 scripts/moon/ancillary/propose_factor1_ko_validation.py \
  --output=/tmp/factor1_ko_alternating_review
```

Repeat the actual-count audit with
`Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R --output=NEW_DIR`
and supply its `gene_expression.csv` with `--expression`. Direction is functional
loss (−1) times the signed-path product under a local monotone branch assumption.
Validation checks alternation, unique sets, shared pair endpoints, absence of
serial co-targets, frozen path signs/layers and separate unscored provenance
for corrected TF prefixes. Ranking remains explicit biological judgement,
not an estimated success probability.
