# Ranked knockout and knockdown targets for TF to lipid validation

The strongest first experiments are **Pld1 loss, Pld1/Pld2 double loss, and
Pla2g4a loss**. They test relatively short, interpretable chemical branches
with measured lipid readouts. The ten sets below cover individual enzymes,
alternative routes and pathway mediation. They are overlapping experimental
designs, not ten independent, validated TF mechanisms.

This proposal assumes primary mouse hepatocytes, consistent with the previous
validation experiments; the upcoming experimental system has not been
confirmed. Ranking combines chemical specificity, primary literature,
expression in the project data, interpretability and the additional information
from combinations. It is not a ranking by absolute MOON score. Directions are
conditional, local pathway predictions, not simulated knockout effects or
guaranteed steady-state concentration changes.

## Important source audit before choosing experiments

The 138-edge conservative network passed computational sign, score-layer and
cofactor checks. A source-level biological audit now identifies several
relationships that should not underpin an experiment-ready TF chain:

- **Ppara→Plpp1:** its sole cited paper studies **PPARγ**, not PPARα. The
  family identity does not support the selected Ppara-specific edge.
  [Huang et al., 2014](https://pubmed.ncbi.nlm.nih.gov/24955089/).
- **Cebpa→Pld1:** its sole citation supports **C/EBPβ**, and explicitly reports
  that C/EBPα did not activate or bind the tested PLD1 promoter elements.
  [Gao et al., 2010](https://pubmed.ncbi.nlm.nih.gov/19950202/).
- **Gata6→Mapk1:** its sole citation describes ERK phosphorylating GATA6,
  followed by GATA6 activating NOX1, rather than GATA6 transcriptionally
  activating MAPK1. [Adachi et al., 2008](https://pubmed.ncbi.nlm.nih.gov/18454176/).
- **Enpp2 as diacyl-PC→PA catalyst:** verified autotaxin biochemistry supports
  lysophosphatidylcholine→lysophosphatidic acid. It does not establish the
  selected diacyl-PC→diacyl-PA assignment.
  [Hausmann et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3064516/).
- **Rxra→Pld1:** one supporting reference is retracted. The unretracted
  reference supports vitamin-D/VDR–RXRα-dependent promoter regulation in human
  keratinocytes, so hepatocyte transfer remains conditional.
  [Unretracted primary study](https://pubmed.ncbi.nlm.nih.gov/17433303/),
  [retraction notice](https://pubmed.ncbi.nlm.nih.gov/39888011/).

These issues change the interpretation of the earlier attractive Ppara,
Gata6 and Rora/Cebpa routes. The original fitted graphs and scores remain frozen
for audit. The current proposal excludes unsupported edges and Enpp2 catalyst
instances from its path calculations; **it does not claim to be a newly refitted
MOON network**. Its numerical scores were obtained before this source audit.
The experimental case therefore rests chiefly on verified chemistry and
explicit candidate branches, rather than on those scores.

## Primary lipid panels

Use these exact measured features as prespecified endpoints:

| Panel | Features | Interpretation |
|---|---|---|
| **PA3** | PA(16:0_18:2), PA(16:1_18:2), PA(18:1_18:2) | Candidate PC→PA products |
| **DG2** | DG(16:0_22:6), DG(18:0_20:4) | Candidate PA→DG products of Plpp1 |
| **DG3** | DG2 plus DG(18:0_22:6) | Candidate PLC-associated DG readouts |
| **LPC18** | LPC(18:0), LPC(18:0/0:0) | Two assay annotations that may represent overlapping chemical identity; not independent replication |

LPC18 is a better primary cPLA2 readout than the full seven-LPC network panel:
release of sn-2 arachidonate from an 18:0/20:4 PC leaves an 18:0 LPC. In
contrast, many generated LPC(20:4), LPC(20:5) and LPC(22:4) hypotheses retain
the PUFA while releasing the other chain. Their structural compatibility does
not establish PLA2G4A specificity. Treat them as exploratory tests, not required
pass/fail endpoints. [Primary molecular-specificity study](https://pubmed.ncbi.nlm.nih.gov/7496992/).

## Ten sets ranked by relevance

Arrows denote the local pathway expectation after functional target loss.
Every direction remains conditional on the branch contributing to the measured
pool. Combination effects are not assumed additive or synergistic.

| Rank | KO/KD set | Primary directional expectation | Why it merits this position |
|---:|---|---|---|
| **1** | **Pld1** | **PA3 ↓** | Best direct chemical anchor: short PC→PA branch, expressed in historical cultures, and hepatocyte-specific knockout precedent supports reduced PA. |
| **2** | **Pld1 + Pld2** | **PA3 ↓** | Tests an alternative PLD isoform that could mask Pld1 loss. Pld2 is present in the six-layer sensitivity and historical cultures; compensation or a larger effect is not guaranteed. |
| **3** | **Pla2g4a** | **LPC18 ↓** | Tests the expanded PC-remodelling hypothesis using the most chemically credible LPC product. Other LPC changes are secondary specificity questions. |
| **4** | **Pld1 + Pla2g4a** | **PA3 ↓ and LPC18 ↓** | Tests two distinct PC-cleavage outputs together. This should be interpreted against both singles, without assuming synergy. |
| **5** | **Plpp1** | **DG2 ↓**; matched PA accumulation is secondary | Direct test of the proposed PA→DG conversion, but the enzyme's extracellular/luminal active site and the PKN's unknown compartments limit confidence for a bulk intracellular lipid phenotype. |
| **6** | **Mapk1 + Pla2g4a** | **LPC18 ↓** | Tests signalling-to-lipase mediation through the independently supported ERK–cPLA2 link. Compare Mapk1 alone, Pla2g4a alone and the pair; ERK1 compensation and broad kinase effects limit interpretation. |
| **7** | **Plcb1 + Plcg1** | **DG3 ↓**, conditional on the relevant signalling pool | Tests two alternative PLC routes. Their usual physiological substrate is a phosphoinositide signalling pool, whereas the resource instantiates unphosphorylated PI species. Exact bulk PI→DG specificity is consequently less secure. |
| **8** | **Pld1 + Pld2 + Pla2g4a** | **PA3 ↓ and LPC18 ↓** | Later-stage integration of the PLD redundancy test with the remodelling branch. This broadens pathway blockade but cannot establish a three-way interaction without the requisite factorial controls. |
| **9** | **Rxra** | **PA3 ↓ only if RXRA loss first reduces PLD1** | The remaining identity-supported TF/intermediate-to-enzyme test. It is a conditional hepatocyte transfer of a ligand/VDR-dependent mechanism, not an established liver prediction. |
| **10** | **Rxra + Pld1** | PA3 ↓ within the branch; test loss of the RXRA-associated effect after PLD1 loss | A conditional mediation experiment following rank 9. The pair should be compared with both singles; similar effects could support pathway occlusion but are not sufficient proof of epistasis. |

Ranks 1–4 are the most straightforward chemical-anchor experiments. Ranks 5–8
test increasingly conditional pool accessibility, signalling or combined
pathway hypotheses. Ranks 9–10 should be reserved until the RXRA→PLD1 response
is demonstrated in the intended hepatocytes. Historical cultured controls have
low Vdr expression; a lipidomic phenotype alone would not establish this
specific TF mechanism.

## Interpretation of each module

**PLD module, ranks 1–2.** Hepatocyte-specific Pld1 loss has reduced PA/LPA in
published mouse liver lipidomics, while either PLD1 or PLD2 overexpression can
raise PA in primary mouse hepatocytes. This supports the class-level direction;
the three exact acyl compositions are our hypotheses. Their matched PC
precursors are useful secondary readouts, with possible accumulation after
reduced hydrolysis. Neither a PC increase nor a stronger double-loss phenotype
is guaranteed: other PA-producing enzymes, compartmental pools and turnover
remain. [Hepatocyte Pld1 study](https://pubmed.ncbi.nlm.nih.gov/37138676/),
[primary-hepatocyte PLD1/2 study](https://pubmed.ncbi.nlm.nih.gov/25512376/).

**PLA2 module, ranks 3–4 and 6.** The primary expectation is reduced production
of the LPC18 pool from arachidonoyl PCs. Parent PCs, liberated arachidonate and
relevant downstream products can distinguish this mechanism if the platform
captures them; their changes are secondary chemical expectations rather than
MOON-measured predictions. MAP kinase phosphorylation can activate cPLA2, but
MAPK1 perturbation can affect many other functions and be bypassed by MAPK3.
A Mapk1/Pla2g4a pair should be used for a mediation question, rather than merely
to obtain a larger lipid effect. [cPLA2 specificity](https://pubmed.ncbi.nlm.nih.gov/7496992/),
[MAP kinase activation of cPLA2](https://pubmed.ncbi.nlm.nih.gov/8381049/),
[ERK2-specific Reactome reaction](https://reactome.org/content/detail/R-HSA-111898).

**PLPP1 module, rank 5.** A local decrease of PA dephosphorylation predicts DG2
decrease and potentially matched PA increase. This chemical inference adds a
substrate-consumption expectation that the all-positive PKN itself does not
encode. It requires access to the relevant PA pool. Hepatocyte PLPP1 literature
more directly supports LPA regulation, and PLPP enzymes have extracellular or
luminal catalytic topology. LPA is therefore a useful context readout, while a
negative DG result cannot invalidate all PLPP1 biology.
[PLPP1 hepatocyte study](https://pubmed.ncbi.nlm.nih.gov/28273928/),
[catalytic topology study](https://pubmed.ncbi.nlm.nih.gov/18215144/).

**PLC module, rank 7.** PI hydrolysis is chemically possible in purified-enzyme
assays, but physiological PLC signalling commonly uses PIP2. An observed
same-chain DG decrease would support a lipid-production role; it would not
identify the unmeasured PIP2 precursor or establish that bulk PI was consumed.
Do not require accumulation of the model's measured PI species as a validation
criterion. Other PLC isoforms can contribute.
[Purified PLCγ1 study](https://pubmed.ncbi.nlm.nih.gov/9915811/),
[PLCβ1 substrate study](https://pubs.acs.org/doi/10.1021/bi961606w).

**RXRA module, ranks 9–10.** RXRA is a TF intermediary in the selected graph,
rather than one of its eight supplied upstream TF activities. Its reliable
PLD1 transcription evidence involves VDR and vitamin D in HaCaT cells. Before
assigning the PA prediction in hepatocytes, establish the receptor context and
whether RXRA loss reduces PLD1 expression or functional availability. A broad
RXRA lipid phenotype without this mediator response would not validate the
proposed route. The 2020 cervical-cancer paper cited in the PKN was retracted
in 2025 and is excluded from this proposal's supporting evidence.
[Unretracted promoter study](https://pubmed.ncbi.nlm.nih.gov/17433303/),
[retraction notice](https://pubmed.ncbi.nlm.nih.gov/39888011/).

## Expression and cell-context evidence

Actual integer RNA counts were used to check expression, separately from
Factor 1 loadings and MOON scores. The audit covers 43 FCG hepatocyte
pseudobulks and four historical negative-control cultured-hepatocyte RNA
libraries. Approximate mean culture CPM values are:

| Gene | Mean culture CPM | Implication |
|---|---:|---|
| Pld1 | 38.5 | Detected in the prior culture system |
| Pld2 | 46.6 | Supports feasibility of the isoform-pair test |
| Pla2g4a | 8.1 | Detected; activity still depends on activation context |
| Plpp1 | 14.8 | Detected; compartment accessibility remains unresolved |
| Mapk1 | 194 | Detected; broad signalling perturbation |
| Plcb1 | 15.9 | Detected; relevant receptor/signalling context required |
| Plcg1 | 93.1 | Detected; relevant signalling pool required |
| Rxra | 125 | Detected; does not establish ligand-dependent activity |
| Vdr | 0.46 | Important limitation for ranks 9–10 |

These are RNA library measurements, not protein/activity estimates, single-cell
detection fractions or proof of independence between biological donors. Confirm
expression and target engagement in the new experimental preparations.

Several longer TF routes have additional context limitations: Btk is low in
historical cultures and immune-enriched, Kdr is endothelial-enriched, and Ntrk2
is low in the prior culture system with unresolved functional isoform. Nr2f1,
Tbx2 and Hoxa10 are weak or absent in those culture libraries. This disfavors
prioritizing them solely because they occur in a network scored without an
additional expression filter.

Ppara remains an interesting independent hepatic biology control, but the
present Ppara→Plpp1 chain is not source-validated. Its expression was much
higher in fresh FCG hepatocyte nuclear pseudobulk than in historical cultures
(approximately 1,557 versus 3.6 CPM). These different RNA preparations do not
provide an absolute fold-change or protein-activity estimate. Broad hepatic knockout lipid effects can
also differ from a small proposed local DG branch.
[Hepatocyte Ppara knockout study](https://pubmed.ncbi.nlm.nih.gov/32300166/).

## Combinations and experimental comparisons

The ten ranked sets require **four additional constituent-single controls**
if every combination is to be interpretable: Pld2, Mapk1, Plcb1 and Plcg1.
Together with the four ranked singles, this gives eight single-target arms and
ten ranked designs, or **14 distinct target-containing conditions**, before
non-targeting controls, rescue conditions or replication. A complete test of
three-way interaction at rank 8 also needs the Pld2/Pla2g4a pair, which is not
in the ranked ten. That would increase the target-containing conditions to 15.

A practical first wave is Pld1, Pld2, their pair, and Pla2g4a with matched
controls. It tests the strongest chemistry and the main compensation question.
Plpp1 and the signalling/TF-mediated arms can follow once the corresponding
pool or mediator response is established. The ranking orders hypotheses;
biological dependencies can require running a lower-ranked constituent single
before its higher-ranked combination.

For serial pairs, compare the pair with both singles and measure the mediator.
An effect similar to enzyme loss alone is compatible with pathway occlusion,
but cannot establish mediation by itself. For parallel isoforms or branches,
an additional effect can support residual pathway capacity, while no additional
effect can reflect inactive branches, efficient single blockade or compensation.
The PKN does not predict additivity, effect sizes or synergy.

Treat KO and partial knockdown as different functional perturbations. Transcript
loss alone does not establish loss of protein or enzyme activity. Interpret
lipidomics alongside target engagement, hepatocyte identity and cell-health
evidence, using matched perturbation burden across combination and single arms.
Use KO/KD-versus-control effects within the experimental system; do not infer
their direction from whether a feature has a positive or negative Factor 1
loading. A terminal product change alone supports the local branch, while a
TF→enzyme→lipid claim also requires the intervening regulatory response.

Prespecify the readout groups and primary versus exploratory endpoints above.
Retain the full lipidomic profile, but distinguish planned direction tests from
post hoc discoveries and use appropriate multiplicity control. Bulk lipid
abundance is a net pool: synthesis, consumption, transport, feedback and
reacylation can reverse or mask the local expectation. A disagreement is useful
evidence against the proposed branch dominating that condition; it is not
automatically a refutation of the underlying enzyme reaction.

## Provenance and reproduction

The supporting directory is
`results/causal_modelling/factor1_ko_validation_targets/`. It contains:

- `prioritization/ranked_target_sets.csv`: the ten reviewed sets and rationale.
- `prioritization/expected_lipid_directions.csv`: one row per set/readout,
  including evidence layer, direction and caveat.
- `prioritization/audited_paths.csv` and `path_audit.json`: exact signed paths
  after the declared source exclusions, including reaction-instance IDs.
- `prioritization/excluded_path_edges.csv`, `review_snapshot.json` and
  `validation_summary.json`: biological-screen decisions and program checks.
- `expression_audit/`: expression, cell-type/sex and individual culture-library
  tables, actual-count checks, source hashes and R session.

The reviewed choices and primary-source ledger are
`support/cosmos_mouse/factor1_ko_validation_review.json` and
`support/cosmos_mouse/ko_validation_literature_sources.json`.

From the repository root, using fresh output directories:

```bash
Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R --validate-only
Rscript scripts/moon/ancillary/audit_factor1_ko_expression.R \
  --output=/tmp/factor1_ko_review/expression
python3 scripts/moon/ancillary/propose_factor1_ko_validation.py \
  --expression=/tmp/factor1_ko_review/expression/gene_expression.csv \
  --output=/tmp/factor1_ko_review/prioritization
```

The direction calculation is **functional loss (−1) × the product of signed
path edges**, independent of the sign of the fitted MOON score. The scripts
validate unique biological target sets, exact node identities, consecutive
score layers, sign agreement, candidate-readout correspondence and supported
directions. Ranking is explicit expert judgement in the configuration, not an
automatically estimated probability of success. No new MOON fitting,
intervention simulation or PKN mutation was performed for this proposal.
