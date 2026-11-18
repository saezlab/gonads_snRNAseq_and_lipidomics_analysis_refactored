# Workflow 05c - Lipid overlap audit and pooled Figure 3 conservation

**Compatibility phase documentation — 9 October 2026:** the active complete
work package is [12](12_wdsd_network_comparison.md), entered through
`scripts/workflows/12_wdsd_network_comparison.R`. This former entry is a thin wrapper for its
pooled/annotation phase; orchestration now lives in `scripts/core/wp2_workflow.R`.
Inputs are repository-local, with byte-preserving import provenance in
`data/wp1_wp2_input_migration_manifest.csv`. The recorded production results
below remain unchanged. The combined WP2 report now presents pooled results
first, then zonation/diet, as authorized on 9 October. Earlier separate-report
instructions below document the original deliveries.

**Entry point:** `scripts/workflows/05c_wdsd_lipid_overlap_and_pooled_network.R`.
**Scope:** The 7 October 2026 user follow-up to WP2: explain the 752/1,254
strict molecular overlap and compare the non-zonated, sex-pooled Figure 3
gene-lipid network between the original FCG cohort and WDSD standard diet.
Preserve all completed WP1/WP2 and original FCG outputs.

Read Workflow 04, 05b and the cross-workflow overview before modifying the
entry point. Write new results only to
`results/wdsd_network_comparison/followup/` or a fresh user-specified destination.
Save a separate Markdown/PDF addendum in Meeting_Summaries after numeric and
visual QA. Do not silently replace the prior WP2 report or its strict matches.

## Inputs and Figure 3 definition

- Original whole-hepatocyte scaled RNA from Workflow 01's three CSV parts,
  original metadata/group labels and log2 lipid intensities.
- Figure 3 species-level archive from Workflow 04, its upstream sparse
  correlation matrix and original selected gene/lipid identities.
- Current `docs/reference_figures/FCGL_Fig3.ai` is PDF-readable and identifies the
  sparse species-level/top-three-gene branch. Its "FULL SPECIES-LEVEL NETWORK"
  panel is not the separate permissive top-33%-tail network.
- Updated WDSD whole-hepatocyte pseudobulks, CF453 lipid export and the
  validated 28-mouse intake. Pooled SD contains 14 distinct animals, seven
  recorded females and seven males. No zonation profiles enter this comparison.
- Raw original annotation metadata in
  `data/metabolomics/fcg/lipidomics_fcg.csv`, when available, for
  annotation aliases, adduct/mass/retention and structural-identifier audits.

Reconstruct and independently replay all original pooled correlations before
comparison. Retain the 43 original FCG mouse pairings and fixed MOFA feature
selection; do not refit MOFA or select pairs on new observations. Extract only
gene-lipid edges from the displayed species branch; do not compare compressed
lipid-class edges with individual species. Use the full original sparse
gene-by-lipid coefficient matrix as a clearly labeled complementary universe.
The much larger permissive Workflow 04 full network is an inventory reference.

## Lipid overlap audit

Reproduce 752 strict matches against 1,254 original network measurements and
separately count all supplied versus retained WDSD rows. Classify every missing
label with complete candidate ledgers, class-level counts and examples.
Distinguish unique identity-preserving changes from structural-resolution
compatibility, alias ambiguity and no recorded compatible candidate. A missing
exact annotation is not proof that the lipid was experimentally undetected.

Preserve slash positions, ether/plasmalogen/oxygenation, chain composition and
CL-specific pairing. Underscore sorting is not globally sufficient for CL.
Raw pipe aliases can be sum compositions shared by different molecular species.
Full InChIKey/SMILES equality cannot override contradictory measured-label
structures or unreliable library provenance. Check bidirectional uniqueness,
including collisions with already strict matched rows. Never silently merge or
sum multiple measured/adduct features to inflate overlap.

Keep strict identity primary. Any compatible-alias mapping is an explicitly
conditional sensitivity at its reported resolution, with mass/adduct support
and collision exclusions; it does not upgrade candidates to confirmed identity.
Use primary nomenclature literature when a chemical matching rule is uncertain.

## Pooled conservation and sex structure

Primary summaries: original-versus-SD coefficient-vector Pearson r, native
intercept ULM slope t and nominal ULM p. Displayed pairs remain frozen, including
weak or reversing SD correlations. Complementary sparse-matrix and predefined
Cyp2b9/Cyp2b13/Cyp2c38/Cyp2d9 gene summaries remain separately labeled.
Eligibility requires finite variable paired measurements; save all exclusions.

Pooled correlations include between-sex mean differences. Report two distinct
animal-level null checks so pooled replication is not confused with evidence
beyond shared sex effects:

1. Unrestricted 14-mouse RNA-lipid pairing permutations and independent
   empirical-marginal null bootstrap test broad absence of pairing information.
   Shared sex differences can alone yield positive reference alignment.
2. Sex-stratified pairings and independent marginal resampling retain the
   seven female/seven male margins, assessing alignment beyond pairing expected
   from that shared sex structure. The null can be centered above zero; report
   upper-tail probabilities for positive excess alignment explicitly as one-sided.

Use 19,999 Monte Carlo pairing draws by default (7! x 7! stratified pairings
are not exhaustively enumerated), identity saved separately as draw zero,
and `(1+exceedances)/(1+valid draws)`. Use 4,999 null-bootstrap draws, 999
stratified paired stability draws and original FCG group-stratified reference
resampling. Save upper/lower/absolute/equal-tail diagnostics with their definitions.
All genes and lipid features share each whole-mouse draw. Undefined edges
invalidate their frozen-vector summary; no feature deletion within draws.

A within-sex-centered correlation-network sensitivity removes each cohort's
sex means before calculating correlations, on the same fixed pairs. It is a
descriptive residual-network comparison, not a replay of Figure 3 coefficients;
report r/ULM/nominal p without assigning unimplemented residual-null tests.

Nominal edge p-values are descriptive under dependent edges. Neither pooled
replication nor conditional tests identify causal regulation. Preserve
uncertainty about age/background, normalization/imputation and nucleus QC.
No multiplicity-adjusted discovery, confidence intervals or reversed signed
score is requested. Do not assume the hoped-for pooled conservation outcome.

## Execution and verification

Source `scripts/support_functions.R` first and use the existing local R library.
New helpers cover molecular audit, pooled reference and pooled statistics;
run every new/modified R script immediately and test actual inputs. The entry
point records seeds/settings, mouse identities, frozen edge and molecular
ledgers, all resampling indices/draws, reference replay, code/input SHA256 and
session versions. Scientific figures and the self-contained report must report
native ULM and both clearly distinguished animal-level nulls.

Tests independently verify exact archived Pearson values, lipid-resolution and
ambiguity examples, animal pairing/sex strata, ULM/OLS equivalence, identity and
Monte Carlo tails, shared complete feature draws and centered-network behavior.
Update this plan with actual counts, outcomes and unresolved limitations after
the completed run; register the extension in the repository query map/meta.

## Completed run: 7 October 2026

Default settings completed: seed 20261007, 19,999 random pairings per null,
4,999 independent-marginal bootstrap draws per null, and 999 paired SD/reference
stability draws. Four analyses use identical mouse index banks: strict,
conditional alias compatibility, QC-CV <=30%, and supplied Fill >=0.5.

### Annotation audit

- Original raw export 1,386 rows; original network universe 1,254. New export
  1,240 rows, 1,177 retained after 14 standards and 49 excluded-class rows.
- The original pre-existing cleaning omits 132 raw rows: 111 canonical classes
  in the repository exclusion list and 21 SM/Cer ontology rows. Their historical
  filtering reasons are not inferred from this descriptive inventory; complete
  row and class/ontology ledgers and a reproducible audit are under `QA/`.
- Strict overlap 752 against either all supplied or retained rows. Exclusions
  account for none of the 502 unmatched original labels.
- Literal recorded-alias coverage 959 (=752+207); cumulative candidate coverage
  after chain-order/positional checks 963. Broader coverage is not identity.
- The 502 partition into 202 coarse-alias-only, five unresolved CL chain-order,
  four positional-resolution, seven library-ID-only mismatch, four same-class
  formula-only, and 280 without a recorded annotation/formula candidate.
- Conditional compatibility map 854 (=752+102): shared coarse alias, same
  class/adduct, <=5 ppm and globally bidirectionally unique; strict feature
  locks prevent reuse. No candidate is upgraded to confirmed identity.
- All strict pairs share formulas; inferred neutral mass median difference
  1.40 ppm, maximum 7.79 ppm. Seven strict matches fail the QC-CV restriction
  globally (745 remain), but all 36 strict Figure 3 lipids pass it.
- O/P notation may represent differing structural resolution, including the
  plasmalogen vinyl-bond counting convention; label disagreement alone cannot
  prove distinct chemical identity. CL constituent-chain mismatches and feature
  collisions remain unresolved. Raw spectra/peak confirmation is unavailable.

### Frozen Figure 3 and actual conservation

192 displayed gene-lipid pairs comprise 25 genes/64 lipids; sparse precursor
is 139x64=8,896. The 43 original mice comprise 21 ovaries/22 testes. Maximum
replay discrepancy is 9.99e-16. Strict SD comparison contains 36/64 lipids,
108/192 displayed pairs and 15/25 displayed genes; precursor contains 5,004
pairs/139 genes. Conditional compatibility adds three lipids (39), nine
displayed pairs (117), and one displayed gene (16); precursor 5,421 pairs.

| Universe | Vector r | Native ULM | Nominal p | Overall permutation / bootstrap p | Within-sex upper permutation / bootstrap p |
|---|---:|---:|---:|---:|---:|
| Strict displayed | 0.973266 | 43.6275 | 1.50982e-69 | 0.0104 / 0.0104 | 0.4256 / 0.4480 |
| Strict precursor | 0.855150 | 116.671 | underflow; <2.23e-308 | 0.00005 / 0.0002 | 0.0633 / 0.0960 |

Sex-centered descriptive vector r is 0.269196 displayed (ULM2.87777, p0.00484425)
and 0.0816838 precursor (ULM5.79644, p7.18886e-9). These nominal p-values are
not animal-calibrated residual-network tests. Conditional alias mapping gives
displayed r0.975008, ULM47.06196 and p5.62921e-77; overall bootstrap p0.0098,
within-sex upper bootstrap p0.4316. Fill restriction removes one Figure 3 lipid:
35 lipids/105 displayed pairs and 4,865 precursor pairs; interpretation persists.

The pooled pattern reproduces strongly; shared sex structure explains much of
its concordance. Extra positive alignment beyond sex is not established for
the displayed network. Within-sex bootstrap means vary; unlike permutations,
the bootstrap does not retain observed sex means exactly. Its empirical
generating distributions and group membership are retained. This addendum does
not revise the earlier within-sex zonation result or establish causal effects.
An independent descriptive check keeps only SD between-sex covariance divided
by its actual pooled variation: coefficient-vector alignment is 0.974513
displayed and 0.861048 precursor. This is `cor(reference, E[SD coefficient])`,
not `E[cor(reference, SD coefficient)]`, and supplies no extra p-value.

### Outputs and QA

`followup/README.md` and `WP2_lipid_overlap_and_pooled_Figure3.pdf` are separate
from the preceding WP2 report. Six scientific figure pages cover annotation
classes/partition, pooled/centered coefficient conservation, null distributions,
four predefined Cyp summaries and originally selected individual-mouse pairs.
CSV/RDS outputs retain all candidate graphs, complete original-edge eligibility,
strict/conditional maps, frozen contexts, reference replay, all shared mouse
indices/null/stability draws and native ULM/tails. `pooled_contexts.rds` preserves
four complete input contexts; `run_manifest.json`, SHA256 and sessionInfo record
provenance. The stored RNA data layer affine replay is a scaling check only.

Independent actual QA: 116 scalar OLS/Pearson checks (maximum difference1.34e-12),
160 null-tail summaries, 720 independently reconstructed scoped draw statistics,
and complete mouse/sex/FCG-group index ledgers passed for all four analyses.
See `tests/qa_wdsd_pooled_actual.R` and `followup/QA/`. Input matching, archival
reference and synthetic conditional-null tests also passed. PDF text/figure
layout is rendered and checked before the separate Meeting_Summaries delivery.

Final PDF has 12 pages (six report, six vector figures). All 223 Markdown
paragraph/table parity checks, including 159 table cells, and page-character
bounds passed. Every page was rendered and visually inspected; nine inputs
shared with the preceding WP2 provenance retain identical SHA256. Separate
Markdown and PDF were delivered as
`../Meeting_Summaries/2026-10-07_WP2_lipid_overlap_and_pooled_Figure3.{md,pdf}`.
The copied PDF is byte-identical to the checked report; delivery hashes and
QA evidence are under `followup/QA/`. Build the PDF with
`scripts/reports/build_wdsd_pooled_followup_pdf.py` using the bundled Python
runtime (reportlab/pypdf); no original workflow outputs are replaced.
