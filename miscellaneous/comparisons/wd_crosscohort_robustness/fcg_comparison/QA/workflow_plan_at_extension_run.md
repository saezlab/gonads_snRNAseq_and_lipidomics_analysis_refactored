# Workflow 13 - Western-diet cross-cohort robustness (WP3)

**Entry:** `scripts/workflows/13_wd_crosscohort_robustness.R`.
**Scope fixed before outcome calculation, 14 October 2026:** meeting WP3
initially compares male CF459 Western-diet networks (four mice) with male CF453 WD
networks (seven mice), separately Periportal and Pericentral. No CF459 SD arm
exists. Female analyses, TF-gene replication, MOFA retraining and nucleus-level
composition analyses are outside this work package.

Read the original meeting plan WP3 section, Workflow 12 plan and meta overview
before modification. Existing WP1/WP2/FCG outputs and unrelated work remain fixed.

## Sections and dependencies

1. Source support functions and optional existing local R library; parse CLI.
2. Read repository-local CF459 inputs (`wd_intake.R`) and CF453 intake
   (`wdsd_intake.R`). Reconstruct frozen Figure 5 FCG references with
   `wdsd_reference.R`, checking archived log2 coefficients before matching.
3. `wd_crosscohort_setup.R`: strict molecular crosswalk, explicit mouse order,
   fixed full/displayed eligibility and independent CF453 WD replay against WP2.
4. `wd_crosscohort_statistics.R`: observed cross-cohort coefficients, native
   ULM/Pearson, exact pairing/null bootstrap and paired/omission stability.
5. Figures, CSV/RDS/index ledgers, Markdown report, input/code SHA256 and
   session/CLI provenance. PDF rendering is a separate reproducible Python step.

Required local CF459 files: `data/metabolomics/wd/lipidomics_wd.csv`,
`data/RNA/wd/wd_pseudobulks.rds`, `data/RNA/wd/wd_zonation_pseudobulks.rds`.
Copy bytes unchanged from shared original_datasets archives, preserve sources,
and record hashes in `data/wp3_input_migration_manifest.csv`.
CF453 inputs remain under `data/metabolomics/wdsd/` and `data/RNA/wdsd/`.
Frozen FCG references are Workflow 05 Figure 5 networks and original matrices;
Workflow 12 outputs supply a numerical replay check, not feature selection.

## Intake and fixed analysis choices

CF459 CSV uses semicolons, decimal commas and UTF-8 BOM. It has 1,463 supplied
feature rows and eight animals. Metadata f-ctrl/m-ctrl denotes female/male WD
according to meeting provenance, not a control diet. IDs1-4 are females and
13-16 males. Validate this independently across lipid and RNA records. Each
mouse supplies one whole-liver lipid vector and two RNA zone profiles; zones
are not independent animals. Supplied objects contain five cell types or two
zones and 20,123 genes. Use named Seurat v5 layers, no new normalization or
imputation. Main RNA is scale.data; data-layer correlations are a sensitivity
and count/library/marker checks are aggregate diagnostics only.

Retain exact existing formatting-only molecular matching (last pipe alias and
class parentheses), preserving chain order, positional, ether and oxygenation
resolution. Never join by Alignment ID alone, sum features or force ambiguous
matches. Keep complete exclusions/duplicates/identity ledgers. Use existing
repository excluded classes and standard/isotope handling. CF459 lacks analytical
QC columns, adduct/mass/RT metadata and Fill values; do not invent QC-CV or
additional chemical confirmation. CF453 QC-CV<=30% and Fill>=0.5 restrictions
are one-cohort sensitivities, explicitly not two-cohort QC validation.

Freeze original FCG male context-specific genes and lipid universe. Primary
quantitative comparison uses all eligible original full gene-lipid pairs;
originally displayed |r|>=0.65 pairs are a secondary fixed subset. Eligibility
requires unique strict measurement identities, measured exact genes and finite,
varying features in both male WD cohorts. No observed-correlation threshold or
outcome-based outlier removal. Save missing/excluded/invariant original edges.
Main axes: CF453 WD Pearson versus CF459 WD Pearson, using log2 lipids. FCG
coefficients establish the initial selection/provenance; the dated extension
below also uses them as a measured standard-diet comparison anchor.

## Inference and stability

Report coefficient-vector Pearson r, native intercept ULM slope t and nominal
ULM p. Edges share genes/lipids and mice; nominal edge-vector p-values are
descriptive and do not establish independent biological replication.

Primary animal null: exhaustive 4!=24 CF459 RNA-lipid pairings anchored to the
observed CF453 WD coefficient vector. Complement: exhaustive 7!=5,040 CF453
pairings anchored to observed CF459. These are distinct conditional no-pairing
hypotheses; do not combine or select the smaller result. Use two-sided absolute
ULM tails including observed identity; all features and both zones share indices.
CF459 exact resolution is 1/24 (before ties/symmetries). Independent empirical-
marginal CF459 null bootstrap (4,999 draws; fixed CF453) is a secondary approximate
no-association null with shared whole-mouse RNA and lipid indices; Monte Carlo
p uses +1 correction. No both-cohort joint-null test is required.

Independent paired whole-mouse resampling in both cohorts (999 draws) describes
stability and supplies no null p-value. At n=4 many draws have few distinct mice;
save unique-mouse counts and undefined frozen-vector draws. Undefined edges
invalidate the frozen summary, never trigger within-draw feature deletion.
Single-cohort omissions (7+4) and all 7x4=28 joint omissions are descriptive
influence checks; CF459 then has only three mice. No claim of calibrated absence
of conservation, formal differential correlation or diet-specific rewiring.

Mtmr11 is predefined where present; species-level example selection uses frozen
FCG membership/ranking, not new cross-cohort responses. Additional highlighted
changes must be labelled descriptive selection. Compare age/background/diet
exposure only if source metadata establishes them; otherwise record unavailable.

## CLI and outputs (implementation target)

Defaults: seed 20261014, bootstrap 4,999, paired 999; exact permutations fixed 24/5040.
Fresh output root `results/wd_crosscohort_robustness/`; reject nonempty analysis
outputs. `--validate-only` checks actual local inputs/reference/matching without
writing. `--output`, `--seed`, `--bootstraps`, `--paired-bootstraps` control runs.
`source()` defines functions; Rscript executes the complete workflow.

Outputs: sample/QC/input and identity/edge ledgers; coverage/replay tables;
full and displayed comparison inferences; all species coefficients and SIFs;
complete shared null/paired/permutation index banks; LOCO tables; sensitivities;
plots and a self-contained Markdown/PDF report, delivered to Meeting_Summaries.
Preserve native ULM/nominal and explicit bootstrap p in technical report; focus
interpretation on correlation magnitude and animal-level uncertainty. No FDR,
reversed score or confidence-interval reporting is requested here.

## Verification and limitations

Run each new/modified R script immediately and meaningful tests for locale,
mouse/zone pairing, strict resolution, frozen membership, direct Pearson/OLS
and decoupleR agreement, exact enumeration/identity, shared indices, invalid
constant draws, seeded replay, overwrite protection and CLI validation.
Independently audit actual output coefficients, resampling tails, mouse indices
and original WP2 coefficient replay before conclusions. Render and inspect all
PDF pages; verify table/text/figure and linked-file integrity before delivery.
Update this plan with implemented APIs/counts/tests/results after completion.

Main limits fixed in advance: four male CF459 animals, unknown upstream
processing/batch/litter/nucleus contributions, absent CF459 analytical QC and
chemical metadata, frozen selection, dependent edges and confounded cohort
age/background/exposure. Positive agreement concerns those cohorts only;
negative/inconclusive results are not proof against lipid regulation.

## Requested original-cohort comparison - amendment of 14 October 2026

The user explicitly expands WP3 to ask whether CF459 is closer to the original
paper male standard-diet network than to CF453 male Western diet. Compare both
anchors with CF459 separately within each hepatocyte zone. The original
network cohort is **22 testes-bearing FCG mice (11 XXT and 11 XYT)**, not the
five separate wild-type male controls and not only the XYT subgroup.

**User clarification, 14 October 2026:** CF459 is **3-month WD**, original
FCG is **3-month SD**, and CF453 is **10-month WD**. The user corrected the
initially transposed ages; CF453 and FCG now agree with the meeting record.
These are supplied cohort ages, not fields present in the RNA/lipid exports.
The explicit question is whether the young WD cohort is more similar to young
SD than to older WD. Diet duration, exact background, batch/litter and other
cohort differences remain unestablished. The comparison cannot identify causal
age and diet effects separately.

**Primary comparison specification, fixed before this extension's resampling:**
use exactly the same strict three-cohort gene-lipid pairs for both CF459–FCG and
CF459–CF453. Retain 610 original lipid identities and fixed context genes:
15,250 periportal and 18,300 pericentral full pairs (primary), with 69 and 202
original displayed pairs (secondary). FCG displayed pairs were selected in the
original data; their comparison remains conditional on that selection.

Reconstruct original FCG coefficients from paired original matrices and require
exact replay against the archive. Compute all three pairwise coefficient-vector
correlations and define relative closeness as
`delta_r = cor(CF459, FCG) - cor(CF459, CF453)` on identical edges. Positive
values mean CF459 has a more positive coefficient association with FCG; this
alone does not imply strong absolute conservation. Native ULM/nominal p and
CF459 exact 24/empirical-marginal null-bootstrap checks accompany the new FCG
anchor using the same CF459 draw banks as the original WP3 comparison.

Assess delta_r descriptively with 999 paired three-cohort resamples: reuse the
same CF459 mouse draw simultaneously against both anchors, independently
resample CF453 paired mice, and independently resample FCG paired RNA/lipid
mice within 11 XXT/11 XYT strata. Save the FCG seed stream/index ledger and all
valid/undefined counts. A paired-bootstrap delta-positive fraction is stability,
not a null-test p-value. Do not use an independent-correlations Fisher test.
Preserve frozen edges in every draw, invalidating summaries if any needed edge
is undefined. Check each cohort's individual-mouse omissions and existing
CF453/CF459 joint omissions. No genotype/feature reselection is performed.

Store extension data/provenance under `fcg_comparison/`, preserving original
WP3 numeric outputs. Fresh Workflow 13 runs include this comparison; an explicit
extension-only mode can add it to completed original WP3 outputs. Update the
original research plan itself and both Markdown/PDF results, showing matched
comparison/scatter and relative-closeness stability. Preserve prior document
versions. Conclude only about measured network similarity across these cohorts;
this is not an identified causal age-versus-diet contrast.

## Initial implementation and production verification - 14 October 2026

The numerical results in this section are the historical first WP3 run. The
original-cohort extension and its independent audit add an exact numerical-zero
null-bootstrap correction without changing observed or paired/omission results.
Corrected null denominators and p-values are exported separately and explained
in the current report; historical files and the initial audit record remain
unchanged.

The implemented entry and `scripts/core/wd_crosscohort_workflow.R` run all five
sections above. `wp3_run()` accepts `--output`, `--seed`, `--bootstraps`,
`--paired-bootstraps`, `--validate-only`, `--report-only` and
`--fcg-comparison-only` (optional fresh `--comparison-output`). Counts must be
integers (bootstrap/paired >=19); default exact banks stay 24/5,040.
Unknown/duplicate options and nonempty full-run destinations are rejected.
Source loading defines functions; validation writes no files or saved WP3
results; report-only regenerates Markdown from `report_bundle.rds`. Scientific
plot/PDF rendering is separate from numerical refitting.

Three local inputs total 13,118,772 bytes with verified source/destination hashes.
CF459 retains 1,395 of 1,463 lipid rows after 68 class exclusions, no standards
or additional isotope rows. There are 244 duplicate Alignment-ID excess rows
and zero canonical-name duplicates. All eight RNA/lipid identities and 16 zone
profiles align; four males are IDs 13-16. Counts/data/scale.data contain 20,123
genes; scaling replay errors are <=1.25e-14. CF459 QC-CV/Fill/nucleus counts
remain unavailable, not zero or a substitute pass criterion.

The strict three-way crosswalk retains 610 of 1,254 frozen original lipids.
Periportal: 25 genes, 15,250/31,350 full pairs and 69/155 originally displayed pairs
(the displayed common set has 18 genes/62 lipids). Pericentral: 30 genes,
18,300/37,620 full and 202/387 displayed pairs (15 genes/146 displayed lipids).
No additional finite/variable-feature exclusions occur in observed data.
All 33,550 CF453 WD overlapping coefficients replay WP2 within 5.55e-16.

| Context/universe | Coefficient-vector r | Native ULM | Nominal ULM p | Exact CF459 p | Complementary CF453 p | Conditional null-bootstrap p |
|---|---:|---:|---:|---:|---:|---:|
| Periportal full (primary) | -0.06309949 | -7.807259 | 6.222888e-15 | 0.1250 | 0.0422619 | 0.06784232 |
| Pericentral full (primary) | -0.01857083 | -2.512511 | 0.01199607 | 0.7916667 | 0.5966270 | 0.76555246 |
| Periportal displayed | 0.03167365 | 0.259390 | 0.7961296 | 0.8750 | 0.7944444 | 0.83609959 |
| Pericentral displayed | -0.07563094 | -1.072655 | 0.2847184 | 0.6666667 | 0.5513889 | 0.54392220 |

Full-count production and an independent full-count pilot were exactly identical
for the complete analysis object, including all banks/indices/omissions. All
three sensitivities (CF453_QC30, CF453_Fill50, RNA_data) rerun full counts using
identical mouse-index banks. Near-zero agreement persists. Native ULM agrees
with decoupleR; the independent actual verifier passes 500 checks, maximum
absolute discrepancy 6.40e-14, including every 24/5,040 bank coefficient and tail,
all observed edges/OLS and all 39 omission configurations. Focused intake,
setup/statistics and CLI guards pass. All 239 prior WP1/WP2 files remain
byte-identical; see `QA/prior_WP1_WP2_result_preservation.json`.

All exact pairings remain defined. Null bootstrap is defined 4,819/4,999 for
both periportal universes, 2,153/4,999 pericentral full and 3,084/4,999 pericentral
displayed. Full paired stability is defined 972/999 periportal and 391/999
pericentral. `feature_bootstrap_degeneracy.csv` records why: pericentral CF459
Shisa7 and Krt4 each have only two expression levels with one distinguishing
mouse, becoming constant 1,611 and 1,565 null RNA draws respectively; Pitx3 is
constant 390 draws. No post-hoc gene deletion is used. Undefined-mass completion
bounds are mathematical bookkeeping, not confidence intervals or calibration.

The data do not demonstrate robust positive cross-cohort coefficient-vector
conservation in these measured male WD features. Negative full ULM values are
small negative coefficient associations amplified by many dependent edges,
not large biological effects. Four CF459 animals, sparse expression support and
unknown cohort differences prevent conclusions of biological absence or an
identified age/background/diet effect. Mtmr11 remains a predefined periportal
example (610 full/6 displayed pairs); selected mouse examples use frozen FCG
ranking only, not observed agreement.

Primary exports include `contexts.rds`, `analysis.rds`, `crosscohort_inference.csv`,
`network_edges.csv`, full crosswalk/eligibility/QC tables, every null/paired index
and draw table,39-omission summaries/matrices, three complete sensitivity fits,
feature-degeneracy diagnostics, eight SIFs, input/code SHA256 and session/CLI.
`wd_crosscohort_presentation.R` creates six vector figure PDFs, fixed example
selection, `report_bundle.rds` and README. `build_wp3_report_pdf.py` renders six
narrative pages plus the six scientific figures; final table/text/vector/page
QA is under `QA/`. Report source and numerical-run code hashes are separately
preserved when presentation layout is polished; `workflow_plan_at_run.md` is
the pre-execution plan, while this document includes completed results.
