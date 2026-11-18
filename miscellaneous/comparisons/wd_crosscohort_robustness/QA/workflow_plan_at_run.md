# Workflow 13 - Western-diet cross-cohort robustness (WP3)

**Entry:** `scripts/workflows/13_wd_crosscohort_robustness.R`.
**Scope fixed before outcome calculation, 14 October 2026:** meeting WP3
compares male CF459 Western-diet networks (four mice) with male CF453 WD
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
coefficients only establish selection/provenance and do not become the main axis.

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

Defaults: seed20261014, bootstrap4,999, paired999; exact permutations fixed24/5040.
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
