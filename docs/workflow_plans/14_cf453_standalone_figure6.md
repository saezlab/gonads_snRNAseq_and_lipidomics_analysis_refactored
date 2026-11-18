# Workflow 14: Standalone CF453 Figure 6

## Current entry and dependencies (18 November 2026)

The current pipeline uses this entry's `--models-only` stage. Default output is
`analysis_current("cf453_models")` under the approved expanded-marker release;
new fits require a fresh `--output` outside that approved root. Fresh full and
`--models-only` modes call `analysis_assert_not_approved_output()` before any
directory creation or model work, rejecting the approved root/descendants even
when generated caches are absent (including symlink aliases). Input defaults are
`current_wdsd_zonation_transport`, `wdsd_hepatocytes_transport` and
`data/metabolomics/wdsd/lipidomics_cf453.csv`. Both RNA transports are mandatory:
shared `wp2_load_inputs()` validates all 28 whole-hepatocyte profiles as well as
56 paired zone profiles, even though the primary RNA fit uses only the zones.
All counts/data/scale layers and metadata/QC are loaded. These inputs are supplied
in the public tree; no source Seurat object, original FCG result or MOFA model is
needed for a fresh CF453 fit.

New model caches record each actual manifest/payload hash. Frozen approved
caches retain their original source-RDS identities; `analysis_verify_inputs()`
can verify an expected original SHA against a complete lossless CSV transport
without deserializing the RDS. It does not substitute a CSV hash for an RDS hash.

Modes (there is no `--validate-only`):

- Default: fit/save models, run the full legacy association/density/bootstrap
  analysis and export its 30-panel presentation.
- `--models-only`: fit/save `models.rds`, model/design/eligibility tables,
  provenance and session information, then stop. This is the upstream stage for
  current Workflows15/16; it does not compute the legacy all-edge networks.
- `--finish`: read `<output>/models.rds`, reject an existing `analysis.rds`, then
  run the legacy network/panel remainder. It is not an input-only validation.
- `--plot-only`: read `<output>/analysis.rds` and rewrite its panels/reports.
  The approved current `cf453_models/` contains only the model stage, so it has
  no full-analysis object for this plotting mode.

`--input-zonation=PATH` overrides only the zone source for a fresh fit; shared
whole-hepatocyte intake remains required. It is rejected with `--finish` or
`--plot-only`. Other values are `--output`, `--bootstraps` (999) and `--seed`
(20261030). Saved modes require generated RDS caches, which are ignored in Git;
a public clone must regenerate them into fresh directories. Source-loading the
entry defines functions without executing analysis. `--finish` and `--plot-only`
remain explicitly mutating saved-result modes and skip that new-fit guard; use
a separate/copied run destination to preserve the approved release.

Production parent reuse checks the recorded input-data identities, not every
mutable code/plan row in a historical provenance file. New exports record current
code/docs hashes. Historical strict export-QA can separately compare all recorded
source/plan hashes and therefore report source-snapshot drift after documentation
edits. Preserve approved manifests/results; do not refresh them to hide that drift.

Approved by the user on 30 October 2026. This workflow implements the accepted
standalone analysis plan, preserving all previous cohort-comparison outputs.

## Scope and dependencies

Use only local CF453 paired RNA/lipid inputs documented in `data/WP1_WP2_INPUTS.md`.
Reuse `scripts/core/wdsd_intake.R` for named-layer access and explicit mouse
alignment, not its frozen-reference workflow. No original FCG, CF459, MOFA,
original selected genes/thresholds, cached enrichment or cross-cohort lipid
overlap enters the statistical models or panel selection. Canonical lipid
formatting is the pure helper in `scripts/core/cyp2c38_reference.R`.

Entry: `scripts/workflows/14_cf453_standalone_figure6.R`.
Statistics: `scripts/core/cf453_standalone.R`.
Plots: `scripts/core/cf453_figure6_plots.R`.
Default output: `analysis_current("cf453_models")`; use a fresh `--output` for any new execution.

## Design and primary analysis

28 mice, seven per female/male x SD/WD group; 56 paired hepatocyte-zone
profiles. Whole-liver lipids are measured once per mouse. RNA model: eight
sex/diet/zone cell means, limma empirical Bayes with mouse duplicateCorrelation.
Count provenance must be audited before choosing a count-derived model. In
the absence of verified raw-count provenance, use supplied normalized RNA.
The intake audit establishes exact ln(1 + 10000 * counts / profile-library)
replay. Primary model uses this unscaled layer, mean-variance trend and robust
empirical Bayes; label effects differences in log-normalized expression, never
log2 fold changes. The initial standardized-expression fit is preserved as a
separate sensitivity under QA/standardized_sensitivity/.
Filter expression using a prespecified expression/variation rule independent
of lipid outcomes; save eligibility and the full design/contrast matrices.

Primary contrasts: male-minus-female difference in PP-minus-PC within SD,
the corresponding WD contrast, and WD-minus-SD of those interactions.
BH FDR across all eligible genes x three primary tests. Secondary simple
zonation effects in four sex/diet groups receive a separate declared family.
Sex/diet simple effects and zone-specific diet differences are secondary.
Gene context membership requires diet-specific primary interaction q<=0.05
and within-sex/diet zonation q<=0.05. Positive = PP, negative = PC; shared
memberships are retained. No weak other-sex t gate or entropy optimization.

Lipid model: sex x diet on 28 mouse log2 supplied intensities. Exclude supplied
standards/isotopes and documented unsupported repository classes, with a full
ledger. Analytical QC CV<=30% is the primary standalone quality set; all
otherwise eligible features and low-Fill exclusions are sensitivity universes.
Do not impose original-cohort molecular overlap. Denominators use all eligible
measured CF453 features within class and sensitivity universe.

## Association analysis and selection

Eight sex/diet/zone contexts, n=7 each. RNA-selected genes only; use the union
of candidates for each sex/zone across diets so SD-WD coefficients are evaluated
on identical pairs. Preserve full coefficients, nominal Pearson probabilities,
and BH q across all evaluated contexts; these probabilities are conditional on
RNA selection. |r|>=0.65 is an exploratory plotting rule, not significance.
Save leave-one-mouse-out estimates and 999 paired mouse bootstrap draws within
sex/diet. The same draw indices carry the lipid vector and both RNA zones;
report undefined correlations explicitly. Resampling is conditional on fixed
selection, not full-pipeline selection validation. Delta r is WD minus SD;
intervals/stability are descriptive, not a null test of diet rewiring.

Use deterministic CF453-only representative ranking: number of displayed
associations stable in sign across all seven omissions, then interaction
evidence, then gene symbol. Select at most one gene per sex/zone, using the
same representative in both diets. Up to six stable partners ranked by
cross-diet maximum absolute correlation; retain both-diet values and never
fabricate partners. Save the complete selection ledger. If no candidates pass,
export continuous effect plots and honest empty network panels.

## Panels and delivery

A: female x male zonation-statistic scatter, separate SD/WD.
B: CF453 design/matching schematic.
C: within-class proportions with unique/shared/unassigned context memberships,
separate SD/WD; do not equate lipid sex abundance with RNA association context.
The observed any-edge rule saturates (all primary lipids shared in both diets).
Keep that result and add C_pair_density_SD/WD: fraction of all primary-universe
gene-lipid pairs with |r|>=0.65, using the fixed candidate union within sex/zone
in both diets. This is a descriptive companion, not an alternative significance
test or new threshold. Export continuous absolute correlations and adjusted-p
fractions alongside the display fractions.
D: independent representative subnetworks, each sex/zone/diet.
E: diet modification of zonation effect sizes with three-way interaction q.
F: individual-mouse RNA-lipid examples, separate files.
Supplementary: input/marker diagnostics, lipid factorial effects, display
threshold and lipid-quality sensitivities. No assembled manuscript figure.
Every scientific panel: PDF/SVG/PNG plus source data, manifest and legend.

## Verification and limitations

Run every new/modified R script immediately. Independently replay design signs,
paired-model contrasts, BH families, every network coefficient and omission,
class denominators, representative selection, and absence of reference inputs.
Render final PDFs and inspect all panel previews for clipping/legibility.
Export package/session versions, seed, input/code SHA256 and resampling indices.
Results describe sex associations and whole-liver lipid relationships, not
isolated hormone/chromosome effects, causal regulation or spatial lipid levels.
Missing nucleus counts, raw-processing and batch/litter provenance remain
limitations. Update this plan with actual results and any justified deviations.

## Historical two-marker full-analysis results and verification (30 October 2026)

Primary modeled RNA: 18,486 genes; detection in >=7 of 56 profiles and nonzero
variation. Supplied normalization replays ln(1+10000*counts/library) to <1e-15.
Mouse duplicate correlation is 0.7160454. Sex-by-zone interaction FDR<=0.05:
536 SD genes and 607 WD genes; three-way sex-by-diet-by-zone: 194 genes.
These share one 55,458-test primary family. The preserved standardized model
has 474/554/191 detections, exactly concordant coefficient signs, and primary
t-vector correlations 0.9911/0.9917/0.9930. That sensitivity also differs in
moderation settings and is a robustness comparison, not an isolated scale test.

Lipid primary universe: 1,129 analytical QC-CV<=30% features; broad endogenous
universe 1,177. No input changes or model dependence on another cohort.
The fixed candidate unions are 244 female PP, 277 female PC, 301 male PP,
396 male PC genes (genes can occur in more than one context). Broad paired
coefficients: 2,867,172; primary coefficients: 2,750,244 (1,375,122 per diet).
Zero primary correlation tests pass global BH q<=0.05; all association displays
remain exploratory. All 1,129 primary lipids have shared membership in both
diets at |r|>=0.65. Fixed-universe class/context density ranges 0.02657807 to
0.31046931. These fractions are descriptive and dependent across features.

Representatives selected by the prespecified ranking: female PC St3gal6,
female PP Auts2, male PC Unc5c, male PP Hspa12a. Each uses the same six partners
in SD and WD; below-threshold or not-currently-nominated results remain clearly
marked. The example genes are not asserted to be causal lipid regulators.
Save all 999 representative correlation draws, 23,976 matched independent-arm
WD-minus-SD difference draws and their descriptive bounds in addition to
the all-edge per-diet bootstrap summaries/indices. Bootstrap uncertainty is
conditional on fixed candidate and example selection.

Independent base-GLS/OLS, input reconstruction, contrast/BH, every broad
coefficient and all 20,070,204 omission values, membership/sensitivity,
representative selection and 143,856 bootstrap replay values passed. The
extended audit also reconstructs all representative/delta draw banks and seed.
The bootstrap variance-guard audit found zero exact varying/constant support
misclassifications in the saved banks; the minimum genuine weighted variance
is 3.439965e-12, above the implementation's 1e-12 guard. Audit CSVs and code
remain under QA/ and tests/. Empty-candidate delivery was separately smoke
tested and preserves unassigned membership and continuous panels.

The initial run used separate model and network stages while provenance and
rendering were developed. Editing the running entry file caused a trailing
Rscript parser error after numerical exports and panels were already saved;
the initial log is preserved. Current sources parse cleanly and the subsequent
actual --plot-only run completes successfully. Presentation/representative
supplements are reproducible through the current full entry. Final code hashes
and QA distinguish completed statistical artifacts from that initial CLI exit.

Final delivery contains 30 independent panel groups: 30 PDF, 30 vector-only
SVG and 30 300-dpi PNG files, plus 30 source-data CSVs. Every PDF was rendered
and all five review contact sheets inspected; bounds, searchable text, glyphs,
SVG raster-content checks and PNG dimensions pass. The plotting fallback reads
markers directly from saved intake QC when no separately generated audit exists.
No figure assembly is performed. The compact panel/source/legend ZIP is built
by `python3 scripts/reports/package_cf453_figure6.py`; full tables remain outside
that bundle. CRC and per-member SHA256 verification pass.

A fresh complete CLI run with 19 bootstrap draws exits successfully; the final
saved-input marker fallback was then verified by its --plot-only mode. All
observed models, full correlations/omissions, selections, memberships, density
and all 30 panel source tables match production in 44 replay checks (bootstrap
summary columns intentionally excluded because draw counts differ). Production
keeps 999 draws. Logs/replay checks are preserved in QA/fresh_cli; temporary
smoke outputs are removed after verification. Full independent actual checks:
720 groups, maximum absolute difference 1.509903e-13; density: 2,400 checks.

## Expanded-marker model-stage production record (16 November 2026)

Workflow 14 now accepts an optional `--input-zonation=PATH` for a fresh fit;
the default is now the current expanded-marker CSV transport. `wp2_load_inputs()` and
`cf6_models()` forward the explicit override and preserve its actual path.
Models-only runs now export a complete parent cache: `models.rds`, model/design
tables, input/model/code SHA256 and session provenance, without all-edge analysis
or the historical 30-panel set. Saved-result modes reject a new input override.

The marker-score input is staged under
`results/zonation_marker_update_20261116/inputs/wdsd/`; the new cache is
`cf453_models/` under that version root. The unchanged eligibility rule retains
18,604 genes; mouse correlation is 0.7859836. The existing three-test family
identifies 570 SD sex-by-zone, 621 WD sex-by-zone and 185 three-way genes.
Mouse pairing and all 1,129 primary lipid values/results are unchanged.
Historical inputs, caches, results and panels are preserved. Old/new eligibility,
contrast, named-gene and selection comparisons are exported separately under
the version root's `comparisons/cf453_*.csv`; see Workflow 17 for the overall
update audit. These new counts replace no historical records.
