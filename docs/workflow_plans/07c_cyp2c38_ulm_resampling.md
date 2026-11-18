# Workflow 07c — Cyp2c38 ULM reporting and resampling audit

**Compatibility phase documentation — 9 October 2026:** the active complete
work package is [11](11_cyp2c38_in_vivo_validation.md), entered through
`scripts/workflows/11_cyp2c38_in_vivo_validation.R`. This former entry is a thin wrapper for its
ULM/resampling phase; orchestration now lives in `scripts/core/wp1_workflow.R`.
Inputs are repository-local, with byte-preserving import provenance in
`data/wp1_wp2_input_migration_manifest.csv`. The recorded production results
below remain unchanged. The complete WP1 entry passes fresh intake/results between phases; standalone
phase execution retains its documented dependencies.

**Entry point:** `scripts/workflows/07c_cyp2c38_ulm_resampling.R`.
**Scope:** October 2026 WP1 follow-up: focus reporting on ULM scores, nominal ULM p-values and explicit bootstrap p-values; regenerate separate ovarian/testicular Cyp2c38–lipid networks; check calibration versus low power. Execution and reporting dated 2026-10-07.

Read Workflow 07b's plan and the cross-workflow overview before editing. This extension reads the completed 07b KO matrix/metadata and reconstructs the original paired FCG reference. It writes a separate `results/validation/cyp2c38_in_vivo/ulm_followup/` directory, preserving the original analysis.

## Specification

- Recompute all 1,254 original lipidomics features against Cyp2c38 in the paired original FCG hepatocytes, separately for ovaries (21 mice), testes (22 mice), and pooled (43). Export all coefficients without a correlation/significance cutoff. Also export the fixed 828-feature historical network and verify its exact replay. Match the retained 07b KO features by strict molecular identity: 856 full-universe and 552 historical matches. Report both results so the expanded feature universe is explicit.
- Retain the 07b log2, four-group, sex-precision-weighted robust/trend limma model. Use matched female KO-minus-WT versus ovarian coefficients, male KO-minus-WT versus testicular coefficients, and equal-sex pooled contrast versus pooled coefficients.
- Compute the original decoupleR ULM with intercept; negative score denotes concordant response to reduction of Cyp2c38. Validate any fast resampling implementation against decoupleR scores and nominal p-values. Display scores and p-values, omitting the earlier sign-reversed correlation score, intervals and FDR from the new report.
- Distinguish resampling methods: the previous ordinary combined bootstrap estimated uncertainty, while the previous p-values were from genotype permutations. The main wild bootstrap imposes the single tested mean contrast equal to zero for every lipid; remaining group-mean dimensions stay unrestricted in the constrained weighted fit. This null is stronger than zero ULM alignment with arbitrary nonaligned KO effects. Use full four-group residuals with HC2 leverage correction and one shared Rademacher sign per whole mouse vector. Refit sex weights, robust/trend moderation and ULM per draw. Fix the reference in the primary test; independently resample paired FCG animals within XXO=10, XYO=11, XXT=11, XYT=11 for the combined sensitivity. Null-restricted residuals provide an additional sensitivity. Every empirical p is two-sided with the +1 correction.
- Compare bootstrap results with whole-mouse genotype permutation. Examine null calibration under a fitted covariance model and known aligned-effect injection curves at the actual sample sizes. Low power and excessive conservatism are different properties; simulations diagnose the chosen fitted scenario and do not establish universal calibration or post-hoc power for the unknown true effect.
- Save complete network, identity and ULM tables, all null draws, calibration/power diagnostics, software/input provenance and a self-contained Markdown/PDF report. Retain primary-hepatocyte siRNA context separately.

## Tests and limitations

Run every new/modified R script immediately and run meaningful tests for ULM equivalence, contrast sign, model/mouse ordering, null centering, shared mouse multipliers across lipids, reference pairing, strict one-to-one joins and RNG reproducibility. Calibration requires finite valid draws and reported simulation replicate counts/seeds. Bootstrap null assumptions and uncertainty due to missing batch/litter/processing metadata remain explicit. Correlations do not establish substrates or flux. The original 07b output remains available for audit.

## Implementation and saved artifacts

Helpers in `scripts/core/`: `cyp2c38_sex_networks.R` reconstructs and exports
the six reference networks; `cyp2c38_ulm_bootstrap.R` implements fast validated
ULM, weighted contrast-null projection, HC2 and paired reference draws;
`cyp2c38_ulm_calibration.R` runs conditional nested calibration and known-effect
injections; `cyp2c38_ulm_presentation.R` validates and writes figures and report.
The unchanged 07b reference/model helpers are reused.

Default run: seed 20261007, 4,999 main null draws, 1,999 restricted-residual
draws, 1,999 within-sex genotype permutations, three workers, 120 null datasets
with 199 inner draws each, 399 injected datasets for each of strengths
0/0.1/0.2/0.4/0.8 and each contrast. CLI controls also accept a fresh `--output`
directory and `--validate-only`. Per-draw seeds make worker scheduling
irrelevant to RNG. Ten comparison rows include three primary and two secondary
contrasts in each of two feature universes.

Core outputs: `ULM_summary.csv`, `ULM_matched_lipid_edges.csv`, six `networks/`
SIF CSVs and `network_manifest.csv`, full identity/mouse ledgers, four null draw
banks, distribution and previous-bootstrap variance diagnostics, nested null
and injection tables in `diagnostics/`, input/code SHA256, provenance/session,
three-page `plots/ULM_concordance.pdf` and self-contained `README.md`.
`scripts/reports/build_cyp2c38_ulm_pdf.py` renders that Markdown with the vector
figure appendix. It needs reportlab/pypdf and DejaVu fonts; the bundled Codex
Python runtime is suitable. Report copies are saved in sibling Meeting_Summaries.

Calibration simulates global absence of KO effects using fitted whole-mouse
HC2 vectors with sign-symmetric errors, then performs a fresh contrast-specific
test for each dataset. Exact binomial checks against 5% are unadjusted internal
diagnostics. With 120 datasets, moderate conservatism cannot be ruled out.
Injection detection uses fixed observed-data null 95th-percentile absolute ULM
thresholds; these are conditional sensitivity curves for idealized full-network
aligned effects, not prospective power or power for the unknown true effect.

Focused tests: `tests/test_cyp2c38_ulm_bootstrap.R` and
`tests/test_cyp2c38_sex_networks.R`; the full entry point independently compares
all observed results with decoupleR and the saved 07b t-statistics. Presentation
and calibration helpers also passed synthetic/actual-data smoke tests, invalid
input checks and worker-invariant RNG checks during development.

## Verified full execution, 2026-10-07

All 4,999 primary and combined bootstrap draws were finite for all ten rows;
the full entry point completed successfully with the defaults above. Original
reference replay differed by at most 5.56e-16; unchanged KO moderated t-values
replayed within 5.19e-13. All observed ULM values and p-values agreed with
decoupleR within 1e-10. Both permanent focused tests passed.

Full-network primary pooled/female/male native ULM scores are
-1.930037 / +0.331177 / -5.706102, nominal p-values
0.0539332 / 0.7405921 / 1.593682e-08. KO-only null-bootstrap p-values are
0.8178 / 0.9528 / 0.4828; KO-plus-reference p-values
0.8200 / 0.9452 / 0.4488. Within-sex genotype permutation gives
0.8060 / 0.9455 / 0.4590. These support the same broad inference while the
nominal male association has the predicted direction.

Full-network conditional null rejections at p<=0.05 are 5/120, 10/120 and
9/120 (4.2%, 8.3%, 7.5%); exact binomial diagnostics do not reject a 5% rate.
The historical female subset rejects 14/120 (11.7%; unadjusted diagnostic
p=0.0028), indicating imperfect calibration rather than uniform conservatism.
No claim of universal or exact calibration is justified. In the full network,
fixed-threshold sensitivity at injection strengths 0.1/0.2/0.4 is
27.1/68.9/100.0% pooled, 26.6/71.4/99.7% female and
14.0/36.6/90.5% male. These are fitted, idealized aligned-effect scenarios,
not power for the unknown true effect. The ULM report records all settings,
physical injected log2 scales, counts and methodological qualifications.
