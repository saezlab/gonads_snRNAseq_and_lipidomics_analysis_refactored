# Workflow 15: Direct diet impact on CF453 sex and zonation differences

## CSV transport update - 18 November 2026

Parent verification uses `analysis_verify_inputs()` before model use. For frozen
source-RDS identities it verifies the expected original SHA against the CSV export
manifest, all CSV file hashes and exact numeric-content hashes without reading the
RDS. `input_verification` records the opened manifest separately from the original
identity. New CSV parent inputs are checked as ordinary actual files. This is an
input representation change only; no contrasts, estimates or selection rules
changed and no new model fit was needed for validation.

---

## Current model entry - 16 November 2026

Default parent is the registry's current `cf453_models/models.rds`; default output
is current `cf453_direct_diet/`. Workflow17 invokes `--models-only` here before16.
The historical four-panel presentation remains under `historical_cf453_direct`.
New fits require fresh output. Raw-source hashes are checked against the original
recorded manifest identities while opening the actual relocated files; fresh
provenance records retain those identities with a separate `resolved_path`.
No statistical contrast, family or selection rule changed during path promotion.

---

Requested 2 November 2026 to replace the 30-panel presentation with four
independent panel groups focused on direct differences between diets.
Preserve Workflow 14 and its historical outputs. This workflow contains no
cross-cohort/network conservation analysis or separate-diet network discovery.

## Inputs and model

Use the validated `results/cf453_standalone_figure6/models.rds` containing
18,486 eligible genes, 56 profiles from 28 mice, named eight-cell RNA design,
mouse blocking correlation, and the four-cell/28-mouse lipid model. Verify the
three source-input SHA256 hashes against Workflow 14 provenance before use.
This is a model extension on the same supplied normalized RNA and log2 lipid
data, not a new normalization or change to animal/feature eligibility.

Refit limma on saved modeled_rna and its design, with saved mouse rho, and
retain mean-variance-trend/robust empirical Bayes moderation. Positive sex
direction is male minus female; positive zone is PP minus PC; diet is WD minus
SD. Equal averaging defines the marginal two-way interactions.

Let DS_PC and DS_PP be the WD-minus-SD change in male-minus-female expression
within each zone. Let DZ_female and DZ_male be the WD-minus-SD change in
PP-minus-PC expression within each sex.

Three primary direct interaction contrasts are:
- diet_sex = (DS_PC + DS_PP)/2;
- diet_zone = (DZ_female + DZ_male)/2;
- diet_sex_zone = DS_PP - DS_PC = DZ_male - DZ_female.

BH correction pools all genes x these three contrasts. A separate three-df
moderated F test of exactly these independent contrasts asks whether diet
changes any sex/zone difference, with BH across genes. Its covariance must
include correlations between contrasts. Use only this three-column fit for
the omnibus F, not a fit containing redundant conditional contrasts.

Secondary BH family: all genes x four conditional interactions. Identities:
DS_PC = diet_sex - diet_sex_zone/2;
DS_PP = diet_sex + diet_sex_zone/2;
DZ_female = diet_zone - diet_sex_zone/2;
DZ_male = diet_zone + diet_sex_zone/2.

The new primary family differs from Workflow 14's SD/WD sex-zone family;
do not transfer its FDR detection counts, including the old 194 count.
Raw three-way estimates/p values must replay the prior model.

Lipid inference retains the direct sex x diet contrast on 1,129 primary
features (one measurement per mouse, no lipid zone replicates), BH across
that interaction's species. Plot female and male WD-minus-SD effects; the
direct interaction tests their difference. Do not compare separate-diet
significance or interpret a uniform diet main effect as changed sex difference.

## Four compact panel groups

A. One three-facet overview of the primary RNA diet interactions, showing
effect estimates and adjusted evidence. Label only a few directly ranked genes.
B. One heatmap of the top up-to-20 omnibus-FDR genes, ranked q, p, gene symbol,
with the four conditional interaction estimates and secondary-FDR markers.
C. One four-example plot: largest absolute diet_sex estimate among primary
FDR<=0.05 genes (mean across zones per mouse); largest distinct diet_zone
estimate (PP-PC per mouse); two further distinct genes ranked by absolute
diet_sex_zone estimate (PP-PC). Ties break by q, p, gene symbol. These are large
estimated changes selected from the same data, not independently validated
or necessarily the most biologically important genes. Leave unsupported groups
empty rather than relax FDR criteria.
Show independent mouse points, model group means and nominal 95% intervals.
Lines connect SD/WD group means; they never imply the same animals changed diet.
D. One whole-liver lipid plot: female WD-SD versus male WD-SD, identity line
and direct diet x sex FDR highlighting. Labels rank by interaction q/p/ID.

Four panel files in PDF, plus alternative SVG/PNG formats (12 graphic exports),
four matching source CSVs, concise methods/legends and a portable ZIP. No extra
main or supplementary plots and no assembled manuscript figure. Store under
`results/cf453_diet_impact_figure6/`.

## Implementation and verification

Entry: `scripts/workflows/15_cf453_diet_impact_figure6.R`.
Core: `scripts/core/cf453_diet_impact.R`; plots and report in matching modules.
Run every new/modified R script immediately. Verify contrasts with independent
base GLS, marginal half-weights, four decomposition identities, 3df F with full
covariance, all BH families, old raw three-way replay, model-derived group
intervals, example selection and exact plotted mouse summaries. Verify lipid
interaction as a difference of diet effects and same 28 animals. Render and
inspect all four PDFs and check vector SVGs and source-data completeness.
Record input/model/code hashes, package versions and interpretation limits.

RNA effects are natural-log normalized expression differences, not log2FC.
Lipid effects are log2 abundance differences. Same limitations on original
count aggregation, nuclei/batches/litter and randomization provenance remain.
Diet-modified observational associations do not isolate hormones or establish
causality. Selected examples and nominal intervals are not independent validation.

## Completed results and verification (2 November 2026)

The pooled 55,458-test direct-diet RNA family identifies 166 diet x sex,
350 diet x zone and 137 three-way genes at BH q<=0.05. The separate omnibus
3-df test identifies 605 genes. The four conditional contrasts have their own
73,944-test BH family. All raw three-way effects and p values replay Workflow
14 exactly; 137 versus its previous 194 is a change of testing family, not a
change of input expression or raw three-way fit.

The 20 omnibus-ranked heatmap genes decompose the supported changes.
Magnitude-ranked significant examples are Egfr (diet x sex, -2.1065775),
Limch1 (diet x zone, +0.8022413), Meg3 (three-way, +0.8781388) and Cyp3a11
(three-way, -0.7334390), in natural-log normalized expression units. All 28
mouse values per example and all 16 fitted group means/nominal intervals
are retained. No original network ranking determines these examples.

Direct lipid sex x diet interaction: 267/1,129 primary species (228 positive,
39 negative). This is a test of differential mean diet response, distinct from
the exploratory gene-lipid correlation tests in Workflow 14.

Independent actual-data verification passes 170 check groups, maximum absolute
difference 1.239009e-13, covering all model covariances, F/t probabilities, BH
families, decomposition identities, example selection and mouse summaries,
lipid inference and input hashes. A fresh complete default CLI run passes 17
exact replay checks for the numerical results and all four source tables.
Input/model/current-code SHA256 records are written by the workflow.
All four PDFs were rendered
and individually inspected; PDF bounds/text, vector-only SVG and 300-dpi PNG
checks pass. Exactly four PDF/SVG/PNG panel groups and four CSV sidecars are
delivered, without an assembled figure or extra supplementary plots.

Reproduce with `Rscript scripts/workflows/15_cf453_diet_impact_figure6.R
--output=results/cf453_diet_impact_figure6_reproduction` (one shell line).
`--plot-only` uses saved `diet_analysis.rds`; `--models-only` fits and exports
without plotting. New full runs require a fresh destination. The validated
Workflow 14 models.rds is the input; if absent, first run Workflow 14 with
`--models-only` rather than its network analysis. Portable packaging uses
`python3 scripts/reports/package_cf453_diet_impact.py`.

## Versioned score-zonation update (16 November 2026)

The parent provenance manifest now defaults to `provenance_sha256.csv` beside
the supplied `--parent` cache, with an optional `--parent-provenance=PATH`.
The historical parent/default behaviour is preserved. Gene eligibility is
validated against the parent ledger rather than a historical 18,486-gene count;
the 28-mouse/56-profile balanced design and all contrast/FDR rules remain fixed.

The versioned `cf453_direct_diet/` stage under
`results/zonation_marker_update_20261116/` uses the newly fitted CF453 model
cache and runs models-only: no older four-panel alternative is regenerated.
Its 18,604 genes yield 163 diet-by-sex, 333 diet-by-zone and 140 three-way
discoveries in the declared pooled family, and 593 omnibus discoveries.
The 267 lipid interactions and four magnitude-ranked example genes are unchanged.
All raw three-way estimates/p values replay the matching updated Workflow 14
fit; this does not assert equality with the historical two-marker-input fit.
Independent actual-data verification passes 170 groups. Cohort-reference
guards distinguish the actual historical `results/zonation/` directory from
the unrelated version-root name `results/zonation_marker_update_20261116/`.
