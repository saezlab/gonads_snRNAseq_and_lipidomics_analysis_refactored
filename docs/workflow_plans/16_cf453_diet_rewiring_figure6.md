# Workflow 16: Diet shifts in selected sex-specific zonation programmes

## Current parents, modes and layout entry (18 November 2026)

Defaults are `analysis_current("cf453_models/models.rds")`,
`analysis_current("cf453_direct_diet/diet_analysis.rds")` and output
`analysis_current("cf453_rewiring")`. Both generated parent caches are mandatory;
the selected Workflow14 parent's adjacent `provenance_sha256.csv` is also read,
unless `--parent-provenance` overrides it. The two parents must contain identical
metadata and modeled RNA. Supplying `--parent` alone does not redirect `--direct`
to a sibling fresh run: specify both explicitly when reproducing elsewhere.
The code does not run14 or15 when a cache is absent.

CLI accepts `--parent`, `--direct`, `--parent-provenance`, `--output` and either
`--models-only` or `--plot-only`. There is no `--validate-only`, bootstrap-count
or seed CLI option. Default/full execution computes all anchor/edge/class/LOO and
999-bootstrap statistics before plotting/reporting. `--models-only` saves this
complete numerical analysis as `rewiring_analysis.rds` plus tables/provenance;
it is not a cheap input/model check. `--plot-only` reads that object from
`--output` and rewrites exports, without reading another supplied parent.
New analyses reject nonempty destinations except a lone `QA/` directory.

Public RNA inputs are CSV transports; expected frozen RDS hashes are verified
against their exact exports with original identity and opened-file provenance
kept separate. Generated model/analysis caches are ignored. A public clone
must run14 models-only, then15 models-only from that parent, then16 with both
fresh parent paths. It does not need the historical30-panel network output.

The separate `scripts/reports/cf453_figure6a_layouts.R` accepts only
`--parent=<rewiring_analysis.rds>` and `--output`; defaults are current rewiring
and `analysis_current("cf453_layouts")`. It reads saved anchors/effects, preserves
selection and creates two independent layout alternatives. It does not fit or
resample. It calls `analysis_assert_not_approved_output()` before writes, so the approved
root/descendants are rejected even when model caches are absent. The declared
default output is an approved lookup location and rendering needs an explicit
fresh `--output`. Outside the protected root there is no general nonempty-output
guard, so existing review exports can be overwritten. Retained-original
hashes are checked through the relocation/CSV identity resolver. PNG/PDF/SVG
packaging/review scripts are separate tools, not implicit statistical stages.

Production parent reuse checks the recorded input-data identities, not every
mutable code/plan row in a historical provenance file. New exports record current
code/docs hashes. Historical strict export-QA can separately compare all recorded
source/plan hashes and therefore report source-snapshot drift after documentation
edits. Preserve approved manifests/results; do not refresh them to hide that drift.

## Current manuscript 6B and 6E

`scripts/reports/cf453_manuscript_panels.R` renders the approved A4 panels from
the Workflow15 `tables/lipid_results.csv` and Workflow16
`source_data/Fig6C_network_gain_loss.csv`. The defaults read the shipped
tables. For a fresh run, pass both `--lipids=` and `--edges=`. Supply an explicit
fresh `--output=` outside the published and legacy roots. Workflow17 supplies
the two paths from its own 15/16 stages. Rendering uses ggplot2 ≥ 3.5, scales,
ragg, digest and ggrepel for labels on alternate selections.

Manuscript 6B uses female and male moderated WD-minus-SD t statistics with
asinh(t/2) display axes. The fixed Figure 5A display rule requires absolute
diet-by-sex t above 4.4 and absolute t below 1 in the opposing sex. The
published categories contain 9 female decreases, 2 female increases and
5 male increases among 1,129 eligible lipids. The low-t band is descriptive.
The saved interaction BH q values remain available in the source table.

Manuscript 6E uses the saved representative and up to six saved partners in
each context. Each pair uses the same lipid ID for the within-sex abundance
effect and both correlations. Each arrow connects open SD to filled WD at a
fixed log2(WD/SD) on square panels with shared scales. The current genes are
Fmo3, Cyp3a41a, Col27a1 and Osbpl3. The source edge selection and its inference
limits are documented below. For other valid selections the renderer lays out
the saved partners, with empty contexts shown explicitly.

The published manuscript exports are under `analysis_current("cf453_manuscript_panels")`.
The native anchor, class and network filenames below retain their Workflow16
numbers. Run `Rscript tests/test_cf453_manuscript_panels.R` to verify selection
counts, arrow geometry and source-table replay.

Requested 2 November 2026 as three independent Figure 6 panel groups closer
to Figure 5. Preserve workflows 14/15 and their outputs. This is a standalone
analysis of 28 approximately 10-month-old wild-type mice: seven per sex and
standard-diet (SD)/Western-diet (WD) group. No original FCG data enter the
analysis and no cross-cohort conservation claim is intended.

## Inputs and RNA selection

Read validated Workflow 14 models.rds and Workflow 15 diet_analysis.rds.
Verify every recorded parent input against Workflow14 provenance; retain its
data-driven eligible-gene/lipid ledgers and the balanced 56-zone-profile roster.
Reuse the eight-cell mouse-blocked limma fit and declared FDR families.

For each selection diet and sex, require sex x zone q<=.05 (Workflow 14
primary family), own-sex zonation q<=.05 (secondary zonation family), and
other-sex absolute zonation t<1. The latter is a descriptive specificity
gate, not evidence of no zonation/equivalence in the other sex. Positive
PP-minus-PC effect defines PP, negative defines PC. Select five genes per
sex/zone, ranking decreasing absolute own-sex zonation effect, then
interaction q, own zonation q, gene symbol. Selection uses only anchor-diet
RNA, not the diet shift or lipid correlations. Never relax gates to fill a
panel. Track the same genes in the other diet regardless of its significance.

## Three panel groups

A: two scatter facets, SD-selected SD->WD and WD-selected WD->SD, each with
up to20 selected genes (underfilled bins retain all eligibility gates). Female PP-PC t is x, male PP-PC t is y, matching Figure5A.
Open endpoints are SD; filled endpoints WD. Arrow starts in selection diet.
Colour encodes anchor sex/zone. Solid arrows have direct three-way diet x
sex x zone q<=.05 from Workflow 15's genes x three direct contrasts family;
other arrows are dashed. Export raw effects and direct interaction evidence.
In this balanced model each gene has exactly the same moderated SE across
the four simple zonation contrasts: its t trajectories are effect trajectories
divided by that common gene-specific SE. Across genes the SE scaling differs.
Selection can induce regression to weaker destination estimates. Arrows
themselves are not tests of diet differences.

B: four-context heatmap (female/male x PP/PC), using the union of A-selected
genes in that context, fixed in both diets. Correlate matched per-zone RNA
and whole-liver lipids across seven mice per sex/diet. No lipid spatial
measurement or duplicate lipid replication is implied. At |r|>=.65, define
E = class pair-hit rate minus pair-hit rate across all other lipid classes.
Display E_WD-E_SD in percentage points. Export both class/rest rates, counts
and denominators, raw density changes and continuous mean-|r| excess change.
This describes relative class association in selected programmes, not a
genome-wide enrichment test. Do not treat dependent edges as independent.
Recalculate nominal Pearson p and global BH over the explicitly selected
context/diet/gene/lipid family for transparency, not selection.

Use 999 independent-arm mouse bootstrap draws (seed 20261102), resampling
matched RNA/lipid vectors together. The same indices apply to both zones
and all features in a sex/diet group. If any modelled gene or lipid is constant
in a context/draw, mark the entire context/draw invalid, preserving fixed
denominators. Export valid/undefined counts, conditional percentile intervals
and directional sign fractions. Optional markers indicate >=.90 bootstrap
directional stability, not significance. Table-only threshold sensitivities
.60/.65/.70/.80 and observed broader all-strict-eligible programme sensitivity
do not create extra figure panels.

C: up to four gene-centred example networks, one per sex/zone. Genes must
come from the same A-selected union. Gained edge: |r_SD|<=.4, |r_WD|>=.8,
and change in |r|>=.5; lost is the reverse. Require all seven high-diet
leave-one-mouse-out correlations |r|>=.65 with the full-data sign and all
seven low-diet omissions |r|<=.65. Export min high-minus-max low strength
(equivalent to minimum across all49 independent-arm omission combinations).
Require >=90% of 999 bootstrap draws defined and >=.90 directional stability
for the change in |r|. Examples remain exploratory, selected from the same
data, without calibrated rewiring p values or causal interpretation.
Rank genes by eligible edge count, maximum absolute strength change, best
anchor interaction q, then symbol. Select up to three gained and three lost
partners, fill to six from remaining eligible edges. Show identical partners'
SD and WD r values, gain/loss colours, gene-above-lipid-boxes motif. Leave an
unsupported context empty rather than lower thresholds.

## Implementation, outputs and verification

Entry: scripts/workflows/16_cf453_diet_rewiring_figure6.R. Core, plots and
report: scripts/core/cf453_diet_rewiring*.R. Output directory:
`analysis_current("cf453_rewiring")` (use a fresh output for reproduction). Save rewiring_analysis.rds, bootstrap
indices, eligibility/selection/edge/class ledgers, three source CSVs, three
PDF/SVG/300-dpi PNG groups, manifest, legends/methods and portable ZIP.
No assembled manuscript figure or extra supplementary plots.

Run every new/modified R script immediately. Independently verify anchor
eligibility/ranking, exact r against stats::cor, matched mouse order, fixed
denominators, class/rest definitions, bootstrap invalidity and pairing,
LOO gates and all49 separation, source-table joins, direct q annotations,
input/code hashes and CLI replay. Render and inspect all PDFs, check bounds,
vector SVG and PNG resolution. Report small n, conditional selection and
resampling limits. Supplied normalized RNA is not raw counts; aggregation,
nuclei/batch/litter and randomisation provenance remain unestablished.

## Historical two-marker three-panel results (2 November 2026)

Both anchor diets support all four five-gene bins: 40 selections, 34 distinct
genes. Fixed context unions are female PP10/PC7 and male PP8/PC9 genes.
Six SD-selected and three WD-selected arrows have the direct three-way
q<=.05 annotation. All 1,129 primary-QC lipids vary in both diets.

The explicit selected family contains 76,772 correlations; none passes BH
q<=.05. Class ΔE spans -33.41 to +72.43 percentage points. All 120 conditional
95% bootstrap ranges include zero. Two cells have >=.90 directional fractions;
the largest shift is female-PC LPG, represented by only one lipid species.
Class-summary paired-diet draws are defined for 993/999 female and 991/999
male draws, using the entire-context gate. Edge-specific valid counts are
separate and can be higher. No class shift is presented as confirmed enrichment.

LOO and bootstrap gates retain 228/136/232/255 display-eligible edges in
female PP/PC and male PP/PC respectively. Ranked examples are Fmo3,
Cyp3a41a, Col27a1 and Osbpl3, with six partners each (Fmo3 five gained/one
lost, others three/three). Their percentile ranges and full conditional
stability ledgers remain available; threshold-defined gain/loss is exploratory.

Independent numerical audit passes 677 check groups, maximum discrepancy
1.67e-15, including all 76,772 correlations, 537,404 omission coefficients,
2,033,964 candidate-edge bootstrap coefficients, complete class bootstrap
arrays/invalidity, all49 omission separation and deterministic selection.
A fresh CLI model run reproduces all 16 saved top-level components exactly.
All three PDFs were rendered and individually inspected; page bounds,
vector-only SVGs and 300-dpi PNG dimensions pass. Source and CLI-export
checks and bundle verification are recorded under QA/ and the completion
manifest. Exact historical reproduction needs matching historical14/15 parent
paths; defaults now use the expanded-marker release. With explicit parents and
a fresh --output,
--models-only exports numerical results, --plot-only regenerates panels from
saved rewiring_analysis.rds. Package with
`python3 scripts/reports/package_cf453_diet_rewiring.py`.

## Historical two-marker Figure6A layout production, 4 November 2026

The user requested two review candidates without replacing the existing
figure: a four-state heatmap and row-wise SD/WD paired dots. Read the saved
Workflow 16 anchors only; no model fitting, reselection, RNA resampling or
network changes. Export under `results/cf453_diet_rewiring_figure6/`
`layout_comparison_20261104/`, using a separate report entry
`scripts/reports/cf453_figure6a_layouts.R` and matching core modules.

Deduplicate the 40 anchor selections to 34 genes, preserving each gene's
sex/zone programme and its SD/WD/both selection membership. Within each
programme, order by decreasing absolute direct three-way effect, then gene;
this is display ordering only. Both views use this same order, four fitted
female/male x SD/WD PP-minus-PC effects (natural-log normalized-expression
differences, not log2FC), and the same direct three-way q<=.05 marker.
Nine of34 genes carry that marker; it tests sex-dependent diet modification
of zonation and is distinct from anchor selection.

Option1: 34 aligned gene rows, Female SD/WD then Male SD/WD columns; one
common symmetric signed effect scale, no row standardization or clipping.
Use a zone-direction palette distinct from female/male colours; show programme
row blocks, SD/WD/both selection annotation and direct-q marks separately.
Option2: four programme blocks in2x2 layout; each gene has separate female
and male tracks, SD open and WD filled, joined by horizontal segments. The
same signed effect axis applies to every block, with zero and PC/PP direction
labels; names are aligned outside the point field. All movement reads SD->WD.

Export each candidate as independent PDF/SVG/300dpiPNG with source CSV,
plus common34-gene ledger, manifest, short legends, provenance and ZIP.
Retain original panels/B/C, analysis object and main bundle unchanged.
Verify exact saved-effect/q joins, deduplication and membership, all136
plotted values in each candidate, common scales, no row transformation,
nine test markers, graphic bounds/vector/dimensions, and inspect both PDFs.

## Versioned score-zonation update (16 November 2026)

Workflow 16 derives parent provenance beside its supplied model cache (or uses
`--parent-provenance=PATH`). All statistical, FDR, anchor and resampling gates
remain fixed. Anchor-bin availability/selected/requested counts are exported;
underfilled bins retain their eligibility gates rather than forcing five genes.
Gene and lipid counts and both 6A layout counts are data-driven; historical
34-gene/9-marker assertions and labels no longer constrain new fitted inputs.
Current defaults now select the approved expanded-marker parents/output;
explicit historical parents remain supported.

The new full result is `results/current/06_western_diet_and_gene_lipid_circuitry/cf453_rewiring/`.
All eight bins still support five anchors: 40 selections, 35 distinct genes,
fixed female PP9/PC7 and male PP9/PC10 programmes. Four SD-selected and six
WD-selected arrows have direct three-way q<=0.05. There are 79,030 selected-family
correlation tests, none BH-significant. Class Delta-E ranges -27.13 to +72.96
percentage points; all 120 conditional 95% ranges include zero and two cells
meet the descriptive directional marker. Class-valid draws remain 993/999 for
female and 991/999 for male contexts.

Full omission/bootstrap gates retain 169/138/242/216 eligible display edges
in female PP/PC and male PP/PC. Representative genes remain Fmo3, Cyp3a41a,
Col27a1 and Osbpl3, with six partners each; individual partner lists are compared
without assuming that unchanged representative names imply unchanged networks.
Independent full numerical replay passes 678 groups. Both 6A alternatives are
under `cf453_layouts/` in the version root, retaining all 140 fitted effects and
ten direct-q marks; independent layout/export checks pass 70 groups. All five
PDF panels were rendered and visually inspected, with vector/bounds/resolution
checks passing. The original panels and statistical artifacts remain intact.

`scripts/reports/compare_cf453_zonation_update.R` writes 15 old/new audit CSVs,
including RNA/FDR, anchor membership, class denominators and network examples.
Matched-programme class comparisons additionally use identical intersection
genes/pairs, separating coefficient changes from programme-reselection effects;
these extra comparisons are observed tables, without new significance tests.

### Historical two-marker layout geometry verification (4 November 2026)

Those original candidates were exported at6.8x9.5in (heatmap) and9x8in (paired dots),
each with a136-row source table. The gene ledger retains34 unique genes,
six selected in both diets and nine direct three-way q marks. The heatmap
uses a common symmetric colour domain +/-2.58409055449904. Dots use one
common asymmetric effect axis covering the complete raw range
[-2.58409055449904,1.30430984547253] with padding and a separate label gutter;
no per-facet numeric rescaling is used. Selection badges and gene order agree.
Both PDFs were rendered and individually inspected, with no off-page text,
vector-only SVGs and exact300dpiPNG dimensions. Independent source/geometry
checks pass70 check groups, recorded in `QA/` within the comparison folder.
Eight original analysis/panel/source/bundle files have retained SHA256 records.
Use `Rscript scripts/reports/cf453_figure6a_layouts.R` to regenerate only these
candidate layouts, followed by `scripts/reports/review_cf453_figure6a_layouts.py`
with the bundled Python runtime for rendering and optional `--package`.
