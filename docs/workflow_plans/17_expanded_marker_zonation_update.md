# Workflow 17 Expanded marker zonation update

## CSV-only public input contract - 18 November 2026

Public scientific defaults now use the seven registry-selected transports under
`data/RNA/pseudobulk_csv/`; no RNA/Seurat source RDS is required or included in Git.
Each transport is a small JSON manifest, file-index CSV, complete metadata/identity
CSVs and split counts/data/scale.data CSVs, each below50,000,000 bytes. Import
verifies all file hashes and exact IEEE754 content hashes and preserves the full
original gene/profile order and metadata types. No renormalization occurs.

Default `--source-root` is `pseudobulk_csv_root`. Staging copies verified manifests
and payloads into fresh `inputs/{fcg,gdx,wdsd}/`, then reproduces the existing FCG/GDX
scaled CSV conventions exactly. Workflow14 receives the staged WDSD manifest;
GDX's unchanged fit reads the staged plain matrices. An explicit local source root
may still contain original RDS files, but a CSV manifest always takes precedence.
Completed release inputs/results are protected against staging/fit overwrites.

`--validate-only` verifies the current CSV transports without source-RDS reads.
The approved cache bridge retains original RDS identities, verifies their source
SHA against the complete CSV export, and records actual opened CSV hashes
separately. Historical input identities are not rewritten. Actual validation
passed85 consumer/provenance checks under a source-RDS rejection guard, including
all current/historical CF453 and CF459 intake values and23 approved scientific
invariants. FCG/GDX staging tables were identical to approved tables; Workflow01's
targeted CSV conversion also matched. No scientific model or plot was regenerated
for this representation-only change. QA is local under
`miscellaneous/local/pseudobulk_csv_migration/`; export roundtrip proof is recorded
in `docs/provenance/pseudobulk_csv_roundtrip.json`.

Unbound provenance for an existing explicit local RDS first hashes its actual
bytes. A changed same-name file cannot inherit an older CSV export identity;
unmatched files retain direct-file provenance. Expected-SHA frozen-cache checks
remain CSV-first and source-RDS-free. A temporary changed-basename fixture in
`tests/test_analysis_provenance_identity.R` verifies both behaviors without
altering source objects or scientific results.

---

## Default release contract - 16 November 2026

`scripts/config/analysis_release.json` selects the authoritative expanded-marker
release at `results/zonation_marker_update_20261116/`. Its approved files stay
immutable. Current inputs and preprocessing references are repository-local under
`data/RNA/zonation/current/`; no external original-datasets folder is needed for
normal scientific execution. Entry validation performs no output writes.

A normal reproduction command is:
`Rscript scripts/workflows/17_expanded_marker_zonation_update.R --output=results/runs/expanded_markers_reproduction`.
The default branch stages current inputs, reruns current FCG05/zonal07 and paired
GDX, builds CF45314/15 model caches, runs16 and both6A layouts, and runs current
scientific invariant QA. It skips old/new comparison stages and historical fits.
Original05/07 engines are frozen under `scripts/core/engines/`; validated adapters
select the appropriate scientific branches. Independent KD/global validation can
be regenerated separately with Workflow07 `--global-only --output=<fresh>`.

`--with-manuscript-audit` additionally requests the retrospective baselines,
private source manuscript inventory, report and bundle. `--report-only` explicitly
refreshes saved audit outputs. Those modes require the separately documented
archived baseline objects and private manuscript. Python is selected with
`--python`, `ZONATION_PYTHON`, then `Sys.which("python3")`; no Codex runtime is
hardcoded. Available Python stages the13 formal table copies with `--tables-only`
in science mode. All scientific CSV sources exist independently of that optional
formal-copy step. Approved reports are not refreshed by default lookup/validation.

Promotion verification included an actual isolated current-only FCG/GDX rerun:
48 scientific/source equivalence checks pass against the approved output;
23 current-release invariant checks pass for both approved and fresh paths.
Both current and historical CF453 cached parent input hashes pass without a new
bootstrap run. Separate historical fixture tests explicitly retain their old
marker-order assertions and data role. Original section/output descriptions below
are retained as the approved16November audit's provenance, not current defaults.

---

Authorized 16 November2026 after the read-only input and dependency audit.
The user requests the affected reruns and a final report tracing every affected
figure panel and manuscript numerical replacement. Preserve all original
inputs, analyses, figure bundles and the manuscript itself.

## Scope and inputs

Output root: `results/zonation_marker_update_20261116/`.
Updated source files reside in the user-specified original-datasets
snRNAseq_seurat_objects folder: fcg/gdx/wdsd_zonation_pseudobulks_score.rds.
Copy bytes to versioned inputs/{fcg,gdx,wdsd}/, record SHA256 and source paths.
Capture historical result/input hashes before any production execution.
FCG86, GDX44 and CF45356 paired profiles have20123 genes and unchanged
profile/mouse/group identities; counts and supplied normalized layers differ.
FCG/GDX extraction uses existing scale.data conventions; CF453 retains the
validated normalized data-layer model. Do not introduce TMM/voom changes
into this assignment-only comparison. Keep the models and FDR definitions
from workflows05/14/15/16; GDX's repaired branch uses paired mouse blocking
and explicit four-group contrasts, separately labelled from legacy hormone
replacement notebooks that do not match the supplied cohort.

## Execution branches

1. Stage inputs and audit RNA layers, identities, paired zones and normalization.
   Generate only the updated FCG zone CSV pair, not unrelated WF01 outputs.
2. FCG: isolated execution of WF05 with versioned input/output paths and fresh
   pathway ULM rather than the old cache. Recompute seven contrasts, TFs,
   threshold optimization, candidate categories, four gene-lipid networks and
   class associations. Preserve old-threshold/old-candidate comparisons.
3. Validate refreshed zonation networks against unchanged knockdown results
   (WF07 sections8/9 and relevant structural consistency). Explicitly preserve
   prespecified male-PP Mtmr11 alongside any adaptive context sensitivity.
4. CF453: WF14 models-only on updated zones; WF15 models-only from that cache;
   full WF16 from the new models/evidence; both revised6A layout candidates.
   Route provenance with each parent and derive counts from data, not old
   observed18486/34/nine counts. Keep selection gates and999 paired-mouse
   bootstraps; permit underfilled anchor bins/unsupported network contexts.
5. GDX: a safe standalone paired analysis of the22 mice's sex x surgery x zone
   groups and old/new descriptive comparison. Never source the legacy converter
   that writes FCG filenames or the legacy hormone-replacement notebook.
6. Independent QA of new estimates, contrast directions, q families, source
   joins, paired matching, fixed association denominators and network selection.
   Render generated PDFs and inspect main-panel exports. Recheck old hashes.
7. Final readable report plus machine-readable panel manifest, manuscript
   numerical replacement ledger, old/new comparisons, current-source hashes,
   changed-panel bundles and explicitly unresolved deliverables.

## Panel and manuscript tracing

- Figure4A and S4B/C are explicitly deferred to James by the user during this
  run. Do not generate these panels or infer nuclei/denominators from expression.
  The report must assign James the nucleus counts/proportions/composition tests
  and their linked manuscript replacements; this does not block the pseudobulk work.
- Figure4B: both axes are zonation-model contrasts; verify signs against legend.
- Figure5A/B: refit/selection/threshold/counts and schematic labels.
- Figure5C/D: class association rates and representative network coefficients.
- Figure5F: refreshed male-PP correlation versus unchanged KD lipid effects;
  verify actual sex-specific n and ULM/Pearson statistics. Figure5E unchanged.
- S5A and TablesS10/S13/S14/S15 are recomputed where their source is available.
- Figure6A/B/C and both6A alternatives are regenerated from the revised
  CF453 cache. Panel letters remain current export labels before manual assembly.
- Independent spatial experiments, whole-hepatocyte DE, MOFA, lipid-only
  results and WF11 KO remain unchanged. WF12 zonal rebase is conditional on
  inclusion of those comparisons, and must update reference and RNA together;
  pooled WF12 remains unchanged. WF13 expanded-marker rebase is on hold
  without equivalent CF459 zone inputs. Figure3 editorial revisions are separate.

## Audit rules

Extract exact manuscript paragraphs/captions/table references, their historical
numbers and source analysis. Every flagged value must have a new measured
replacement, an explicit unchanged status, or a stated missing-input reason.
Differentiate scientific changes from pre-existing sign/sample-size errors,
rounding-only differences and newly generated supporting analyses with no
existing manuscript claim. The preliminary99.3% statement is not an acceptance
criterion. Preserve Methods/report caveats on observational correlations and
the role of scoring-panel genes in zone definition.

The GDX supporting branch uses scale.data, genes detected in>=5profiles and
nonzero variation, mouse duplicateCorrelation and default limma eBayes.
Both old/new inputs are refitted identically; eleven explicit contrasts have
separate four-simple-zone, four-conditional-GDX and three-interaction BH families.
This is supporting analysis of the available four vehicle groups, not a replay
of the mismatched legacy hormone-replacement notebook or a new main-figure claim.

Run every modified/new R script immediately. Inputs and results outside the
versioned run are read-only. Any transform of production code is saved with
an explicit patch manifest and code hash. Do not claim blocked composition
panels have been regenerated without nucleus-level data.

## Observed results and execution record

FCG historical model replay reproduces all seven archived contrasts (maximum
t discrepancy4.1e-12). The repaired zero-category threshold implementation
also reproduces the archived5.4 optimum. Updated reoptimization selects4.4,
giving173 genes (93TPC/39TPP/17OPC/24OPP), versus78 originally and82 when
the updated data retain the historical5.4 cutoff. Associated species are580
of1254 (46.3%), versus442 originally and452 at the fixed historical cutoff.
Figure4B's actual FDR definition reproduces the archived4276/2783 baseline;
the revised values are4449 sex-FDR genes with2901 also zone-FDR (1539PC,
1362PP). Do not substitute the unrelated whole-hepatocyte WF02 denominator.

Prespecified male-PP Mtmr11 validation replays the archived nativeULM exactly
on all661 KD lipid targets (zero prediction weights for unmatched targets),
-7.211103/1.5276e-12 -> -6.791136/2.4875e-11. The plotted scatter has345
matched lipid species and22 reference animals; these are three distinct units.
Experimental KD t-values remain unchanged. Gene-programme class fractions
and the historical butterfly's lipid-abundance-sex fractions are both exported
with explicit, different denominators and a legend clarification.

CF45314/15/16 and both6A alternatives completed:18604 eligible genes,
direct diet x sex163, diet x zone333, three-way140, omnibus593; unchanged
lipid interaction267. Forty selections give35 unique genes and ten direct
marks. Four network representatives persist, with18/24 selected edges retained.
All120 class bootstrap ranges still include zero and none of79030 selected
correlation tests passes BH. Independent suites170/678/102/70, exact999-draw
CLI replay and all five main PDF visual checks pass.

The repaired GDX branch completes for old/new inputs,18176 ->18410 eligible
genes; no updated primary GDX-zone interaction passes its declared BH family.
It remains supporting analysis with no current manuscript-panel replacement.

Main report outputs are report/zonation_update_report.pdf and.html, the
89-claim manuscript_numeric_inventory.csv and new_replacements.csv with
additional24 historical pictured-edge coefficients. Panel/source navigation
uses report/panel_trace.csv and fcg_panel_manifest.csv. Source tablesS10,
S13,S14,S15 are copied under tables/ with hashes; S11/S12 remain unchanged.
Figure4A/S4B-C and associated composition statements are explicitly deferred
to James, as are final retained-marker-list Methods provenance details.

Entry `scripts/workflows/17_expanded_marker_zonation_update.R --validate-only`
checks setup without writing. A full run requires a fresh --output and invokes
staging, the isolated FCG/GDX branches and CF45314-16/layouts before audit and
reporting. --report-only refreshes independent checks, inventory and report
from saved results. The output root retains104 historical file SHA256 records;
bundle creation rechecks them. The original manuscript is never edited.

Before inventory/report generation, `scripts/reports/prepare_zonation_update_report_sources.py`
stages the 13 Table S10/S13/S14/S15 CSV copies and writes `tables/table_manifest.csv`.
It also copies and SHA256-verifies the three supplied James preprocessing
scripts, candidate marker list and source manuscript under
`inputs/preprocessing_reference/`, using the recorded provenance origins and
frozen hashes. Derived table copies may refresh; differing frozen reference
bytes are refused. The helper is idempotent, accepts `--output`, and accepts
an optional `--source-root` for the snRNAseq source folder. Existing per-run
reference provenance remains authoritative.
