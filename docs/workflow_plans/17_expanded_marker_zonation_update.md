# Workflow 17: Current expanded-marker zonation and diet pipeline

## Maintained entry and release

Entry: `scripts/workflows/17_expanded_marker_zonation_update.R`.
`analysis_path("current_root")` selects the published `results/current/` layout. The previous complete run resides in the ignored archive. A fresh scientific run must supply
`--output=results/runs/NAME`. Before any staging or writes, full scientific and
full manuscript-audit modes call `analysis_assert_not_approved_output()`, rejecting
the published root, legacy archive and their descendants even when ignored caches are absent;
home-directory and symlink aliases are resolved; unresolved dangling symlinks
are refused. Existing `cf453_models/models.rds` or
`fcg05/zonation_ttop_list.RData` additionally trigger the completed-run guard.
Staging and internal FCG/GDX components retain their checkpoint/cache guards.
There is no general resume mode for a partial release.

`stage_zonation_marker_update.R`, the GDX driver and `zonation_update_*` helpers
are internal components of this protected chain, not independent public rerun
entries with the same approved-output guarantee. Their local cache guards are
not a substitute for the main17 guard. The maintained standalone FCG driver
`run_zonation_update_fcg.R` also has the shared full-run guard; use the public
05/17 entries and documented modes for ordinary reruns.

Source loading defines `zu17_run()` without executing it. Run from the repository
root. The public scientific default does not require the private manuscript,
source Seurat objects or a Codex runtime. Python is chosen by `--python`, then
`ZONATION_PYTHON`, then system `python3`; `.local/R-library` is only an optional
local R-library preference, not part of the distributed environment.

## Actual upstream prerequisites

| Input/artifact | Producer and role | Public-clone status |
|---|---|---|
| `current_fcg_zonation_transport`, `current_gdx_zonation_transport`, `current_wdsd_zonation_transport` | Supplied expanded-marker counts/data/scale CSV transports | Shipped manifests and all payloads under `data/RNA/pseudobulk_csv/`; original RDS identities retained only as provenance |
| `wdsd_hepatocytes_transport` and `data/metabolomics/wdsd/lipidomics_cf453.csv` | Mandatory additional Workflow14 intake, even in models-only mode | Shipped; not replaced by `--source-root` or the staged WDSD zones |
| `data/metabolomics/Lipid_data_file_JPC_cleaned.csv` | Original FCG whole-liver lipids used by the current FCG gene-lipid calculations | Shipped; loaded directly, with the established log2 transform |
| `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData` | Workflow03 gonadal lipid effects used by class/sex displays | Required full-run upstream output; shipped current artifact |
| `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/enrichment/lipid_class_enrichment_scores.RData` | Workflow03 lipid-class enrichment display source | Required full-run upstream output; shipped current artifact |
| `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/structural_model/lipids_with_structure.csv` | Workflow06 structural reference used by the zonal07 structural-consistency plot ranges | Required when that structural branch has matched data; supplied current artifact; regenerate06 if absent |
| `results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv`, `knockdown_ttop_list.RData` | Workflow07 independent KD fit artifacts, reused without refitting raw experiments | Required: CSV supplies KD t statistics; RData is also required by validation/provenance checks |
| `support/TF_regulons.RData`, `support/Gonadal_Hormone_Genes.csv`, `support/m5.all.v2023.2.Mm.symbols.gmt` | Supplied TF/pathway resources used by FCG05 | Shipped; pathway ULM is recalculated, not read from the old05 cache |
| `engine_fcg05`, `engine_validation07` | Frozen internal implementation engines, with versioned adapters | Shipped code; do not execute as unreviewed public legacy entries |

Workflow17 does not execute03,06 or the independent KD model producer. A
successful entry `--validate-only` is not proof that all these upstream result
files or plotting dependencies are present: its checks are deliberately narrower.
The actual full-run engine reads the03 results and06 structural reference later.
No original whole-hepatocyte MOFA fit or historical FCG zonation result is needed
for the default current scientific branch.

## Modes, validation and overwrite behavior

Accepted flags are `--validate-only`, `--report-only`, `--with-manuscript-audit`;
values are `--output`, `--source-root`, `--python`. Unknown/duplicate options and
simultaneous validation/report-only are rejected.

- **Default science:** stage three current zone inputs; run current FCG05 and
  zonal07 adapters plus paired current GDX; run14 `--models-only`,15
  `--models-only`,16 full and the two6A layouts; run current scientific QA.
  If Python is available, copy the formal source tables with `--tables-only`.
  It does not run the old/new comparisons, private manuscript inventory or
  report/bundle generation.
- **`--validate-only`:** verify the three current CSV transports/layers and
  balanced zone counts, check the two engine files and saved KD CSV/RData.
  It bypasses the approved-output write guard, writes nothing and does not fit models, enumerate bootstrap draws or
  validate every downstream03/06/resource prerequisite.
- **`--validate-only --with-manuscript-audit`:** additionally check Python and
  selected historical FCG/CF453 cache files. It still does not prove that every
  private manuscript/preprocessing/report prerequisite is restorable.
- **`--with-manuscript-audit`:** execute the scientific stages with old/new
  FCG/GDX/CF453 comparison branches, then the retrospective QA, private inventory,
  author report and delivery bundle. Needs historical matrices/networks/model
  caches and the private DOCX/preprocessing reference sources, many excluded
  from Git. It is not the public-clone default. The same `--source-root` is forwarded to
  preprocessing-reference staging, which looks for `james_scripts_20261104/` and
  `Zonation_Genes.csv`. A fresh audit therefore needs a root containing both
  acceptable zone inputs and those references, or exact pre-staged snapshots;
  the default public CSV root alone does not contain the reference sources.
- **`--report-only`:** implies manuscript-audit mode and skips model/staging
  stages. It reads saved FCG/GDX/14/15/16 caches and comparison tables, runs
  scientific/retrospective numerical QA, stages source tables/references and
  rewrites inventory/report/bundle. This explicit mode bypasses the approved-output write guard. It is not a read-only command and cannot
  reconstruct missing comparisons/model caches. Main model fits are not rerun.

`--source-root` defaults to registry `pseudobulk_csv_root`. For each cohort,
staging prefers `current_<cohort>_zonation/manifest.json`, then
`<cohort>/manifest.json`, then an explicit local `<cohort>_zonation_pseudobulks_score.rds`.
It copies and SHA-verifies all transport payloads into `<output>/inputs/<cohort>/`.
FCG/GDX scale-layer CSVs reproduce the established sample-name conventions;
Workflow14 receives the staged WDSD manifest. New provenance hashes actual CSV
files separately from the original source-RDS identity. No normalization changes.

## Outputs and parent graph

Fresh run components under `--output` are `inputs`, `fcg05`, `fcg07`, `gdx_zonation`,
`cf453_models`, `cf453_direct_diet`, `cf453_rewiring`, `cf453_layouts`, `runtime`,
`logs` and QA/manifests. `cf453_models/models.rds` feeds15 and16;
`cf453_direct_diet/diet_analysis.rds` additionally feeds16;
`cf453_rewiring/rewiring_analysis.rds` feeds the separate6A layout entry.
No stage infers a sibling fresh parent unless the orchestrator passes it.

Curated panels, source tables, methods/legends and public manifests are supplied.
Generated RDS/RData model, bootstrap and runtime caches are local/ignored and
must be regenerated into a fresh run or restored before saved-result modes.
The original report and numerical replacement ledger are kept in the local archive. `--report-only` requires `--output` pointing to a complete writable saved run outside the published tree and legacy archive.
Current source tables exist independently of the optional formal-copy step.

Figure4A/S4B–C composition/nucleus QC and exact executed retained-marker
confirmation remain James's tasks. Pseudobulk expression cannot recover those
inputs. Current diet association panels remain exploratory, with paired RNA
zones and whole-liver lipid measurements.

## Historical authorization and approved production procedure (16 November 2026)

The following records the original approved run and its original audit/report
procedure. Its original-datasets paths and RDS copy procedure are historical;
the maintained runtime/dependency contract is specified above.

Authorized 16 November2026 after the read-only input and dependency audit.
The user requests the affected reruns and a final report tracing every affected
figure panel and manuscript numerical replacement. Preserve all original
inputs, analyses, figure bundles and the manuscript itself.

### Scope and inputs

Output root: `results/legacy/pre_reorganization_baseline/zonation_marker_update_20261116/`.
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

### Execution branches

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

### Panel and manuscript tracing

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

### Audit rules

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

### Observed results and execution record

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

The original production entry validated setup without writing. That approved
full run staged inputs and executed the FCG/GDX/CF453 stages before its audit
and report. The maintained CLI now separates scientific and manuscript-audit
modes as described above. The output root retains104 historical file SHA256 records;
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
