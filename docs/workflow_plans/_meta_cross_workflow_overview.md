# Cross-workflow dependency overview

Audited against execution code on **18 November 2026**. The current release is `expanded_markers_20261116`; its executable selection is [the release registry](../../scripts/config/analysis_release.json). This overview replaces the layered historical overview, byte-preserved in [the pre-audit snapshot](../../miscellaneous/documentation/dependency_audit_20261118/README.md).

## Which documentation to use

- [PIPELINE.md](../../PIPELINE.md): scientific graph, orchestration graph, fresh-checkout starting points and commands.
- [Complete workflow reference](../WORKFLOW_REFERENCE.md): all numbered entries, compatibility phases,6A layouts and optional PKN/MOON branches; exact inputs, producers, modes, outputs and validation scope.
- [Dependency CSV](../workflow_dependencies.csv): filter by consumer, producer or availability. [Workflow catalog](../analysis_catalog.csv): one row per workflow branch.
- [JSON contract](../workflow_dependencies.json): editable source for both generated indices and the complete reference. Keep it aligned with source and these per-workflow plans.
- [Result lookup](../AGENT_GUIDE.md) and [panel navigation](../../results/current_panels.md): inspection of approved outputs without fitting.
- [Archive/provenance](../REPOSITORY_ORGANIZATION.md) and [causal map](../CAUSAL_ANALYSIS_MAP.md): frozen dependencies, restoration and optional extension details.

## Main handoffs

| Consumer | Data/artifact producers | Important boundary |
| --- | --- | --- |
| 01 preparation | Externally supplied original cell-type Seurat objects | 00 is a legacy HPC reference; an explicit import/name handoff is needed. Full01 requires excluded local originals. |
| 02 RNA;03 lipid;04 integration | Supplied01 cell-type CSV/metadata plus branch-specific lipids/resources | Public downstream runs can use shipped CSVs.02 DE is not an input to04's fit or current05. |
| 06 lipid structural | 03 lipid DE / curated lipid annotations | Produces `lipids_with_structure.csv`, needed by07 global and current zonal validation. |
| 07 global validation | 04 global networks,06 structure, independent KD measurements | Independent KD fit artifacts are preserved and reused by current zonal07; these are not rerun by17. |
| 07 Figure 3F–G replay | Saved 07 knockdown t-table and 04 global network | Produces matched-only scores and panels from the 210 unique lipids per gene. The full global validation also uses this plotting branch. |
| 05 FCG zonation | Current FCG zone transport, 03 lipid DE/class scores, original FCG lipids, TF/pathway resources | Produces the FCG models, networks and panels. |
| 07 zonal validation | Workflow 05 checkpoint and adjacent networks, saved KD table/fit and 06 structure | Uses the FCG result in the same root or the checkpoint passed with `--parent=`. |
| 05 paired GDX | Current GDX zone transport | `--cohort=gdx` selects the paired supporting analysis. |
| 14 standalone diet | Current diet zone transport, diet whole-hepatocyte transport, diet whole-liver lipids | Whole-hepatocyte intake is mandatory even in models-only mode. No original-cohort statistical dependency. |
| 15 direct diet tests | 14 `models.rds` and adjacent provenance | Manual runs must explicitly select a new parent. |
| 16 current diet presentation | 14 model cache/provenance and15 `diet_analysis.rds` | Neither parent is inferred from the latest output directory. Layouts use16 `rewiring_analysis.rds`. |
| Manuscript Figure 6B and 6E | 15 `tables/lipid_results.csv` and 16 `source_data/Fig6C_network_gain_loss.csv` | 17 passes the two saved tables from its own run. The report can also render from the shipped tables. |
| 08 variance | Raw lipid measurements and metadata | Standalone; not03's matched-RNA lipid subset. |
| 11 Cyp2c38 in vivo | Original whole-hepatocyte inputs,02 checked groups,04 fixed network, independent KO lipids | No zonation dependency.07b/07c are compatibility phase wrappers. |
| 12/13 historical comparisons | Explicit historical two-marker inputs and networks |05b/05c are compatibility wrappers. Preserve both sides of the historical definition;13 has no matching expanded-marker CF459 input. |
| 09 human;10 Xenium | External human/signature resources; independent raw imaging |09's signature is not emitted by current02.10 has documented execution blockers and is not supplied by02's panel-restricted volcano. |
| Optional PKN/MOON | Frozen factor inputs, priors, fitted/extracted states and branch-specific references | Follow the causal map; some archived pilots remain live frozen dependencies. |

The current FCG/GDX and diet branches run individually. Optional Workflow 17
coordinates the same entries and runs the anchor layouts and manuscript 6B/6E
exports. The diet cohort is analysed independently because of age differences.
RNA zones are paired within mouse and whole-liver lipids are measured once per mouse.

## Reproducibility limits and checks

Shipped tables and plots can be inspected without ignored caches. A shipped derived prerequisite can be consumed without rerunning its producer. Saved-result rendering and downstream chaining may require ignored `.rds` caches/provenance; regenerate into a fresh destination or restore the exact artifacts. A successful narrow validation command does not establish that every later plotting or private-audit input exists.

Workflows01–04/06/08 execute on source and use fixed output roots. The [parallel manuscript replay](../MANUSCRIPT_RESULTS_REPLAY.md) routes those paths within a separate workspace. It selects Workflow04's model-reuse-only mode and verifies all eleven retained HDF5 files. It routes the new Workflow02 PCA and Workflow04 network to Workflow11, the new Workflow03 lipid DE to Workflow06, the new Workflow04 network and Workflow06 annotations to global Workflow07, and the new lipid and KD artifacts to Workflow17. Workflow11 validation fits the observed model in memory, whereas Workflow17 validation checks inputs. Explicit report refresh modes can write to existing trees.

Historical source-snapshot QA compares recorded code/plan hashes with the current working tree and may flag intentional documentation drift. Preserve those hashes. Repository organization checks separately verify approved scientific bytes. The [audit record](../provenance/workflow_documentation_audit_20261118.json) records this distinction, documentation coverage and approved-output guard tests.

Figure4A/S4B–C, nucleus composition/counts and final executed marker-list confirmation remain James's tasks. The pseudobulk pipelines do not reconstruct those missing inputs.

After changing a workflow, update its plan and JSON contract and run `python3 scripts/maintenance/build_workflow_documentation.py`, then `--check`, organization verification and strict GitHub readiness. The generator checks declared contracts and coverage; it does not infer every possible read from code or execute scientific analyses.
