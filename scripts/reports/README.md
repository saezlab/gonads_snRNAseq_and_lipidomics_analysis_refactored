# Panel and report tools

Run from the repository root. Main scientific execution is Workflow17, not a report script. Existing panels are separate files for manual figure assembly.

| Group | Tools | Current behavior |
| --- | --- | --- |
| Current scientific stages | `stage_zonation_marker_update.R`, `run_zonation_update_fcg.R`, `run_zonation_update_gdx.R`, `cf453_figure6a_layouts.R` | Called by17 with explicit versioned inputs/output |
| Current tables | `prepare_zonation_update_report_sources.py --tables-only` | Public source-table copies; no manuscript required |
| Current Figure6 delivery / visual review | `package_cf453_diet_rewiring.py`, `render_cf453_diet_rewiring.py`, `review_cf453_figure6a_layouts.py` | Defaults follow current release components; recorded manifest paths are resolved externally |
| Optional author-facing audit | `zonation_update_manuscript_audit.py`, `build_zonation_update_report.py`, `package_zonation_marker_update.py`, `compare_cf453_zonation_update.R` | Requires historical references and, for manuscript steps, the private DOCX; explicit17 audit mode |
| Earlier Figure6 presentations | `package_cf453_figure6.py`, `package_cf453_diet_impact.py` | Require an explicit full Workflow14/15 output directory; these are not the current main Figure6 |
| Optional historical comparison PDFs | `build_wdsd_wp2_pdf.py`, `build_wp2_report_pdf.py`, `build_wdsd_pooled_followup_pdf.py`, `build_wp3_report_pdf.py` | Resolve saved WF12/13 inputs through historical release keys |
| Other optional authoring | `build_cyp2c38_ulm_pdf.py`, `build_research_plan_pdf.py`, `render_cosmos_ko_validation_pdf.py` | Separate validation/planning/causal reports, not main scientific model entry points |

Review CLI help before regeneration. Some older PDF builders require an explicit `--font-dir` on a different host; default scientific R execution does not use those tools. Python report rendering requires its documented libraries, whereas repository maintenance uses the standard library.

Do not run a packaging/report command on an immutable release just to validate paths: it may replace a derived ZIP/report. Use a fresh output copy or the read-only organization/path checks. Original embedded paths/hashes are retained as provenance; the registry resolver finds moved files without rewriting those records.
