# Directionality-audited STITCH v5 exploratory overlay

This directory is a separate exploratory COSMOS PKN overlay. It removes every
orientation-unknown STITCH record from the frozen base exploratory PKN and
replaces them with raw STITCH v5 mouse records only when the chemical is raw
partner A, `a_is_acting=t`, the raw action is explicit activation/inhibition,
the action score is at least 700, and detailed-link text mining is below 700.
Signs are derived from the raw `action`, never `mode`.

`pkn_exploratory_stitch_directional.csv` and `.rds` are the signed graph for a
future explicitly named sensitivity analysis. This overlay was not used in the
reported TF→lipid or lipid→TF MOON fits. It is not a substitute for
`pkn_lipid_curated` or `pkn_lipid_rules`.

Read `stitch_directional_summary.json`, `stitch_directional_manifest.json`,
`stitch_directional_validation.json`, and `r_validation.json` before use.
`stitch_directional_evidence.csv`, the mapping review and exclusions retain the
raw records and all admission decisions. STITCH v5 action/link files were last
modified on 16 August 2016; source-record directionality does not establish
currentness, measurement identity, tissue context or experimental causality.
