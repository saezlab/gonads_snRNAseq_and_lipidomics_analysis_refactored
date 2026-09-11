# Manuscript causal-model audit

Run from the repository root: `Rscript scripts/moon/audit_manuscript_causal_model.R`.

All checks passed. The relocated PKN and input RDS hashes match the archived run manifest; the original-ID final graph projects exactly to the compressed scored graph. Re-extracting the selected network from its frozen scores reproduces all 388 signed edges and 153 nodes, with four TFs and three measured DGs. Input signs, score-retention counts, threshold-grid counts and the saved stronger-compression comparison also agree.

This audit does not refit MOFA, recompute MOON scores, modify the archived inputs, or demonstrate end-to-end replay of historical runners. Some runners and RDS metadata retain pre-move paths and require migration before such a replay. `audit_summary.json` and `artifact_manifest.csv` record the audited values and files.
