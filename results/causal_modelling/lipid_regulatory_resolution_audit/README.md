# Regulatory chemical resolution audit

Run `python3 scripts/moon/ancillary/audit_regulatory_resolution_candidates.py --output NEW_EMPTY_DIRECTORY`. This audit uses explicit ChEBI subclass relations to detect specific regulatory molecules beneath measured composition-level terms. It does not change chemical identities, measurements, networks or scores.

`specific_regulatory_candidates.csv` retains the exact subclass statements, original mapping status, formula checks, resource tier and evidence IDs. `resource_summary.csv` counts unique feature opportunities, excluding formula-incompatible rows. Empty SwissLipids abbreviations explain some candidates missed by the earlier abbreviation-only search. These candidates require explicit identity assumptions or analytical confirmation before use.
