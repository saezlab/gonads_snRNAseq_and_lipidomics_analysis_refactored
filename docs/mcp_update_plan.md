# MCP Server Update Plan

**Date**: 2026-02-19 (created) / 2026-03-04 (completed)
**File**: `mcp/fcgl_server.py`
**Status**: ✅ COMPLETED — All Priority 1 and Priority 2 items implemented.

---

## Completed: Priority 1 — New Tools Added

### Tool 1: `query_sex_specific_zonation_genes` ✅
- Data file: `results/zonation/sex_specific_zonation_genes_summary.csv` (138 genes)
- Parameters: gonad (Testes/Ovaries), zone (pericentral/periportal), gene, top_n
- Filters by category column, sorts by interaction t-value

### Tool 2: `query_variance_analysis` ✅
- Data files: `variance_analysis_results.csv`, `variance_analysis_by_class.csv`, `variance_analysis_summary.csv`
- Parameters: lipid, level (summary/class/lipid), significant_only, top_n
- Supports three analysis levels with proper formatting

### Tool 3: `query_lipid_species_zonation` ✅
- Data file: `results/zonation/lipid_species_zonation_association.csv` (390 lipids)
- Parameters: lipid, zone (periportal/pericentral/all), gonad_bias, top_n
- Shows TPC/TPP/OPC/OPP correlations with top correlated gene

### Tool 4: `query_lipid_network` ✅ (new, not in original plan)
- Data files: `significant_lipids_network_SIF.csv` (6014 edges) + `_ATT.csv` (794 nodes)
- Parameters: lipid, top_n
- Shows LINEX2 enzymatic conversion relationships between gonad-dependent lipids

### Tool 5: `query_tf_gene_lipid_network` ✅ (new, not in original plan)
- Data files: `hepatocyte_TF_gene_lipid_network_compressed_SIF.csv` (140 edges) + `_ATT.csv` (57 nodes)
- Parameters: query (gene/TF/lipid), node_type filter
- Shows TF→Gene→Lipid regulatory cascade from MOFA integration

---

## Completed: Priority 2 — Enhancements

### Enhancement A: Expanded `get_zonation_info` ✅
- Added sex-specific zonation gene counts (138 genes by category)
- Added top 10 zonation-associated pathways from `rna_zonation/pathway_scores_df.csv`

### Enhancement B: Expanded `get_dataset_summary` ✅
- Added precise finding summaries for all workflows (02-08)
- Added complete list of all 17 available tools organized by category
- Added paper title and sample group sizes

### Enhancement C: Expanded `get_validation_results` ✅
- Added per-gene scores from `global_network_validation_scores.csv` (3895 scores)
- Added per-gene scores from `zonation_network_validation_scores.csv` (690 scores)
- Shows top hits per knockdown gene with p-values

---

## Priority 3 (Optional / Not Implemented)
These remain for future consideration if needed:
- `query_gdx_validation` — GDX hormone replacement validation
- `query_zonation_networks` — Zone×gonad-specific gene-lipid correlation networks
- `query_lipid_structural_model_lipid_only` — Lipid-only structural effects

---

## Final State

| Metric | Before | After |
|--------|--------|-------|
| Total tools | 12 | 17 |
| CSV files loaded | 8 | 21 |
| Workflow coverage | WF02-04 (partial WF05-07) | WF02-08 (comprehensive) |
| Regex safety | No escaping | `re.escape()` for lipid names with parentheses |

---

## Bug Fixes Applied
- Fixed regex injection in lipid name searches — `PE(19:0_20:3)` now searches correctly
- Fixed `false` → `False` Python boolean in tool schema
