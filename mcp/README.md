# FCGL Analysis MCP Server

A Model Context Protocol (MCP) server that provides tools for querying the Four Core Genotypes Liver (FCGL) multi-omic analysis results.

## Overview

This MCP server allows Claude (or other MCP-compatible AI assistants) to dynamically query:

- **Differential Expression Results**: RNA-seq and lipidomics DE results across contrasts and cell types
- **Lipid Class Enrichment**: Enrichment of lipid classes by gonadal condition
- **Lipid Conversion Network**: LINEX2 enzymatic conversion relationships between significant lipids
- **MOFA Integration**: Gene-lipid correlation networks and TF→Gene→Lipid regulatory cascades
- **Lipid Structural Preferences**: Chain length and saturation effects on gene-lipid associations
- **Hepatocyte Zonation**: Periportal vs pericentral spatial differences, sex-specific zonation genes
- **Validation**: siRNA knockdown experimental validation of network predictions
- **Variance Analysis**: Lipid variance differences between gonadal groups

## Installation

### 1. Install Python Dependencies

```bash
cd mcp
pip install -r requirements.txt
```

### 2. Configure for Claude Code

The repository includes a `.mcp.json` configuration file. Claude Code will automatically detect it.

Alternatively, add manually:

```bash
# From the repository root directory
claude mcp add --transport stdio fcgl-analysis -- python mcp/fcgl_server.py
```

## Available Tools (19 total)

### Gene Expression (WF02)

| Tool | Description |
|------|-------------|
| `query_gene` | DE stats for a gene across all contrasts and cell types |
| `get_top_genes` | Top DE genes for a specific cell type and contrast |
| `compare_gene_celltypes` | Compare a gene's effect across all 5 liver cell types |
| `count_significant_features` | Count genes/lipids passing significance + FC thresholds (works for both) |

### Lipidomics (WF03)

| Tool | Description |
|------|-------------|
| `get_lipid_de_summary` | Comprehensive lipid DE summary (counts, enrichment, GDX correlation) |
| `query_lipid` | DE stats for lipids (gonadal effect) |
| `query_lipid_class_enrichment` | Lipid class enrichment scores (male vs female bias) |
| `query_lipid_network` | LINEX2 lipid-lipid conversion network |

### Multi-Omic Integration (WF04)

| Tool | Description |
|------|-------------|
| `query_gene_lipid_network` | MOFA Factor 1 gene-lipid correlations |
| `query_tf_gene_lipid_network` | TF→Gene→Lipid regulatory network |

### Lipid Structure (WF06)

| Tool | Description |
|------|-------------|
| `get_structural_associations` | Gene-lipid structural preferences (chain length, saturation) |

### Zonation (WF05)

| Tool | Description |
|------|-------------|
| `get_zonation_info` | Zonation analysis overview with pathway enrichment |
| `query_gene_zonation` | Gene zonation bias (periportal vs pericentral) |
| `query_lipid_zonation` | Lipid class zonation associations |
| `query_sex_specific_zonation_genes` | 137 genes with sex-specific zonation |
| `query_lipid_species_zonation` | Per-lipid species zonation associations |

### Validation (WF07)

| Tool | Description |
|------|-------------|
| `get_validation_results` | siRNA knockdown validation with per-gene scores |

### Variance Analysis (WF08)

| Tool | Description |
|------|-------------|
| `query_variance_analysis` | Lipid variance differences (summary, class, or lipid level) |

### Overview

| Tool | Description |
|------|-------------|
| `get_dataset_summary` | Full dataset overview with key findings and available queries |

## Example Queries

Once connected to Claude, you can ask questions like:

- "What is the differential expression of Cyp2b9 in hepatocytes?"
- "Show me the top 20 upregulated genes in hepatocytes for the gonadal effect"
- "Which lipid classes are enriched in males vs females?"
- "What genes are correlated with triglycerides in the MOFA network?"
- "Compare Esr1 expression across all liver cell types"
- "Which TFs regulate Cyp2d9 in the regulatory network?"
- "Show me sex-specific zonation genes in testes pericentral hepatocytes"
- "What lipids are connected to PE(19:0_20:3) in the conversion network?"
- "Are testes or ovaries lipids more variable?"
- "What validation scores did Cyp2b9 knockdown achieve?"
- "How many genes with FDR <= 0.05 and |logFC| >= 0.4 in hepatocytes?"
- "Count significant lipids with logFC above 0.5"

## Data Files Used

The server reads CSV files from the `results/` directory:

| Workflow | File | Purpose |
|----------|------|---------|
| WF02 | `results/pseudobulk/tval_full_df.csv` | RNA t-statistics |
| WF02 | `results/pseudobulk/ttop_full_long.csv` | Full RNA DE results (logFC, P.Value, adj.P.Val) |
| WF03 | `results/lipidomic/gonade_type_ttop.csv` | Lipid DE results |
| WF03 | `results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv` | Lipid class enrichment |
| WF03 | `results/lipidomic/network/significant_lipids_network_SIF.csv` | LINEX2 network edges |
| WF03 | `results/lipidomic/network/significant_lipids_network_ATT.csv` | LINEX2 node attributes |
| WF04 | `results/mofa/networks/gene_lipid_correlations_factor1.csv` | Gene-lipid correlations |
| WF04 | `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv` | TF network edges |
| WF04 | `results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_ATT.csv` | TF network attributes |
| WF05 | `results/zonation/summary_statistics.csv` | Zonation summary |
| WF05 | `results/zonation/t_table_zonation.csv` | Gene zonation t-values |
| WF05 | `results/zonation/lipid_class_association_summary.csv` | Lipid class zonation |
| WF05 | `results/zonation/sex_specific_zonation_genes_summary.csv` | Sex-specific zonation genes |
| WF05 | `results/zonation/lipid_species_zonation_association.csv` | Per-lipid zonation |
| WF05 | `results/pseudobulk/rna_zonation/pathway_scores_df.csv` | Zonation pathway scores |
| WF06 | `results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv` | Structural associations |
| WF07 | `results/validation/validation_summary.csv` | Validation overview |
| WF07 | `results/validation/global_network_validation_scores.csv` | Global validation scores |
| WF07 | `results/validation/zonation_network_validation_scores.csv` | Zonation validation scores |
| WF08 | `results/lipidomic/variance_analysis/variance_analysis_summary.csv` | Variance overview |
| WF08 | `results/lipidomic/variance_analysis/variance_analysis_results.csv` | Per-lipid variance |
| WF08 | `results/lipidomic/variance_analysis/variance_analysis_by_class.csv` | Per-class variance |

## Environment Variables

- `FCGL_DATA_PATH`: Path to the repository root directory (optional, defaults to parent of script)

## Troubleshooting

### "Data file not found" errors
Ensure you have run the analysis workflows (01-08) to generate the required CSV files.

### Testing the server
Run the server directly to check for errors:
```bash
python mcp/fcgl_server.py
```

## Development

To extend the server with new tools:

1. Add a new data loader method in `FCGLDataLoader`
2. Add a new tool definition in `list_tools()`
3. Add the handler function (e.g., `handle_new_tool()`)
4. Add the routing in `call_tool()`

## License

Same as the main repository (GNU GPL v3.0).
