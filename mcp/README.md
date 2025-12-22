# FCGL Analysis MCP Server

A Model Context Protocol (MCP) server that provides tools for querying the Four Core Genotypes Liver (FCGL) multi-omic analysis results.

## Overview

This MCP server allows Claude (or other MCP-compatible AI assistants) to dynamically query:

- **Differential Expression Results**: RNA-seq and lipidomics DE results across contrasts and cell types
- **Pathway Enrichment**: GO Biological Process pathway scores
- **TF Activity**: Transcription factor activity inferred from target genes
- **Lipid Class Enrichment**: Enrichment of lipid classes by condition
- **MOFA Integration**: Factor weights and gene-lipid correlation networks

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

### 3. Configure for Claude Desktop

Add to your Claude Desktop configuration (`claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "fcgl-analysis": {
      "command": "python",
      "args": ["/path/to/gonads_snRNAseq_and_lipidomics_analysis_refactored/mcp/fcgl_server.py"],
      "env": {
        "FCGL_DATA_PATH": "/path/to/gonads_snRNAseq_and_lipidomics_analysis_refactored"
      }
    }
  }
}
```

## Available Tools

### 1. `get_dataset_summary`
Get an overview of the FCGL dataset including sample counts, cell types, contrasts analyzed, and key findings.

### 2. `query_gene_de`
Query differential expression results for a specific gene.

**Parameters:**
- `gene` (required): Gene symbol (e.g., "Cyp2b9", "Esr1")
- `cell_type`: Cell type to query (default: "Hepatocytes")
- `contrast`: Contrast to query (default: "testes_vs_ova")

**Example:**
```
Query the differential expression of Cyp2b9 in hepatocytes for the gonadal effect
```

### 3. `query_lipid_de`
Query differential abundance results for a lipid or lipid class.

**Parameters:**
- `lipid` (required): Lipid name or class (e.g., "PC", "TG", "Cer")
- `contrast`: Contrast to query
- `top_n`: Number of results to return

### 4. `get_top_de_genes`
Get the top differentially expressed genes for a contrast and cell type.

**Parameters:**
- `cell_type` (required): Cell type to query
- `contrast`: Contrast to query
- `direction`: "up", "down", or "both"
- `top_n`: Number of genes to return
- `fdr_threshold`: Significance threshold

### 5. `query_pathway_enrichment`
Query pathway enrichment scores for a specific pathway or keyword.

**Parameters:**
- `pathway` (required): Pathway name or keyword (e.g., "lipid", "hormone")
- `cell_type`: Filter by cell type
- `top_n`: Number of results

### 6. `query_tf_activity`
Query transcription factor activity scores.

**Parameters:**
- `tf` (required): TF name (e.g., "Esr1", "Ar", "Stat5a")
- `cell_type`: Filter by cell type

### 7. `query_lipid_class_enrichment`
Query lipid class enrichment results.

**Parameters:**
- `lipid_class`: Lipid class to query (e.g., "TG", "PC", "Cer")
- `contrast`: Filter by contrast

### 8. `query_mofa_gene_lipid_network`
Query gene-lipid correlations from MOFA Factor 1.

**Parameters:**
- `query` (required): Gene or lipid to search for
- `correlation_threshold`: Minimum |r| (default: 0.5)
- `top_n`: Number of correlations to return

### 9. `query_mofa_factor_weights`
Query MOFA factor weights for specific views and factors.

**Parameters:**
- `factor` (required): Factor number (1-12)
- `view`: View to query (e.g., "Hepatocytes", "Lipidomic")
- `top_n`: Number of features to return

### 10. `compare_gene_across_celltypes`
Compare a gene's differential expression across all cell types.

**Parameters:**
- `gene` (required): Gene symbol to compare
- `contrast`: Contrast to query

## Example Queries

Once connected to Claude, you can ask questions like:

- "What is the differential expression of Cyp2b9 in hepatocytes?"
- "Show me the top 20 upregulated genes in hepatocytes for the gonadal effect"
- "Which lipid classes are enriched in males vs females?"
- "What genes are correlated with triglycerides in the MOFA network?"
- "Compare Esr1 expression across all liver cell types"
- "What pathways are enriched for lipid metabolism?"
- "Show me the transcription factor activity for Stat5a"

## Data Files Used

The server reads from these RData files in the `results/` directory:

- `results/pseudobulk/ttop_full_list.RData` - RNA DE results
- `results/lipidomic/ttop_full_list.RData` - Lipid DE results
- `results/mofa/main_model_RNA_lipid_weight.RData` - MOFA weights
- `results/mofa/networks/gene_lipid_correlations_factor1.RData` - Gene-lipid correlations
- `results/lipidomic/enrichment/lipid_class_enrichment_scores.RData` - Lipid enrichment
- `results/pseudobulk/pathways/pathway_scores_df_all.RData` - Pathway enrichment
- `results/pseudobulk/pathways/tf_scores_df_all.RData` - TF enrichment

## Environment Variables

- `FCGL_DATA_PATH`: Path to the repository root directory (optional, defaults to parent of script)

## Troubleshooting

### "Data file not found" errors
Ensure you have run the analysis workflows (01-06) to generate the required RData files.

### pyreadr import errors
Install pyreadr: `pip install pyreadr`

On Windows, you may need to install the R runtime or use conda:
```bash
conda install -c conda-forge r-base pyreadr
```

### Testing the server
Run the server directly to check for errors:
```bash
python mcp/fcgl_server.py
```

## Development

To extend the server with new tools:

1. Add a new tool definition in `list_tools()`
2. Add the handler function (e.g., `handle_new_tool()`)
3. Add the routing in `call_tool()`
4. Add a new data loader method if needed in `FCGLDataLoader`

## License

Same as the main repository (GNU GPL v3.0).
