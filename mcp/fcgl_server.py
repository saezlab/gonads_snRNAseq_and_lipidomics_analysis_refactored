#!/usr/bin/env python3
"""
FCGL Analysis MCP Server

This MCP server provides tools for querying the Four Core Genotypes Liver (FCGL)
multi-omic analysis results using CSV files for reliable data access.

Usage:
    python fcgl_server.py
"""

import os
import sys
import json
import logging
from pathlib import Path
from typing import Optional
import asyncio

# MCP SDK imports
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

# Data handling
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, stream=sys.stderr)
logger = logging.getLogger("fcgl-mcp")


class FCGLDataLoader:
    """Handles loading and caching of FCGL analysis data files."""

    def __init__(self, data_path: str):
        self.data_path = Path(data_path)
        self.cache = {}
        logger.info(f"Data path: {self.data_path}")

    def _load_csv(self, rel_path: str) -> pd.DataFrame:
        """Load a CSV file with caching."""
        if rel_path in self.cache:
            return self.cache[rel_path]

        full_path = self.data_path / rel_path
        if not full_path.exists():
            raise FileNotFoundError(f"File not found: {full_path}")

        logger.info(f"Loading: {rel_path}")
        df = pd.read_csv(full_path)
        self.cache[rel_path] = df
        return df

    def get_rna_tvalues(self) -> pd.DataFrame:
        """Get RNA t-statistics (genes x contrasts_celltypes)."""
        return self._load_csv("results/pseudobulk/tval_full_df.csv")

    def get_lipid_de(self) -> pd.DataFrame:
        """Get lipid DE results for gonadal type comparison."""
        return self._load_csv("results/lipidomic/gonade_type_ttop.csv")

    def get_gene_lipid_correlations(self) -> pd.DataFrame:
        """Get gene-lipid correlation network from MOFA Factor 1."""
        return self._load_csv("results/mofa/networks/gene_lipid_correlations_factor1.csv")

    def get_gene_lipid_network(self) -> pd.DataFrame:
        """Get the full gene-lipid network edges."""
        return self._load_csv("results/mofa/networks/gene_lipid_network_edges.csv")

    def get_structural_associations(self) -> pd.DataFrame:
        """Get significant gene-lipid class structural associations."""
        return self._load_csv("results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv")

    def get_zonation_summary(self) -> pd.DataFrame:
        """Get zonation analysis summary."""
        return self._load_csv("results/zonation/summary_statistics.csv")

    def get_zonation_genes(self) -> pd.DataFrame:
        """Get gene-level zonation t-statistics."""
        return self._load_csv("results/zonation/t_table_zonation.csv")

    def get_lipid_zonation(self) -> pd.DataFrame:
        """Get lipid class zonation associations."""
        return self._load_csv("results/zonation/lipid_class_association_summary.csv")

    def get_validation_scores(self) -> pd.DataFrame:
        """Get validation experiment scores."""
        return self._load_csv("results/validation/validation_summary.csv")

    def get_lipid_class_enrichment(self) -> pd.DataFrame:
        """Get lipid class enrichment scores from gonadal comparison.

        First tries the new file with p-values, falls back to legacy file.
        """
        try:
            # Try new file with p-values first
            return self._load_csv("results/lipidomic/enrichment/lipid_class_enrichment_gonadal.csv")
        except FileNotFoundError:
            # Fall back to legacy file (average scores only, no p-values)
            df = self._load_csv("results/lipidomic/lipid_weights_classes.csv")
            # Rename columns to match new format
            df = df.rename(columns={
                'Lipid Class': 'lipid_class',
                'average_value': 'enrichment_score'
            })
            return df


# Initialize server
server = Server("fcgl-analysis")
data_loader: Optional[FCGLDataLoader] = None


def get_loader() -> FCGLDataLoader:
    global data_loader
    if data_loader is None:
        data_path = os.environ.get("FCGL_DATA_PATH", str(Path(__file__).parent.parent))
        data_loader = FCGLDataLoader(data_path)
    return data_loader


@server.list_tools()
async def list_tools():
    """List all available tools."""
    return [
        Tool(
            name="get_dataset_summary",
            description="Get an overview of the FCGL dataset: sample counts, cell types, contrasts, and key findings about gonadal hormone effects on liver.",
            inputSchema={"type": "object", "properties": {}, "required": []}
        ),
        Tool(
            name="query_gene",
            description="Query differential expression statistics for a gene across all cell types and contrasts. Returns t-values showing effect direction and strength.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (e.g., 'Cyp2b9', 'Esr1', 'Ar', 'Cyp2d9')"
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="query_lipid",
            description="Query differential abundance for lipids matching a pattern. Returns logFC, p-values for gonadal effect (Testes vs Ovaries).",
            inputSchema={
                "type": "object",
                "properties": {
                    "lipid": {
                        "type": "string",
                        "description": "Lipid name or class to search (e.g., 'PC', 'TG', 'Cer', 'sphingo')"
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of results to return (default: 20)",
                        "default": 20
                    }
                },
                "required": ["lipid"]
            }
        ),
        Tool(
            name="get_top_genes",
            description="Get top differentially expressed genes for a specific cell type and contrast.",
            inputSchema={
                "type": "object",
                "properties": {
                    "cell_type": {
                        "type": "string",
                        "description": "Cell type",
                        "enum": ["Hepatocytes", "Endothelial_cells", "Immune_cells", "Stellate_cells", "Cholangiocytes"]
                    },
                    "contrast": {
                        "type": "string",
                        "description": "Contrast (default: testes_vs_ova)",
                        "enum": ["testes_vs_ova", "XY_vs_XX", "XXT_vs_XXO", "XYT_vs_XYO", "XYO_vs_XXO", "XYT_vs_XXT"]
                    },
                    "direction": {
                        "type": "string",
                        "description": "Filter: 'up' (positive), 'down' (negative), or 'both'",
                        "enum": ["up", "down", "both"]
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of genes (default: 20)"
                    }
                },
                "required": ["cell_type"]
            }
        ),
        Tool(
            name="query_gene_lipid_network",
            description="Find genes correlated with lipids (or lipids correlated with genes) from MOFA Factor 1 integration.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Gene name or lipid pattern to search for correlations"
                    },
                    "min_correlation": {
                        "type": "number",
                        "description": "Minimum absolute correlation (default: 0.5)"
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of results (default: 20)"
                    }
                },
                "required": ["query"]
            }
        ),
        Tool(
            name="get_structural_associations",
            description="Get genes with significant lipid structural preferences (chain length, saturation effects).",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Filter by gene name (optional)"
                    },
                    "lipid_class": {
                        "type": "string",
                        "description": "Filter by lipid class (optional)"
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="compare_gene_celltypes",
            description="Compare a gene's effect across all 5 liver cell types for the gonadal comparison.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol to compare"
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="get_zonation_info",
            description="Get hepatocyte zonation analysis results (periportal vs pericentral differences).",
            inputSchema={"type": "object", "properties": {}, "required": []}
        ),
        Tool(
            name="query_gene_zonation",
            description="Check if a gene is periportal or pericentral biased in hepatocytes. Positive portal_vs_Central = periportal-biased.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (e.g., 'Cyp2b9', 'Cyp2c38')"
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="query_lipid_zonation",
            description="Get lipid classes associated with periportal or pericentral zonation, including male/female bias.",
            inputSchema={
                "type": "object",
                "properties": {
                    "zonation": {
                        "type": "string",
                        "description": "Filter by zonation: 'periportal', 'pericentral', or 'all'",
                        "enum": ["periportal", "pericentral", "all"]
                    },
                    "gonad_bias": {
                        "type": "string",
                        "description": "Filter by gonad bias: 'Testes' (male), 'Ovaries' (female), or 'all'",
                        "enum": ["Testes", "Ovaries", "all"]
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="get_validation_results",
            description="Get siRNA knockdown validation results comparing predictions to experiments.",
            inputSchema={"type": "object", "properties": {}, "required": []}
        ),
        Tool(
            name="query_lipid_class_enrichment",
            description="Get lipid class enrichment scores for the gonadal comparison (Testes vs Ovaries). Shows which lipid classes are male-biased (positive) or female-biased (negative). Lipid classes include: TG (Triglycerides), DG (Diglycerides), PI (Phosphatidylinositols), PE (Phosphatidylethanolamines), LPE (Lysophosphatidylethanolamines), PC (Phosphatidylcholines), SM (Sphingomyelins), Cer (Ceramides), HexCer (Hexosylceramides), and more.",
            inputSchema={
                "type": "object",
                "properties": {
                    "lipid_classes": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of lipid class abbreviations to query (e.g., ['TG', 'DG', 'SM']). Leave empty for all classes."
                    }
                },
                "required": []
            }
        )
    ]


@server.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    """Handle tool calls."""
    try:
        loader = get_loader()

        if name == "get_dataset_summary":
            return handle_summary()
        elif name == "query_gene":
            return handle_query_gene(loader, arguments)
        elif name == "query_lipid":
            return handle_query_lipid(loader, arguments)
        elif name == "get_top_genes":
            return handle_top_genes(loader, arguments)
        elif name == "query_gene_lipid_network":
            return handle_network(loader, arguments)
        elif name == "get_structural_associations":
            return handle_structural(loader, arguments)
        elif name == "compare_gene_celltypes":
            return handle_compare_celltypes(loader, arguments)
        elif name == "get_zonation_info":
            return handle_zonation(loader)
        elif name == "query_gene_zonation":
            return handle_gene_zonation(loader, arguments)
        elif name == "query_lipid_zonation":
            return handle_lipid_zonation(loader, arguments)
        elif name == "get_validation_results":
            return handle_validation(loader)
        elif name == "query_lipid_class_enrichment":
            return handle_lipid_class_enrichment(loader, arguments)
        else:
            return [TextContent(type="text", text=f"Unknown tool: {name}")]

    except Exception as e:
        logger.exception(f"Error in {name}")
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_summary() -> list[TextContent]:
    text = """# FCGL Dataset Summary

## Overview
Multi-omic study of sex-specific liver physiology using the Four Core Genotypes (FCG) mouse model.

## Key Statistics
- **snRNA-seq**: 466,986 nuclei from 48 liver samples (43 FCG + 5 WT males)
- **Lipidomics**: 1,386 lipid species across 42 classes
- **Cell types**: Hepatocytes, Endothelial cells, Immune cells, Stellate cells, Cholangiocytes

## FCG Groups
| Group | Chromosomes | Gonads |
|-------|-------------|--------|
| XXO | XX | Ovaries |
| XXT | XX | Testes |
| XYO | XY | Ovaries |
| XYT | XY | Testes |

## Key Findings
- **~4,000 gonad-dependent genes** in hepatocytes (FDR ≤ 0.05)
- **733/1,386 lipids (53%)** are gonad-dependent
- **NO lipids** show sex chromosome effects
- **Male-biased**: TG, DG, glycerophospholipids (PI, PE, PC)
- **Female-biased**: Sphingolipids (SM, Cer, HexCer)
- **Key genes**: Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9, Mtmr11
- **Key TFs**: Esr1, Ar, Stat5a/b, Nfic, Dbp

## Available Queries
- `query_gene`: DE stats for any gene
- `query_lipid`: DE stats for lipids
- `query_gene_lipid_network`: MOFA correlations
- `get_structural_associations`: Lipid structural preferences
- `compare_gene_celltypes`: Cross-cell type comparison
"""
    return [TextContent(type="text", text=text)]


def handle_query_gene(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    gene = args.get("gene", "").strip()
    if not gene:
        return [TextContent(type="text", text="Error: gene name required")]

    try:
        df = loader.get_rna_tvalues()

        # Find gene column (ID or gene)
        gene_col = 'ID' if 'ID' in df.columns else ('gene' if 'gene' in df.columns else df.columns[0])

        # Search case-insensitive
        mask = df[gene_col].str.lower() == gene.lower()
        result = df[mask]

        if len(result) == 0:
            # Try partial match
            mask = df[gene_col].str.lower().str.contains(gene.lower(), na=False)
            result = df[mask].head(5)
            if len(result) == 0:
                return [TextContent(type="text", text=f"Gene '{gene}' not found")]
            text = f"No exact match for '{gene}'. Partial matches:\n\n"
            text += result.to_string()
            return [TextContent(type="text", text=text)]

        row = result.iloc[0]

        text = f"# Differential Expression: {gene}\n\n"
        text += "## T-statistics by Contrast and Cell Type\n"
        text += "(Positive = higher in first condition, e.g., Testes > Ovaries)\n\n"
        text += "| Contrast | Cell Type | t-value |\n"
        text += "|----------|-----------|--------|\n"

        for col in df.columns:
            if col != gene_col and '_tval_' in col:
                val = row[col]
                if pd.notna(val):
                    # Parse column name: "Testes vs Ovaries_tval_Hepatocytes"
                    parts = col.split('_tval_')
                    contrast = parts[0] if len(parts) > 0 else col
                    cell_type = parts[1] if len(parts) > 1 else ""
                    text += f"| {contrast} | {cell_type} | {val:.3f} |\n"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_query_lipid(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    lipid = args.get("lipid", "").strip()
    top_n = args.get("top_n", 20)

    if not lipid:
        return [TextContent(type="text", text="Error: lipid name/class required")]

    try:
        df = loader.get_lipid_de()

        # Find lipid column (ID or lipid)
        lipid_col = 'ID' if 'ID' in df.columns else ('lipid' if 'lipid' in df.columns else df.columns[0])

        # Search
        mask = df[lipid_col].str.contains(lipid, case=False, na=False)
        result = df[mask].copy()

        if len(result) == 0:
            return [TextContent(type="text", text=f"No lipids matching '{lipid}'")]

        # Sort by absolute logFC if available
        if 'logFC' in result.columns:
            result = result.sort_values('logFC', key=abs, ascending=False)

        result = result.head(top_n)

        text = f"# Lipid DE Results: {lipid}\n\n"
        text += f"**Contrast**: Testes vs Ovaries (gonadal effect)\n"
        text += f"Found {len(result)} matches\n\n"
        text += result.to_string(index=False)

        # Summary
        if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
            all_matches = df[df[lipid_col].str.contains(lipid, case=False, na=False)]
            sig_up = len(all_matches[(all_matches['adj.P.Val'] < 0.05) & (all_matches['logFC'] > 0)])
            sig_down = len(all_matches[(all_matches['adj.P.Val'] < 0.05) & (all_matches['logFC'] < 0)])
            text += f"\n\n**Summary**: {sig_up} male-biased (logFC > 0), {sig_down} female-biased (logFC < 0) at FDR < 0.05"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_top_genes(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    cell_type = args.get("cell_type", "Hepatocytes")
    contrast = args.get("contrast", "testes_vs_ova")
    direction = args.get("direction", "both")
    top_n = args.get("top_n", 20)

    try:
        df = loader.get_rna_tvalues()
        gene_col = 'ID' if 'ID' in df.columns else ('gene' if 'gene' in df.columns else df.columns[0])

        # Map user-friendly contrast names to actual column patterns
        contrast_map = {
            "testes_vs_ova": "Testes vs Ovaries",
            "XY_vs_XX": "XY vs XX",
            "XXT_vs_XXO": "XXT vs XXO",
            "XYT_vs_XYO": "XYT vs XYO",
            "XYO_vs_XXO": "XYO vs XXO",
            "XYT_vs_XXT": "XYT vs XXT"
        }
        contrast_pattern = contrast_map.get(contrast, contrast)

        # Map cell type names (handle underscores vs spaces)
        cell_type_clean = cell_type.replace("_", " ")

        # Find matching column
        target_col = None
        for col in df.columns:
            if contrast_pattern in col and cell_type_clean in col:
                target_col = col
                break

        if target_col is None:
            # Try with original cell_type
            for col in df.columns:
                if contrast_pattern in col and cell_type in col:
                    target_col = col
                    break

        if target_col is None:
            available = [c for c in df.columns if '_tval_' in c or 'tval' in c.lower()][:5]
            return [TextContent(type="text",
                text=f"Column not found for '{cell_type_clean}' + '{contrast_pattern}'.\nAvailable columns: {available}")]

        # Filter and sort
        result = df[[gene_col, target_col]].dropna()

        if direction == "up":
            result = result[result[target_col] > 0]
        elif direction == "down":
            result = result[result[target_col] < 0]

        result = result.sort_values(target_col, key=abs, ascending=False).head(top_n)

        # Better formatting
        dir_label = "male-biased (t > 0)" if direction == "up" else "female-biased (t < 0)" if direction == "down" else "all"
        text = f"# Top DE Genes: {cell_type_clean}\n\n"
        text += f"**Contrast**: {contrast_pattern}\n"
        text += f"**Direction**: {dir_label}\n"
        text += f"**Top {top_n}**\n\n"
        text += "| Gene | t-value |\n"
        text += "|------|--------|\n"
        for _, row in result.iterrows():
            text += f"| {row[gene_col]} | {row[target_col]:.3f} |\n"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_network(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    query = args.get("query", "").strip()
    min_cor = args.get("min_correlation", 0.5)
    top_n = args.get("top_n", 20)

    if not query:
        return [TextContent(type="text", text="Error: query required")]

    try:
        df = loader.get_gene_lipid_correlations()

        # The file is in wide format: genes as rows (Unnamed: 0), lipids as columns
        # First column contains gene names
        gene_col = df.columns[0]  # Usually "Unnamed: 0"
        lipid_cols = df.columns[1:]  # All other columns are lipids

        # Check if query matches a gene
        gene_match = df[df[gene_col].str.lower() == query.lower()]

        if len(gene_match) > 0:
            # Query is a gene - find correlated lipids
            row = gene_match.iloc[0]
            gene_name = row[gene_col]

            # Get correlations with all lipids
            correlations = []
            for lipid in lipid_cols:
                cor_val = row[lipid]
                if pd.notna(cor_val) and abs(cor_val) >= min_cor:
                    correlations.append({
                        'lipid': lipid,
                        'correlation': cor_val
                    })

            if not correlations:
                return [TextContent(type="text",
                    text=f"No lipid correlations >= {min_cor} found for gene '{gene_name}'")]

            # Sort by absolute correlation
            correlations.sort(key=lambda x: abs(x['correlation']), reverse=True)
            correlations = correlations[:top_n]

            text = f"# Lipids Correlated with {gene_name}\n\n"
            text += f"**Min |correlation|**: {min_cor}\n"
            text += f"**Found**: {len(correlations)} lipids\n\n"
            text += "| Lipid | Correlation | Direction |\n"
            text += "|-------|-------------|----------|\n"
            for c in correlations:
                direction = "Positive" if c['correlation'] > 0 else "Negative"
                text += f"| {c['lipid']} | {c['correlation']:.3f} | {direction} |\n"
            text += "\n*From MOFA Factor 1 (strongest gonadal association)*"

            return [TextContent(type="text", text=text)]

        # Check if query matches a lipid (column name)
        matching_lipids = [col for col in lipid_cols if query.lower() in col.lower()]

        if matching_lipids:
            # Query is a lipid - find correlated genes
            all_correlations = []
            for lipid in matching_lipids:
                for idx, row in df.iterrows():
                    cor_val = row[lipid]
                    if pd.notna(cor_val) and abs(cor_val) >= min_cor:
                        all_correlations.append({
                            'gene': row[gene_col],
                            'lipid': lipid,
                            'correlation': cor_val
                        })

            if not all_correlations:
                return [TextContent(type="text",
                    text=f"No gene correlations >= {min_cor} found for lipid pattern '{query}'")]

            # Sort by absolute correlation
            all_correlations.sort(key=lambda x: abs(x['correlation']), reverse=True)
            all_correlations = all_correlations[:top_n]

            text = f"# Genes Correlated with Lipids Matching '{query}'\n\n"
            text += f"**Min |correlation|**: {min_cor}\n"
            text += f"**Found**: {len(all_correlations)} associations\n\n"
            text += "| Gene | Lipid | Correlation |\n"
            text += "|------|-------|-------------|\n"
            for c in all_correlations:
                text += f"| {c['gene']} | {c['lipid']} | {c['correlation']:.3f} |\n"
            text += "\n*From MOFA Factor 1 (strongest gonadal association)*"

            return [TextContent(type="text", text=text)]

        return [TextContent(type="text",
            text=f"No gene or lipid matching '{query}' found in the correlation network")]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_structural(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    gene = args.get("gene", "")
    lipid_class = args.get("lipid_class", "")

    try:
        df = loader.get_structural_associations()

        if gene:
            df = df[df['gene'].str.contains(gene, case=False, na=False)]
        if lipid_class:
            df = df[df['lipid_class'].str.contains(lipid_class, case=False, na=False)]

        if len(df) == 0:
            return [TextContent(type="text", text="No structural associations found")]

        # Add interpretation column
        def interpret_term(row):
            term = row['term']
            est = row['estimate']
            if 'chain_length' in term:
                if est > 0:
                    return "longer chains"
                else:
                    return "shorter chains"
            elif 'db_position' in term:
                if est > 0:
                    return "more unsaturated"
                else:
                    return "more saturated"
            elif term == 'logfc':
                if est > 0:
                    return "higher abundance"
                else:
                    return "lower abundance"
            return ""

        df = df.copy()
        df['interpretation'] = df.apply(interpret_term, axis=1)

        text = "# Lipid Structural Associations\n\n"
        text += "Genes whose correlation with lipids depends on structural features.\n\n"
        text += "**Interpretation guide:**\n"
        text += "- chain_length: positive estimate = stronger correlation with longer chain lipids\n"
        text += "- db_position: positive estimate = stronger correlation with more unsaturated lipids\n"
        text += "- logfc:X (interaction): effect of X depends on lipid abundance\n\n"

        # Group by gene if searching for a specific gene
        if gene and len(df) > 0:
            gene_name = df['gene'].iloc[0]
            text += f"## {gene_name} Structural Preferences\n\n"
            text += "| Lipid Class | Term | Estimate | p-value | Preference |\n"
            text += "|-------------|------|----------|---------|------------|\n"
            for _, row in df.head(20).iterrows():
                lc = row['lipid_class'].replace('_with_name', '')
                text += f"| {lc} | {row['term']} | {row['estimate']:.4f} | {row['p_value']:.4f} | {row['interpretation']} |\n"

            # Summary
            chain_rows = df[df['term'].str.contains('chain_length', na=False) & ~df['term'].str.contains(':', na=False)]
            if len(chain_rows) > 0:
                avg_chain = chain_rows['estimate'].mean()
                text += f"\n**Chain length summary**: Average estimate = {avg_chain:.4f} "
                text += f"({'prefers longer chains' if avg_chain > 0 else 'prefers shorter chains'})\n"

            db_rows = df[df['term'].str.contains('db_position', na=False) & ~df['term'].str.contains(':', na=False)]
            if len(db_rows) > 0:
                avg_db = db_rows['estimate'].mean()
                text += f"**Saturation summary**: Average estimate = {avg_db:.4f} "
                text += f"({'prefers unsaturated' if avg_db > 0 else 'prefers saturated'})\n"
        else:
            text += df.head(30).to_string(index=False)

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_compare_celltypes(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    gene = args.get("gene", "").strip()

    if not gene:
        return [TextContent(type="text", text="Error: gene required")]

    try:
        df = loader.get_rna_tvalues()
        gene_col = 'gene' if 'gene' in df.columns else df.columns[0]

        mask = df[gene_col].str.lower() == gene.lower()
        result = df[mask]

        if len(result) == 0:
            return [TextContent(type="text", text=f"Gene '{gene}' not found")]

        row = result.iloc[0]

        text = f"# {gene} Across Cell Types\n\n"
        text += "## Gonadal Effect (Testes vs Ovaries)\n\n"
        text += "| Cell Type | t-value | Interpretation |\n"
        text += "|-----------|---------|----------------|\n"

        cell_types = ["Hepatocytes", "Endothelial_cells", "Immune_cells", "Stellate_cells", "Cholangiocytes"]

        for ct in cell_types:
            col = f"testes_vs_ova_{ct}"
            if col in row.index:
                val = row[col]
                if pd.notna(val):
                    interp = "Male-biased" if val > 0 else "Female-biased"
                    if abs(val) < 2:
                        interp = "Weak effect"
                    text += f"| {ct.replace('_', ' ')} | {val:.3f} | {interp} |\n"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_zonation(loader: FCGLDataLoader) -> list[TextContent]:
    try:
        df = loader.get_zonation_summary()
        text = "# Hepatocyte Zonation Analysis Summary\n\n"
        text += "Periportal (Cyp2f2+) vs Pericentral (Cyp2e1+) hepatocytes\n\n"
        text += "## DE Statistics by Comparison\n\n"
        text += df.to_string(index=False)

        # Also show lipid zonation summary
        try:
            lipid_df = loader.get_lipid_zonation()
            text += "\n\n## Lipid Class Zonation Associations\n\n"
            text += "TPC=Testes Pericentral, TPP=Testes Periportal, OPC=Ovaries Pericentral, OPP=Ovaries Periportal\n\n"
            text += lipid_df.to_string(index=False)
        except:
            pass

        return [TextContent(type="text", text=text)]
    except Exception as e:
        return [TextContent(type="text", text=f"Error loading zonation data: {str(e)}")]


def handle_gene_zonation(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    gene = args.get("gene", "").strip()

    if not gene:
        return [TextContent(type="text", text="Error: gene name required")]

    try:
        df = loader.get_zonation_genes()

        # Search case-insensitive
        mask = df['gene'].str.lower() == gene.lower()
        result = df[mask]

        if len(result) == 0:
            # Try partial match
            mask = df['gene'].str.lower().str.contains(gene.lower(), na=False)
            result = df[mask].head(5)
            if len(result) == 0:
                return [TextContent(type="text", text=f"Gene '{gene}' not found in zonation data")]
            text = f"No exact match for '{gene}'. Partial matches:\n\n"
            text += result.to_string(index=False)
            return [TextContent(type="text", text=text)]

        row = result.iloc[0]
        gene_name = row['gene']

        # Interpret zonation bias
        portal_t = row['portal_vs_Central']
        zonation_bias = "periportal-biased" if portal_t > 0 else "pericentral-biased"
        zonation_strength = "strongly " if abs(portal_t) > 5 else "moderately " if abs(portal_t) > 2 else "weakly "

        # Interpret gonad bias
        gonad_t = row['testes_vs_ova']
        gonad_bias = "male-biased (testes)" if gonad_t > 0 else "female-biased (ovaries)"

        text = f"# Zonation Analysis: {gene_name}\n\n"
        text += f"## Summary\n"
        text += f"- **Zonation**: {zonation_strength}{zonation_bias} (t = {portal_t:.2f})\n"
        text += f"- **Gonad effect**: {gonad_bias} (t = {gonad_t:.2f})\n\n"

        text += "## All Contrasts\n\n"
        text += "| Contrast | t-value | Interpretation |\n"
        text += "|----------|---------|----------------|\n"

        contrasts = {
            'portal_vs_Central': 'Portal vs Central (overall zonation)',
            'testes_vs_ova': 'Testes vs Ovaries (overall gonad effect)',
            'portalova_vs_centralova': 'Portal vs Central in Ovaries',
            'portaltestes_vs_centraltestes': 'Portal vs Central in Testes',
            'centraltestes_vs_centralova': 'Testes vs Ovaries in Central',
            'portaltestes_vs_portalova': 'Testes vs Ovaries in Portal',
            'testes_effect_zonation': 'Testes effect on zonation (interaction)'
        }

        for col, desc in contrasts.items():
            if col in row.index:
                val = row[col]
                if pd.notna(val):
                    direction = "UP" if val > 0 else "DOWN"
                    text += f"| {desc} | {val:.3f} | {direction} |\n"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_lipid_zonation(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    zonation = args.get("zonation", "all")
    gonad_bias = args.get("gonad_bias", "all")

    try:
        df = loader.get_lipid_zonation()

        # Filter by gonad bias
        if gonad_bias != "all":
            df = df[df['gonad_bias'] == gonad_bias]

        # Filter by zonation (TPC, TPP, OPC, OPP columns indicate association)
        # TPP/OPP = periportal, TPC/OPC = pericentral
        if zonation == "periportal":
            df = df[(df['TPP'] == 1) | (df['OPP'] == 1)]
        elif zonation == "pericentral":
            df = df[(df['TPC'] == 1) | (df['OPC'] == 1)]

        if len(df) == 0:
            return [TextContent(type="text", text="No lipid classes match the criteria")]

        text = f"# Lipid Class Zonation Associations\n\n"
        text += f"**Filters**: zonation={zonation}, gonad_bias={gonad_bias}\n\n"
        text += "Columns: TPC=Testes-Pericentral, TPP=Testes-Periportal, OPC=Ovaries-Pericentral, OPP=Ovaries-Periportal\n"
        text += "(1 = associated with that condition)\n\n"
        text += df.to_string(index=False)

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_validation(loader: FCGLDataLoader) -> list[TextContent]:
    try:
        df = loader.get_validation_scores()
        text = "# Validation Results\n\n"
        text += "siRNA knockdown experiments validating gene-lipid predictions:\n\n"
        text += df.to_string(index=False)
        return [TextContent(type="text", text=text)]
    except Exception as e:
        return [TextContent(type="text", text=f"Error loading validation data: {str(e)}")]


def handle_lipid_class_enrichment(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Query lipid class enrichment scores from gonadal comparison."""
    lipid_classes = args.get("lipid_classes", [])

    try:
        df = loader.get_lipid_class_enrichment()

        # Check if we have the new format with p-values
        has_pvalues = 'p_value' in df.columns

        # Rename columns for clarity (handle both old and new formats)
        if 'Lipid Class' in df.columns:
            df = df.rename(columns={
                'Lipid Class': 'lipid_class',
                'average_value': 'enrichment_score'
            })

        # Map abbreviations to full names
        lipid_class_names = {
            'Cer': 'Ceramides',
            'DG': 'Diglycerides',
            'FA': 'Fatty Acids',
            'HexCer': 'Hexosylceramides',
            'LPE': 'Lysophosphatidylethanolamines',
            'LPC': 'Lysophosphatidylcholines',
            'PA': 'Phosphatidic Acids',
            'PC': 'Phosphatidylcholines',
            'PE': 'Phosphatidylethanolamines',
            'PG': 'Phosphatidylglycerols',
            'PI': 'Phosphatidylinositols',
            'PS': 'Phosphatidylserines',
            'SHexCer': 'Sulfatides',
            'SM': 'Sphingomyelins',
            'ST': 'Sterols',
            'TG': 'Triglycerides'
        }

        # Add full names column
        df['full_name'] = df['lipid_class'].map(lipid_class_names).fillna(df['lipid_class'])

        # Filter by requested classes if specified
        if lipid_classes and len(lipid_classes) > 0:
            # Case-insensitive matching
            classes_upper = [c.upper() for c in lipid_classes]
            df = df[df['lipid_class'].str.upper().isin(classes_upper)]

            if len(df) == 0:
                available = loader.get_lipid_class_enrichment()['Lipid Class'].tolist()
                return [TextContent(type="text",
                    text=f"No matching lipid classes found. Available classes: {', '.join(available)}")]

        # Sort by enrichment score (most male-biased first)
        df = df.sort_values('enrichment_score', ascending=False)

        # Add interpretation
        def interpret_score(score):
            if score > 1.5:
                return "Strongly male-biased (Testes)"
            elif score > 0.5:
                return "Moderately male-biased (Testes)"
            elif score > 0:
                return "Weakly male-biased (Testes)"
            elif score > -0.5:
                return "Weakly female-biased (Ovaries)"
            elif score > -1.5:
                return "Moderately female-biased (Ovaries)"
            else:
                return "Strongly female-biased (Ovaries)"

        df['sex_bias'] = df['enrichment_score'].apply(interpret_score)

        text = "# Lipid Class Enrichment Analysis\n\n"
        text += "## Gonadal Comparison: Testes vs Ovaries\n\n"
        text += "**Interpretation:**\n"
        text += "- Positive enrichment score = Male-biased (higher in Testes)\n"
        text += "- Negative enrichment score = Female-biased (higher in Ovaries)\n\n"

        if has_pvalues:
            text += "| Lipid Class | Full Name | Enrichment Score | p-value | Sex Bias |\n"
            text += "|-------------|-----------|-----------------|---------|----------|\n"

            for _, row in df.iterrows():
                pval = row['p_value']
                pval_str = f"{pval:.2e}" if pval < 0.001 else f"{pval:.4f}"
                sig_marker = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
                text += f"| {row['lipid_class']} | {row['full_name']} | {row['enrichment_score']:.3f} | {pval_str}{sig_marker} | {row['sex_bias']} |\n"
        else:
            text += "| Lipid Class | Full Name | Enrichment Score | Sex Bias |\n"
            text += "|-------------|-----------|-----------------|----------|\n"

            for _, row in df.iterrows():
                text += f"| {row['lipid_class']} | {row['full_name']} | {row['enrichment_score']:.3f} | {row['sex_bias']} |\n"

        text += "\n## Summary\n\n"

        male_biased = df[df['enrichment_score'] > 0]
        female_biased = df[df['enrichment_score'] < 0]

        if len(male_biased) > 0:
            text += "**Male-biased lipid classes (Testes):**\n"
            for _, row in male_biased.iterrows():
                if has_pvalues:
                    pval = row['p_value']
                    pval_str = f"p = {pval:.2e}" if pval < 0.001 else f"p = {pval:.4f}"
                    text += f"- {row['full_name']} ({row['lipid_class']}): ES = {row['enrichment_score']:.3f}, {pval_str}\n"
                else:
                    text += f"- {row['full_name']} ({row['lipid_class']}): ES = {row['enrichment_score']:.3f}\n"

        if len(female_biased) > 0:
            text += "\n**Female-biased lipid classes (Ovaries):**\n"
            for _, row in female_biased.iterrows():
                if has_pvalues:
                    pval = row['p_value']
                    pval_str = f"p = {pval:.2e}" if pval < 0.001 else f"p = {pval:.4f}"
                    text += f"- {row['full_name']} ({row['lipid_class']}): ES = {row['enrichment_score']:.3f}, {pval_str}\n"
                else:
                    text += f"- {row['full_name']} ({row['lipid_class']}): ES = {row['enrichment_score']:.3f}\n"

        text += "\n*Enrichment scores from decoupleR run_ulm() on lipid t-statistics (workflow 03, step 5)*"
        if has_pvalues:
            text += "\n*Significance: *** p<0.001, ** p<0.01, * p<0.05*"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


async def main():
    """Run the MCP server."""
    logger.info("Starting FCGL Analysis MCP Server")

    global data_loader
    data_path = os.environ.get("FCGL_DATA_PATH", str(Path(__file__).parent.parent))
    data_loader = FCGLDataLoader(data_path)

    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream, server.create_initialization_options())


if __name__ == "__main__":
    asyncio.run(main())
