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
import re
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

    # --- WF03: Lipid network ---
    def get_lipid_network_sif(self) -> pd.DataFrame:
        """Get LINEX2 lipid conversion network edges."""
        return self._load_csv("results/lipidomic/network/significant_lipids_network_SIF.csv")

    def get_lipid_network_att(self) -> pd.DataFrame:
        """Get LINEX2 lipid network node attributes."""
        return self._load_csv("results/lipidomic/network/significant_lipids_network_ATT.csv")

    # --- WF04: TF-Gene-Lipid network ---
    def get_tf_network_sif(self) -> pd.DataFrame:
        """Get TF-Gene-Lipid regulatory network edges."""
        return self._load_csv("results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_SIF.csv")

    def get_tf_network_att(self) -> pd.DataFrame:
        """Get TF-Gene-Lipid regulatory network node attributes."""
        return self._load_csv("results/mofa/networks/hepatocyte_TF_gene_lipid_network_compressed_ATT.csv")

    # --- WF05: Zonation ---
    def get_zonation_ttop_long(self) -> pd.DataFrame:
        """Get full zonation DE results (long format: gene, contrast, logFC, t, P.Value, adj.P.Val)."""
        return self._load_csv("results/zonation/ttop_full_long.csv")

    def get_sex_specific_zonation(self) -> pd.DataFrame:
        """Get sex-specific zonation genes (137 genes)."""
        return self._load_csv("results/zonation/sex_specific_zonation_genes_summary.csv")

    def get_lipid_species_zonation(self) -> pd.DataFrame:
        """Get per-lipid species zonation associations."""
        return self._load_csv("results/zonation/lipid_species_zonation_association.csv")

    def get_zonation_pathways(self) -> pd.DataFrame:
        """Get zonation pathway enrichment scores."""
        return self._load_csv("results/pseudobulk/rna_zonation/pathway_scores_df.csv")

    # --- WF06: Structural analysis ---
    def get_structural_summary(self) -> pd.DataFrame:
        """Get structural analysis summary metrics."""
        return self._load_csv("results/multi-omic/structural_model/analysis_summary.csv")

    # --- WF07: Validation ---
    def get_global_validation(self) -> pd.DataFrame:
        """Get global network validation scores."""
        return self._load_csv("results/validation/global_network_validation_scores.csv")

    def get_zonation_validation(self) -> pd.DataFrame:
        """Get zonation network validation scores."""
        return self._load_csv("results/validation/zonation_network_validation_scores.csv")

    # --- WF08: Variance analysis ---
    def get_variance_results(self) -> pd.DataFrame:
        """Get per-lipid variance analysis results."""
        return self._load_csv("results/lipidomic/variance_analysis/variance_analysis_results.csv")

    def get_variance_by_class(self) -> pd.DataFrame:
        """Get per-class variance analysis summary."""
        return self._load_csv("results/lipidomic/variance_analysis/variance_analysis_by_class.csv")

    def get_variance_summary(self) -> pd.DataFrame:
        """Get variance analysis overview metrics."""
        return self._load_csv("results/lipidomic/variance_analysis/variance_analysis_summary.csv")

    # --- GDX validation ---
    def get_gdx_lipid_tvalues(self) -> pd.DataFrame:
        """Get GDX lipid t-values for correlation with FCG."""
        return self._load_csv("results/lipidomic/gdxhr_t_table.csv")

    def get_gdx_rna_tvalues(self) -> pd.DataFrame:
        """Get GDX RNA t-values for correlation with FCG."""
        return self._load_csv("results/pseudobulk/gdxrh/gdhxr_groups_tval_df.csv")

    # --- Full ttop for counting ---
    def get_rna_ttop_full(self) -> pd.DataFrame:
        """Get full RNA DE results (gene, contrast, cell_type, logFC, P.Value, adj.P.Val, B)."""
        return self._load_csv("results/pseudobulk/ttop_full_long.csv")


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
        ),
        Tool(
            name="get_lipid_de_summary",
            description="Get a comprehensive summary of gonad-dependent lipid differential expression: total counts, male/female split, per-class breakdown with counts and significance, and GDX validation correlation. Essential for writing the lipidomics results section.",
            inputSchema={"type": "object", "properties": {}, "required": []}
        ),
        Tool(
            name="query_lipid_network",
            description="Query the LINEX2 lipid-lipid conversion network for significant gonad-dependent lipids. Shows enzymatic conversion relationships between lipids, with fold changes and significance. Useful for understanding how male/female-biased lipids are metabolically connected.",
            inputSchema={
                "type": "object",
                "properties": {
                    "lipid": {
                        "type": "string",
                        "description": "Lipid name or class to search (e.g., 'PE(19:0_20:3)', 'TG', 'Cer')"
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
            name="query_tf_gene_lipid_network",
            description="Query the TF->Gene->Lipid regulatory network from MOFA integration. Shows how transcription factors (Nfic, Dbp, Esr1) regulate lipid metabolism genes (Cyp2b9, Cyp2d9, etc.) which in turn correlate with lipid classes. Central to the paper's regulatory cascade findings.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Gene, TF, or lipid class name to search (e.g., 'Cyp2d9', 'Nfic', 'Cer')"
                    },
                    "node_type": {
                        "type": "string",
                        "description": "Filter results by node type (optional)",
                        "enum": ["TF", "enzyme", "lipid_class", "all"]
                    }
                },
                "required": ["query"]
            }
        ),
        Tool(
            name="query_sex_specific_zonation_genes",
            description="Query the 137 genes that show zonation bias (periportal vs pericentral) only in one gonadal context. These genes have significant zonation in testes OR ovaries but not both, representing sex-specific spatial regulation.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gonad": {
                        "type": "string",
                        "description": "Filter by gonadal context: 'Testes' or 'Ovaries'",
                        "enum": ["Testes", "Ovaries"]
                    },
                    "zone": {
                        "type": "string",
                        "description": "Filter by zonation: 'pericentral' or 'periportal'",
                        "enum": ["pericentral", "periportal"]
                    },
                    "gene": {
                        "type": "string",
                        "description": "Search for a specific gene"
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of results (default: 20)",
                        "default": 20
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="query_lipid_species_zonation",
            description="Query individual lipid species for zonation associations (periportal vs pericentral) with gonad-specific bias. Shows which lipids are associated with specific spatial zones and which genes correlate with them.",
            inputSchema={
                "type": "object",
                "properties": {
                    "lipid": {
                        "type": "string",
                        "description": "Lipid name or class to search (e.g., 'TG', 'Cer', 'PC(16:0')"
                    },
                    "zone": {
                        "type": "string",
                        "description": "Filter by zonation: 'periportal', 'pericentral', or 'all'",
                        "enum": ["periportal", "pericentral", "all"]
                    },
                    "gonad_bias": {
                        "type": "string",
                        "description": "Filter by gonad bias: 'Testes', 'Ovaries', or 'all'",
                        "enum": ["Testes", "Ovaries", "all"]
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of results (default: 20)",
                        "default": 20
                    }
                },
                "required": ["lipid"]
            }
        ),
        Tool(
            name="query_variance_analysis",
            description="Query lipid variance differences between gonadal groups. Tests whether lipids show different variability (not just mean) in testes vs ovaries samples. Key finding: testes samples are significantly less variable overall (Wilcoxon p = 2.2e-07).",
            inputSchema={
                "type": "object",
                "properties": {
                    "lipid": {
                        "type": "string",
                        "description": "Lipid name or class to search (optional, returns overview if empty)"
                    },
                    "level": {
                        "type": "string",
                        "description": "Analysis level: 'lipid' for individual species, 'class' for lipid classes, 'summary' for overview",
                        "enum": ["lipid", "class", "summary"],
                        "default": "summary"
                    },
                    "significant_only": {
                        "type": "boolean",
                        "description": "Only return significant results (Levene's p < 0.05)",
                        "default": False
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of results (default: 20)",
                        "default": 20
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="count_significant_features",
            description="Count genes or lipids that pass significance and fold-change thresholds for a given contrast. Returns total counts, up/down split, and optionally a per-class or per-cell-type breakdown. Use this to get precise numbers for the paper (e.g., 'how many genes with FDR <= 0.05 and |logFC| > 0.4 in hepatocytes for testes vs ovaries?').",
            inputSchema={
                "type": "object",
                "properties": {
                    "feature_type": {
                        "type": "string",
                        "description": "Type of features to count",
                        "enum": ["genes", "lipids"]
                    },
                    "contrast": {
                        "type": "string",
                        "description": "Contrast to query. For genes: 'testes_vs_ova', 'XY_vs_XX', 'XXT_vs_XXO', 'XYT_vs_XYO', 'XYO_vs_XXO', 'XYT_vs_XXT'. For lipids: 'testes_vs_ova' (only contrast available).",
                        "default": "testes_vs_ova"
                    },
                    "cell_type": {
                        "type": "string",
                        "description": "Cell type (genes only). Default: 'Hepatocytes'.",
                        "enum": ["Hepatocytes", "Endothelial_cells", "Immune_cells", "Stellate_cells", "Cholangiocytes", "all"],
                        "default": "Hepatocytes"
                    },
                    "fdr_threshold": {
                        "type": "number",
                        "description": "FDR (adj.P.Val) threshold (default: 0.05)",
                        "default": 0.05
                    },
                    "pvalue_threshold": {
                        "type": "number",
                        "description": "Raw P.Value threshold (optional, use instead of fdr_threshold if set)"
                    },
                    "logfc_threshold": {
                        "type": "number",
                        "description": "Minimum absolute logFC threshold (default: 0, meaning no FC filter)",
                        "default": 0
                    },
                    "breakdown": {
                        "type": "boolean",
                        "description": "If true, also show per-cell-type breakdown (genes) or per-class breakdown (lipids)",
                        "default": False
                    }
                },
                "required": ["feature_type"]
            }
        ),
        Tool(
            name="count_zonation_categories",
            description="Count genes in zonation-vs-sex-bias categories using FDR and logFC thresholds. Classifies genes into 9 categories (Testes Pericentral, Testes Periportal, Testes only, Ovaries Pericentral, Ovaries Periportal, Ovaries only, Pericentral only, Periportal only, None) based on significance in the zonation (portal_vs_Central) and sex bias (testes_vs_ova) contrasts. Use this to answer: 'how many genes are sex-biased and zonation-biased at FDR 5% and |logFC| > 0.4?'",
            inputSchema={
                "type": "object",
                "properties": {
                    "fdr_threshold": {
                        "type": "number",
                        "description": "FDR (adj.P.Val) threshold (default: 0.05)",
                        "default": 0.05
                    },
                    "pvalue_threshold": {
                        "type": "number",
                        "description": "Raw P.Value threshold (optional, use instead of fdr_threshold if set)"
                    },
                    "logfc_threshold": {
                        "type": "number",
                        "description": "Minimum absolute logFC threshold (default: 0, meaning no FC filter)",
                        "default": 0
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="compare_fcg_gdx",
            description="Compare FCG (Four Core Genotypes) and GDX (gonadectomy hormone replacement) results. Computes correlations between FCG gonadal effects and GDX gonadectomy effects for genes or lipids. Can also look up specific genes/lipids to compare their FCG and GDX t-values side by side.",
            inputSchema={
                "type": "object",
                "properties": {
                    "feature_type": {
                        "type": "string",
                        "description": "Type of features to compare",
                        "enum": ["genes", "lipids"]
                    },
                    "cell_type": {
                        "type": "string",
                        "description": "Cell type for gene comparison (default: Hepatocytes). Only used for genes.",
                        "enum": ["Hepatocytes", "Endothelial_cells", "Immune_cells", "Stellate_cells", "Cholangiocytes"],
                        "default": "Hepatocytes"
                    },
                    "contrast": {
                        "type": "string",
                        "description": "FCG contrast to compare against GDX (default: testes_vs_ova)",
                        "default": "testes_vs_ova"
                    },
                    "query": {
                        "type": "string",
                        "description": "Specific gene or lipid name(s) to look up (comma-separated). If omitted, returns the global correlation summary."
                    },
                    "gdx_contrast": {
                        "type": "string",
                        "description": "GDX contrast: 'testes_gonadosectomy' (castration) or 'ova_gonadosectomy' (ovariectomy). Default: testes_gonadosectomy.",
                        "enum": ["testes_gonadosectomy", "ova_gonadosectomy", "both"],
                        "default": "both"
                    }
                },
                "required": ["feature_type"]
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
        elif name == "get_lipid_de_summary":
            return handle_lipid_de_summary(loader)
        elif name == "query_lipid_network":
            return handle_lipid_network(loader, arguments)
        elif name == "query_tf_gene_lipid_network":
            return handle_tf_network(loader, arguments)
        elif name == "query_sex_specific_zonation_genes":
            return handle_sex_specific_zonation(loader, arguments)
        elif name == "query_lipid_species_zonation":
            return handle_lipid_species_zonation(loader, arguments)
        elif name == "query_variance_analysis":
            return handle_variance_analysis(loader, arguments)
        elif name == "count_significant_features":
            return handle_count_significant(loader, arguments)
        elif name == "count_zonation_categories":
            return handle_count_zonation_categories(loader, arguments)
        elif name == "compare_fcg_gdx":
            return handle_compare_fcg_gdx(loader, arguments)
        else:
            return [TextContent(type="text", text=f"Unknown tool: {name}")]

    except Exception as e:
        logger.exception(f"Error in {name}")
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_summary() -> list[TextContent]:
    text = """# FCGL Dataset Summary

## Overview
Multi-omic study of sex-specific liver physiology using the Four Core Genotypes (FCG) mouse model.
**Paper**: "Gonadal hormone coordination of sex-specific adult liver physiology"

## Key Statistics
- **snRNA-seq**: 466,986 nuclei from 48 liver samples (43 FCG + 5 WT males)
- **Lipidomics**: 1,254 lipid species analyzed (after QC filtering) across 20 lipid classes from 48 samples (21 ovaries, 27 testes)
- **Cell types**: Hepatocytes, Endothelial cells, Immune cells, Stellate cells, Cholangiocytes

## FCG Groups
| Group | Chromosomes | Gonads | N |
|-------|-------------|--------|---|
| XXO | XX | Ovaries | ~10 |
| XXT | XX | Testes | ~11 |
| XYO | XY | Ovaries | ~11 |
| XYT | XY | Testes | ~11 |
| XYT-WT | XY | Testes (wild-type) | 5 |

## Key Findings

### RNA-seq (Workflow 02)
- **~4,000 gonad-dependent genes** in hepatocytes (FDR ≤ 0.05)
- Minimal sex chromosome effects (mainly Y-linked or X-escape genes)
- Hepatocytes show strongest gonad-dependent gene expression among all cell types

### Lipidomics (Workflow 03)
- **773/1,254 lipids (61.6%)** are gonad-dependent (FDR ≤ 0.05)
- 488 male-biased (logFC > 0), 285 female-biased (logFC < 0)
- **NO lipids** show sex chromosome effects
- **Male-biased lipid classes** (enrichment scores): PI (4.68), TG (4.67), PE (3.34), DG (2.40), LPE (1.96)
- **Female-biased lipid classes**: SM (-8.67), Cer (-8.55), HexCer (-6.95), ST (-4.71)
- Top male/female-biased lipids form interconnected LINEX2 network
- GDX validation: FCG lipid effects correlate with testes gonadectomy (r = -0.721, p = 6.5e-118)

### Multi-Omic Integration (Workflow 04)
- MOFA Factor 1 captures strongest gonadal association
- **Key regulatory cascade**: TFs (Nfic, Dbp, Esr1) → Cyp genes → Lipid classes
- **Key genes**: Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9
- **Key TFs**: Esr1, Ar, Stat5a/b, Nfic, Dbp

### Zonation (Workflow 05)
- 43% of gonad-dependent genes are pericentral-biased
- 18% are periportal-biased
- **137 genes** show sex-specific zonation (bias only in one gonadal context)

### Lipid Structural Preferences (Workflow 06)
- Gene-lipid correlations depend on lipid chain length and saturation
- Multiple genes show structural substrate preferences

### Validation (Workflow 07)
- siRNA knockdowns validated for Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9, Mtmr11

### Variance Analysis (Workflow 08)
- Testes samples are significantly less variable than ovaries (Wilcoxon p = 2.2e-07)
- SHexCer shows significantly different variance by class

## Available Queries

### Gene Expression
- `query_gene`: DE stats for any gene across all contrasts and cell types
- `get_top_genes`: Top DE genes for a specific cell type and contrast
- `compare_gene_celltypes`: Cross-cell type comparison for a gene
- `count_significant_features`: Count genes/lipids passing significance + FC thresholds

### Lipidomics
- `get_lipid_de_summary`: Comprehensive lipid DE summary (counts, per-class breakdown, enrichment, GDX correlation)
- `query_lipid`: DE stats for individual lipids (gonadal effect)
- `query_lipid_class_enrichment`: Lipid class enrichment scores
- `query_lipid_network`: LINEX2 lipid conversion network
- `query_variance_analysis`: Lipid variance differences between gonadal groups

### Multi-Omic Integration
- `query_gene_lipid_network`: MOFA Factor 1 gene-lipid correlations
- `query_tf_gene_lipid_network`: TF→Gene→Lipid regulatory network
- `get_structural_associations`: Lipid structural preferences (chain length, saturation)

### Zonation
- `get_zonation_info`: Zonation analysis overview
- `query_gene_zonation`: Gene zonation bias (periportal vs pericentral)
- `query_lipid_zonation`: Lipid class zonation associations
- `query_sex_specific_zonation_genes`: 137 sex-specific zonation genes
- `query_lipid_species_zonation`: Per-lipid species zonation associations

### Validation
- `get_validation_results`: siRNA knockdown validation results
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

        # Search (escape regex special chars like parentheses in lipid names)
        lipid_escaped = re.escape(lipid)
        mask = df[lipid_col].str.contains(lipid_escaped, case=False, na=False)
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
            all_matches = df[df[lipid_col].str.contains(lipid_escaped, case=False, na=False)]
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

        # Sex-specific zonation genes
        try:
            ssz = loader.get_sex_specific_zonation()
            cat_counts = ssz['category'].value_counts()
            text += f"\n\n## Sex-Specific Zonation Genes\n\n"
            text += f"**{len(ssz)} genes** show zonation bias only in one gonadal context:\n"
            for cat, count in cat_counts.items():
                text += f"- {cat}: {count} genes\n"
            text += "\nUse `query_sex_specific_zonation_genes` for details.\n"
        except:
            pass

        # Lipid zonation summary
        try:
            lipid_df = loader.get_lipid_zonation()
            text += "\n## Lipid Class Zonation Associations\n\n"
            text += "TPC=Testes Pericentral, TPP=Testes Periportal, OPC=Ovaries Pericentral, OPP=Ovaries Periportal\n\n"
            text += lipid_df.to_string(index=False)
        except:
            pass

        # Top zonation pathways
        try:
            pathways = loader.get_zonation_pathways()
            # Get top pathways by absolute t-value for portal_vs_Central contrast
            if 'portal_vs_Central' in pathways.columns:
                top_pw = pathways.copy()
                top_pw['abs_t'] = top_pw['portal_vs_Central'].abs()
                top_pw = top_pw.sort_values('abs_t', ascending=False).head(10)
                text += "\n\n## Top Zonation-Associated Pathways (periportal vs pericentral)\n\n"
                text += "| Pathway | t-value | Bias |\n"
                text += "|---------|---------|------|\n"
                for _, row in top_pw.iterrows():
                    val = row['portal_vs_Central']
                    bias = "periportal" if val > 0 else "pericentral"
                    pw_name = row['source'].replace('GOBP_', '').replace('_', ' ').title()
                    text += f"| {pw_name} | {val:.2f} | {bias} |\n"
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
        text += "## Overview\n\n"
        text += df.to_string(index=False)

        # Global network validation scores
        try:
            global_val = loader.get_global_validation()
            # Show per-gene scores (each gene knockdown tested against all network genes)
            # Group by condition (KD gene) and show top scores
            text += "\n\n## Global Network Validation Scores\n\n"
            text += "DecoupleR ULM scores testing whether gene knockdown effects match predicted network relationships.\n\n"

            for kd_gene in global_val['condition'].unique():
                gene_data = global_val[global_val['condition'] == kd_gene].copy()
                gene_data = gene_data.sort_values('p_value')
                sig = gene_data[gene_data['p_value'] < 0.05]
                text += f"\n### {kd_gene} knockdown\n"
                text += f"- Tested against {len(gene_data)} network genes\n"
                text += f"- {len(sig)} significant (p < 0.05)\n"

                # Show top 5 hits
                top = gene_data.head(5)
                text += "| Source Gene | Score | p-value |\n"
                text += "|------------|-------|---------|\n"
                for _, row in top.iterrows():
                    pval = row['p_value']
                    pval_str = f"{pval:.2e}" if pval < 0.01 else f"{pval:.4f}"
                    text += f"| {row['source']} | {row['score']:.3f} | {pval_str} |\n"
        except:
            pass

        # Zonation network validation scores
        try:
            zon_val = loader.get_zonation_validation()
            text += "\n\n## Zonation Network Validation Scores\n\n"
            text += "Validation against zonation-specific gene-lipid networks.\n\n"

            for kd_gene in zon_val['condition'].unique():
                gene_data = zon_val[zon_val['condition'] == kd_gene].copy()
                gene_data = gene_data.sort_values('p_value')
                sig = gene_data[gene_data['p_value'] < 0.05]

                # Get zonation context if available
                zon_context = gene_data['gonade_zonation'].iloc[0] if 'gonade_zonation' in gene_data.columns else ""
                text += f"\n### {kd_gene} knockdown"
                if zon_context:
                    text += f" ({zon_context})"
                text += f"\n- {len(sig)}/{len(gene_data)} significant (p < 0.05)\n"

                top = gene_data.head(3)
                text += "| Source Gene | Score | p-value |\n"
                text += "|------------|-------|---------|\n"
                for _, row in top.iterrows():
                    pval = row['p_value']
                    pval_str = f"{pval:.2e}" if pval < 0.01 else f"{pval:.4f}"
                    text += f"| {row['source']} | {row['score']:.3f} | {pval_str} |\n"
        except:
            pass

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


def handle_lipid_de_summary(loader: FCGLDataLoader) -> list[TextContent]:
    """Comprehensive summary of gonad-dependent lipid DE results."""
    try:
        from scipy import stats as scipy_stats

        df = loader.get_lipid_de()
        lipid_col = 'ID' if 'ID' in df.columns else df.columns[0]

        # Overall counts
        total = len(df)
        sig = df[df['adj.P.Val'] < 0.05]
        n_sig = len(sig)
        male_biased = sig[sig['logFC'] > 0]
        female_biased = sig[sig['logFC'] < 0]

        text = "# Gonad-Dependent Lipid DE Summary\n\n"
        text += "## Overall Counts (Testes vs Ovaries, FDR ≤ 0.05)\n\n"
        text += f"- **Total lipids analyzed**: {total}\n"
        text += f"- **Significant**: {n_sig} ({100*n_sig/total:.1f}%)\n"
        text += f"- **Male-biased** (logFC > 0, higher in testes): {len(male_biased)}\n"
        text += f"- **Female-biased** (logFC < 0, higher in ovaries): {len(female_biased)}\n"

        # Per-class breakdown
        def extract_class(name):
            # Extract lipid class from name like "PC(16:0/18:1)" -> "PC"
            if '(' in name:
                return name.split('(')[0]
            return name.split(' ')[0] if ' ' in name else name

        df['lipid_class'] = df[lipid_col].apply(extract_class)
        sig_with_class = df[df['adj.P.Val'] < 0.05].copy()

        text += "\n## Per-Class Breakdown (FDR ≤ 0.05)\n\n"
        text += "| Lipid Class | Total | Male-biased | Female-biased | % Significant |\n"
        text += "|-------------|-------|-------------|---------------|---------------|\n"

        class_stats = []
        for cls in sorted(df['lipid_class'].unique()):
            total_cls = len(df[df['lipid_class'] == cls])
            sig_cls = sig_with_class[sig_with_class['lipid_class'] == cls]
            n_male = len(sig_cls[sig_cls['logFC'] > 0])
            n_female = len(sig_cls[sig_cls['logFC'] < 0])
            n_total_sig = n_male + n_female
            pct = 100 * n_total_sig / total_cls if total_cls > 0 else 0
            if n_total_sig > 0:
                class_stats.append((cls, total_cls, n_total_sig, n_male, n_female, pct))

        # Sort by total significant (descending)
        class_stats.sort(key=lambda x: x[2], reverse=True)
        for cls, total_cls, n_total_sig, n_male, n_female, pct in class_stats:
            text += f"| {cls} | {total_cls} ({n_total_sig} sig) | {n_male} | {n_female} | {pct:.0f}% |\n"

        # Enrichment scores
        try:
            enr = loader.get_lipid_class_enrichment()
            enr = enr.sort_values('enrichment_score', ascending=False)

            text += "\n## Lipid Class Enrichment Scores (decoupleR ULM)\n\n"
            text += "Positive = male-biased (testes), Negative = female-biased (ovaries)\n\n"
            text += "| Lipid Class | Enrichment Score | p-value | Significance |\n"
            text += "|-------------|-----------------|---------|-------------|\n"

            for _, row in enr.iterrows():
                pval = row['p_value']
                pval_str = f"{pval:.2e}" if pval < 0.001 else f"{pval:.4f}"
                sig_marker = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
                text += f"| {row['lipid_class']} | {row['enrichment_score']:.3f} | {pval_str} | {sig_marker} |\n"
        except Exception:
            pass

        # GDX validation correlation
        try:
            gdx = loader.get_gdx_lipid_tvalues()
            merged = pd.merge(df[['ID', 't']], gdx[['ID', 'testes_gonadosectomy', 'ova_gonadosectomy']], on='ID')

            r_testes, p_testes = scipy_stats.pearsonr(merged['t'], merged['testes_gonadosectomy'])
            r_ova, p_ova = scipy_stats.pearsonr(merged['t'], merged['ova_gonadosectomy'])

            text += f"\n## GDX Validation (Lipid-Level Correlation)\n\n"
            text += f"Correlation between FCG gonadal effect t-values and GDX gonadectomy t-values:\n\n"
            text += f"- **Testes gonadectomy**: r = {r_testes:.3f}, p = {p_testes:.2e} (n = {len(merged)} matched lipids)\n"
            text += f"- **Ovaries gonadectomy**: r = {r_ova:.3f}, p = {p_ova:.2e}\n"
            text += f"\nInterpretation: Removing testes reverses the testes-biased lipid pattern (strong negative correlation),\n"
            text += f"while removing ovaries has no significant effect — confirming testes as the primary driver.\n"
        except ImportError:
            text += "\n*Note: scipy not installed — GDX correlation not computed.*\n"
        except Exception:
            pass

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_lipid_network(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Query the LINEX2 lipid conversion network."""
    lipid = args.get("lipid", "").strip()
    top_n = args.get("top_n", 20)

    if not lipid:
        return [TextContent(type="text", text="Error: lipid name/class required")]

    try:
        sif = loader.get_lipid_network_sif()
        att = loader.get_lipid_network_att()

        # Find edges involving the queried lipid (escape regex special chars like parentheses)
        lipid_escaped = re.escape(lipid)
        mask = (sif['source'].str.contains(lipid_escaped, case=False, na=False) |
                sif['target'].str.contains(lipid_escaped, case=False, na=False))
        edges = sif[mask].copy()

        if len(edges) == 0:
            return [TextContent(type="text",
                text=f"No network edges found for '{lipid}'. This lipid may not be in the significant gonad-dependent network.")]

        # Get all nodes involved
        all_nodes = set(edges['source'].tolist() + edges['target'].tolist())

        # Merge with node attributes
        node_info = att[att['ID'].isin(all_nodes)].copy() if 'ID' in att.columns else att[att.iloc[:, 0].isin(all_nodes)].copy()

        text = f"# LINEX2 Lipid Conversion Network: {lipid}\n\n"
        text += f"**Found {len(edges)} conversion relationships**\n\n"

        text += "## Network Edges (enzymatic conversions)\n\n"
        text += "| Source | Interaction | Target |\n"
        text += "|--------|-------------|--------|\n"
        for _, row in edges.head(top_n).iterrows():
            interaction = row.get('interaction', 'lipid_conversion')
            text += f"| {row['source']} | {interaction} | {row['target']} |\n"

        if len(node_info) > 0:
            text += "\n## Node Attributes\n\n"
            text += node_info.head(top_n).to_string(index=False)

        text += "\n\n*Network from LINEX2 lipid-enzyme conversion database (workflow 03)*"
        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_tf_network(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Query the TF-Gene-Lipid regulatory network."""
    query = args.get("query", "").strip()
    node_type_filter = args.get("node_type", "all")

    if not query:
        return [TextContent(type="text", text="Error: query required")]

    try:
        sif = loader.get_tf_network_sif()
        att = loader.get_tf_network_att()

        # Find edges involving the queried node
        mask = (sif['node1'].str.contains(query, case=False, na=False) |
                sif['node2'].str.contains(query, case=False, na=False))
        edges = sif[mask].copy()

        if len(edges) == 0:
            # Show available nodes
            all_nodes = set(sif['node1'].tolist() + sif['node2'].tolist())
            partial = [n for n in all_nodes if query.lower() in n.lower()]
            if partial:
                return [TextContent(type="text",
                    text=f"No exact match for '{query}'. Partial matches: {', '.join(partial[:10])}")]
            return [TextContent(type="text",
                text=f"No network edges found for '{query}' in the TF-Gene-Lipid network.")]

        # Filter by node type if requested
        if node_type_filter != "all" and 'type' in att.columns:
            typed_nodes = set(att[att['type'] == node_type_filter]['node'].tolist())
            edges = edges[
                edges['node1'].isin(typed_nodes) | edges['node2'].isin(typed_nodes)
            ]

        # Get node attributes for involved nodes
        all_nodes = set(edges['node1'].tolist() + edges['node2'].tolist())
        node_info = att[att['node'].isin(all_nodes)].copy()

        text = f"# TF-Gene-Lipid Regulatory Network: {query}\n\n"
        text += "Regulatory cascade: TF → Enzyme/Gene → Lipid Class\n\n"

        text += "## Network Edges\n\n"
        text += "| Node 1 | Node 2 | Correlation | |Correlation| |\n"
        text += "|--------|--------|-------------|---------------|\n"
        for _, row in edges.iterrows():
            text += f"| {row['node1']} | {row['node2']} | {row['mor']:.3f} | {row['abs_mor']:.3f} |\n"

        if len(node_info) > 0:
            text += "\n## Node Attributes\n\n"
            text += "| Node | Type | Weight |\n"
            text += "|------|------|--------|\n"
            for _, row in node_info.iterrows():
                text += f"| {row['node']} | {row['type']} | {row['weight']:.4f} |\n"

        # Summarize by type
        if 'type' in node_info.columns:
            type_counts = node_info['type'].value_counts()
            text += "\n## Node Type Summary\n"
            for t, c in type_counts.items():
                text += f"- **{t}**: {c} nodes\n"

        text += "\n*From MOFA Factor 1 integration (workflow 04). Compressed network showing top correlated features.*"
        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_sex_specific_zonation(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Query sex-specific zonation genes."""
    gonad = args.get("gonad", "")
    zone = args.get("zone", "")
    gene = args.get("gene", "")
    top_n = args.get("top_n", 20)

    try:
        df = loader.get_sex_specific_zonation()

        # Filter by gene if specified
        if gene:
            mask = df['gene'].str.contains(gene, case=False, na=False)
            df = df[mask]
            if len(df) == 0:
                return [TextContent(type="text", text=f"Gene '{gene}' not found among sex-specific zonation genes")]

        # Filter by gonad and zone using category column
        # Categories are like "Testes pericentral", "Ovaries periportal", etc.
        if gonad:
            df = df[df['category'].str.contains(gonad, case=False, na=False)]
        if zone:
            df = df[df['category'].str.contains(zone, case=False, na=False)]

        if len(df) == 0:
            return [TextContent(type="text", text="No genes match the criteria")]

        # Sort by absolute interaction t-value (strongest sex-specific zonation first)
        df = df.sort_values('interaction_t', key=abs, ascending=False)

        text = f"# Sex-Specific Zonation Genes\n\n"
        text += f"Genes with zonation bias (periportal vs pericentral) only in one gonadal context.\n\n"

        # Show category counts
        cat_counts = df['category'].value_counts()
        text += "## Category Counts\n"
        for cat, count in cat_counts.items():
            text += f"- **{cat}**: {count} genes\n"
        text += f"- **Total**: {len(df)} genes\n\n"

        text += "## Gene Details\n\n"
        text += "| Gene | Category | Zonation (Ovaries) t | Zonation (Testes) t | Interaction t | Gonad Effect t |\n"
        text += "|------|----------|---------------------|--------------------|--------------|--------------|\n"

        for _, row in df.head(top_n).iterrows():
            text += f"| {row['gene']} | {row['category']} | {row['zonation_ovaries_t']:.2f} | {row['zonation_testes_t']:.2f} | {row['interaction_t']:.2f} | {row['gonad_effect_t']:.2f} |\n"

        text += "\n**Interpretation:**\n"
        text += "- **Zonation (Ovaries) t**: periportal vs pericentral t-statistic in ovaries mice\n"
        text += "- **Zonation (Testes) t**: periportal vs pericentral t-statistic in testes mice\n"
        text += "- **Interaction t**: gonad × zonation interaction (significant = different zonation pattern by gonad)\n"
        text += "- **Gonad Effect t**: overall testes vs ovaries effect (positive = male-biased)\n"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_lipid_species_zonation(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Query per-lipid species zonation associations."""
    lipid = args.get("lipid", "").strip()
    zone = args.get("zone", "all")
    gonad_bias = args.get("gonad_bias", "all")
    top_n = args.get("top_n", 20)

    if not lipid:
        return [TextContent(type="text", text="Error: lipid name/class required")]

    try:
        df = loader.get_lipid_species_zonation()

        # Filter by lipid pattern (escape regex special chars)
        lipid_escaped = re.escape(lipid)
        mask = df['lipid'].str.contains(lipid_escaped, case=False, na=False)
        df = df[mask].copy()

        if len(df) == 0:
            return [TextContent(type="text", text=f"No lipid species matching '{lipid}' in zonation data")]

        # Filter by zone
        # TPC/OPC = pericentral, TPP/OPP = periportal
        if zone == "pericentral":
            df = df[(df['TPC'] != 0) | (df['OPC'] != 0)]
        elif zone == "periportal":
            df = df[(df['TPP'] != 0) | (df['OPP'] != 0)]

        # Filter by gonad bias
        # TPC/TPP = Testes associations, OPC/OPP = Ovaries associations
        if gonad_bias == "Testes":
            df = df[(df['TPC'] != 0) | (df['TPP'] != 0)]
        elif gonad_bias == "Ovaries":
            df = df[(df['OPC'] != 0) | (df['OPP'] != 0)]

        if len(df) == 0:
            return [TextContent(type="text", text="No lipids match the criteria")]

        # Sort by max correlation
        df = df.sort_values('max_cor', ascending=False)

        text = f"# Lipid Species Zonation: {lipid}\n\n"
        text += f"**Found {len(df)} lipids**\n\n"
        text += "TPC=Testes Pericentral, TPP=Testes Periportal, OPC=Ovaries Pericentral, OPP=Ovaries Periportal\n"
        text += "Values represent correlation with zonation-specific gene expression.\n\n"

        text += "| Lipid | TPC | TPP | OPC | OPP | Top Gene | Max Cor |\n"
        text += "|-------|-----|-----|-----|-----|----------|--------|\n"

        for _, row in df.head(top_n).iterrows():
            text += f"| {row['lipid']} | {row['TPC']:.3f} | {row['TPP']:.3f} | {row['OPC']:.3f} | {row['OPP']:.3f} | {row['gene']} | {row['max_cor']:.3f} |\n"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_variance_analysis(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Query lipid variance differences between gonadal groups."""
    lipid = args.get("lipid", "")
    level = args.get("level", "summary")
    significant_only = args.get("significant_only", False)
    top_n = args.get("top_n", 20)

    try:
        if level == "summary":
            # Return overview
            df = loader.get_variance_summary()
            text = "# Lipid Variance Analysis Summary\n\n"
            text += "Comparing lipid abundance variability between testes and ovaries samples.\n\n"
            for _, row in df.iterrows():
                text += f"- **{row['Metric']}**: {row['Value']}\n"

            # Also add class-level highlights
            try:
                class_df = loader.get_variance_by_class()
                sig_classes = class_df[class_df.get('sig_levene', class_df.get('sig_nominal', pd.Series(dtype=bool))) == True]
                if len(sig_classes) > 0:
                    text += "\n## Significant Lipid Classes (Levene's test)\n\n"
                    text += sig_classes.to_string(index=False)
            except:
                pass

            return [TextContent(type="text", text=text)]

        elif level == "class":
            df = loader.get_variance_by_class()

            if lipid:
                lipid_escaped = re.escape(lipid)
                df = df[df['lipid_class'].str.contains(lipid_escaped, case=False, na=False)]

            if significant_only:
                if 'sig_levene' in df.columns:
                    df = df[df['sig_levene'] == True]
                elif 'sig_nominal' in df.columns:
                    df = df[df['sig_nominal'] == True]

            if len(df) == 0:
                return [TextContent(type="text", text="No lipid classes match the criteria")]

            df = df.sort_values('levene_p')

            text = f"# Lipid Variance by Class\n\n"
            text += "| Class | N Lipids | CV Ovaries | CV Testes | Levene p | More Variable In |\n"
            text += "|-------|----------|------------|-----------|----------|------------------|\n"

            for _, row in df.head(top_n).iterrows():
                pval = row['levene_p']
                pval_str = f"{pval:.2e}" if pval < 0.001 else f"{pval:.4f}"
                more_var = row.get('more_var_in', row.get('label', 'NS'))
                text += f"| {row['lipid_class']} | {row['n_lipids']} | {row['cv_ovaries']:.3f} | {row['cv_testes']:.3f} | {pval_str} | {more_var} |\n"

            return [TextContent(type="text", text=text)]

        else:  # level == "lipid"
            df = loader.get_variance_results()

            if lipid:
                lipid_escaped = re.escape(lipid)
                df = df[df['lipid'].str.contains(lipid_escaped, case=False, na=False)]

            if significant_only:
                if 'sig_levene' in df.columns:
                    df = df[df['sig_levene'] == True]
                elif 'levene_padj' in df.columns:
                    df = df[df['levene_padj'] < 0.05]

            if len(df) == 0:
                return [TextContent(type="text", text="No lipids match the criteria")]

            df = df.sort_values('levene_p')

            text = f"# Lipid Variance Analysis (Individual Species)\n\n"
            text += f"**Found {len(df)} lipids**\n\n"
            text += "| Lipid | CV Ovaries | CV Testes | Levene p | Levene padj | More Variable In |\n"
            text += "|-------|------------|-----------|----------|-------------|------------------|\n"

            for _, row in df.head(top_n).iterrows():
                pval = row['levene_p']
                pval_str = f"{pval:.2e}" if pval < 0.001 else f"{pval:.4f}"
                padj = row.get('levene_padj', float('nan'))
                padj_str = f"{padj:.2e}" if pd.notna(padj) and padj < 0.001 else f"{padj:.4f}" if pd.notna(padj) else "NA"
                more_var = row.get('more_var_in', 'NS')
                text += f"| {row['lipid']} | {row['cv_ovaries']:.3f} | {row['cv_testes']:.3f} | {pval_str} | {padj_str} | {more_var} |\n"

            return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_count_significant(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Count genes or lipids passing significance and fold-change thresholds."""
    feature_type = args.get("feature_type", "genes")
    contrast_key = args.get("contrast", "testes_vs_ova")
    cell_type = args.get("cell_type", "Hepatocytes")
    fdr_threshold = args.get("fdr_threshold", 0.05)
    pvalue_threshold = args.get("pvalue_threshold", None)
    logfc_threshold = args.get("logfc_threshold", 0)
    breakdown = args.get("breakdown", False)

    # Map contrast keys to actual column values
    contrast_map = {
        "testes_vs_ova": "Testes vs Ovaries",
        "XY_vs_XX": "XY vs XX",
        "XXT_vs_XXO": "XXT vs XXO",
        "XYT_vs_XYO": "XYT vs XYO",
        "XYO_vs_XXO": "XYO vs XXO",
        "XYT_vs_XXT": "XYT vs XXT"
    }

    try:
        if feature_type == "genes":
            return _count_genes(loader, contrast_map, contrast_key, cell_type,
                                fdr_threshold, pvalue_threshold, logfc_threshold, breakdown)
        else:
            return _count_lipids(loader, fdr_threshold, pvalue_threshold,
                                 logfc_threshold, breakdown)
    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def _count_genes(loader, contrast_map, contrast_key, cell_type,
                 fdr_threshold, pvalue_threshold, logfc_threshold, breakdown):
    """Count significant genes."""
    df = loader.get_rna_ttop_full()
    contrast_name = contrast_map.get(contrast_key, contrast_key)

    # Determine significance column and threshold
    if pvalue_threshold is not None:
        sig_col = "P.Value"
        sig_thresh = pvalue_threshold
        sig_label = f"P.Value ≤ {pvalue_threshold}"
    else:
        sig_col = "adj.P.Val"
        sig_thresh = fdr_threshold
        sig_label = f"FDR ≤ {fdr_threshold}"

    fc_label = f"|logFC| ≥ {logfc_threshold}" if logfc_threshold > 0 else "no FC filter"

    text = f"# Significant Gene Counts\n\n"
    text += f"**Thresholds**: {sig_label}, {fc_label}\n"
    text += f"**Contrast**: {contrast_name}\n\n"

    def count_in_subset(sub):
        """Apply thresholds and return counts."""
        total = len(sub)
        mask = sub[sig_col] <= sig_thresh
        if logfc_threshold > 0:
            mask = mask & (sub['logFC'].abs() >= logfc_threshold)
        sig = sub[mask]
        up = sig[sig['logFC'] > 0]
        down = sig[sig['logFC'] < 0]
        return total, len(sig), len(up), len(down)

    if cell_type == "all" or breakdown:
        # Per-cell-type breakdown
        filtered = df[df['contrast'] == contrast_name]
        cell_types = ["Hepatocytes", "Endothelial cells", "Immune cells",
                      "Stellate cells", "Cholangiocytes"]

        text += "| Cell Type | Total Genes | Significant | Up (male-biased) | Down (female-biased) | % Significant |\n"
        text += "|-----------|-------------|-------------|-------------------|----------------------|---------------|\n"

        for ct in cell_types:
            sub = filtered[filtered['cell_type'] == ct]
            if len(sub) == 0:
                continue
            total, n_sig, n_up, n_down = count_in_subset(sub)
            pct = 100 * n_sig / total if total > 0 else 0
            text += f"| {ct} | {total} | {n_sig} | {n_up} | {n_down} | {pct:.1f}% |\n"

        if cell_type != "all":
            # Also show the specific cell type summary prominently
            ct_clean = cell_type.replace("_", " ")
            sub = filtered[filtered['cell_type'] == ct_clean]
            total, n_sig, n_up, n_down = count_in_subset(sub)
            text += f"\n**{ct_clean}**: {n_sig}/{total} significant ({n_up} up, {n_down} down)\n"
    else:
        # Single cell type
        ct_clean = cell_type.replace("_", " ")
        sub = df[(df['contrast'] == contrast_name) & (df['cell_type'] == ct_clean)]
        if len(sub) == 0:
            return [TextContent(type="text",
                text=f"No data found for contrast '{contrast_name}' and cell type '{ct_clean}'")]

        total, n_sig, n_up, n_down = count_in_subset(sub)
        pct = 100 * n_sig / total if total > 0 else 0

        text += f"**Cell type**: {ct_clean}\n\n"
        text += f"| Metric | Count |\n"
        text += f"|--------|-------|\n"
        text += f"| Total genes tested | {total} |\n"
        text += f"| Significant | {n_sig} ({pct:.1f}%) |\n"
        text += f"| Up-regulated (male-biased, logFC > 0) | {n_up} |\n"
        text += f"| Down-regulated (female-biased, logFC < 0) | {n_down} |\n"

    return [TextContent(type="text", text=text)]


def _count_lipids(loader, fdr_threshold, pvalue_threshold, logfc_threshold, breakdown):
    """Count significant lipids."""
    df = loader.get_lipid_de()
    lipid_col = 'ID' if 'ID' in df.columns else df.columns[0]

    # Determine significance column and threshold
    if pvalue_threshold is not None:
        sig_col = "P.Value"
        sig_thresh = pvalue_threshold
        sig_label = f"P.Value ≤ {pvalue_threshold}"
    else:
        sig_col = "adj.P.Val"
        sig_thresh = fdr_threshold
        sig_label = f"FDR ≤ {fdr_threshold}"

    fc_label = f"|logFC| ≥ {logfc_threshold}" if logfc_threshold > 0 else "no FC filter"

    # Apply thresholds
    total = len(df)
    mask = df[sig_col] <= sig_thresh
    if logfc_threshold > 0:
        mask = mask & (df['logFC'].abs() >= logfc_threshold)
    sig = df[mask]
    up = sig[sig['logFC'] > 0]
    down = sig[sig['logFC'] < 0]

    text = f"# Significant Lipid Counts\n\n"
    text += f"**Thresholds**: {sig_label}, {fc_label}\n"
    text += f"**Contrast**: Testes vs Ovaries (gonadal effect)\n\n"
    text += f"| Metric | Count |\n"
    text += f"|--------|-------|\n"
    text += f"| Total lipids tested | {total} |\n"
    text += f"| Significant | {len(sig)} ({100*len(sig)/total:.1f}%) |\n"
    text += f"| Male-biased (logFC > 0) | {len(up)} |\n"
    text += f"| Female-biased (logFC < 0) | {len(down)} |\n"

    if breakdown:
        # Per-class breakdown
        def extract_class(name):
            if '(' in name:
                return name.split('(')[0]
            return name.split(' ')[0] if ' ' in name else name

        df_with_class = df.copy()
        df_with_class['lipid_class'] = df_with_class[lipid_col].apply(extract_class)
        sig_with_class = df_with_class[mask]

        text += "\n## Per-Class Breakdown\n\n"
        text += "| Lipid Class | Total | Significant | Male-biased | Female-biased | % Sig |\n"
        text += "|-------------|-------|-------------|-------------|---------------|-------|\n"

        class_stats = []
        for cls in sorted(df_with_class['lipid_class'].unique()):
            total_cls = len(df_with_class[df_with_class['lipid_class'] == cls])
            sig_cls = sig_with_class[sig_with_class['lipid_class'] == cls]
            n_male = len(sig_cls[sig_cls['logFC'] > 0])
            n_female = len(sig_cls[sig_cls['logFC'] < 0])
            n_total_sig = n_male + n_female
            pct = 100 * n_total_sig / total_cls if total_cls > 0 else 0
            class_stats.append((cls, total_cls, n_total_sig, n_male, n_female, pct))

        class_stats.sort(key=lambda x: x[2], reverse=True)
        for cls, total_cls, n_total_sig, n_male, n_female, pct in class_stats:
            text += f"| {cls} | {total_cls} | {n_total_sig} | {n_male} | {n_female} | {pct:.0f}% |\n"

    return [TextContent(type="text", text=text)]


def handle_count_zonation_categories(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Count genes in zonation-vs-sex-bias categories using FDR/p-value and logFC thresholds."""
    fdr_threshold = args.get("fdr_threshold", 0.05)
    pvalue_threshold = args.get("pvalue_threshold", None)
    logfc_threshold = args.get("logfc_threshold", 0)

    try:
        df = loader.get_zonation_ttop_long()

        # Determine significance column and threshold
        if pvalue_threshold is not None:
            sig_col = "P.Value"
            sig_thresh = pvalue_threshold
            sig_label = f"P.Value ≤ {pvalue_threshold}"
        else:
            sig_col = "adj.P.Val"
            sig_thresh = fdr_threshold
            sig_label = f"FDR ≤ {fdr_threshold}"

        fc_label = f"|logFC| ≥ {logfc_threshold}" if logfc_threshold > 0 else "no FC filter"

        # Get zonation contrast (portal_vs_Central) and sex contrast (testes_vs_ova)
        zon_df = df[df['contrast'] == 'portal_vs_Central'].copy()
        sex_df = df[df['contrast'] == 'testes_vs_ova'].copy()

        if len(zon_df) == 0 or len(sex_df) == 0:
            return [TextContent(type="text", text="Error: zonation or sex contrast data not found")]

        # Apply significance + logFC filter
        zon_sig_mask = zon_df[sig_col] <= sig_thresh
        if logfc_threshold > 0:
            zon_sig_mask = zon_sig_mask & (zon_df['logFC'].abs() >= logfc_threshold)

        sex_sig_mask = sex_df[sig_col] <= sig_thresh
        if logfc_threshold > 0:
            sex_sig_mask = sex_sig_mask & (sex_df['logFC'].abs() >= logfc_threshold)

        # Build gene-level classification
        zon_sig_genes = set(zon_df.loc[zon_sig_mask, 'gene'])
        sex_sig_genes = set(sex_df.loc[sex_sig_mask, 'gene'])

        # Get logFC direction for each gene
        zon_lfc = zon_df.set_index('gene')['logFC']  # positive = periportal
        sex_lfc = sex_df.set_index('gene')['logFC']  # positive = testes-biased

        all_genes = set(zon_df['gene'])
        total = len(all_genes)

        # Categorize
        categories = {
            "Testes Pericentral": 0,
            "Testes Periportal": 0,
            "Testes (sex only)": 0,
            "Ovaries Pericentral": 0,
            "Ovaries Periportal": 0,
            "Ovaries (sex only)": 0,
            "Pericentral only": 0,
            "Periportal only": 0,
            "None": 0
        }

        for gene in all_genes:
            is_zon = gene in zon_sig_genes
            is_sex = gene in sex_sig_genes
            zon_dir = zon_lfc.get(gene, 0)  # positive = periportal, negative = pericentral
            sex_dir = sex_lfc.get(gene, 0)  # positive = testes, negative = ovaries

            if is_zon and is_sex:
                if sex_dir > 0 and zon_dir <= 0:
                    categories["Testes Pericentral"] += 1
                elif sex_dir > 0 and zon_dir > 0:
                    categories["Testes Periportal"] += 1
                elif sex_dir <= 0 and zon_dir <= 0:
                    categories["Ovaries Pericentral"] += 1
                else:
                    categories["Ovaries Periportal"] += 1
            elif is_sex and not is_zon:
                if sex_dir > 0:
                    categories["Testes (sex only)"] += 1
                else:
                    categories["Ovaries (sex only)"] += 1
            elif is_zon and not is_sex:
                if zon_dir <= 0:
                    categories["Pericentral only"] += 1
                else:
                    categories["Periportal only"] += 1
            else:
                categories["None"] += 1

        # Summary groupings
        n_zon = len(zon_sig_genes)
        n_sex = len(sex_sig_genes)
        n_both = len(zon_sig_genes & sex_sig_genes)
        n_zon_only = len(zon_sig_genes - sex_sig_genes)
        n_sex_only = len(sex_sig_genes - zon_sig_genes)
        n_neither = total - len(zon_sig_genes | sex_sig_genes)

        text = f"# Zonation × Sex Bias Gene Categories\n\n"
        text += f"**Thresholds**: {sig_label}, {fc_label}\n"
        text += f"**Total genes**: {total}\n\n"

        text += "## Quadrant Counts (both zonation AND sex biased)\n\n"
        text += "| | Pericentral | Periportal |\n"
        text += "|---|---|---|\n"
        text += f"| **Testes** | {categories['Testes Pericentral']} | {categories['Testes Periportal']} |\n"
        text += f"| **Ovaries** | {categories['Ovaries Pericentral']} | {categories['Ovaries Periportal']} |\n"

        text += f"\n## Single-Bias Categories\n\n"
        text += f"| Category | Count |\n"
        text += f"|----------|-------|\n"
        text += f"| Testes (sex only) | {categories['Testes (sex only)']} |\n"
        text += f"| Ovaries (sex only) | {categories['Ovaries (sex only)']} |\n"
        text += f"| Pericentral only | {categories['Pericentral only']} |\n"
        text += f"| Periportal only | {categories['Periportal only']} |\n"
        text += f"| None | {categories['None']} |\n"

        text += f"\n## Summary\n\n"
        text += f"| Grouping | Count |\n"
        text += f"|----------|-------|\n"
        text += f"| Zonation-biased (any) | {n_zon} |\n"
        text += f"| Sex-biased (any) | {n_sex} |\n"
        text += f"| Both zonation AND sex biased | {n_both} |\n"
        text += f"| Zonation-biased but NOT sex-biased | {n_zon_only} |\n"
        text += f"| Sex-biased but NOT zonation-biased | {n_sex_only} |\n"
        text += f"| Neither | {n_neither} |\n"

        return [TextContent(type="text", text=text)]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


def handle_compare_fcg_gdx(loader: FCGLDataLoader, args: dict) -> list[TextContent]:
    """Compare FCG and GDX results: global correlation and per-feature lookup."""
    from scipy import stats as scipy_stats

    feature_type = args.get("feature_type", "genes")
    cell_type = args.get("cell_type", "Hepatocytes").replace("_", " ")
    contrast_key = args.get("contrast", "testes_vs_ova")
    query = args.get("query", "")
    gdx_contrast = args.get("gdx_contrast", "both")

    contrast_map = {
        "testes_vs_ova": "Testes vs Ovaries",
        "XY_vs_XX": "XY vs XX",
        "XXT_vs_XXO": "XXT vs XXO",
        "XYT_vs_XYO": "XYT vs XYO",
        "XYO_vs_XXO": "XYO vs XXO",
        "XYT_vs_XXT": "XYT vs XXT"
    }
    contrast_name = contrast_map.get(contrast_key, contrast_key)

    try:
        if feature_type == "genes":
            fcg_df = loader.get_rna_tvalues()
            gdx_df = loader.get_gdx_rna_tvalues()

            # Find the FCG column for the given contrast + cell type
            fcg_col = None
            for col in fcg_df.columns:
                if contrast_name in col and cell_type in col:
                    fcg_col = col
                    break
            if fcg_col is None:
                return [TextContent(type="text",
                    text=f"No FCG column found for '{contrast_name}' + '{cell_type}'")]

            id_col = 'ID' if 'ID' in fcg_df.columns else fcg_df.columns[0]
            merged = pd.merge(
                fcg_df[[id_col, fcg_col]].rename(columns={id_col: 'ID', fcg_col: 'fcg_t'}),
                gdx_df,
                on='ID'
            ).dropna()

            feature_label = "genes"
            data_label = f"FCG {contrast_name} ({cell_type}) vs GDX"

        else:  # lipids
            fcg_df = loader.get_lipid_de()
            gdx_df = loader.get_gdx_lipid_tvalues()

            id_col = 'ID' if 'ID' in fcg_df.columns else fcg_df.columns[0]
            merged = pd.merge(
                fcg_df[[id_col, 't']].rename(columns={id_col: 'ID', 't': 'fcg_t'}),
                gdx_df,
                on='ID'
            ).dropna()

            feature_label = "lipids"
            data_label = f"FCG Testes vs Ovaries vs GDX"

        text = f"# FCG vs GDX Comparison ({feature_label.title()})\n\n"
        text += f"**{data_label}**\n"
        text += f"**Matched {feature_label}**: {len(merged)}\n\n"

        # Global correlation
        gdx_contrasts = []
        if gdx_contrast in ("testes_gonadosectomy", "both"):
            gdx_contrasts.append(("testes_gonadosectomy", "Testes gonadectomy (castration)"))
        if gdx_contrast in ("ova_gonadosectomy", "both"):
            gdx_contrasts.append(("ova_gonadosectomy", "Ovaries gonadectomy (ovariectomy)"))

        text += "## Global Correlation\n\n"
        text += "| GDX Contrast | Pearson r | p-value | n |\n"
        text += "|--------------|-----------|---------|---|\n"

        for gdx_col, gdx_label in gdx_contrasts:
            r, p = scipy_stats.pearsonr(merged['fcg_t'], merged[gdx_col])
            p_str = f"{p:.2e}" if p < 0.001 else f"{p:.4f}"
            text += f"| {gdx_label} | {r:.3f} | {p_str} | {len(merged)} |\n"

        text += "\n*Negative correlation with castration = testes-driven effect confirmed by removal.*\n"

        # Per-feature lookup
        if query:
            queries = [q.strip() for q in query.split(",")]
            text += f"\n## Specific {feature_label.title()} Lookup\n\n"
            text += f"| {feature_label[:-1].title()} | FCG t-value |"
            for gdx_col, gdx_label in gdx_contrasts:
                text += f" GDX {gdx_label.split('(')[0].strip()} t |"
            text += " Interpretation |\n"
            text += "|" + "---|" * (3 + len(gdx_contrasts)) + "\n"

            for q in queries:
                match = merged[merged['ID'].str.lower() == q.lower()]
                if len(match) == 0:
                    # Try partial match
                    q_escaped = re.escape(q)
                    match = merged[merged['ID'].str.contains(q_escaped, case=False, na=False)]
                    if len(match) == 0:
                        text += f"| {q} | *not found* |"
                        text += " - |" * len(gdx_contrasts)
                        text += " - |\n"
                        continue
                    match = match.head(5)

                for _, row in match.iterrows():
                    fcg_t = row['fcg_t']
                    text += f"| {row['ID']} | {fcg_t:.2f} |"

                    interps = []
                    for gdx_col, gdx_label in gdx_contrasts:
                        gdx_t = row[gdx_col]
                        text += f" {gdx_t:.2f} |"

                        # Interpret concordance
                        if "testes" in gdx_col:
                            if fcg_t < 0 and gdx_t > 0:
                                interps.append("FCG female-biased, reversed by castration")
                            elif fcg_t > 0 and gdx_t < 0:
                                interps.append("FCG male-biased, reversed by castration")
                            elif abs(gdx_t) < 2:
                                interps.append("weak GDX castration effect")
                        elif "ova" in gdx_col:
                            if abs(gdx_t) < 2:
                                interps.append("no ovariectomy effect")
                            else:
                                interps.append(f"ovariectomy effect (t={gdx_t:.1f})")

                    text += f" {'; '.join(interps)} |\n"
        else:
            # Show top concordant and discordant features
            text += f"\n## Top Concordant {feature_label.title()} (FCG effect reversed by castration)\n\n"

            if 'testes_gonadosectomy' in merged.columns:
                merged['concordance'] = -merged['fcg_t'] * merged['testes_gonadosectomy']
                top_concordant = merged.sort_values('concordance', ascending=False).head(10)

                text += f"| {feature_label[:-1].title()} | FCG t | GDX castration t | Product |\n"
                text += "|---|---|---|---|\n"
                for _, row in top_concordant.iterrows():
                    text += f"| {row['ID']} | {row['fcg_t']:.2f} | {row['testes_gonadosectomy']:.2f} | {row['concordance']:.1f} |\n"

                text += f"\n## Top Discordant {feature_label.title()} (FCG effect NOT reversed by castration)\n\n"
                top_discord = merged.sort_values('concordance', ascending=True).head(10)
                text += f"| {feature_label[:-1].title()} | FCG t | GDX castration t | Product |\n"
                text += "|---|---|---|---|\n"
                for _, row in top_discord.iterrows():
                    text += f"| {row['ID']} | {row['fcg_t']:.2f} | {row['testes_gonadosectomy']:.2f} | {row['concordance']:.1f} |\n"

        return [TextContent(type="text", text=text)]

    except ImportError:
        return [TextContent(type="text", text="Error: scipy is required for correlation computation. Install with: pip install scipy")]
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
