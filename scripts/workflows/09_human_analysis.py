# Workflow 09: Human liver comparison
#
# Purpose:
#   Compare the mouse FCG hepatocyte gonadal signature against sex-associated
#   gene expression in human bulk liver RNA-seq.
#
# External human data:
#   This workflow uses the Zhou et al. human liver cohort from GEO accession
#   GSE179250 ("A resource for integrated genomic analysis of the human liver").
#   GEO series page:
#       https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179250
#   The large GEO-derived files are not committed to this repository. Before
#   running, place the expected files under DATA_DIR with this layout:
#
#       bulk/zhou/GSE179250_gene_reads.txt
#       bulk/zhou/GSE179250_gene_tpm_normalized.txt
#       bulk/zhou/GSE179250_family.soft.gz
#
# Mouse signature input:
#   results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/hepa_gonad_signature.csv is the repository copy of the
#   mouse hepatocyte gonadal t-value signature. The script currently reads a
#   copy from DATA_DIR / "hepa_gonad_signature.csv". The later mouse-human
#   ranked-gene dot plots also expect DATA_DIR / "male_vs_female_DA_res.csv".
#
# Main analysis steps:
#   1. Load human bulk liver count and TPM matrices plus GEO metadata.
#   2. Build an AnnData object, normalize counts, and plot PCA by sex.
#   3. Translate the mouse signature to human ortholog symbols with pypath.
#   4. Score the translated signature across human samples with decoupler ULM.
#   5. Run PyDESeq2 male-vs-female DE for all ages and for ages 12-49.
#   6. Compare concordant mouse and human sex-biased genes.
#   7. Test young-vs-old signal among female human samples.
#
# Outputs:
#   Use --data-dir and --output to run outside the collaborator's HPC paths.

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Dict, Iterable, Tuple

import anndata as ad
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import mannwhitneyu
from adjustText import adjust_text

plt.rcParams.update({'font.size': 14})
sns.set_style("whitegrid", {"grid.color":'#EEE9E9'})

DATA_DIR = Path("/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/data/")
RESULTS_DIR = Path("/home/hd/hd_hd/hd_lx260/sds-hd/sd22b002/projects/sexdifferences_humanliver/results/bulk_manuscript")
GENDER_PALETTE = {'female': '#FF0000', 'male': '#5771CF'}

VOLCANO_PADJ_THR = 0.05
VOLCANO_LOG2FC_THR = 1.0

print(f"decoupler version: {dc.__version__}")


def save_fig(fig: plt.Figure, filename: str, bbox_inches = "tight") -> Path:
    """Save a matplotlib figure as PDF in the results directory."""
    path = RESULTS_DIR / filename
    fig.tight_layout()
    fig.savefig(path, format="pdf", bbox_inches=bbox_inches)
    plt.close(fig)
    return path


def save_clustergrid(grid, filename: str) -> Path:
    """Save a seaborn ClusterGrid without squeezing its layout."""
    path = RESULTS_DIR / filename
    grid.savefig(path, format="pdf")
    plt.close(grid.fig)
    return path



def load_metadata(soft_path: Path) -> pd.DataFrame:
    """Parse GEO soft file for sample metadata (title, sex, age, etc.)."""
    metadata: Dict[str, Dict[str, str]] = {}
    current = None
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("^SAMPLE"):
                current = line.split("=")[1].strip()
                metadata[current] = {}
            elif current and line.startswith("!Sample_title"):
                metadata[current]["title"] = line.split("=", 1)[1].strip()
            elif current and line.startswith("!Sample_characteristics_ch1"):
                value = line.split("=", 1)[1].strip()
                if ":" in value:
                    key, val = [x.strip() for x in value.split(":", 1)]
                    metadata[current][key.lower()] = val
    return pd.DataFrame([{"sample_id": sid, **attrs} for sid, attrs in metadata.items()])


def prepare_adata(counts: pd.DataFrame, counts_tpm: pd.DataFrame, meta_df: pd.DataFrame) -> ad.AnnData:
    """Build AnnData, keep raw counts, and normalize."""
    if not counts.index.equals(counts_tpm.index) or not counts.columns.equals(counts_tpm.columns):
        raise ValueError("GEO count and TPM tables must have identical gene and sample order")
    adata = ad.AnnData(
        X=counts.T.values,
        obs=meta_df.copy(),
        var=pd.DataFrame(index=counts.index),
    )
    adata.layers["tpm_normed"] = counts_tpm.T.values
    adata.var_names_make_unique()
    adata.raw = adata.copy()
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    sc.pp.filter_genes(adata, min_counts=10)
    return adata


def translate_signature(signature: pd.DataFrame, mapping: Dict[str, Iterable[str]], genes_present: Iterable[str]) -> pd.DataFrame:
    """Translate mouse signature genes to human and drop duplicates."""
    signature = signature.copy()
    signature["human_gene_orth"] = [
        sorted(mapping[g])[0] if g in mapping else g.upper() for g in signature["ID"]
    ]
    present_both = set(signature["human_gene_orth"]) & set(genes_present)
    print(f"{len(signature) - len(present_both)} genes cannot be translated with ortho & caps")

    signature = signature[signature["human_gene_orth"].isin(present_both)].copy()
    signature["source"] = "sex_signature"
    signature = signature.rename(columns={"t_Hepatocytes": "value"})

    duplicates = {
        key for key, group in signature.groupby("human_gene_orth") if len(group) > 1
    }
    signature = signature[~signature["human_gene_orth"].isin(duplicates)]
    return signature


def plot_gender_boxplot(ax: plt.Axes, data: pd.DataFrame, palette: Dict[str, str], y_col: str = "sex_signature", x_col: str = "sex") -> float | None:
    """Box/strip plot with Mann-Whitney p-value annotation."""
    data_plot = data.dropna(subset=[x_col, y_col])
    if data_plot.empty:
        ax.set_visible(False)
        return None
    order = [g for g in palette if g in data_plot[x_col].unique()]

    sns.boxplot(
        x=x_col,
        y=y_col,
        data=data_plot,
        order=order,
        palette=palette,
        ax=ax,
        showfliers=True,
        linewidth=0.8
    )
    sns.stripplot(
        x=x_col,
        y=y_col,
        data=data_plot,
        order=order,
        palette=palette,
        ax=ax,
        dodge=False,
        jitter=True,
        alpha=0.5,
        size=5,
        linewidth=0.8,
    )

    ax.set_xlabel("")
    ax.set_ylabel("enrichment")
    ax.set_xticklabels([label.capitalize() for label in order])

    cats = data_plot[x_col].unique()
    if len(cats) >= 2:
        vals_1 = data_plot[data_plot[x_col] == order[0]][y_col]
        vals_2 = data_plot[data_plot[x_col] == order[1]][y_col]
        if len(vals_1) > 0 and len(vals_2) > 0:
            p_value = mannwhitneyu(vals_1, vals_2, alternative="less").pvalue
            y_max = data_plot[y_col].max()
            y_min = data_plot[y_col].min()
            y_range = y_max - y_min or 1.0
            h = 0.08 * y_range
            line_y = y_max + h
            ax.plot([0, 0, 1, 1], [line_y, line_y + h, line_y + h, line_y], color="black", linewidth=1)
            ax.text(0.5, line_y + h * 1.1, f"p = {p_value:.3f}", ha="center", va="bottom", fontsize=10)
            ax.set_ylim(top=line_y + h * 2)
            ax.margins(y=0.1)
            return float(p_value)
    ax.margins(y=0.1)
    sns.despine(bottom = True, left = True)
    return None


def run_ulm_scores(adata: ad.AnnData, signature_df: pd.DataFrame) -> pd.DataFrame:
    """Run decoupler ULM and return the enrichment estimates."""
    dc.run_ulm(
        adata,
        net=signature_df,
        source="source",
        target="human_gene_orth",
        min_n=1,
        verbose=True,
        weight="value",
        use_raw=False,
    )
    estimate = adata.obsm["ulm_estimate"].copy()
    estimate["sample"] = adata.obs["sample_id"].values
    estimate["sex"] = adata.obs["sex"].values
    estimate["age"] = adata.obs["age"].astype(int).values
    return estimate


def calc_deg(adata_to_use: ad.AnnData, cat: str, cat1: str, cat2: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run DESeq2 differential expression."""
    adata_deg = adata_to_use.copy()
    adata_deg.X = adata_deg.layers["counts"]

    genders = adata_deg.obs[cat].astype(str).str.lower()
    gender_groups = genders.unique().tolist()
    if not {cat1, cat2}.issubset(set(gender_groups)):
        raise ValueError("Missing pseudobulks for DESeq2.")

    dds = DeseqDataSet(
        adata=adata_deg,
        design_factors=cat,
        ref_level=[cat, cat1],
        refit_cooks=True,
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=[cat, cat2, cat1])
    stat_res.summary()

    results_df = stat_res.results_df
    mat = results_df[["stat"]].T.rename(index={"stat": "Hepatocyte"})
    return results_df, mat


def prepare_volcano(results_df: pd.DataFrame) -> pd.DataFrame:
    """Add helper volcano columns."""
    results_df = results_df.copy()
    #results_df = results_df[results_df['padj']!= 0]
    results_df["padj"] = results_df["padj"].replace(0, np.nextafter(0, 1))
    results_df["neg_log10_padj"] = -np.log10(results_df["padj"])
    results_df["sig"] = "not_sig"
    results_df.loc[
        (results_df["padj"] < VOLCANO_PADJ_THR) & (results_df["log2FoldChange"] > VOLCANO_LOG2FC_THR),
        "sig",
    ] = "up"
    results_df.loc[
        (results_df["padj"] < VOLCANO_PADJ_THR) & (results_df["log2FoldChange"] < -VOLCANO_LOG2FC_THR),
        "sig",
    ] = "down"

    
    return results_df

def volcano_nice(df,
                 hAss=0.05,
                 FCIndex=None,
                 pValIndex=None,
                 IDIndex=None,
                 vAss=None,
                 label=False,
                 straight=False,
                 nlabels=None,
                 manual_labels=None,
                 col_high_neut_low=("#FF0000", "#424242", "#5771CF"),
                 ax=None, labels = []):
    """
    Python translation of the R function volcano_nice.

    Parameters
    ----------
    df : pandas.DataFrame
    hAss : float
        p-value / FDR threshold.
    FCIndex, pValIndex, IDIndex :
        Column *names* or integer positions for logFC, p-value and ID.
    vAss : float or None
        Fold-change threshold (if None, set from data range).
    label : bool
        If True, draw labels for selected points.
    straight : bool
        If True, use straight thresholds (vAss & hAss), otherwise curved lines.
    manual_labels : list-like or None
        IDs to keep as labels.
    col_high_neut_low : tuple(str, str, str)
        (high, neutral, low) colours = (red, grey, blue).
    ax : matplotlib Axes or None
        If None, a new figure/axes are created.
    """

    # --- prep data -----------------------------------------------------------
    df = df.dropna().copy()
    df.columns = list(df.columns)
    #df.columns.values[0] = "X1"

    hAssOri = hAss
    hAss = -np.log10(hAss)

    # rename "adj.P.Val" -> "FDR" anywhere in column names
    df.columns = [c.replace("adj.P.Val", "FDR") for c in df.columns]

    # resolve column indices / names
    def _get_col(col_or_idx):
        if isinstance(col_or_idx, int):
            return df.columns[col_or_idx]
        return col_or_idx

    fc_col = _get_col(FCIndex)
    p_col = _get_col(pValIndex)
    id_col = _get_col(IDIndex)

    df = df.rename(columns={fc_col: "logFC", p_col: "adj.P.Val"})
    df = df[df["adj.P.Val"] > 0].copy()
    df["neg_log10_adjP"] = -np.log10(df["adj.P.Val"])
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)

    # ranges
    fc_max = np.max(np.abs(df["logFC"]))
    p_max = np.nanmax(np.abs(df["neg_log10_adjP"].values))
    if np.isnan(p_max):
        p_max = 0.0

    print(f"pmax is {p_max}")
    print(f"fc_max is {fc_max}")

    xlimAbs = fc_max
    ylimAbs = p_max

    if vAss is None:
        vAss = xlimAbs / 10.0

    # curved threshold functions
    def xneg(x):
        return np.abs(hAss - 1 + x / (x + vAss))

    def xpos(x):
        return np.abs(hAss - 1 + x / (x - vAss))

    def test(x, y):
        # x: logFC, y: adj.P.Val
        if x < -vAss:
            return "1" if xneg(x) < -np.log10(y) else "0"
        elif x > vAss:
            return "1" if xpos(x) < -np.log10(y) else "0"
        else:
            return "0"

    # --- classification (couleur) ------------------------------------------
    if straight:
        df["couleur"] = np.where(
            (np.abs(df["logFC"]) >= vAss) & (df["adj.P.Val"] <= hAssOri),
            "1",
            "0",
        )
    else:
        df["couleur"] = [
            test(x, y) for x, y in zip(df["logFC"].values, df["adj.P.Val"].values)
        ]

    # order by p-value and set labels
    df = df.sort_values("adj.P.Val", ascending=True)

    print(df.head())
    df["condLabel"] = df[id_col]

    if manual_labels is None:
        manual_labels = []
    labels_to_keep = set(manual_labels)
    df.loc[~df["condLabel"].isin(labels_to_keep), "condLabel"] = np.nan

    # split "up" vs "down" for color
    df.loc[(df["couleur"] == "1") & (df["logFC"] < 0), "couleur"] = "2"

    print(df["couleur"].value_counts())

    # --- plotting -----------------------------------------------------------
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 5), dpi=150)

    # choose mapping 0/1/2 -> colours
    high, neut, low = col_high_neut_low
    if label:
        # c(col_high_neut_low[2], col_high_neut_low[3], col_high_neut_low[1])
        color_map = {"0": neut, "1": low, "2": high}
    else:
        # c(col_high_neut_low[2], col_high_neut_low[1], col_high_neut_low[3])
        color_map = {"0": neut, "1": high, "2": low}

    colours = df["couleur"].map(color_map)

    ax.scatter(
        df["logFC"],
        df["neg_log10_adjP"],
        c=colours,
        alpha=0.5,
        s=15,
        linewidths=0,
    )

    # thresholds
    ax.set_ylim(0, ylimAbs + 10)
    ax.set_xlim(-0.1 -xlimAbs, xlimAbs + 0.1)

    if straight:
        # vertical & horizontal straight lines
        ax.axvline(-vAss, color=low)
        ax.axvline(vAss, color=low)
        ax.axhline(hAss, color=high)
    else:
        # curved lines
        x_neg = np.linspace(-xlimAbs, -vAss, 3000)
        x_pos = np.linspace(vAss, xlimAbs, 3000)
        ax.plot(x_neg, xneg(x_neg), color="black", alpha=0.7, linewidth=0.5)
        ax.plot(x_pos, xpos(x_pos), color="black", alpha=0.7, linewidth=0.5)

    # labels if requested (no repel, just simple boxes)
    if label:
        labeled = df[~df["condLabel"].isna()]
        texts = []
        xs = labeled["logFC"].values
        ys = labeled["neg_log10_adjP"].values

        # Create text objects first
        for _, row in labeled.iterrows():
            if row['logFC'] < -1:
                color = col_high_neut_low[0]   # red (down)
            elif row['logFC'] > 1:
                color = col_high_neut_low[2]   # blue (up)
            else:
                color = col_high_neut_low[1]   # neutral

            t = ax.text(
                row["logFC"],
                row["neg_log10_adjP"],
                row["condLabel"],
                fontsize=10,
                color=color,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=color, lw=0.5),
            )
            texts.append(t)

        # Let adjust_text move the labels slightly and draw small arrows
        adjust_text(
            texts,
            x=xs,
            y=ys,
            ax=ax,
            arrowprops=dict(
                arrowstyle="-",
                color="black",
                lw=0.4,
            ),
            force_text=0.3,      # ← prevents large jumps
            force_points=0.1,
            expand_points=(1.5, 1.5),
            expand_text=(1.02, 1.02)
        )

    # minimal theme
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("logFC")
    ax.set_ylabel("-log10(adj. p-value)")
    ax.legend().set_visible(False)

    xmin, xmax = ax.get_xlim()
    plt.text(xmin, -0.8, labels[1], color=col_high_neut_low[0], fontsize=10, ha="left", va="top")
    plt.text(xmax, -0.8, labels[0], color=col_high_neut_low[2], fontsize=10, ha="right", va="top")


    return fig



def plot_dotplot(data,path,size_range=(0, 3), size_title="-log10(adj. pval)", color_title="log2(FC)", title= 'test'):
    plot_data = pd.DataFrame({
        'ID': list(data.index) * 2,
        'dataset': ['mouse'] * len(data) + ['human'] * len(data),
        'effect': pd.concat([data['logFC'], data['log2FoldChange']], ignore_index=True),
        'pval': pd.concat([-np.log10(data['adj.P.Val']), -np.log10(data['padj'])], ignore_index=True),
    })
    plot_data['ID'] = pd.Categorical(plot_data['ID'], categories=data.index.tolist()[::-1], ordered=True)
    print((plot_data.shape[0] * 0.5) + 1)
    plt.figure(figsize=(3, 5))
    color_norm = plt.Normalize(vmin=-1, vmax=1)
    size_range = (0, 3)  # fixed pval range for sizing
    size_norm = plt.Normalize(vmin=size_range[0], vmax=size_range[1])
    size_min, size_max = 10, 310
    ax = sns.scatterplot(
        data=plot_data,
        x='dataset',
        y='ID',
        hue='effect',
        size='pval',
        palette='coolwarm_r',
        hue_norm=color_norm,
        size_norm=size_norm,
        sizes=(size_min, size_max),
        legend=False,
    )

    ax.set_ylim([-0.5, 5.5])
    ax.set_xlim([-0.5, 1.5])

    sm = plt.cm.ScalarMappable(cmap='coolwarm_r', norm=color_norm)
    sm.set_array([])

    legend_vals = np.linspace(size_range[0], size_range[1], num=4)
    legend_sizes = [size_min + size_norm(val) * (size_max - size_min) for val in legend_vals]
    size_handles = [plt.scatter([], [], s=s, color='gray') for s in legend_sizes]
    size_labels = [f'{val:.1f}' for val in legend_vals]
    ax.legend(
        size_handles,
        size_labels,
        loc='center left',
        bbox_to_anchor=(1.8, 0.5),
        title=size_title,
    )
    sub_ax = plt.axes([1, 0.35, 0.1, 0.25])
    ax.figure.colorbar(sm, label=color_title, cax=sub_ax)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    labels = ax.get_xticklabels()
    ax.set_xticklabels(labels, rotation=45, ha='right', rotation_mode='anchor')
    ax.set_title(title)
    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.savefig(path, format="pdf", bbox_inches="tight")

def main(de_only: bool = False) -> None:
    # Load bulk counts and TPM data.
    bulk_data = pd.read_csv(
        DATA_DIR / "bulk/zhou/GSE179250_gene_reads.txt",
        sep="\t",
        dtype={"CHR": str},
        low_memory=False,
    )
    tpm_bulk_data = pd.read_csv(
        DATA_DIR / "bulk/zhou/GSE179250_gene_tpm_normalized.txt",
        sep="\t",
        dtype={"CHR": str},
        low_memory=False,
    )
    # Parse GEO metadata from the SOFT file and align to sample columns.
    meta_df = load_metadata(DATA_DIR / "bulk/zhou/GSE179250_family.soft.gz")
    counts = bulk_data.set_index("Description").iloc[:, 3:]
    counts_tpm = tpm_bulk_data.set_index("Description").iloc[:, 3:]

    # Map numeric column headers (SJLB IDs) to GEO samples via the title field (e.g., SJLB1 -> column "1").
    meta_df["title_num"] = meta_df["title"].str.extract(r"([0-9]+)")[0]
    sample_ids = counts.columns.astype(str)
    missing = set(sample_ids) - set(meta_df["title_num"])
    if missing:
        raise ValueError(f"Missing metadata for columns: {sorted(missing)}")
    meta_df = meta_df.set_index("title_num").loc[sample_ids].reset_index()
    meta_df["age"] = meta_df["age"].astype(int)
    meta_df["age_group"] = np.where(meta_df["age"] < 50, "below_50", "above_50")
    meta_df["age_sex"] = meta_df["age_group"].astype(str) + meta_df["sex"].astype(str)
    meta_df["sample_sex"] = meta_df["sample_id"].astype(str) + meta_df["sex"].astype(str)

    # Assemble AnnData and keep both raw counts and TPM-normalized layers.
    adata = prepare_adata(counts, counts_tpm, meta_df)

    # QC/visualization: PCA by sex.
    adata_pca = adata.copy()
    sc.pp.scale(adata_pca, max_value=10)
    sc.tl.pca(adata_pca, svd_solver='arpack')
    fig_pca = sc.pl.pca(adata_pca, color="sex", palette=GENDER_PALETTE, return_fig=True, show=False)
    save_fig(fig_pca, "pca_by_sex.pdf")

    if not de_only:
        from pypath.utils import homology

        signature = pd.read_csv(DATA_DIR / "hepa_gonad_signature.csv")
        homology_map = homology.get_dict(source="mouse", target="human", id_type="genesymbol")
        if not homology_map:
            raise RuntimeError("pypath returned an empty mouse-to-human orthology mapping")
        # Translate signature genes to human symbols and drop duplicates.
        signature = translate_signature(signature, homology_map, adata.var_names)

        # Run decoupler ULM to score the signature.
        estimate = run_ulm_scores(adata, signature)

        adata.write(RESULTS_DIR / "adata_processed.h5ad")
        estimate.to_csv(RESULTS_DIR / "ulm_sex_signature_estimates.csv", index=False)

        estimate["age_group_bins"] = pd.cut(
            estimate["age"],
            bins=np.arange(0, estimate["age"].max() + 15, 15),
            labels=None,
            include_lowest=True,
            right=False,
        )

        fig, ax = plt.subplots(figsize=(3, 4))
        young_cutoff, old_cutoff = 12, 49
        estimate_young_range = estimate[(estimate["age"] >= young_cutoff) & (estimate["age"] <= old_cutoff)]
        plot_gender_boxplot(ax, estimate_young_range, palette=GENDER_PALETTE)
        ax.set_title(f"ages {young_cutoff}-{old_cutoff}, n = {estimate_young_range.shape[0]}", fontsize = 16)
        save_fig(fig, "ulm_enrichment_12_49.pdf")

    # Differential expression: female vs male across all ages.
    results_df, mat = calc_deg(adata, "sex", "female", "male")
    mat.to_csv(RESULTS_DIR / "deg_zhou.csv")
    results_df.to_csv(RESULTS_DIR / "full_results_deg_zhou_malevsfemale.csv")

    plt.rcParams.update({'font.size': 11})

    volcano = prepare_volcano(results_df)
    print(volcano.head())
    significant = volcano[volcano["padj"] < 0.05]
    top3 = significant.nlargest(3, "log2FoldChange").index
    bottom3 = significant.nsmallest(3, "log2FoldChange").index
    genes = list(dict.fromkeys(list(top3) + list(bottom3)))
    print(genes)
    fig = volcano_nice(volcano.reset_index(),
                    hAss=0.05,
                    FCIndex="log2FoldChange",
                    pValIndex="padj",
                    IDIndex="Description",
                    manual_labels=genes,
                    straight=False,
                    label=True,
                    labels = ['male', 'female'])
    save_fig(fig, "volcano_male_vs_female_highlight.pdf")
    

    adata_sub = adata[(adata.obs["age"] >= 12) & (adata.obs["age"] <= 49)].copy()
    results_df_sub, mat_sub = calc_deg(adata_sub, "sex", "female", "male")
    mat_sub.to_csv(RESULTS_DIR / "12_49_malevsfemale_deg_zhou.csv")
    results_df_sub.to_csv(RESULTS_DIR / "12_49_malevsfemale_full_results_deg_zhou.csv")


    volcano = prepare_volcano(results_df_sub)
    print(volcano.head())
    significant = volcano[volcano["padj"] < 0.05]
    top3 = significant.nlargest(3, "log2FoldChange").index
    bottom3 = significant.nsmallest(3, "log2FoldChange").index
    genes = list(dict.fromkeys(list(top3) + list(bottom3)))
    print(genes)
    fig = volcano_nice(volcano.reset_index(),
                    hAss=0.05,
                    FCIndex="log2FoldChange",
                    pValIndex="padj",
                    IDIndex="Description",
                    manual_labels=genes,
                    straight=False,
                    label=True,
                    labels = ['male', 'female'])
    save_fig(fig, "volcano_male_vs_female_highlight_1249.pdf")


    if not de_only:
        full_results = pd.read_csv(DATA_DIR / "male_vs_female_DA_res.csv")
        full_results["human_gene_orth"] = [
            list(homology_map[g])[0] if g in homology_map else g.upper() for g in full_results["ID"]
        ]
        present_both = set(full_results["human_gene_orth"]) & set(adata.var_names)
        print(f"{len(full_results) - len(present_both)} genes cannot be translated with ortho & caps")

        all_results = full_results.merge(results_df_sub, left_on="human_gene_orth", right_index=True)
        all_results["rank_mouse_fem"] = all_results["t_Hepatocytes"].rank()
        all_results["rank_human_fem"] = all_results["stat"].rank()
        all_results["rank_mouse_male"] = all_results["t_Hepatocytes"].rank(ascending=False)
        all_results["rank_human_male"] = all_results["stat"].rank(ascending=False)
        all_results["avg_male_rank"] = (all_results["rank_mouse_male"] + all_results["rank_human_male"]) / 2
        all_results["avg_female_rank"] = (all_results["rank_mouse_fem"] + all_results["rank_human_fem"]) / 2
        all_results["con_malemouse_rank"] = (all_results["rank_mouse_male"] + all_results["rank_human_fem"]) / 2
        all_results["con_femalemouse_rank"] = (all_results["rank_mouse_fem"] + all_results["rank_human_male"]) / 2
        all_results = all_results.set_index("ID")


        sig = all_results[(all_results['adj.P.Val'] < 0.1) & (all_results['padj'] < 0.1)]
        # female
        data = sig[(sig['logFC'] < 0) & (sig['log2FoldChange'] < 0)].sort_values(by = 'avg_female_rank').iloc[:5]
        plot_dotplot(data, title = 'top female', path = RESULTS_DIR / 'top_ranked_genes_female_dot.pdf')
        # male
        data = sig[(sig['logFC'] > 0) & (sig['log2FoldChange'] > 0)].sort_values(by = 'avg_male_rank').iloc[:5]
        plot_dotplot(data, title = 'top male', path = RESULTS_DIR / 'top_ranked_genes_male_dot.pdf')

    adata_sub_female = adata[(adata.obs["age"] >= 12) & (adata.obs["sex"] == "female")].copy()
    adata_sub_female.obs["agecat"] = "young"
    adata_sub_female.obs.loc[adata_sub_female.obs["age"] > 49, "agecat"] = "old"

    # Young vs old contrast within females.
    results_df_yo, mat_yo = calc_deg(adata_sub_female, "agecat", "young", "old")
    results_df_yo.to_csv(RESULTS_DIR / "olderthan11female_oldvsyoung_full_results_deg_zhou.csv")
    mat_yo.to_csv(RESULTS_DIR / "olderthan11female_oldvsyoung_deg_statmatrix_zhou.csv")

    if not de_only:
        fig, ax = plt.subplots(figsize=(3, 4))
        estimate_sub = estimate[(estimate["age"] >= 12) & (estimate["sex"] == "female")].copy()
        estimate_sub["agecat"] = "young"
        estimate_sub.loc[estimate_sub["age"] > 49, "agecat"] = "old"
        plot_gender_boxplot(ax, estimate_sub, x_col="agecat", palette={"young": "#ED9FA3", "old": "#FF0000"})
        ax.set_title(f"female, n = {estimate_sub.shape[0]}", fontsize = 16)
        save_fig(fig, "ulm_enrichment_female__old_vs_young.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare mouse and human liver sex signatures.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--de-only", action="store_true",
                        help="Run GEO differential expression without mouse orthology or ULM scoring.")
    args = parser.parse_args()
    DATA_DIR = args.data_dir
    if args.output is not None:
        RESULTS_DIR = args.output
        if RESULTS_DIR.is_dir() and any(RESULTS_DIR.iterdir()):
            parser.error("--output requires an empty destination")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    main(de_only=args.de_only)
