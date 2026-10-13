#!/usr/bin/env python3
"""
=============================================================================
Step 5 v4: visual summaries for tree-comparison metrics
=============================================================================

v4 change: the kingdom-monophyly heatmap status colors are now taken from
tanglegram_colors.STATUS_COLORS, and kingdom y-axis labels are tinted with
the manual KINGDOM_COLORS so the heatmap matches the strip tanglegrams.

Creates additional figures requested for comparing RF and other metrics.

Outputs
-------
output/figures/rooting_metric_comparison.png/pdf
output/figures/mantel_correlation_by_rooting.png/pdf
output/figures/monophyly_kingdom_heatmap.png/pdf
output/figures/tanglegram_crossings_by_rooting.png/pdf
=============================================================================
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import tanglegram_colors as tc

METRICS = "output/tree_comparison_metrics_by_rooting.tsv"
KINGDOM_MONO = "output/per_kingdom_monophyly_by_rooting.tsv"
CROSSINGS = "output/tanglegram_crossing_summary.tsv"
OUTDIR = "output/figures"

ROOT_ORDER = [
    "archaea_outgroup",
    "unrooted",
    "selected_root_GCA_021654395_1",
    "midpoint",
]

ROOT_LABELS = {
    "archaea_outgroup": "Archaea\noutgroup",
    "unrooted": "Unrooted /\ninput root",
    "selected_root_GCA_021654395_1": "Selected root\nGCA_021654395.1",
    "midpoint": "Midpoint",
}


def ensure_outdir():
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs("output/tables", exist_ok=True)


def plot_metric_comparison(metrics):
    """Compare topology metrics across rooting strategies."""
    df = metrics.set_index("rooting_method").loc[ROOT_ORDER].reset_index()

    x = np.arange(len(df))
    width = 0.22

    fig, ax = plt.subplots(figsize=(10.5, 5.5))

    ax.bar(x - width, df["normalised_RF"], width, label="Normalized RF")
    ax.bar(x, df["normalised_rooted_cluster_distance"], width, label="Normalized rooted cluster distance")
    ax.bar(x + width, df["normalised_quartet_distance_all"], width, label="Normalized quartet distance")

    ax.set_xticks(x)
    ax.set_xticklabels([ROOT_LABELS[m] for m in df["rooting_method"]])
    ax.set_ylabel("Normalized distance")
    ax.set_ylim(0, 1.05)
    ax.set_title("Tree topology disagreement across rooting choices")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=0.25)

    # Add RF labels above first bar group.
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(i - width, row["normalised_RF"] + 0.025, f"RF={int(row['RF_distance_unrooted'])}", ha="center", va="bottom", fontsize=8, rotation=90)

    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/rooting_metric_comparison.png", dpi=300, facecolor="white")
    fig.savefig(f"{OUTDIR}/rooting_metric_comparison.pdf", facecolor="white")
    plt.close(fig)


def plot_mantel(metrics):
    """Plot Mantel Pearson and Spearman correlations across rooting strategies."""
    df = metrics.set_index("rooting_method").loc[ROOT_ORDER].reset_index()
    x = np.arange(len(df))
    width = 0.32

    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.bar(x - width/2, df["mantel_pearson_r"], width, label="Pearson r")
    ax.bar(x + width/2, df["mantel_spearman_rho"], width, label="Spearman ρ")

    ax.set_xticks(x)
    ax.set_xticklabels([ROOT_LABELS[m] for m in df["rooting_method"]])
    ax.set_ylabel("Mantel correlation")
    ax.set_ylim(0, 1.05)
    ax.set_title("Branch-length-aware agreement: patristic distance correlations")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="y", alpha=0.25)

    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(i, 0.04, f"p={row['mantel_pearson_p']:.3f}", ha="center", va="bottom", fontsize=8, color="#555555")

    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/mantel_correlation_by_rooting.png", dpi=300, facecolor="white")
    fig.savefig(f"{OUTDIR}/mantel_correlation_by_rooting.pdf", facecolor="white")
    plt.close(fig)


def plot_monophyly_heatmap():
    """Create a simple kingdom monophyly heatmap for root-sensitive interpretation."""
    mono = pd.read_csv(KINGDOM_MONO, sep="\t")
    multi = mono[
        (~mono["monophyletic_16S"].astype(str).str.startswith("N/A")) &
        (~mono["group"].astype(str).isin(["Unknown", "nan", "NaN", "None", ""]))
    ].copy()

    # Status codes:
    # 0 = neither
    # 1 = 16S only
    # 2 = FastAAI only
    # 3 = both
    def code(row):
        a = row["monophyletic_16S"] == "yes"
        b = row["monophyletic_FastAAI"] == "yes"
        if a and b:
            return 3
        if a:
            return 1
        if b:
            return 2
        return 0

    multi["group"] = multi["group"].astype(str)
    multi["status_code"] = multi.apply(code, axis=1)

    groups = sorted(multi["group"].dropna().astype(str).unique())
    methods = ROOT_ORDER

    mat = np.full((len(groups), len(methods)), np.nan)
    for i, g in enumerate(groups):
        for j, m in enumerate(methods):
            hit = multi[(multi["group"] == g) & (multi["rooting_method"] == m)]
            if len(hit):
                mat[i, j] = hit.iloc[0]["status_code"]

    # Custom colors via ListedColormap, ordered to match status codes 0..3
    # (neither, 16S only, FastAAI only, both). Defined manually in tanglegram_colors.
    from matplotlib.colors import ListedColormap, BoundaryNorm
    status_order = ["neither", "16S", "FastAAI", "both"]
    cmap = ListedColormap([tc.STATUS_COLORS[s] for s in status_order])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    fig, ax = plt.subplots(figsize=(9.5, max(4, 0.45 * len(groups))))
    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(np.arange(len(methods)))
    ax.set_xticklabels([ROOT_LABELS[m] for m in methods])
    ax.set_yticks(np.arange(len(groups)))
    ax.set_yticklabels(groups)
    # Tint each kingdom label with its manual color so this figure matches
    # the kingdom strip in the step4 tanglegrams.
    for lab in ax.get_yticklabels():
        lab.set_color(tc.resolve(lab.get_text(), tc.KINGDOM_COLORS))
        lab.set_fontweight("bold")
    ax.set_title("Kingdom monophyly by rooting method")
    ax.set_xlabel("Rooting method")
    ax.set_ylabel("Kingdom")

    # Annotate cells
    labels = {0: "neither", 1: "16S", 2: "AAI", 3: "both"}
    for i in range(len(groups)):
        for j in range(len(methods)):
            val = mat[i, j]
            if not np.isnan(val):
                ax.text(j, i, labels[int(val)], ha="center", va="center", fontsize=8)

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=tc.STATUS_COLORS["neither"], label="Neither tree"),
        plt.Rectangle((0, 0), 1, 1, color=tc.STATUS_COLORS["16S"], label="16S only"),
        plt.Rectangle((0, 0), 1, 1, color=tc.STATUS_COLORS["FastAAI"], label="FastAAI only"),
        plt.Rectangle((0, 0), 1, 1, color=tc.STATUS_COLORS["both"], label="Both trees"),
    ]
    ax.legend(handles=handles, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/monophyly_kingdom_heatmap.png", dpi=300, facecolor="white")
    fig.savefig(f"{OUTDIR}/monophyly_kingdom_heatmap.pdf", facecolor="white")
    plt.close(fig)


def plot_crossings():
    if not os.path.exists(CROSSINGS):
        return

    df = pd.read_csv(CROSSINGS, sep="\t").set_index("rooting_method").loc[ROOT_ORDER].reset_index()
    x = np.arange(len(df))
    width = 0.32

    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.bar(x - width/2, df["crossings_before"], width, label="Before untangling")
    ax.bar(x + width/2, df["crossings_after"], width, label="After untangling")

    ax.set_xticks(x)
    ax.set_xticklabels([ROOT_LABELS[m] for m in df["rooting_method"]])
    ax.set_ylabel("Line crossings")
    ax.set_title("Tanglegram line crossings across rooting choices")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/tanglegram_crossings_by_rooting.png", dpi=300, facecolor="white")
    fig.savefig(f"{OUTDIR}/tanglegram_crossings_by_rooting.pdf", facecolor="white")
    plt.close(fig)


def write_interpretation_notes(metrics):
    rows = [
        {
            "topic": "RF distance",
            "interpretation": "Standard RF is unrooted; it should generally stay constant across rerooting because rerooting does not change the unrooted topology.",
        },
        {
            "topic": "Rooted cluster distance",
            "interpretation": "This is root-sensitive and is the better number for asking how clade interpretation changes after Archaea, selected-sample, midpoint, or input-root rooting.",
        },
        {
            "topic": "Quartet distance",
            "interpretation": "Quartet distance is a robust topology summary based on four-taxon relationships and is less sensitive to single-tip perturbations than RF.",
        },
        {
            "topic": "Mantel patristic correlation",
            "interpretation": "Mantel correlation compares pairwise tree path-length matrices; high values indicate similar global distance structure even if exact branching differs.",
        },
        {
            "topic": "Monophyly heatmap",
            "interpretation": "The monophyly heatmap shows whether each kingdom is monophyletic in neither tree, only 16S, only FastAAI, or both under each rooting.",
        },
    ]
    pd.DataFrame(rows).to_csv("output/tables/metric_interpretation_notes.tsv", sep="\t", index=False)


def main():
    ensure_outdir()
    metrics = pd.read_csv(METRICS, sep="\t")

    plot_metric_comparison(metrics)
    plot_mantel(metrics)
    plot_monophyly_heatmap()
    plot_crossings()
    write_interpretation_notes(metrics)

    print("[done] wrote metric comparison figures to output/figures")


if __name__ == "__main__":
    main()
