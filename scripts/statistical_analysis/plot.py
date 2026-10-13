#!/usr/bin/env python3
"""Plot the manuscript ROC and precision-recall panels."""

from argparse import ArgumentParser
from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_IN = ROOT / "outputs" / "statistical_analysis"
METHODS = ("16S rRNA", "FastANI", "Mash", "FastAAI")


def slug(text):
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def main():
    p = ArgumentParser()
    p.add_argument("--input-dir", type=Path, default=DEFAULT_IN)
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()

    out = args.out or args.input_dir / "roc_pr_curves"
    summary = pd.read_csv(args.input_dir / "roc_pr_summary.tsv", sep="\t").set_index("method")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for method in METHODS:
        if method not in summary.index:
            continue
        roc = pd.read_csv(args.input_dir / "curve_points" / f"{slug(method)}_roc.tsv", sep="\t")
        pr = pd.read_csv(args.input_dir / "curve_points" / f"{slug(method)}_pr.tsv", sep="\t")
        axes[0].plot(roc["fpr"], roc["tpr"], label=f"{method} (AUC={summary.loc[method, 'roc_auc']:.4f})")
        axes[1].plot(pr["recall"], pr["precision"], label=f"{method} (AP={summary.loc[method, 'average_precision']:.3f})")

    axes[0].plot([0, 1], [0, 1], linestyle="--", linewidth=1)
    axes[0].set(xlabel="False positive rate", ylabel="True positive rate", title="A. ROC")
    axes[1].set(xlabel="Recall", ylabel="Precision", title="B. Precision–recall")
    for ax in axes:
        ax.legend(fontsize=8)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)


if __name__ == "__main__":
    main()
