#!/usr/bin/env python3
"""ROC and precision-recall analysis for species discrimination."""

from argparse import ArgumentParser
from pathlib import Path
import re

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "data" / "species_discrimination_scores.tsv"
DEFAULT_OUT = ROOT / "outputs" / "statistical_analysis"
METHODS = ("16S rRNA", "FastANI", "Mash", "FastAAI")


def slug(text):
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def eligible_rows(df):
    """Retain numeric SS/SGDS scores from genera containing both classes."""
    df = df[df["class"].isin(["SS", "SGDS"])].copy()
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score", "genus"])
    keep = (
        df.groupby("genus")["class"]
        .agg(lambda x: {"SS", "SGDS"}.issubset(set(x)))
    )
    return df[df["genus"].isin(keep[keep].index)].copy()


def main():
    p = ArgumentParser()
    p.add_argument("input", nargs="?", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--reference", default="NCBI")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    required = {"method", "reference", "genus", "class", "score"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {', '.join(sorted(missing))}")

    df = df[df["reference"].astype(str).str.upper() == args.reference.upper()]
    args.out.mkdir(parents=True, exist_ok=True)
    curve_dir = args.out / "curve_points"
    curve_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for method in METHODS:
        d = eligible_rows(df[df["method"] == method])
        if d.empty:
            continue

        y = (d["class"] == "SS").astype(int).to_numpy()
        score = d["score"].to_numpy(float)
        analysis_score = -score if method == "Mash" else score

        fpr, tpr, roc_thr = roc_curve(y, analysis_score)
        precision, recall, pr_thr = precision_recall_curve(y, analysis_score)

        pd.DataFrame(
            {"fpr": fpr, "tpr": tpr, "threshold": roc_thr}
        ).to_csv(curve_dir / f"{slug(method)}_roc.tsv", sep="\t", index=False)

        pd.DataFrame(
            {
                "recall": recall,
                "precision": precision,
                "threshold": np.r_[pr_thr, np.nan],
            }
        ).to_csv(curve_dir / f"{slug(method)}_pr.tsv", sep="\t", index=False)

        rows.append(
            {
                "method": method,
                "reference": args.reference.upper(),
                "eligible_ss": int((y == 1).sum()),
                "eligible_sgds": int((y == 0).sum()),
                "roc_auc": roc_auc_score(y, analysis_score),
                "average_precision": average_precision_score(y, analysis_score),
            }
        )

    if not rows:
        raise SystemExit("No eligible method/reference rows found.")

    pd.DataFrame(rows).to_csv(
        args.out / "roc_pr_summary.tsv", sep="\t", index=False, float_format="%.6f"
    )


if __name__ == "__main__":
    main()
