#!/usr/bin/env python3
"""Literature-threshold classification performance."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "data" / "species_discrimination_scores.tsv"
DEFAULT_OUT = ROOT / "outputs" / "statistical_analysis" / "literature_thresholds.tsv"

THRESHOLDS = {
    "16S rRNA": ("≥98.65% identity", 98.65, "higher"),
    "FastANI": ("≥95% ANI", 95.0, "higher"),
    "Mash": ("≤0.05 distance", 0.05, "lower"),
    "FastAAI": ("No established universal cutoff", None, None),
}


def eligible_rows(df):
    df = df[df["class"].isin(["SS", "SGDS"])].copy()
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score", "genus"])
    keep = (
        df.groupby("genus")["class"]
        .agg(lambda x: {"SS", "SGDS"}.issubset(set(x)))
    )
    return df[df["genus"].isin(keep[keep].index)].copy()


def ratio(a, b):
    return a / b if b else np.nan


def main():
    p = ArgumentParser()
    p.add_argument("input", nargs="?", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    required = {"method", "reference", "genus", "class", "score"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {', '.join(sorted(missing))}")

    rows = []
    references = list(dict.fromkeys(df["reference"].dropna().astype(str)))
    for reference in references:
        ref = df[df["reference"].astype(str) == reference]
        for method, (criterion, cutoff, direction) in THRESHOLDS.items():
            d = eligible_rows(ref[ref["method"] == method])
            if d.empty:
                continue

            ss = d["class"].eq("SS")
            sgds = d["class"].eq("SGDS")
            row = {
                "method": method,
                "reference": reference,
                "criterion": criterion,
                "eligible_ss": int(ss.sum()),
                "eligible_sgds": int(sgds.sum()),
            }

            if cutoff is None:
                row.update({k: np.nan for k in (
                    "TP", "FN", "TN", "FP", "sensitivity", "specificity",
                    "precision", "F1", "balanced_accuracy"
                )})
            else:
                predicted_positive = (
                    d["score"].ge(cutoff) if direction == "higher"
                    else d["score"].le(cutoff)
                )
                tp = int((ss & predicted_positive).sum())
                fn = int((ss & ~predicted_positive).sum())
                tn = int((sgds & ~predicted_positive).sum())
                fp = int((sgds & predicted_positive).sum())
                sensitivity = ratio(tp, tp + fn)
                specificity = ratio(tn, tn + fp)
                precision = ratio(tp, tp + fp)
                f1 = ratio(2 * precision * sensitivity, precision + sensitivity)
                row.update(
                    TP=tp, FN=fn, TN=tn, FP=fp,
                    sensitivity=sensitivity,
                    specificity=specificity,
                    precision=precision,
                    F1=f1,
                    balanced_accuracy=(sensitivity + specificity) / 2,
                )
            rows.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(
        args.out, sep="\t", index=False, float_format="%.4f", na_rep=""
    )


if __name__ == "__main__":
    main()
