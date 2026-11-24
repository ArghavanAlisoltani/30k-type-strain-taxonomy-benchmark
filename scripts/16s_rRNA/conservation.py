#!/usr/bin/env python3
"""Site-wise 16S conservation from an aligned FASTA, following the supplied entropy plots."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from Bio import AlignIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def compute(alignment, max_missing=0.50, window=12, normalize_u=False):
    pos, fraction, info = [], [], []
    num = len(alignment)
    if num == 0:
        raise ValueError("Empty alignment")
    for i in range(alignment.get_alignment_length()):
        # Default matches supplied source scripts (T and U counted separately).
        column = alignment[:, i].upper()
        if normalize_u:
            column = column.replace("U", "T")
        alphabet = "ACGT" if normalize_u else "ACGTU"
        bases = [b for b in column if b in alphabet]
        missing = 1 - len(bases)/num
        if missing > max_missing:
            continue
        freq = np.array([bases.count(b) for b in alphabet], dtype=float)
        if freq.sum() == 0:
            continue
        probs = freq[freq > 0] / freq.sum()
        entropy = -float(np.sum(probs * np.log2(probs)))
        pos.append(i + 1); fraction.append(missing); info.append(max(0, 2 - entropy))
    if not info:
        raise ValueError("All columns excluded")
    raw = np.asarray(info)
    smoothed = np.convolve(raw, np.ones(window)/window, mode="same") if len(raw) >= window else raw.copy()
    return pos, fraction, raw, smoothed


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--alignment", type=Path, required=True)
    p.add_argument("--output-prefix", type=Path, required=True)
    p.add_argument("--max-missing", type=float, default=0.50)
    p.add_argument("--window", type=int, default=12)
    p.add_argument("--normalize-u", action="store_true", help="Optional four-state analysis (not historical default)")
    args = p.parse_args()
    if not (0 <= args.max_missing < 1) or args.window < 1:
        p.error("Invalid filtering/smoothing parameters")
    aln = AlignIO.read(str(args.alignment), "fasta")
    pos, gaps, raw, smooth = compute(aln, args.max_missing, args.window, args.normalize_u)
    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if prefix.with_suffix(".tsv").exists() or prefix.with_suffix(".pdf").exists():
        p.error("Output exists; refusing to overwrite")
    with prefix.with_suffix(".tsv").open("x", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["alignment_position", "missing_fraction", "information_2_minus_H_bits", "smoothed_information_bits"])
        w.writerows(zip(pos, gaps, raw, smooth))
    fig, ax = plt.subplots(figsize=(11, 3))
    ax.plot(np.arange(1, len(raw) + 1), smooth, linewidth=0.75)
    ax.set(xlabel="Gap-filtered alignment column", ylabel="Information (2 − H, bits)", ylim=(0, 2.05), xlim=(1, len(raw)))
    fig.tight_layout()
    fig.savefig(prefix.with_suffix(".pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
