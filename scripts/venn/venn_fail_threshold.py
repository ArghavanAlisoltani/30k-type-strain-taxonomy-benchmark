#!/usr/bin/env python3
"""
Four-way Venn diagram of species-boundary DISAGREEMENTS across 4 methods.

Input : exact_pair_status_table.csv
        columns: genome1, genome2, pair_id,
                 S16_Vsearch_status, FastANI_status,
                 Mash_status, FastAAI_Jaccard_status

For each method we take the pairs flagged FAIL_THRESHOLD (i.e. the tool ran
successfully but the score fell below the species-level cutoff) and compute
the *exact* pair-level intersections between the four methods.

Sanity-check totals (must match):
    16S              :   349
    Mash             : 1,297
    FastANI          : 2,026
    FastAAI/Jaccard  : 1,271

Outputs (written next to the script):
    venn_fail_threshold.png     high-resolution figure (300 dpi)
    venn_fail_threshold.pdf     vector version
    venn_region_counts.csv      all 15 non-empty region sizes
    venn_pair_regions.csv       every pair_id + which region it falls in
"""

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import venn as venn_lib   # `pip install venn`

# ---------------------------------------------------------------------------
# 1. Load data & build the four sets of FAIL_THRESHOLD pair_ids
# ---------------------------------------------------------------------------
CSV_PATH = Path("/mnt/user-data/uploads/exact_pair_status_table.csv")
OUT_DIR  = Path(__file__).resolve().parent

df = pd.read_csv(CSV_PATH)

METHODS = {                       # display name -> column in the CSV
    "16S":             "S16_Vsearch_status",
    "Mash":            "Mash_status",
    "FastANI":         "FastANI_status",
    "FastAAI/Jaccard": "FastAAI_Jaccard_status",
}

sets = {
    name: set(df.loc[df[col] == "FAIL_THRESHOLD", "pair_id"])
    for name, col in METHODS.items()
}

EXPECTED = {"16S": 349, "Mash": 1297, "FastANI": 2026, "FastAAI/Jaccard": 1271}
print("Sanity check (FAIL_THRESHOLD counts):")
all_ok = True
for m, s in sets.items():
    ok = len(s) == EXPECTED[m]
    all_ok &= ok
    tag = "OK" if ok else "MISMATCH"
    print(f"  {m:<18} {len(s):>5}  (expected {EXPECTED[m]:>5})  [{tag}]")
if not all_ok:
    raise SystemExit("Set sizes do not match expected values.")

# ---------------------------------------------------------------------------
# 2. Compute the 15 exclusive regions of the 4-set Venn (exact pair matching)
# ---------------------------------------------------------------------------
names   = list(sets.keys())
regions = {}                                        # frozenset(subset) -> set of pair_ids
for r in range(1, 5):
    for combo in combinations(names, r):
        included = [sets[n] for n in combo]
        excluded = [sets[n] for n in names if n not in combo]
        region   = set.intersection(*included).difference(*excluded)
        regions[frozenset(combo)] = region

union = set.union(*sets.values())
assert sum(len(v) for v in regions.values()) == len(union)

# Region-size table
rows = [{"region": " ∩ ".join(n for n in names if n in k),
         "n_methods": len(k),
         "size": len(v)}
        for k, v in regions.items()]
region_df = pd.DataFrame(rows).sort_values(["n_methods", "size"],
                                           ascending=[True, False])
region_df.to_csv(OUT_DIR / "venn_region_counts.csv", index=False)

# Per-pair table (which region each disagreeing pair belongs to)
pair_rows = []
for subset, pair_ids in regions.items():
    label = " ∩ ".join(n for n in names if n in subset)
    for pid in pair_ids:
        pair_rows.append({"pair_id": pid, "region": label})
pd.DataFrame(pair_rows).sort_values("pair_id")\
    .to_csv(OUT_DIR / "venn_pair_regions.csv", index=False)

print(f"\nUnion of all FAIL_THRESHOLD pairs: {len(union):,}")
print("\nRegion sizes:")
for _, r in region_df.iterrows():
    print(f"  {r['region']:<55} {r['size']:>5}")

# ---------------------------------------------------------------------------
# 3. Draw the 4-way Venn diagram
# ---------------------------------------------------------------------------
# Colourblind-friendly palette (Okabe–Ito-inspired), one colour per method.
colors = [
    (0.00, 0.45, 0.70, 0.45),   # blue        -> 16S
    (0.90, 0.62, 0.00, 0.45),   # orange      -> Mash
    (0.00, 0.62, 0.45, 0.45),   # bluish green-> FastANI
    (0.80, 0.40, 0.47, 0.45),   # reddish     -> FastAAI/Jaccard
]

# The venn library expects sets keyed by display name and returns a proper
# 4-ellipse layout with region labels centred inside each petal.
size_sets = {name: s for name, s in sets.items()}

fig, ax = plt.subplots(figsize=(11, 9), dpi=150)
venn_lib.venn(
    size_sets,
    ax=ax,
    cmap=[c[:3] for c in colors],   # library wants RGB, alpha applied below
    alpha=0.45,
    fontsize=12,
    legend_loc=None,                # no legend box — corner labels do it
)

# Corner labels: colored dot + method name + total, close to each ellipse.
# Positions in Axes coordinates (0-1). The venn library places the four
# ellipses centred in the axes, so these coordinates sit just outside the
# tips of each ellipse.
corner_labels = {
    "16S":             dict(x=0.16, y=0.90, ha="left",  va="center",
                            display="16S rRNA",       color=colors[0][:3]),
    "Mash":            dict(x=0.84, y=0.90, ha="right", va="center",
                            display="Mash",            color=colors[1][:3]),
    "FastANI":         dict(x=0.16, y=0.10, ha="left",  va="center",
                            display="FastANI",         color=colors[2][:3]),
    "FastAAI/Jaccard": dict(x=0.84, y=0.10, ha="right", va="center",
                            display="FastAAI/Jaccard", color=colors[3][:3]),
}
for key, cfg in corner_labels.items():
    # Small filled circle acting as the legend swatch
    ax.scatter([cfg["x"]], [cfg["y"]],
               s=220, c=[cfg["color"]], edgecolors="black", linewidths=0.8,
               transform=ax.transAxes, zorder=5, clip_on=False)
    # Method name + total, offset slightly from the dot
    dx = 0.018 if cfg["ha"] == "left" else -0.018
    ax.text(cfg["x"] + dx, cfg["y"],
            f"{cfg['display']}\nn = {len(sets[key]):,}",
            transform=ax.transAxes,
            ha=cfg["ha"], va=cfg["va"],
            fontsize=13, fontweight="bold", color="black")

fig.tight_layout()
fig.savefig(OUT_DIR / "venn_fail_threshold.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT_DIR / "venn_fail_threshold.pdf",           bbox_inches="tight")

print(f"\nWrote: {OUT_DIR/'venn_fail_threshold.png'}")
print(f"Wrote: {OUT_DIR/'venn_fail_threshold.pdf'}")
print(f"Wrote: {OUT_DIR/'venn_region_counts.csv'}")
print(f"Wrote: {OUT_DIR/'venn_pair_regions.csv'}")
