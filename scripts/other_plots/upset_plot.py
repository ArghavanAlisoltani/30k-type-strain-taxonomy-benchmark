import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Input:
# fail_16S, fail_FastANI, fail_Mash, n, degree, combination
df = pd.read_csv("upset_intersections.tsv", sep="\t")

# Exclude pairs that passed all three methods
df = df[df["degree"] > 0].copy()

# Largest intersections first
df = df.sort_values("n", ascending=False).reset_index(drop=True)

methods = [
    ("fail_16S", "16S rRNA"),
    ("fail_FastANI", "FastANI"),
    ("fail_Mash", "Mash")
]

x = np.arange(len(df))

fig, (ax_bar, ax_mat) = plt.subplots(
    2, 1,
    figsize=(8, 5),
    sharex=True,
    gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05}
)

# Intersection bars
bars = ax_bar.bar(x, df["n"], color="black")

ax_bar.set_ylabel("Intersection size")

# Bar-value labels
for bar, value in zip(bars, df["n"]):
    ax_bar.text(
        bar.get_x() + bar.get_width() / 2,
        value + max(2, value * 0.02),
        str(int(value)),
        ha="center",
        va="bottom",
        fontsize=9
    )

# UpSet membership matrix
for row, (column, label) in enumerate(methods):

    y = len(methods) - 1 - row

    for i, active in enumerate(df[column]):

        ax_mat.scatter(
            i,
            y,
            s=55,
            facecolor="black" if active else "white",
            edgecolor="black"
        )

    ax_mat.text(
        -0.7,
        y,
        label,
        ha="right",
        va="center"
    )

# Connect active dots
for i, row in df.iterrows():

    active_y = [
        len(methods) - 1 - j
        for j, (column, _) in enumerate(methods)
        if row[column]
    ]

    if len(active_y) > 1:
        ax_mat.plot(
            [i, i],
            [min(active_y), max(active_y)],
            color="black",
            linewidth=1.5
        )

ax_mat.set_ylim(-0.5, len(methods) - 0.5)
ax_mat.set_xlim(-0.5, len(df) - 0.5)

ax_mat.set_yticks([])
ax_mat.set_xticks([])

ax_mat.spines[
    ["top", "right", "bottom", "left"]
].set_visible(False)

fig.subplots_adjust(
    left=0.20,
    right=0.98,
    top=0.96,
    bottom=0.08
)

plt.savefig(
    "Figure4_UpSet.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "Figure4_UpSet.pdf",
    bbox_inches="tight"
)

plt.show()
