import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(
    "CompleteGenome_per_GCA_full_length16S.tsv",
    sep="\t"
)

counts = (
    df.groupby(["Full-length 16S count", "Domain"])
      .size()
      .unstack(fill_value=0)
)

# Include all copy-number positions, even those with zero genomes
counts = counts.reindex(
    range(0, int(df["Full-length 16S count"].max()) + 1),
    fill_value=0
)

x = np.arange(len(counts))
width = 0.42

# Broken y-axis
fig, (ax_top, ax_bottom) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(14, 6.5),
    gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05}
)

for ax in (ax_top, ax_bottom):

    ax.bar(
        x - width/2,
        counts["Bacteria"],
        width,
        color="black",
        label="Bacteria"
    )

    ax.bar(
        x + width/2,
        counts["Archaea"],
        width,
        color="orange",
        label="Archaea"
    )

# Break at 15
ax_bottom.set_ylim(0, 15)
ax_top.set_ylim(15, counts.to_numpy().max() * 1.10)

ax_top.spines["bottom"].set_visible(False)
ax_bottom.spines["top"].set_visible(False)

ax_top.tick_params(
    axis="x",
    bottom=False,
    labelbottom=False
)

# Small diagonal break marks
d = 0.006

kwargs = dict(
    transform=ax_top.transAxes,
    clip_on=False,
    color="black"
)

ax_top.plot((-d, +d), (-d, +d), **kwargs)
ax_top.plot((1-d, 1+d), (-d, +d), **kwargs)

kwargs.update(transform=ax_bottom.transAxes)

ax_bottom.plot((-d, +d), (1-d, 1+d), **kwargs)
ax_bottom.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

# Add value to every non-zero bar
for offset, domain in [
    (-width/2, "Bacteria"),
    ( width/2, "Archaea")
]:

    for i, value in enumerate(counts[domain]):

        if value == 0:
            continue

        xpos = x[i] + offset

        if value >= 15:
            ax_top.text(
                xpos,
                value + max(4, value * 0.01),
                str(int(value)),
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=8
            )

        else:
            ax_bottom.text(
                xpos,
                value + 0.35,
                str(int(value)),
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=8
            )

ax_bottom.set_xticks(x)
ax_bottom.set_xticklabels(
    counts.index.astype(str),
    rotation=90
)

ax_bottom.set_xlabel(
    "Full-length 16S rRNA copy number"
)

fig.supylabel(
    "Number of genomes"
)

handles, labels = ax_top.get_legend_handles_labels()

ax_top.legend(
    handles[:2],
    labels[:2],
    frameon=False,
    loc="upper right"
)

fig.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.15,
    top=0.98
)

plt.savefig(
    "PanelB_CompleteGenome_only.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "PanelB_CompleteGenome_only.pdf",
    bbox_inches="tight"
)

plt.show()
