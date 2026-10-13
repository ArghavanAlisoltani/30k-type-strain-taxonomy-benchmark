
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import to_rgba

df = pd.read_csv("upset_intersections.tsv", sep="\t")

def get_n(a, b, c):
    row = df[
        (df["fail_16S"] == a) &
        (df["fail_FastANI"] == b) &
        (df["fail_Mash"] == c)
    ]
    return int(row["n"].iloc[0])

# Mutually exclusive regions
a_only  = get_n(True,  False, False)
b_only  = get_n(False, True,  False)
c_only  = get_n(False, False, True)
ab_only = get_n(True,  True,  False)
ac_only = get_n(True,  False, True)
bc_only = get_n(False, True,  True)
abc     = get_n(True,  True,  True)

# Set totals
a_total = a_only + ab_only + ac_only + abc
b_total = b_only + ab_only + bc_only + abc
c_total = c_only + ac_only + bc_only + abc

# Same blue / orange / green palette used in the original Venn
colors = {
    "16S": "#1f77b4",
    "FastANI": "#ff7f0e",
    "Mash": "#2ca02c",
}

fig, ax = plt.subplots(figsize=(7, 6))

r = 1.20
centers = {
    "16S": (-0.65,  0.38),
    "FastANI": (0.65, 0.38),
    "Mash": (0.00, -0.75),
}

for name, center in centers.items():
    ax.add_patch(
        Circle(
            center,
            r,
            facecolor=to_rgba(colors[name], 0.35),
            edgecolor=colors[name],
            linewidth=1.5
        )
    )

# Set labels
ax.text(-1.45, 1.55, f"16S rRNA\n(n={a_total})",
        ha="center", va="center")
ax.text( 1.45, 1.55, f"FastANI\n(n={b_total})",
        ha="center", va="center")
ax.text( 0.00,-2.05, f"Mash\n(n={c_total})",
        ha="center", va="center")

# Mutually exclusive region counts
ax.text(-1.30,  0.55, str(a_only),  ha="center", va="center", fontsize=12)
ax.text( 1.30,  0.55, str(b_only),  ha="center", va="center", fontsize=12)
ax.text( 0.00, -1.45, str(c_only),  ha="center", va="center", fontsize=12)

ax.text( 0.00,  1.00, str(ab_only), ha="center", va="center", fontsize=12)
ax.text(-0.62, -0.42, str(ac_only), ha="center", va="center", fontsize=12)
ax.text( 0.62, -0.42, str(bc_only), ha="center", va="center", fontsize=12)
ax.text( 0.00,  0.02, str(abc),     ha="center", va="center",
        fontsize=12, fontweight="bold")

ax.set_aspect("equal")
ax.set_xlim(-2.0, 2.0)
ax.set_ylim(-2.25, 1.95)
ax.axis("off")

plt.savefig("Figure4_Venn_symmetric_colored.png",
            dpi=300, bbox_inches="tight")
plt.savefig("Figure4_Venn_symmetric_colored.pdf",
            bbox_inches="tight")
plt.close()

print("16S total:", a_total)
print("FastANI total:", b_total)
print("Mash total:", c_total)
print("Discordant union:",
      a_only + b_only + c_only + ab_only + ac_only + bc_only + abc)
