import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import AlignIO
import math
import pandas as pd

# ==========================================
# 1. CONFIGURATION
# ==========================================
ARCHEAE_ALN  = "aligned_16s.fna"      # Path to your MAFFT-aligned FASTA
WINDOW_SIZE   = 12                      # Sliding window for smoothing (bp)
MAX_GAP_FRAC  = 0.50                   # Drop columns where >50% are gaps
OUTPUT_CSV    = "archeae_conservation_data.csv"
OUTPUT_SVG    = "entropy.svg"
OUTPUT_PNG    = "entropy.png"

# Variable regions — approximate ARCHEAE 16S reference coordinates
V_REGIONS = [
    ( 57,  87, "V1"),
    (150, 205, "V2"),
    (262, 294, "V3"),
    (384, 440, "V4"),
    (760, 812, "V5"),
    (924, 985, "V6"),
    (1050, 1098, "V7"),
    (1189,1282,"V8"),
    (1353,1438,"V9"),
]

# ==========================================
# 2. SHANNON ENTROPY CALCULATION
# ==========================================
def calculate_information_content(alignment_file, max_gap_fraction=0.50):
    """
    Reads a FASTA alignment and calculates Shannon Information Content
    (bits) for every column, skipping gap-dominant columns.

    Returns
    -------
    scores       : np.array of information content values (0–2 bits)
    positions    : list of 1-based alignment column indices that were kept
    gap_fracs    : list of gap fractions for kept columns
    """
    print(f"\nProcessing: {alignment_file}")

    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except FileNotFoundError:
        print(f"  ERROR: File not found — {alignment_file}")
        return np.array([]), [], []

    num_sequences    = len(alignment)
    alignment_length = alignment.get_alignment_length()
    valid_chars      = set("ACGTU")

    information_scores = []
    kept_positions     = []
    gap_fracs          = []

    for i in range(alignment_length):
        column = alignment[:, i].upper()

        # ── Gap filter ────────────────────────────────────────────────────
        gap_count    = sum(1 for base in column if base not in valid_chars)
        gap_fraction = gap_count / num_sequences

        if gap_fraction > max_gap_fraction:
            continue                        # skip insertion-artifact columns

        kept_positions.append(i + 1)       # 1-based
        gap_fracs.append(round(gap_fraction, 4))

        # ── Nucleotide frequency ──────────────────────────────────────────
        counts      = {c: 0 for c in valid_chars}
        total_valid = 0

        for base in column:
            if base in valid_chars:
                counts[base] += 1
                total_valid  += 1

        if total_valid == 0:
            information_scores.append(0.0)
            continue

        # ── Shannon entropy: H = -Σ p·log2(p) ────────────────────────────
        entropy = 0.0
        for base in valid_chars:
            p = counts[base] / total_valid
            if p > 0:
                entropy -= p * math.log2(p)

        # Information content = max_entropy(2.0) − actual_entropy
        info = max(0.0, 2.0 - entropy)
        information_scores.append(info)

    print(f"  Sequences          : {num_sequences}")
    print(f"  Alignment columns  : {alignment_length}")
    print(f"  Kept after filter  : {len(information_scores)}  "
          f"(gap threshold ≤ {max_gap_fraction*100:.0f}%)")

    return np.array(information_scores), kept_positions, gap_fracs


# ==========================================
# 3. SMOOTHING HELPER
# ==========================================
def smooth(scores, window=WINDOW_SIZE):
    if len(scores) < window:
        return scores
    kernel = np.ones(window) / window
    return np.convolve(scores, kernel, mode='same')


# ==========================================
# 4. RUN CALCULATION
# ==========================================
bac_scores, bac_positions, bac_gap_fracs = calculate_information_content(
    ARCHEAE_ALN, max_gap_fraction=MAX_GAP_FRAC
)

if len(bac_scores) == 0:
    raise SystemExit("No data returned — check your alignment file path.")

smoothed_bac = smooth(bac_scores)

# ==========================================
# 5. EXPORT CSV
# ==========================================
df_export = pd.DataFrame({
    "Sequential_Position"      : np.arange(1, len(bac_scores) + 1),
    "Alignment_Position"       : bac_positions,
    "Gap_Fraction"             : bac_gap_fracs,
    "Raw_Information_Bits"     : bac_scores,
    "Smoothed_Information_Bits": smoothed_bac,
})
df_export.to_csv(OUTPUT_CSV, index=False)
print(f"\nCSV saved → {OUTPUT_CSV}  ({len(df_export)} rows)")

# ==========================================
# 6. PLOTTING
# ==========================================
def plot_domain(ax, scores, smoothed, color, label):
    """Plot smoothed line + filled area for one domain."""
    if len(scores) == 0:
        return
    x = np.arange(1, len(scores) + 1)
    ax.plot(x, smoothed, color=color, alpha=0.9, linewidth=1.2, label=label)
    ax.fill_between(x, smoothed, 0, color=color, alpha=0.15)


print("\nGenerating figure...")

COLOR_ARC   = "#1f77b4"   # blue  — Archaea
# COLOR_ARC = "#d62728"   # red   — Archaea  (uncomment when ready)
# COLOR_EUK = "#2ca02c"   # green — Eukarya  (uncomment when ready)
COLOR_VBAND = "#f5a623"   # amber — variable regions

fig, ax = plt.subplots(figsize=(16, 6))

# ── Variable-region shading ───────────────────────────────────────────────
n_positions = len(bac_scores)
for (start, end, label) in V_REGIONS:
    ax.axvspan(start, min(end, n_positions),
               color=COLOR_VBAND, alpha=0.18, zorder=0)
    mid = (start + min(end, n_positions)) / 2
    ax.text(mid, 2.08, label, ha='center', va='bottom',
            fontsize=7.5, color='#b07000', fontweight='bold')

# ── Entropy traces ────────────────────────────────────────────────────────
plot_domain(ax, bac_scores, smoothed_bac, COLOR_ARC, "Archaea (16S)")
# plot_domain(ax, arc_scores, smoothed_arc, COLOR_ARC, "Archaea (16S)")
# plot_domain(ax, euk_scores, smoothed_euk, COLOR_EUK, "Eukarya (18S)")

# ── Reference lines ───────────────────────────────────────────────────────
ax.axhline(2.0, color="grey", linewidth=0.7, linestyle="--", alpha=0.6)
ax.axhline(1.0, color="grey", linewidth=0.7, linestyle=":",  alpha=0.6)

# ── Axes labels & limits ──────────────────────────────────────────────────
ax.set_title(
    "16S rRNA Conservation Across Archeae\n"
    fr"Shannon Information Content (2 − H, bits)  •  {WINDOW_SIZE}-nt sliding window",
    fontsize=14, fontweight="bold", pad=12
)
ax.set_xlabel("Nucleotide Position (gap-filtered alignment)",
              fontsize=12, fontweight="bold")
ax.set_ylabel("Information Content (bits)", fontsize=12, fontweight="bold")
ax.set_xlim(1, n_positions)
ax.set_ylim(0, 2.20)

# ── Aesthetics ────────────────────────────────────────────────────────────
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", linestyle="--", alpha=0.4)
ax.tick_params(axis="both", labelsize=11)

# ── Legend ────────────────────────────────────────────────────────────────
vband_patch = mpatches.Patch(color=COLOR_VBAND, alpha=0.5,
                              label="Variable regions (V1–V9)")
ax.legend(handles=[vband_patch], loc="lower right",
          fontsize=10, frameon=True)

# ==========================================
# 7. SAVE
# ==========================================
plt.tight_layout()
plt.savefig(OUTPUT_SVG, format="svg", dpi=300)
plt.savefig(OUTPUT_PNG, format="png", dpi=200)
print(f"Figures saved → {OUTPUT_SVG}  /  {OUTPUT_PNG}")
plt.show()
