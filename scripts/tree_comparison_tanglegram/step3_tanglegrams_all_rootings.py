#!/usr/bin/env python3
"""
=============================================================================
Step 3 v4: tanglegrams for multiple rooting strategies
=============================================================================

v4 change: all colors now come from the manual definitions in
`tanglegram_colors.py` (no auto-assignment). The connecting lines and tip
labels are colored by the column named in tanglegram_colors.LINE_COLOR_BY
(default "domain"); edit that file to recolor.

Generates four standard tanglegrams:

  1. Archaea-outgroup-rooted
  2. Unrooted/input-root orientation preserved
  3. Selected-root using GCA_021654395.1
  4. Midpoint-rooted

The Archaea-outgroup-rooted version keeps the v2 behavior and output names
for compatibility. Additional rootings are added as extra outputs.

The Domain legend is placed outside the tip-label area to avoid overlap.
=============================================================================
"""
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import tree_lib
import tanglegram_colors as tc

TREE_LEFT_PATH = "output/16S_reduced.treefile"
TREE_RIGHT_PATH = "inputs/fastaai_nj.treefile"
METADATA = "output/metadata_harmonized.tsv"

OUTDIR = "output/tanglegrams"
DPI = 300

ROOTING_METHODS = [
    "archaea_outgroup",
    "unrooted",
    "selected_root_GCA_021654395_1",
    "midpoint",
]

SELECTED_ROOT_ACCESSION = "GCA_021654395.1"
OUTGROUP_DOMAIN_VALUES = {"Archaea", "Archea"}

FIG_SIZE = (17.5, 14)
LEFT_TITLE_BASE = "16S rRNA"
RIGHT_TITLE_BASE = "FastAAI NJ"
TITLE_FS = 12
TIP_FS = 7
TIP_PAD = 0.01

TREE_LW = 0.9
TREE_COLOR = "#222222"
LEFT_X_RANGE = (0.00, 0.40)
RIGHT_X_RANGE = (0.60, 1.00)
CENTER_GAP = (0.40, 0.60)

DOMAIN_COLORS = tc.DOMAIN_COLORS
MISSING_COLOR = tc.MISSING_COLOR
LINE_LW = 0.7
LINE_ALPHA = 0.55
UNTANGLE_PASSES = 4


def gca_from_tip(name):
    m = re.search(r"(GC[AF]_\d+\.\d+)", str(name))
    return m.group(1) if m else None


def root_with_archaea(tree, meta):
    id_col = meta.columns[0]
    gca_to_domain = dict(zip(meta[id_col], meta["domain"]))

    def is_archaea(tip_name):
        g = gca_from_tip(tip_name)
        return g is not None and gca_to_domain.get(g) in OUTGROUP_DOMAIN_VALUES

    return tree_lib.root_outgroup(tree, is_archaea)


def root_with_selected(tree, selected_accession):
    def is_selected(tip_name):
        return tip_name == selected_accession or gca_from_tip(tip_name) == selected_accession

    return tree_lib.root_outgroup(tree, is_selected)


def apply_rooting(tree, method, meta):
    if method == "archaea_outgroup":
        return root_with_archaea(tree, meta), "Archaea-outgroup-rooted"
    if method == "unrooted":
        return tree, "unrooted / input-root orientation"
    if method.startswith("selected_root_"):
        return root_with_selected(tree, SELECTED_ROOT_ACCESSION), f"rooted on {SELECTED_ROOT_ACCESSION}"
    if method == "midpoint":
        return tree_lib.root_midpoint(tree), "midpoint-rooted"
    raise ValueError(method)


def cladogram_y_positions(tree, leaf_order=None):
    terminals = tree.get_terminals()
    if leaf_order is None:
        leaf_order = [t.name for t in terminals]
    name_to_y = {name: i for i, name in enumerate(leaf_order)}
    y = {}
    for t in terminals:
        y[t] = name_to_y[t.name]
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            continue
        ys = [y[c] for c in clade.clades]
        y[clade] = sum(ys) / len(ys)
    return y, leaf_order


def x_positions(tree):
    return tree.depths()


def reorder_leaves_min_crossings(tree, target_y):
    def reorder(node):
        if node.is_terminal():
            return [node.name]

        child_orders = [reorder(c) for c in node.clades]
        means = []
        for order in child_orders:
            ys = [target_y[n] for n in order if n in target_y]
            means.append(sum(ys) / len(ys) if ys else 0)

        idx_sorted = sorted(range(len(node.clades)), key=lambda i: means[i])
        node.clades = [node.clades[i] for i in idx_sorted]

        merged = []
        for i in idx_sorted:
            merged.extend(child_orders[i])
        return merged

    return reorder(tree)


def count_crossings(leaf_order_left, leaf_order_right):
    pos_right = {n: i for i, n in enumerate(leaf_order_right)}
    sequence = [pos_right.get(n, -1) for n in leaf_order_left if n in pos_right]
    n = len(sequence)
    cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            if sequence[i] > sequence[j]:
                cnt += 1
    return cnt


def draw_tree(ax, tree, y_coords, x_coords, x_min, x_max, mirror=False, color=TREE_COLOR, lw=TREE_LW):
    max_d = max(x_coords.values()) or 1.0

    def x_of(clade):
        d = x_coords[clade] / max_d if max_d > 0 else 0
        if mirror:
            return x_max - d * (x_max - x_min)
        return x_min + d * (x_max - x_min)

    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        px = x_of(clade)
        ys = [y_coords[c] for c in clade.clades]
        ax.plot([px, px], [min(ys), max(ys)], color=color, lw=lw, solid_capstyle="butt")
        for ch in clade.clades:
            cy = y_coords[ch]
            cx = x_of(ch)
            ax.plot([px, cx], [cy, cy], color=color, lw=lw, solid_capstyle="butt")


def make_tanglegram(method):
    os.makedirs(OUTDIR, exist_ok=True)

    meta = pd.read_csv(METADATA, sep="\t")
    id_col = meta.columns[0]

    # v4: lines/labels are colored by whichever column LINE_COLOR_BY names,
    # using the manual color map for that column.
    color_col = tc.LINE_COLOR_BY
    color_mapping, _ = tc.color_map_for_column(color_col)
    value_of = dict(zip(meta[id_col], meta[color_col]))

    t_left = tree_lib.parse_newick(TREE_LEFT_PATH)
    t_right = tree_lib.parse_newick(TREE_RIGHT_PATH)

    t_left, root_label_left = apply_rooting(t_left, method, meta)
    t_right, root_label_right = apply_rooting(t_right, method, meta)

    leaves_l = {t.name for t in t_left.get_terminals()}
    leaves_r = {t.name for t in t_right.get_terminals()}
    common = leaves_l & leaves_r

    print(f"[{method}] common leaves: {len(common)}")

    # 16S (left) is the REFERENCE tree: its leaf order is fixed to the natural
    # post-rooting order (the same order the circular branch-length tree uses).
    # We never rotate the 16S tree, so its structure/order stay untouched and
    # the Archaea block remains contiguous, matching the circular tree. Only the
    # FastAAI (right) tree is rotated to minimize crossings.
    order_l = [t.name for t in t_left.get_terminals() if t.name in common]
    order_r = [t.name for t in t_right.get_terminals() if t.name in common]
    cross0 = count_crossings(order_l, order_r)

    target_y_for_right = {n: i for i, n in enumerate(order_l)}
    for k in range(UNTANGLE_PASSES):
        order_r = reorder_leaves_min_crossings(t_right, target_y_for_right)
        order_r = [n for n in order_r if n in common]

    final_cross = count_crossings(order_l, order_r)

    y_l, _ = cladogram_y_positions(t_left, order_l)
    y_r, _ = cladogram_y_positions(t_right, order_r)
    x_l = x_positions(t_left)
    x_r = x_positions(t_right)

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    n = len(common)
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, n)
    ax.axis("off")

    draw_tree(ax, t_left, y_l, x_l, LEFT_X_RANGE[0], LEFT_X_RANGE[1], mirror=False)
    draw_tree(ax, t_right, y_r, x_r, RIGHT_X_RANGE[0], RIGHT_X_RANGE[1], mirror=True)

    name_to_left_terminal = {t.name: t for t in t_left.get_terminals()}
    name_to_right_terminal = {t.name: t for t in t_right.get_terminals()}

    max_dl = max(x_l.values()) or 1.0
    max_dr = max(x_r.values()) or 1.0

    for name in common:
        yl = y_l[name_to_left_terminal[name]]
        yr = y_r[name_to_right_terminal[name]]

        x_left_actual = LEFT_X_RANGE[0] + (x_l[name_to_left_terminal[name]] / max_dl) * (LEFT_X_RANGE[1] - LEFT_X_RANGE[0])
        x_right_actual = RIGHT_X_RANGE[1] - (x_r[name_to_right_terminal[name]] / max_dr) * (RIGHT_X_RANGE[1] - RIGHT_X_RANGE[0])

        dom = value_of.get(name)
        color = tc.resolve(dom, color_mapping)

        ax.plot(
            [x_left_actual, CENTER_GAP[0], CENTER_GAP[1], x_right_actual],
            [yl, yl, yr, yr],
            color=color,
            lw=LINE_LW,
            alpha=LINE_ALPHA,
            solid_capstyle="round",
        )

        ax.text(CENTER_GAP[0] - TIP_PAD, yl, name, ha="right", va="center", fontsize=TIP_FS, color=tc.TIP_LABEL_COLOR)
        ax.text(CENTER_GAP[1] + TIP_PAD, yr, name, ha="left", va="center", fontsize=TIP_FS, color=tc.TIP_LABEL_COLOR)

    pretty = {
        "archaea_outgroup": "Archaea-outgroup-rooted",
        "unrooted": "unrooted / input-root orientation",
        f"selected_root_{SELECTED_ROOT_ACCESSION.replace('.', '_')}": f"rooted on {SELECTED_ROOT_ACCESSION}",
        "midpoint": "midpoint-rooted",
    }[method]

    fig.suptitle(f"Tanglegram: 16S rRNA vs FastAAI — {pretty}", fontsize=TITLE_FS + 2, weight="bold", y=0.99)
    ax.text((LEFT_X_RANGE[0] + LEFT_X_RANGE[1]) / 2, n + 0.5, f"{LEFT_TITLE_BASE}\n{root_label_left}", ha="center", va="bottom", fontsize=TITLE_FS, weight="bold")
    ax.text((RIGHT_X_RANGE[0] + RIGHT_X_RANGE[1]) / 2, n + 0.5, f"{RIGHT_TITLE_BASE}\n{root_label_right}", ha="center", va="bottom", fontsize=TITLE_FS, weight="bold")

    # Legend lists only the categories actually present, in declaration order,
    # using the manual colors. Title reflects the coloring column.
    legend_map = tc.present_color_map(
        [value_of.get(name) for name in common], color_mapping
    )
    handles = [mpatches.Patch(color=c, label=k) for k, c in legend_map.items()]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.015, 0.965),
        fontsize=9,
        frameon=False,
        title=color_col.capitalize(),
        title_fontsize=10,
    )

    ax.text(
        0.5,
        -0.5,
        f"line crossings before/after untangling: {cross0} → {final_cross}",
        ha="center",
        va="top",
        fontsize=9,
        color="#555555",
        style="italic",
    )

    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.03)

    out_png = f"{OUTDIR}/tanglegram_{method}.png"
    out_pdf = f"{OUTDIR}/tanglegram_{method}.pdf"
    fig.savefig(out_png, dpi=DPI, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)

    # v2-compatible copy names for Archaea-rooted tanglegram.
    if method == "archaea_outgroup":
        fig_src_png = out_png
        fig_src_pdf = out_pdf
        import shutil
        shutil.copy2(fig_src_png, "output/tanglegram_archaea_rooted.png")
        shutil.copy2(fig_src_pdf, "output/tanglegram_archaea_rooted.pdf")

    return {"rooting_method": method, "crossings_before": cross0, "crossings_after": final_cross}


def main():
    rows = []
    for method in ROOTING_METHODS:
        rows.append(make_tanglegram(method))
    pd.DataFrame(rows).to_csv("output/tanglegram_crossing_summary.tsv", sep="\t", index=False)
    print("[done] wrote all tanglegrams")


if __name__ == "__main__":
    main()
