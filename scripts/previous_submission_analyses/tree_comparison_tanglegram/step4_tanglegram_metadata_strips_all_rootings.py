#!/usr/bin/env python3
"""
=============================================================================
Step 4 v4: tanglegrams with Domain/Kingdom/Phylum metadata strips
=============================================================================

v4 change: the Domain, Kingdom AND Phylum strip colors now come from the
manual name -> hex maps in `tanglegram_colors.py`. v3 auto-assigned kingdom
and phylum colors by alphabetical index from a palette; that auto-assignment
is removed. Edit tanglegram_colors.py to recolor any strip.

Generates metadata-strip tanglegrams for the same four rooting strategies:

  1. Archaea-outgroup-rooted
  2. Unrooted/input-root orientation preserved
  3. Selected-root using GCA_021654395.1
  4. Midpoint-rooted

Each figure shows the tanglegram plus four vertical metadata strips beside
each tree:
  D = Domain
  K = Kingdom
  P = Phylum (official_phylum)
  W = Carl Woese historical group (carl_woese_historical_group)
=============================================================================
"""
import os
import re
from collections import OrderedDict

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

import tree_lib
import tanglegram_colors as tc

TREE_LEFT_PATH = "output/16S_reduced.treefile"
TREE_RIGHT_PATH = "inputs/fastaai_nj.treefile"
METADATA = "output/metadata_harmonized.tsv"
OUTDIR = "output/tanglegrams_metadata_strips"

ROOTING_METHODS = [
    "archaea_outgroup",
    "unrooted",
    "selected_root_GCA_021654395_1",
    "midpoint",
]
SELECTED_ROOT_ACCESSION = "GCA_021654395.1"
OUTGROUP_DOMAIN_VALUES = {"Archaea", "Archea"}

FIG_SIZE = (21, 14)
DPI = 300
TIP_FS = 6.4
TREE_LW = 0.85
TREE_COLOR = "#222222"

LEFT_X_RANGE = (0.00, 0.31)
RIGHT_X_RANGE = (0.69, 1.00)
CENTER_GAP = (0.43, 0.57)
TIP_PAD = 0.006

# Three strips: D, K, P
# Four strips: D, K, P, W  (Domain, Kingdom, Phylum, carl-Woese group)
LEFT_STRIP_X = [0.330, 0.341, 0.352, 0.363]
RIGHT_STRIP_X = [0.670, 0.659, 0.648, 0.637]
STRIP_W = 0.008
STRIP_H = 0.78

DOMAIN_COLORS = tc.DOMAIN_COLORS
MISSING_COLOR = tc.MISSING_COLOR
KINGDOM_COLORS = tc.KINGDOM_COLORS
PHYLUM_COLORS = tc.PHYLUM_COLORS
WOESE_COLORS = tc.WOESE_COLORS

LINE_LW = 0.65
LINE_ALPHA = 0.45
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
    cnt = 0
    for i in range(len(sequence)):
        for j in range(i + 1, len(sequence)):
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


def add_strip(ax, x, y, color):
    ax.add_patch(Rectangle((x - STRIP_W / 2, y - STRIP_H / 2), STRIP_W, STRIP_H, facecolor=color, edgecolor="none"))


def make_metadata_strip_tanglegram(method):
    os.makedirs(OUTDIR, exist_ok=True)

    meta = pd.read_csv(METADATA, sep="\t")
    id_col = meta.columns[0]
    meta_lookup = meta.set_index(id_col).to_dict("index")

    domain_cmap = DOMAIN_COLORS
    kingdom_cmap = KINGDOM_COLORS
    phylum_cmap = PHYLUM_COLORS
    woese_cmap = WOESE_COLORS

    # Which column colors the connecting line (default domain). v4 honors
    # tanglegram_colors.LINE_COLOR_BY so step3 and step4 stay in sync.
    line_col_name = tc.LINE_COLOR_BY
    line_mapping, _ = tc.color_map_for_column(line_col_name)

    t_left = tree_lib.parse_newick(TREE_LEFT_PATH)
    t_right = tree_lib.parse_newick(TREE_RIGHT_PATH)

    t_left, root_label_left = apply_rooting(t_left, method, meta)
    t_right, root_label_right = apply_rooting(t_right, method, meta)

    leaves_l = {t.name for t in t_left.get_terminals()}
    leaves_r = {t.name for t in t_right.get_terminals()}
    common = leaves_l & leaves_r

    # 16S (left) is the REFERENCE tree: its leaf order is fixed to the natural
    # post-rooting order, which is exactly what the circular branch-length tree
    # uses. We never rotate the 16S tree, so its structure and display order are
    # untouched (e.g. the Archaea block stays contiguous, matching the circular
    # tree). Only the FastAAI (right) tree is rotated to minimize crossings.
    order_l = [t.name for t in t_left.get_terminals() if t.name in common]
    order_r = [t.name for t in t_right.get_terminals() if t.name in common]
    cross0 = count_crossings(order_l, order_r)

    target_y_for_right = {n: i for i, n in enumerate(order_l)}
    for _ in range(UNTANGLE_PASSES):
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

        info = meta_lookup.get(name, {})
        dom = info.get("domain")
        kgd = info.get("kingdom")
        phy = info.get("phylum")
        woe = info.get("woese")

        dom_col = tc.resolve(dom, domain_cmap)
        kgd_col = tc.resolve(kgd, kingdom_cmap)
        phy_col = tc.resolve(phy, phylum_cmap)
        woe_col = tc.resolve(woe, woese_cmap)
        line_col = tc.resolve(info.get(line_col_name), line_mapping)

        ax.plot(
            [x_left_actual, CENTER_GAP[0], CENTER_GAP[1], x_right_actual],
            [yl, yl, yr, yr],
            color=line_col,
            lw=LINE_LW,
            alpha=LINE_ALPHA,
            solid_capstyle="round",
        )

        # Left strips: D, K, P, W
        for x, c in zip(LEFT_STRIP_X, [dom_col, kgd_col, phy_col, woe_col]):
            add_strip(ax, x, yl, c)

        # Right strips: D, K, P, W
        for x, c in zip(RIGHT_STRIP_X, [dom_col, kgd_col, phy_col, woe_col]):
            add_strip(ax, x, yr, c)

        ax.text(CENTER_GAP[0] - TIP_PAD, yl, name, ha="right", va="center", fontsize=TIP_FS, color=tc.TIP_LABEL_COLOR)
        ax.text(CENTER_GAP[1] + TIP_PAD, yr, name, ha="left", va="center", fontsize=TIP_FS, color=tc.TIP_LABEL_COLOR)

    # Strip labels
    for label, x in zip(["D", "K", "P", "W"], LEFT_STRIP_X):
        ax.text(x, n + 0.25, label, ha="center", va="bottom", fontsize=8, weight="bold")
    for label, x in zip(["D", "K", "P", "W"], RIGHT_STRIP_X):
        ax.text(x, n + 0.25, label, ha="center", va="bottom", fontsize=8, weight="bold")

    pretty = {
        "archaea_outgroup": "Archaea-outgroup-rooted",
        "unrooted": "unrooted / input-root orientation",
        f"selected_root_{SELECTED_ROOT_ACCESSION.replace('.', '_')}": f"rooted on {SELECTED_ROOT_ACCESSION}",
        "midpoint": "midpoint-rooted",
    }[method]

    fig.suptitle(f"Tanglegram with metadata strips: 16S rRNA vs FastAAI — {pretty}", fontsize=14, weight="bold", y=0.99)
    ax.text((LEFT_X_RANGE[0] + LEFT_X_RANGE[1]) / 2, n + 0.6, f"16S rRNA\n{root_label_left}", ha="center", va="bottom", fontsize=11, weight="bold")
    ax.text((RIGHT_X_RANGE[0] + RIGHT_X_RANGE[1]) / 2, n + 0.6, f"FastAAI NJ\n{root_label_right}", ha="center", va="bottom", fontsize=11, weight="bold")

    # Compact legends outside tip-label area, listing only categories present
    # in this figure, colored by the manual maps.
    common_domains = [meta_lookup.get(nm, {}).get("domain") for nm in common]
    common_kingdoms = [meta_lookup.get(nm, {}).get("kingdom") for nm in common]
    common_phyla = [meta_lookup.get(nm, {}).get("phylum") for nm in common]
    common_woese = [meta_lookup.get(nm, {}).get("woese") for nm in common]

    domain_legend = tc.present_color_map(common_domains, domain_cmap)
    kingdom_legend = tc.present_color_map(common_kingdoms, kingdom_cmap)
    phylum_legend = tc.present_color_map(common_phyla, phylum_cmap)
    woese_legend = tc.present_color_map(common_woese, woese_cmap)

    # Left legend column: Domain, Kingdom, Woese.
    domain_handles = [mpatches.Patch(color=c, label=k) for k, c in domain_legend.items()]
    leg1 = fig.legend(handles=domain_handles, title="Domain", loc="upper left", bbox_to_anchor=(0.715, 0.965), frameon=False, fontsize=8, title_fontsize=9)

    kingdom_handles = [mpatches.Patch(color=c, label=k) for k, c in kingdom_legend.items()]
    leg2 = fig.legend(handles=kingdom_handles, title="Kingdom", loc="upper left", bbox_to_anchor=(0.715, 0.84), frameon=False, fontsize=7.5, title_fontsize=9)

    woese_handles = [mpatches.Patch(color=c, label=k) for k, c in woese_legend.items()]
    leg3 = fig.legend(handles=woese_handles, title="Carl Woese historical group (W)", loc="upper left", bbox_to_anchor=(0.715, 0.60), frameon=False, fontsize=7, title_fontsize=9)

    # Right legend column: Phylum (tall, two columns).
    phylum_handles = [mpatches.Patch(color=c, label=k) for k, c in phylum_legend.items()]
    fig.legend(handles=phylum_handles, title="Phylum (P)", loc="upper left", bbox_to_anchor=(0.845, 0.965), frameon=False, fontsize=6.2, title_fontsize=9, ncol=2)

    ax.text(
        0.5,
        -0.5,
        f"line crossings before/after untangling: {cross0} → {final_cross}; strips: D=domain, K=kingdom, P=phylum, W=Carl Woese group",
        ha="center",
        va="top",
        fontsize=9,
        color="#555555",
        style="italic",
    )

    fig.subplots_adjust(left=0.01, right=0.70, top=0.92, bottom=0.03)

    out_png = f"{OUTDIR}/tanglegram_metadata_strips_{method}.png"
    out_pdf = f"{OUTDIR}/tanglegram_metadata_strips_{method}.pdf"
    fig.savefig(out_png, dpi=DPI, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)

    if method == "archaea_outgroup":
        import shutil
        shutil.copy2(out_png, "output/tanglegram_with_metadata_strips_archaea_rooted.png")
        shutil.copy2(out_pdf, "output/tanglegram_with_metadata_strips_archaea_rooted.pdf")

    return {"rooting_method": method, "crossings_before": cross0, "crossings_after": final_cross}


def main():
    rows = []
    for method in ROOTING_METHODS:
        rows.append(make_metadata_strip_tanglegram(method))
    pd.DataFrame(rows).to_csv("output/tanglegram_metadata_strips_crossing_summary.tsv", sep="\t", index=False)
    print("[done] wrote all metadata-strip tanglegrams")


if __name__ == "__main__":
    main()
