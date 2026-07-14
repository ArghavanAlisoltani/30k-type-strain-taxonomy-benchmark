#!/usr/bin/env python3
"""
=============================================================================
Circular phylogenetic tree with 3 metadata rings (domain, kingdom, phylum)
=============================================================================

PIPELINE:
  1. Load Newick tree and metadata TSV.
  2. Extract GCA ID (e.g. "GCA_000008545.1") from each tip label and merge
     with the metadata table on its first column ("ID").
  3. Draw the tree in circular layout using ete3 (or fall back to matplotlib).
  4. Render 3 concentric heatmap-style rings around the tree:
       Ring 1 (innermost) = domain
       Ring 2 (middle)    = kingdom
       Ring 3 (outermost) = phylum
  5. Show a scale bar (branch-length units) and three legends.

OUTPUTS:
  output/circular_tree.png   (300 dpi raster)
  output/circular_tree.pdf   (vector)
  output/tip_metadata_merged.tsv  (audit table: which tips matched what)

DEPENDENCIES:
  python >= 3.9
  matplotlib, pandas, numpy
  (No biopython / ete3 / dendropy needed — a small Newick parser is bundled
   below so this runs in any headless environment without extra installs.
   If you have biopython, swap `parse_newick()` for Bio.Phylo.read.)

-----------------------------------------------------------------------------
HOW TO CUSTOMISE (search for the matching tag in the CONFIG block below):
  [SPACING]  ring widths, gap between tree and rings, padding
  [FONTS]    tip label font size / family / weight
  [COLORS]   per-category palette assignments
  [LEGEND]   legend position, columns, title font sizes
  [SCALE]    scale-bar length, position, label
  [LABELS]   show/hide tip labels, truncate long names
=============================================================================
"""

import os
import re
import sys
import math
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ----------------------------------------------------------------------------
# Minimal Newick parser (no external deps)
#
# Builds a tree out of lightweight Clade objects:
#   .name        str | None   — tip name or internal label/support
#   .branch_length float       — distance from parent (0 if absent)
#   .clades      list[Clade]   — direct children
# and exposes the small surface area the renderer uses:
#   tree.find_clades(order=)   tree.get_terminals()   tree.depths()
# ----------------------------------------------------------------------------
class Clade:
    __slots__ = ("name", "branch_length", "clades", "confidence")
    def __init__(self):
        self.name = None
        self.branch_length = 0.0
        self.clades = []
        self.confidence = None

    def is_terminal(self):
        return len(self.clades) == 0

    def get_terminals(self):
        out = []
        stack = [self]
        while stack:
            c = stack.pop()
            if c.is_terminal():
                out.append(c)
            else:
                # right-to-left so left-most tips come out first
                stack.extend(reversed(c.clades))
        return out

    def find_clades(self, order="preorder"):
        if order == "preorder":
            stack = [self]
            while stack:
                c = stack.pop()
                yield c
                stack.extend(reversed(c.clades))
        elif order == "postorder":
            # iterative postorder
            stack = [(self, False)]
            while stack:
                node, visited = stack.pop()
                if visited:
                    yield node
                else:
                    stack.append((node, True))
                    for ch in reversed(node.clades):
                        stack.append((ch, False))
        else:
            raise ValueError(order)

    def depths(self):
        """Return {clade: cumulative branch length from this node}."""
        d = {self: 0.0}
        for clade in self.find_clades(order="preorder"):
            for ch in clade.clades:
                d[ch] = d[clade] + (ch.branch_length or 0.0)
        return d


def parse_newick(path):
    """Parse a Newick file. Supports branch lengths and bootstrap labels."""
    with open(path) as f:
        s = f.read().strip()
    if s.endswith(";"):
        s = s[:-1]

    i = 0
    def parse_clade():
        nonlocal i
        node = Clade()
        if s[i] == "(":
            i += 1
            node.clades.append(parse_clade())
            while s[i] == ",":
                i += 1
                node.clades.append(parse_clade())
            assert s[i] == ")", f"expected ) at {i}: {s[i:i+20]}"
            i += 1
        # name / support
        start = i
        while i < len(s) and s[i] not in ",():":
            i += 1
        label = s[start:i]
        # branch length
        bl = 0.0
        if i < len(s) and s[i] == ":":
            i += 1
            start = i
            while i < len(s) and s[i] not in ",()":
                i += 1
            bl = float(s[start:i])
        node.branch_length = bl
        if label:
            if node.is_terminal():
                node.name = label
            else:
                # internal label: usually a bootstrap/support value
                try:
                    node.confidence = float(label)
                except ValueError:
                    node.name = label
        return node

    root = parse_clade()
    return root


# ----------------------------------------------------------------------------
# Tree rooting helpers
#
# IQ-TREE writes UNROOTED trees as a degree-3 root. The two functions below
# transform that into a rooted, binary tree:
#   - root_midpoint(tree)         — places the root at the midpoint of the
#                                   longest tip-to-tip path. No biology
#                                   required; just geometry.
#   - root_outgroup(tree, fn)     — places the root on the branch leading to
#                                   the MRCA of all tips that satisfy a
#                                   predicate `fn(tip_name) -> bool`. Used
#                                   here to root with archaeal tips as the
#                                   outgroup.
#
# Both work by:
#   1. Building parent/child pointers and computing every node's distance
#      from the current root.
#   2. Identifying the branch on which the new root should sit.
#   3. Bisecting that branch and re-orienting all parent pointers so the new
#      root becomes the top of the tree.
# ----------------------------------------------------------------------------
def _build_parents(tree):
    """Return {child -> parent} for every non-root node."""
    parents = {}
    for clade in tree.find_clades(order="preorder"):
        for ch in clade.clades:
            parents[ch] = clade
    return parents


def _reroot_on_branch(tree, target_child, distance_from_child):
    """
    Insert a new root on the branch above `target_child`, at `distance_from_child`
    along that branch (measured from target_child toward its parent).
    Returns the new root Clade.
    """
    parents = _build_parents(tree)
    if target_child not in parents:
        # already at root — nothing to do
        return tree
    parent = parents[target_child]
    bl = target_child.branch_length or 0.0
    d = max(0.0, min(distance_from_child, bl))   # clamp into branch

    new_root = Clade()
    new_root.branch_length = 0.0

    # First piece: target_child gets a shortened branch
    target_child.branch_length = d

    # Second piece: walk from `parent` up to the OLD root, flipping
    # parent/child relationships so `parent` becomes a child of new_root
    # with branch length (bl - d), and everything above it cascades down.
    # Detach target_child from parent first
    parent.clades = [c for c in parent.clades if c is not target_child]

    # Reverse the chain parent -> old_parent -> ... -> old_root
    chain = []
    node = parent
    while node is not None:
        chain.append(node)
        node = parents.get(node)
    # chain[0] = parent, chain[-1] = old root

    # Re-link in reverse: each node becomes the child of the previous one
    # in the chain, with the branch length of the previous node.
    # The first link uses (bl - d) as the branch from new_root to `parent`.
    prev_bl = bl - d
    new_root.clades = [target_child, parent]
    target_child.branch_length = d
    # Detach successive nodes from their original parents and re-attach
    for i in range(len(chain) - 1):
        cur = chain[i]
        nxt = chain[i + 1]
        # remove `cur` from `nxt.clades` (it was cur's old parent)
        nxt.clades = [c for c in nxt.clades if c is not cur]
        # `nxt` becomes a child of `cur`, taking `cur`'s old branch length
        nxt_old_bl = cur.branch_length  # this is the branch *into* cur from nxt
        cur.clades.append(nxt)
        nxt.branch_length = nxt_old_bl

    # Branch from new_root to `parent`
    parent.branch_length = prev_bl

    # Collapse any unifurcations (degree-2 internal nodes) produced by the flip
    _collapse_unifurcations(new_root)
    return new_root


def _collapse_unifurcations(root):
    """Splice out any internal node with a single child (degree-2 chains)."""
    # Iterate until stable
    changed = True
    while changed:
        changed = False
        for clade in list(root.find_clades(order="preorder")):
            new_children = []
            for ch in clade.clades:
                if (not ch.is_terminal()) and len(ch.clades) == 1:
                    only = ch.clades[0]
                    only.branch_length = (only.branch_length or 0.0) + (ch.branch_length or 0.0)
                    new_children.append(only)
                    changed = True
                else:
                    new_children.append(ch)
            clade.clades = new_children


def root_midpoint(tree):
    """
    Midpoint rooting: place the root at the midpoint of the longest path
    between any two tips. Standard approach when no outgroup is available.
    """
    terminals = tree.get_terminals()

    # 1. Find tip farthest from an arbitrary start tip (BFS-style on tree)
    def farthest_tip_from(start):
        # cumulative distances from `start` to every other tip, walking the
        # undirected tree graph
        parents = _build_parents(tree)
        # adjacency: each node connects to its parent (if any) and its children
        def neighbors(n):
            out = [(c, c.branch_length or 0.0) for c in n.clades]
            if n in parents:
                out.append((parents[n], n.branch_length or 0.0))
            return out
        dist = {start: 0.0}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nb, w in neighbors(cur):
                if nb not in dist:
                    dist[nb] = dist[cur] + w
                    stack.append(nb)
        far_tip = max(terminals, key=lambda t: dist.get(t, -1))
        return far_tip, dist

    a, _ = farthest_tip_from(terminals[0])
    b, dist_from_a = farthest_tip_from(a)
    max_dist = dist_from_a[b]

    # 2. Walk from b back toward a along the tree, accumulating branch length,
    #    until we cross the halfway point. That's where the new root goes.
    parents = _build_parents(tree)
    def path_b_to_a():
        # Find the path b -> a by BFS predecessors
        pred = {a: None}
        stack = [a]
        def neighbors(n):
            out = [(c, c.branch_length or 0.0) for c in n.clades]
            if n in parents:
                out.append((parents[n], n.branch_length or 0.0))
            return out
        while stack:
            cur = stack.pop()
            for nb, _w in neighbors(cur):
                if nb not in pred:
                    pred[nb] = cur
                    stack.append(nb)
        # reconstruct b -> ... -> a, then reverse to a -> ... -> b
        path = [b]
        while pred[path[-1]] is not None:
            path.append(pred[path[-1]])
        return list(reversed(path))   # a, ..., b

    path = path_b_to_a()
    half = max_dist / 2.0

    # Walk from a along the path; find the edge we cross at the halfway mark
    acc = 0.0
    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        # edge length is whichever node is the child of the other
        if parents.get(v) is u:
            edge_len = v.branch_length or 0.0
            child, parent_node = v, u
        else:
            edge_len = u.branch_length or 0.0
            child, parent_node = u, v
        if acc + edge_len >= half:
            d_from_child_to_parent = (acc + edge_len) - half
            # we want to root the branch above `child` at distance
            # `d_from_child_to_parent` from `child`'s end
            # ...but path direction matters: we walked from a toward b, so the
            # "from a" distance is `acc + edge_len` when we reach v.
            # If `child` is v (i.e. v is the child of u), distance from child
            # along the branch toward parent = edge_len - (half - acc)
            if child is v:
                d_from_child = edge_len - (half - acc)
            else:
                d_from_child = half - acc
            return _reroot_on_branch(tree, child, d_from_child)
        acc += edge_len
    return tree   # already at midpoint (shouldn't happen)


def root_outgroup(tree, is_outgroup_fn):
    """
    Outgroup rooting: root on the branch leading to the MRCA of all tips
    where `is_outgroup_fn(tip.name) -> True`. The root is placed at the
    midpoint of that branch (a common convention).
    """
    terminals = tree.get_terminals()
    outgroup_tips = [t for t in terminals if is_outgroup_fn(t.name)]
    if not outgroup_tips:
        print("[root] no outgroup tips matched; tree left unrooted")
        return tree
    if len(outgroup_tips) == len(terminals):
        print("[root] every tip matches the outgroup; tree left unrooted")
        return tree

    # Find MRCA of outgroup tips
    parents = _build_parents(tree)
    def ancestors(t):
        chain = [t]
        while chain[-1] in parents:
            chain.append(parents[chain[-1]])
        return chain
    # Intersect ancestor lists; the first common ancestor walking up is MRCA
    common = set(ancestors(outgroup_tips[0]))
    for t in outgroup_tips[1:]:
        common &= set(ancestors(t))
    # Pick the deepest common ancestor (the one farthest from the current root)
    depths = tree.depths()
    mrca = max(common, key=lambda n: depths.get(n, 0))

    # Verify the MRCA's clade contains exactly the outgroup (a clean split).
    # If it doesn't, the outgroup isn't monophyletic — we still root on the
    # MRCA's branch, but warn the user.
    mrca_tips = set(mrca.get_terminals())
    extra = mrca_tips - set(outgroup_tips)
    missing = set(outgroup_tips) - mrca_tips
    if extra or missing:
        print(f"[root] WARNING: outgroup is non-monophyletic — MRCA clade "
              f"contains {len(extra)} extra tip(s) and is missing "
              f"{len(missing)} outgroup tip(s). Rooting on MRCA's branch anyway.")

    bl = mrca.branch_length or 0.0
    return _reroot_on_branch(tree, mrca, bl / 2.0)

# =============================================================================
# CONFIG  --  edit anything in this block to change appearance
# =============================================================================

# --- Input / output paths ----------------------------------------------------
TREE_FILE     = "output/fastaai_nj.treefile"
METADATA_FILE = "inputs/metadata_with_kingdom.tsv"
OUT_DIR       = "output"

# --- Rooting ---------------------------------------------------------------- [ROOT]
# Choose one of:
#   "none"      — leave tree as-is (IQ-TREE's arbitrary unrooted layout)
#   "midpoint"  — root at midpoint of longest tip-to-tip path
#   "outgroup"  — root on the branch leading to MRCA of OUTGROUP_TIPS_FN tips
# Output filenames automatically get a suffix like _midpoint or _outgroup.
ROOT_METHOD = "midpoint"

# Predicate used by ROOT_METHOD="outgroup". Receives a tip name (e.g.
# "GCA_000008545.1_rRNA_AE000512.1_188975-190520_DIR_") and returns True if
# that tip belongs to the outgroup. The default below marks every ARCHAEAL
# tip as outgroup, using the metadata table loaded in main().
#
# The actual lookup is wired up inside main() (see `outgroup_fn`) — this
# constant just names the domain to treat as the outgroup. Set to "bacteria"
# to flip it, or rewrite `outgroup_fn` for a custom rule.
OUTGROUP_DOMAIN = "Archea"   # note metadata spelling — see metadata file

# --- Figure geometry --------------------------------------------------------- [SPACING]
FIG_SIZE      = (22, 16)      # inches; width includes ~6in legend column
DPI           = 300           # raster export DPI
TREE_RADIUS   = 1.00          # outer radius of the tree itself (unitless)
RING_GAP      = 0.04          # radial gap between tree edge and ring 1
RING_WIDTH    = 0.06          # thickness of each metadata ring
RING_SPACING  = 0.005         # small gap between adjacent rings
LABEL_PAD     = 0.025         # gap between outermost ring and tip-label text
PLOT_LIMIT    = 1.75          # axis half-extent; raise if tip labels clip

# --- Branch-length display mode --------------------------------------------- [SPACING]
# AAI / protein-distance trees often have very long terminal branches and
# very short internal ones, which compresses the topology into a tiny spot
# near the center. To see the topology, this script offers two options:
#
#   "phylogram"  — honest branch lengths; everything is exact but internal
#                  structure may be invisible when terminal branches dominate.
#   "compressed" — terminal branches get a fixed share of the radial extent
#                  (TERMINAL_FRACTION), and the remaining radius is allocated
#                  proportionally to the internal structure. Internal branch
#                  proportions are preserved relative to each other; only the
#                  internal-vs-terminal *balance* is rescaled. Recommended for
#                  AAI trees.
#   "cladogram"  — all tips at the same radius, all internal nodes at evenly
#                  spaced "topology depths" (ignores branch lengths entirely).
BRANCH_DISPLAY     = "compressed"
TERMINAL_FRACTION  = 0.20      # share of the radius reserved for terminal branches
                                # (when BRANCH_DISPLAY = "compressed"). 0.15–0.30 works well.

# --- Tip labels ------------------------------------------------------------- [LABELS] [FONTS]
SHOW_TIP_LABELS  = True
TIP_FONT_SIZE    = 7           # 56 tips -> readable size; bump higher if needed
TIP_FONT_FAMILY  = "DejaVu Sans"
TIP_FONT_WEIGHT  = "normal"
TIP_LABEL_AS_GCA = True        # True -> show "GCA_xxxxxxx.x"; False -> full tip
TIP_LABEL_MAXLEN = 40          # truncate labels longer than this (ellipsis added)

# --- Tree drawing ----------------------------------------------------------- [SPACING]
BRANCH_COLOR  = "#222222"
BRANCH_WIDTH  = 0.6
SHOW_SUPPORT  = False          # IQ-TREE bootstrap labels at internal nodes
SUPPORT_MIN   = 90             # only show support >= this value
SUPPORT_FONT  = 4

# --- Tip-to-ring connector lines -------------------------------------------- [SPACING]
# When tips end at varying depths (typical for AAI / NJ trees), a thin dashed
# line connects each tip to the heatmap ring boundary so the rings still line
# up. Set CONNECTOR_LW = 0 or CONNECTOR_ALPHA = 0 to hide.
CONNECTOR_COLOR = "#888888"
CONNECTOR_LW    = 0.4
CONNECTOR_STYLE = (0, (1, 1.5))    # matplotlib dash pattern: (offset, (on, off))
CONNECTOR_ALPHA = 0.55

# --- Scale bar -------------------------------------------------------------- [SCALE]
SHOW_SCALE_BAR    = True
SCALE_BAR_LENGTH  = 0.05       # in AAI-distance units (d = 1 - AAI/100)
SCALE_BAR_XY      = (-1.65, -1.65)  # bottom-left of plot, in data coords
SCALE_BAR_LW      = 1.5
SCALE_FONT_SIZE   = 9
SCALE_BAR_LABEL   = "0.05 AAI distance units"   # shown below the bar

# --- Color palettes --------------------------------------------------------- [COLORS]
# Use distinct, color-blind-friendly palettes per ring. Categories are
# auto-assigned in alphabetical order from each palette; override below if you
# want a fixed mapping (e.g. DOMAIN_COLORS = {"bacteria": "#1b9e77", ...}).

DOMAIN_PALETTE = [
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
]
KINGDOM_PALETTE = [
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928",
]
PHYLUM_PALETTE = [
    "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
    "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928",
    "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
]
MISSING_COLOR = "#cccccc"      # used when a tip has no metadata entry

# Optional explicit overrides — set keys to override auto-assigned palette colors.
# Requested by user: Archaea ring color = orange, Bacteria ring color = black.
DOMAIN_COLORS  = {
    "Archea":   "#ff7f00",   # orange  (note: metadata spells it "Archea")
    "bacteria": "#000000",   # black
}
KINGDOM_COLORS = {}
PHYLUM_COLORS  = {}

# --- Legend ----------------------------------------------------------------- [LEGEND]
LEGEND_TITLE_FS   = 10
LEGEND_LABEL_FS   = 6.5
LEGEND_MARKER_SZ  = 8
LEGEND_FRAME      = False
# Each legend has its own ncol — phylum has 56 entries so it gets 2 columns
LEGEND_DOMAIN_NCOL  = 1
LEGEND_KINGDOM_NCOL = 1
LEGEND_PHYLUM_NCOL  = 2
# Anchor points (in FIGURE-fraction coords). Adjust if you resize FIG_SIZE.
LEGEND_DOMAIN_POS  = (0.80, 0.96)
LEGEND_KINGDOM_POS = (0.80, 0.86)
LEGEND_PHYLUM_POS  = (0.80, 0.62)

# --- Title ------------------------------------------------------------------ [FONTS]
TITLE_TEXT     = "FastAAI Neighbor-Joining phylogeny"
TITLE_FONT_SZ  = 14
TITLE_WEIGHT   = "bold"

# =============================================================================
# END CONFIG
# =============================================================================


# ----------------------------------------------------------------------------
# Helper: pull "GCA_xxxxxxxxx.x" out of a tip label like
# "GCA_000008545.1_rRNA_AE000512.1_188975-190520_DIR_"
# ----------------------------------------------------------------------------
GCA_RE = re.compile(r"GCA_\d+\.\d+")

def gca_from_tip(name: str) -> str | None:
    m = GCA_RE.match(name)
    return m.group(0) if m else None


# ----------------------------------------------------------------------------
# Helper: build {category -> color} dict from a palette
# ----------------------------------------------------------------------------
def build_color_map(categories, palette, override):
    cats = sorted(c for c in set(categories) if pd.notna(c) and c != "")
    cmap = OrderedDict()
    for i, c in enumerate(cats):
        cmap[c] = override.get(c, palette[i % len(palette)])
    return cmap


# ----------------------------------------------------------------------------
# Helper: compute angular position for each tip in radial layout
#
# We post-order traverse to assign tip angles evenly around the circle
# (skipping a small "gap" wedge so the layout doesn't fully close on itself —
# this gap is purely cosmetic; set GAP_FRACTION=0 for a fully closed circle).
# Internal node angles are the mean of their children's angles.
# ----------------------------------------------------------------------------
GAP_FRACTION = 0.0  # 0.0 = closed circle; e.g. 0.02 leaves a 2% wedge at top

def assign_radial_coords(tree):
    """Returns dict node -> (angle_radians, radial_distance_from_root)."""
    terminals = tree.get_terminals()
    n = len(terminals)
    span = 2 * math.pi * (1 - GAP_FRACTION)
    start = math.pi / 2 + (2 * math.pi - span) / 2   # start at top, rotate CCW

    # angle per tip (evenly spaced)
    tip_angle = {}
    # With GAP_FRACTION=0 (closed circle), the angular slot per tip is span/n
    # — NOT span/(n-1). Dividing by n-1 would place the last tip exactly on
    # top of the first one (angle wraps around 2π). Each tip is centered
    # in its own slot, so we add a half-slot offset.
    slot = span / n
    for i, t in enumerate(terminals):
        tip_angle[t] = start - (i + 0.5) * slot

    # radial distance = cumulative branch length from root
    depths = tree.depths()  # {clade: distance_from_root}

    coords = {}
    # tips first
    for t in terminals:
        coords[t] = (tip_angle[t], depths[t])

    # internal nodes: mean angle of descendant tips, depth from .depths()
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            continue
        child_tips = clade.get_terminals()
        ang = sum(tip_angle[t] for t in child_tips) / len(child_tips)
        coords[clade] = (ang, depths[clade])
    return coords


# ----------------------------------------------------------------------------
# Helper: draw a radial branch
#
# In a circular cladogram each branch is rendered in two parts:
#   (a) an arc at the PARENT's radius spanning from parent_angle to child_angle
#   (b) a straight radial line from (child_angle, parent_radius) to
#       (child_angle, child_radius)
# ----------------------------------------------------------------------------
def draw_radial_branch(ax, parent_xy, child_xy, color, lw):
    pa, pr = parent_xy
    ca, cr = child_xy
    # arc at parent radius
    a0, a1 = sorted([pa, ca])
    arc_theta = np.linspace(a0, a1, 30)
    ax.plot(pr * np.cos(arc_theta), pr * np.sin(arc_theta),
            color=color, lw=lw, solid_capstyle="round")
    # radial spoke from parent radius to child radius at child's angle
    ax.plot([pr * math.cos(ca), cr * math.cos(ca)],
            [pr * math.sin(ca), cr * math.sin(ca)],
            color=color, lw=lw, solid_capstyle="round")


# ----------------------------------------------------------------------------
# Helper: draw a wedge for one tip in one metadata ring
# ----------------------------------------------------------------------------
def draw_ring_wedge(ax, angle, half_width, r_inner, r_outer, color):
    """Draw a filled annular wedge centered on `angle`."""
    a0 = angle - half_width
    a1 = angle + half_width
    theta = np.linspace(a0, a1, 12)
    # outer arc forward, inner arc reverse -> closed polygon
    x = np.concatenate([r_outer * np.cos(theta), r_inner * np.cos(theta[::-1])])
    y = np.concatenate([r_outer * np.sin(theta), r_inner * np.sin(theta[::-1])])
    ax.fill(x, y, color=color, linewidth=0)


# ============================================================================
# MAIN
# ============================================================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Output paths get a suffix matching ROOT_METHOD so different rootings
    # don't overwrite each other. We rename "none" -> "unrooted" for the
    # user-facing label since that's what biologists call it.
    label = "unrooted" if ROOT_METHOD == "none" else ROOT_METHOD
    suffix = f"_{label}"
    out_png    = os.path.join(OUT_DIR, f"circular_tree{suffix}.png")
    out_pdf    = os.path.join(OUT_DIR, f"circular_tree{suffix}.pdf")
    out_merged = os.path.join(OUT_DIR, f"tip_metadata_merged{suffix}.tsv")

    # --- 1. Load metadata ---------------------------------------------------
    meta = pd.read_csv(METADATA_FILE, sep="\t")
    meta.columns = [c.strip() for c in meta.columns]
    id_col = meta.columns[0]   # always first column per spec
    print(f"[meta] {len(meta)} rows; columns = {list(meta.columns)}")

    # --- 2. Load tree -------------------------------------------------------
    tree = parse_newick(TREE_FILE)
    print(f"[tree] loaded with {len(tree.get_terminals())} tips, "
          f"root degree = {len(tree.clades)} "
          f"({'unrooted' if len(tree.clades) > 2 else 'rooted'})")

    # --- 2b. Root the tree --------------------------------------------------
    if ROOT_METHOD == "midpoint":
        tree = root_midpoint(tree)
        print(f"[root] applied midpoint rooting (new root degree = {len(tree.clades)})")
    elif ROOT_METHOD == "outgroup":
        # Build a lookup: GCA -> domain, so we can decide outgroup membership
        # from the tip name without parsing metadata again later.
        gca_to_domain = dict(zip(meta[id_col], meta["domain"]))
        def outgroup_fn(tip_name):
            g = gca_from_tip(tip_name)
            return g is not None and gca_to_domain.get(g) == OUTGROUP_DOMAIN
        n_out = sum(1 for t in tree.get_terminals() if outgroup_fn(t.name))
        print(f"[root] outgroup ({OUTGROUP_DOMAIN}) tips found: {n_out}")
        tree = root_outgroup(tree, outgroup_fn)
        print(f"[root] applied outgroup rooting (new root degree = {len(tree.clades)})")
    elif ROOT_METHOD == "none":
        print("[root] no rooting applied; tree remains as IQ-TREE produced it")
    else:
        raise ValueError(f"unknown ROOT_METHOD: {ROOT_METHOD!r}")

    terminals = tree.get_terminals()

    # Build per-tip metadata table (left-join tip -> meta via GCA prefix)
    rows = []
    for t in terminals:
        gca = gca_from_tip(t.name)
        row = {"tip": t.name, "GCA": gca}
        if gca is not None:
            hit = meta[meta[id_col] == gca]
            if len(hit):
                for c in meta.columns:
                    row[c] = hit.iloc[0][c]
        rows.append(row)
    tip_df = pd.DataFrame(rows)
    tip_df.to_csv(out_merged, sep="\t", index=False)
    print(f"[merge] wrote {out_merged} ({len(tip_df)} tips, "
          f"{tip_df['phylum'].notna().sum()} matched)")

    # --- 3. Color maps ------------------------------------------------------
    domain_cmap  = build_color_map(tip_df["domain"],  DOMAIN_PALETTE,  DOMAIN_COLORS)
    kingdom_cmap = build_color_map(tip_df["kingdom"], KINGDOM_PALETTE, KINGDOM_COLORS)
    phylum_cmap  = build_color_map(tip_df["phylum"],  PHYLUM_PALETTE,  PHYLUM_COLORS)

    # --- 4. Layout ----------------------------------------------------------
    coords = assign_radial_coords(tree)

    # Apply the chosen BRANCH_DISPLAY transform. See CONFIG block above for
    # what each mode does. After this block, `coords` maps every node to its
    # final (angle, radius) in tree-space; the tree is then uniformly scaled
    # so the deepest tip sits at TREE_RADIUS.
    if BRANCH_DISPLAY == "phylogram":
        # No transform; honest branch lengths.
        pass

    elif BRANCH_DISPLAY == "cladogram":
        # Replace each node's radius with its topological depth (max # of
        # internal nodes between root and any descendant tip). Tips end up at
        # the same outer radius; internal nodes are evenly spaced.
        def topo_depth(node):
            if node.is_terminal():
                return 0
            return 1 + max(topo_depth(c) for c in node.clades)
        max_topo = topo_depth(tree)
        for n in coords:
            a, _ = coords[n]
            if n.is_terminal():
                coords[n] = (a, max_topo)
            else:
                coords[n] = (a, max_topo - topo_depth(n))

    elif BRANCH_DISPLAY == "compressed":
        # Find the deepest internal node and the deepest tip. Internal nodes
        # collectively get (1 - TERMINAL_FRACTION) of the radius; terminals
        # get TERMINAL_FRACTION. Each internal node's radius is rescaled
        # linearly within the internal band; each tip is placed
        # TERMINAL_FRACTION beyond its parent, scaled by its own terminal
        # branch length so within-clade differences remain visible.
        terminals = tree.get_terminals()
        all_radii = [r for _, r in coords.values()]
        max_tip_r = max(coords[t][1] for t in terminals)
        internal_radii = [coords[n][1] for n in coords if not n.is_terminal()]
        max_internal_r = max(internal_radii) if internal_radii else 0.0

        if max_internal_r <= 0 or max_tip_r <= max_internal_r:
            # Degenerate cases; fall back to phylogram
            pass
        else:
            internal_band = 1.0 - TERMINAL_FRACTION  # 0..internal_band
            terminal_band = TERMINAL_FRACTION         # internal_band..1.0

            # Rescale internal nodes into [0, internal_band]
            for n in list(coords):
                a, r = coords[n]
                if not n.is_terminal():
                    new_r = (r / max_internal_r) * internal_band
                    coords[n] = (a, new_r)

            # Each tip's terminal-branch length is what's left after subtracting
            # its parent's original (un-rescaled) depth. We need parent depths
            # in the ORIGINAL frame to compute terminal lengths.
            parents = _build_parents(tree)
            depths_orig = tree.depths()
            max_term_branch = max(
                (depths_orig[t] - depths_orig[parents[t]]) for t in terminals
                if t in parents
            ) or 1.0

            for t in terminals:
                a, _ = coords[t]
                p = parents.get(t)
                if p is None:
                    new_r = internal_band + terminal_band
                else:
                    parent_new_r = coords[p][1]
                    term_branch = depths_orig[t] - depths_orig[p]
                    new_r = parent_new_r + (term_branch / max_term_branch) * terminal_band
                coords[t] = (a, new_r)

    else:
        raise ValueError(f"unknown BRANCH_DISPLAY: {BRANCH_DISPLAY!r}")

    # Final uniform scale so deepest point reaches TREE_RADIUS
    max_depth = max(r for (_, r) in coords.values())
    scale = TREE_RADIUS / max_depth if max_depth > 0 else 1.0
    coords = {n: (a, r * scale) for n, (a, r) in coords.items()}
    branch_unit_in_data = scale   # display units per original distance unit

    # --- 5. Figure ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.set_aspect("equal")
    ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
    ax.axis("off")

    # --- 6. Draw branches ---------------------------------------------------
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        pa, pr = coords[clade]
        for child in clade.clades:
            ca, cr = coords[child]
            draw_radial_branch(ax, (pa, pr), (ca, cr), BRANCH_COLOR, BRANCH_WIDTH)

        # optional support-value labels at internal nodes
        if SHOW_SUPPORT and clade.confidence is not None:
            try:
                if float(clade.confidence) >= SUPPORT_MIN:
                    ax.text(pr * math.cos(pa), pr * math.sin(pa),
                            f"{int(clade.confidence)}",
                            fontsize=SUPPORT_FONT, ha="center", va="center",
                            color="#555555")
            except (TypeError, ValueError):
                pass

    # --- 7. Metadata rings + tip labels ------------------------------------
    # angular half-width of each tip's wedge (so they tile without gaps)
    n_tips = len(terminals)
    half_w = (2 * math.pi * (1 - GAP_FRACTION)) / n_tips / 2

    r_tree_edge = TREE_RADIUS + RING_GAP
    r_ring1_in  = r_tree_edge
    r_ring1_out = r_ring1_in + RING_WIDTH
    r_ring2_in  = r_ring1_out + RING_SPACING
    r_ring2_out = r_ring2_in + RING_WIDTH
    r_ring3_in  = r_ring2_out + RING_SPACING
    r_ring3_out = r_ring3_in + RING_WIDTH
    r_label     = r_ring3_out + LABEL_PAD

    # Look-up tip -> metadata row for fast access
    tip_meta = {r["tip"]: r for r in tip_df.to_dict("records")}

    for t in terminals:
        ang, r_tip = coords[t]
        row = tip_meta.get(t.name, {})
        col_dom = domain_cmap.get(row.get("domain"),   MISSING_COLOR)
        col_kgd = kingdom_cmap.get(row.get("kingdom"), MISSING_COLOR)
        col_phy = phylum_cmap.get(row.get("phylum"),   MISSING_COLOR)

        # Connector line: thin gray dashed segment from where the tip actually
        # ends (r_tip) out to the inner edge of ring 1. Keeps rings aligned
        # while preserving honest branch lengths in the tree itself. Edit
        # CONNECTOR_* below to restyle or hide.
        if r_tip < r_ring1_in - 1e-6:
            ax.plot([r_tip * math.cos(ang), r_ring1_in * math.cos(ang)],
                    [r_tip * math.sin(ang), r_ring1_in * math.sin(ang)],
                    color=CONNECTOR_COLOR, lw=CONNECTOR_LW,
                    linestyle=CONNECTOR_STYLE, alpha=CONNECTOR_ALPHA,
                    zorder=0)

        draw_ring_wedge(ax, ang, half_w, r_ring1_in, r_ring1_out, col_dom)
        draw_ring_wedge(ax, ang, half_w, r_ring2_in, r_ring2_out, col_kgd)
        draw_ring_wedge(ax, ang, half_w, r_ring3_in, r_ring3_out, col_phy)

        # tip label
        if SHOW_TIP_LABELS:
            label = row.get("GCA") if TIP_LABEL_AS_GCA and row.get("GCA") else t.name
            if label and len(label) > TIP_LABEL_MAXLEN:
                label = label[: TIP_LABEL_MAXLEN - 1] + "…"
            # rotate text along the radial direction; flip on left half so
            # it reads left->right when viewed from outside the circle
            deg = math.degrees(ang)
            if -90 <= deg <= 90 or 270 <= deg <= 360:
                rot, ha = deg, "left"
            else:
                rot, ha = deg + 180, "right"
            ax.text(r_label * math.cos(ang), r_label * math.sin(ang), label,
                    rotation=rot, rotation_mode="anchor",
                    ha=ha, va="center",
                    fontsize=TIP_FONT_SIZE, family=TIP_FONT_FAMILY,
                    weight=TIP_FONT_WEIGHT)

    # --- 8. Scale bar ------------------------------------------------------ [SCALE]
    # Only meaningful in "phylogram" mode (honest branch lengths). For
    # "compressed" and "cladogram" the radial axis no longer corresponds to
    # raw distance, so a scale bar would be misleading.
    if SHOW_SCALE_BAR and BRANCH_DISPLAY == "phylogram":
        x0, y0 = SCALE_BAR_XY
        bar_len = SCALE_BAR_LENGTH * branch_unit_in_data
        ax.plot([x0, x0 + bar_len], [y0, y0], color="black", lw=SCALE_BAR_LW)
        ax.plot([x0, x0], [y0 - 0.01, y0 + 0.01], color="black", lw=SCALE_BAR_LW)
        ax.plot([x0 + bar_len, x0 + bar_len], [y0 - 0.01, y0 + 0.01],
                color="black", lw=SCALE_BAR_LW)
        ax.text(x0 + bar_len / 2, y0 - 0.04,
                SCALE_BAR_LABEL,
                ha="center", va="top", fontsize=SCALE_FONT_SIZE)
    elif SHOW_SCALE_BAR:
        # Show a small note instead so the reader knows branch lengths are
        # transformed.
        ax.text(SCALE_BAR_XY[0], SCALE_BAR_XY[1],
                f"branch lengths: {BRANCH_DISPLAY} (radial scale not linear)",
                ha="left", va="bottom", fontsize=SCALE_FONT_SIZE,
                style="italic", color="#555555")

    # --- 9. Legends -------------------------------------------------------- [LEGEND]
    def make_handles(cmap):
        return [mpatches.Patch(color=c, label=str(k)) for k, c in cmap.items()]

    leg1 = fig.legend(handles=make_handles(domain_cmap), title="Domain",
                      loc="upper left", bbox_to_anchor=LEGEND_DOMAIN_POS,
                      fontsize=LEGEND_LABEL_FS, title_fontsize=LEGEND_TITLE_FS,
                      frameon=LEGEND_FRAME, ncol=LEGEND_DOMAIN_NCOL)
    leg2 = fig.legend(handles=make_handles(kingdom_cmap), title="Kingdom",
                      loc="upper left", bbox_to_anchor=LEGEND_KINGDOM_POS,
                      fontsize=LEGEND_LABEL_FS, title_fontsize=LEGEND_TITLE_FS,
                      frameon=LEGEND_FRAME, ncol=LEGEND_KINGDOM_NCOL)
    leg3 = fig.legend(handles=make_handles(phylum_cmap), title="Phylum",
                      loc="upper left", bbox_to_anchor=LEGEND_PHYLUM_POS,
                      fontsize=LEGEND_LABEL_FS, title_fontsize=LEGEND_TITLE_FS,
                      frameon=LEGEND_FRAME, ncol=LEGEND_PHYLUM_NCOL)

    # ring labels: tiny annotations placed at the 6 o'clock outer edge of each
    # ring so the viewer can identify which ring is which without consulting
    # the legend. (Change the angle below to relocate; e.g. -math.pi/2 = 6 o'clock,
    # 0 = 3 o'clock, math.pi/2 = 12 o'clock.)
    label_angle = -math.pi / 2
    for r_mid, name in [
        ((r_ring1_in + r_ring1_out) / 2, "domain"),
        ((r_ring2_in + r_ring2_out) / 2, "kingdom"),
        ((r_ring3_in + r_ring3_out) / 2, "phylum"),
    ]:
        ax.text(r_mid * math.cos(label_angle), r_mid * math.sin(label_angle),
                name, fontsize=6, color="#222222",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.18", fc="white",
                          ec="#888888", lw=0.4, alpha=0.95))

    # --- 10. Title --------------------------------------------------------- [FONTS]
    # Title includes the rooting method so the reader can tell figures apart
    title_with_root = (f"{TITLE_TEXT}  —  unrooted" if ROOT_METHOD == "none"
                       else f"{TITLE_TEXT}  —  {ROOT_METHOD}-rooted")
    fig.suptitle(title_with_root, fontsize=TITLE_FONT_SZ, weight=TITLE_WEIGHT, y=0.96)

    # The tree occupies the left ~70% of the figure; legends sit in the right margin.
    # If you change FIG_SIZE, also re-tune these and LEGEND_*_POS together.
    fig.subplots_adjust(left=0.02, right=0.72, top=0.95, bottom=0.05)
    fig.savefig(out_png, dpi=DPI, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    print(f"[done] {out_png}\n[done] {out_pdf}")


if __name__ == "__main__":
    main()
