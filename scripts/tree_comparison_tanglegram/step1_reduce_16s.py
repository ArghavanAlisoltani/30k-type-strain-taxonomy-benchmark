#!/usr/bin/env python3
"""
=============================================================================
Step 1: reduce the 110-tip 16S tree to one tip per genome (56 tips)
=============================================================================

Each genome (GCA ID) is represented by multiple 16S copies in the input tree.
To compare topologies with the FastAAI tree (which has one tip per genome),
we keep one representative 16S tip per genome.

SELECTION RULE
  Default = "longest": pick the 16S copy with the longest sequence in the
  FASTA, since this is the convention most published 16S workflows use.
  Alternatives: "first" (deterministic, alphabetical) and "shortest".

METHOD
  Tree-pruning, not re-inference. We remove non-representative tips from the
  existing IQ-TREE topology and collapse the resulting degree-2 internal
  nodes by summing their branch lengths. This preserves the original
  inference (no need for MAFFT/IQ-TREE binaries) at the cost of inheriting
  any branch-length compromises in the original 110-tip estimate.

FASTA <-> tree name matching
  The tree converted '+' to '_' in tip names (Newick reserved chars), so we
  normalise FASTA IDs the same way before matching. After this normalisation
  all 110 FASTA records match all 110 tree tips by full name.

INPUT
  inputs/16S_full.treefile        — original 110-tip IQ-TREE tree
  inputs/all_rnammer_16S.fasta    — 110 sequences with rich headers

OUTPUT
  output/16S_reduced.treefile     — pruned tree, 56 tips (one per GCA)
  output/16S_representatives.tsv  — audit: which copy was chosen per GCA
  output/16S_reduced_tip_to_gca.tsv  — tip-name -> GCA mapping
=============================================================================
"""
import os
import re
import sys
from collections import defaultdict

# Reuse the Newick parser and tree-manipulation helpers from tree_lib
import tree_lib

# --- CONFIG -----------------------------------------------------------------
TREE_IN     = "inputs/16S_full.treefile"
FASTA_IN    = "inputs/all_rnammer_16S.fasta"
TREE_OUT    = "output/16S_reduced.treefile"
REPS_OUT    = "output/16S_representatives.tsv"
MAP_OUT     = "output/16S_reduced_tip_to_gca.tsv"

SELECTION_RULE = "longest"   # "longest" | "shortest" | "first"
# ----------------------------------------------------------------------------


def normalise_fasta_id(raw_id):
    """Match the transformation IQ-TREE/Newick applied to tip names.
    Specifically, trailing '+' and '-' get replaced with '_' in some
    Newick writers; we mirror that so FASTA IDs align with tree tips."""
    if raw_id.endswith("+"):
        return raw_id[:-1] + "_"
    # '-' is kept as-is in this dataset's Newick; no change needed
    return raw_id


def load_fasta_lengths(path):
    """Return {normalised_id: sequence_length}. Streams the file once."""
    lengths = {}
    cur_id = None
    cur_len = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    lengths[cur_id] = cur_len
                cur_id = normalise_fasta_id(line[1:].split()[0])
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_id is not None:
            lengths[cur_id] = cur_len
    return lengths


def gca_from_tip(name):
    m = re.match(r"GCA_\d+\.\d+", name)
    return m.group(0) if m else None


def choose_representative(tips_for_gca, lengths, rule):
    """Given a list of tip names for one GCA, pick the chosen one.
    `lengths` maps tip name -> sequence length; missing tips fall back to 0."""
    if rule == "first":
        return sorted(tips_for_gca)[0]
    paired = [(t, lengths.get(t, 0)) for t in tips_for_gca]
    if rule == "longest":
        return max(paired, key=lambda x: (x[1], x[0]))[0]
    if rule == "shortest":
        return min(paired, key=lambda x: (x[1], x[0]))[0]
    raise ValueError(f"unknown rule: {rule}")


# ---------------------------------------------------------------------------
# Tree pruning: keep only the `keep` set of tips, splicing out resulting
# degree-2 internal nodes.
#
# Strategy: traverse postorder; for any internal node, drop children whose
# subtree contains no `keep` tips; collapse single-child internal nodes by
# summing branch lengths.
# ---------------------------------------------------------------------------
def prune_tree(root, keep_set):
    def prune(node):
        if node.is_terminal():
            return node if node.name in keep_set else None
        # Recurse and filter children
        kept = []
        for ch in node.clades:
            pruned = prune(ch)
            if pruned is not None:
                kept.append(pruned)
        if not kept:
            return None
        node.clades = kept
        # Collapse degree-2 internal nodes: a node with only one child has no
        # topological signal, so absorb the child's branch into ours.
        if len(node.clades) == 1 and not node.is_terminal():
            only = node.clades[0]
            only.branch_length = (only.branch_length or 0.0) + (node.branch_length or 0.0)
            return only
        return node

    pruned_root = prune(root)
    # Edge case: if the root itself became degree-1 after pruning, walk down
    # until we find a real bifurcation.
    while pruned_root is not None and not pruned_root.is_terminal() \
            and len(pruned_root.clades) == 1:
        only = pruned_root.clades[0]
        only.branch_length = (only.branch_length or 0.0) + (pruned_root.branch_length or 0.0)
        pruned_root = only
    return pruned_root


def write_newick(node):
    """Serialise a Clade tree back to Newick. Branch lengths included; no
    bootstrap labels (the pruned topology's support is no longer trustworthy
    after collapsing nodes, so we omit it rather than mislead)."""
    def serialise(n):
        if n.is_terminal():
            return f"{n.name}:{n.branch_length or 0:.8f}"
        inner = ",".join(serialise(c) for c in n.clades)
        return f"({inner}):{n.branch_length or 0:.8f}"
    # Root has no incoming branch length
    if node.is_terminal():
        return f"{node.name};"
    inner = ",".join(serialise(c) for c in node.clades)
    return f"({inner});"


def main():
    os.makedirs("output", exist_ok=True)

    print(f"[load] tree: {TREE_IN}")
    tree = tree_lib.parse_newick(TREE_IN)
    tips = [t.name for t in tree.get_terminals()]
    print(f"[load] {len(tips)} tips")

    print(f"[load] fasta: {FASTA_IN}")
    lengths = load_fasta_lengths(FASTA_IN)
    print(f"[load] {len(lengths)} sequences")

    # Sanity: how many tree tips have a matching FASTA entry?
    matched = sum(1 for t in tips if t in lengths)
    print(f"[match] tree tips with FASTA length: {matched}/{len(tips)}")

    # Group tips by GCA ID
    gca_to_tips = defaultdict(list)
    for t in tips:
        g = gca_from_tip(t)
        if g is None:
            print(f"[warn] tip without GCA prefix: {t!r}; skipping")
            continue
        gca_to_tips[g].append(t)
    print(f"[group] {len(gca_to_tips)} unique GCA IDs")

    # Choose one representative per GCA
    reps = {}
    rep_rows = []
    for gca, tip_list in sorted(gca_to_tips.items()):
        chosen = choose_representative(tip_list, lengths, SELECTION_RULE)
        reps[gca] = chosen
        rep_rows.append({
            "GCA": gca,
            "n_copies": len(tip_list),
            "chosen_tip": chosen,
            "chosen_length": lengths.get(chosen, 0),
            "all_tips": ";".join(tip_list),
            "all_lengths": ";".join(str(lengths.get(t, 0)) for t in tip_list),
        })

    # Write representatives audit table
    with open(REPS_OUT, "w") as f:
        cols = ["GCA", "n_copies", "chosen_tip", "chosen_length",
                "all_tips", "all_lengths"]
        f.write("\t".join(cols) + "\n")
        for row in rep_rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(f"[write] {REPS_OUT} ({len(rep_rows)} GCA IDs)")

    # Tip -> GCA mapping for downstream relabelling
    with open(MAP_OUT, "w") as f:
        f.write("tip\tGCA\n")
        for gca, tip in reps.items():
            f.write(f"{tip}\t{gca}\n")
    print(f"[write] {MAP_OUT}")

    # Prune the tree, keeping only the chosen representative tips
    keep = set(reps.values())
    pruned = prune_tree(tree, keep)
    pruned_tips = [t.name for t in pruned.get_terminals()]
    print(f"[prune] {len(pruned_tips)} tips remain "
          f"(expected {len(reps)})")

    # Relabel tips from the long rRNA name to just the GCA so the pruned
    # tree shares a leaf-naming convention with the FastAAI tree.
    tip_to_gca = {tip: gca for gca, tip in reps.items()}
    for t in pruned.get_terminals():
        t.name = tip_to_gca.get(t.name, t.name)

    # Write pruned tree
    newick = write_newick(pruned)
    with open(TREE_OUT, "w") as f:
        f.write(newick + "\n")
    print(f"[write] {TREE_OUT}")


if __name__ == "__main__":
    main()
