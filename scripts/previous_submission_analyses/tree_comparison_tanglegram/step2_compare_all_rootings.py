#!/usr/bin/env python3
"""
=============================================================================
Step 2 v3: compare 16S and FastAAI trees under multiple rooting strategies
=============================================================================

This script extends the v2 comparison in four ways:

1. Uses the current metadata file and harmonizes the domain spelling:
      "Archea" -> "Archaea"
   The original metadata is not overwritten; a harmonized copy is written to:
      output/metadata_harmonized.tsv

2. Computes metrics under four rooting/drawing strategies:
      - archaea_outgroup
      - unrooted
      - selected_root_GCA_021654395_1
      - midpoint

3. Keeps the v2 quantitative metrics and adds a root-wise comparison table:
      - unrooted RF distance
      - normalized RF distance
      - rooted cluster distance
      - weighted RF / branch-score distance
      - exact quartet distance
      - Mantel Pearson and Spearman correlations on patristic distances
      - per-domain, per-kingdom, per-phylum monophyly

4. Saves rooted/transformed tree files for all root modes.

Important interpretation
------------------------
- Standard RF is reported as an unrooted split/bipartition comparison.
  It is therefore expected to be identical or nearly identical across rooting
  choices because rerooting does not change the unrooted topology.

- Rooted cluster distance and monophyly checks depend on the selected root.
  These are the values to compare when asking whether root choice changes
  clade interpretation.

Inputs
------
  output/16S_reduced.treefile
  inputs/fastaai_nj.treefile
  inputs/metadata_with_kingdom.tsv

Outputs
-------
  output/tree_comparison_metrics_by_rooting.tsv
  output/tree_comparison_metrics.tsv  [archa­ea_outgroup copy for v2 compatibility]
  output/per_domain_monophyly_by_rooting.tsv
  output/per_kingdom_monophyly_by_rooting.tsv
  output/per_phylum_monophyly_by_rooting.tsv
  output/trees/*.treefile
=============================================================================
"""
import os
import math
import itertools
import copy
from pathlib import Path

import numpy as np
import pandas as pd

import tree_lib

# --- CONFIG -----------------------------------------------------------------
TREE_16S = "output/16S_reduced.treefile"
TREE_AAI = "inputs/fastaai_nj.treefile"
METADATA = "inputs/metadata_with_official_phyla_and_woese_groups.csv"

OUTDIR = "output"
TREE_OUTDIR = "output/trees"

OUTGROUP_DOMAIN_VALUES = {"Archaea", "Archea"}
SELECTED_ROOT_ACCESSION = "GCA_021654395.1"

N_MANTEL_PERMUTATIONS = 999
RANDOM_SEED = 12345

ROOTING_METHODS = [
    "archaea_outgroup",
    "unrooted",
    f"selected_root_{SELECTED_ROOT_ACCESSION.replace('.', '_')}",
    "midpoint",
]
# ----------------------------------------------------------------------------


def write_newick(node):
    """Serialize the lightweight tree_lib.Clade tree to Newick."""
    def serialise(n):
        name = n.name or ""
        bl = n.branch_length or 0.0
        if n.is_terminal():
            return f"{name}:{bl:.10f}"
        inner = ",".join(serialise(c) for c in n.clades)
        label = ""
        if getattr(n, "confidence", None) is not None:
            label = str(n.confidence)
        elif name:
            label = name
        return f"({inner}){label}:{bl:.10f}"

    if node.is_terminal():
        return f"{node.name};"
    return serialise(node).rsplit(":", 1)[0] + ";"


def gca_from_tip(name):
    """Extract a versioned GCA/GCF accession from a tip name."""
    import re
    m = re.search(r"(GC[AF]_\d+\.\d+)", str(name))
    return m.group(1) if m else None


def harmonize_metadata(meta):
    """
    Harmonize taxonomy strings for plotting and analysis.

    The new metadata CSV uses:
      - `official_phylum`              -> aliased to `phylum`
      - `carl_woese_historical_group`  -> aliased to `woese`
      - domain spelled "Bacteria"/"Archaea"

    Domain is standardized to Title case ("Bacteria"/"Archaea"). The new file
    has no missing kingdoms, so no "Unknown" kingdom is produced.
    """
    meta = meta.copy()

    # Alias the new column names to the canonical ones used downstream.
    if "official_phylum" in meta.columns and "phylum" not in meta.columns:
        meta["phylum"] = meta["official_phylum"]
    if "carl_woese_historical_group" in meta.columns and "woese" not in meta.columns:
        meta["woese"] = meta["carl_woese_historical_group"]

    meta = meta.fillna("Unknown")
    for col in meta.columns:
        meta[col] = (
            meta[col]
            .astype(str)
            .str.strip()
            .replace({"": "Unknown", "nan": "Unknown", "NaN": "Unknown", "None": "Unknown"})
        )
    if "domain" in meta.columns:
        meta["domain_original"] = meta["domain"]
        meta["domain"] = meta["domain"].replace({
            "Archea": "Archaea",
            "archaea": "Archaea",
            "bacteria": "Bacteria",
        })
    return meta


def root_with_archaea_outgroup(tree, metadata):
    """Root a tree using all archaeal tips as the outgroup."""
    id_col = metadata.columns[0]
    gca_to_domain = dict(zip(metadata[id_col], metadata["domain"]))

    def is_archaea_tip(tip_name):
        gca = gca_from_tip(tip_name)
        return gca is not None and gca_to_domain.get(gca) in OUTGROUP_DOMAIN_VALUES

    n_out = sum(1 for t in tree.get_terminals() if is_archaea_tip(t.name))
    print(f"[root] archaeal outgroup tips found: {n_out}")
    return tree_lib.root_outgroup(tree, is_archaea_tip)


def root_with_selected_sample(tree, selected_accession):
    """Root a tree on the branch leading to one selected accession/tip."""
    def is_selected_tip(tip_name):
        gca = gca_from_tip(tip_name)
        return tip_name == selected_accession or gca == selected_accession

    n_out = sum(1 for t in tree.get_terminals() if is_selected_tip(t.name))
    print(f"[root] selected-root tips found for {selected_accession}: {n_out}")
    return tree_lib.root_outgroup(tree, is_selected_tip)


def apply_rooting(tree, method, metadata):
    """Apply one of the configured rooting strategies."""
    if method == "archaea_outgroup":
        return root_with_archaea_outgroup(tree, metadata), "Archaea-outgroup-rooted"

    if method == "unrooted":
        # Do not alter the tree. We mark the interpretation in output tables.
        return tree, "input/unrooted orientation preserved"

    if method.startswith("selected_root_"):
        return root_with_selected_sample(tree, SELECTED_ROOT_ACCESSION), f"rooted on {SELECTED_ROOT_ACCESSION}"

    if method == "midpoint":
        return tree_lib.root_midpoint(tree), "midpoint-rooted"

    raise ValueError(f"Unknown rooting method: {method}")


def leaf_names(tree):
    return {t.name for t in tree.get_terminals()}


def canonical_split(desc, all_leaves):
    """Canonicalize an unrooted split by using the smaller side."""
    desc = frozenset(desc)
    comp = frozenset(all_leaves - desc)
    if len(desc) > len(comp):
        return comp
    if len(desc) < len(comp):
        return desc
    return min(desc, comp, key=lambda x: tuple(sorted(x)))


def split_lengths(tree, leaves):
    """
    Return a dictionary: canonical unrooted split -> branch length.

    Trivial splits are excluded.
    """
    all_leaves = frozenset(leaves)
    n = len(all_leaves)
    out = {}

    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            continue

        desc = frozenset(t.name for t in node.get_terminals() if t.name in all_leaves)
        if len(desc) < 2 or len(desc) >= n - 1:
            continue

        key = canonical_split(desc, all_leaves)
        out[key] = node.branch_length or 0.0

    return out


def bipartitions(tree, leaves):
    """Return the set of unrooted non-trivial bipartitions."""
    return set(split_lengths(tree, leaves).keys())


def robinson_foulds(tree_a, tree_b, common_leaves):
    """Return RF, normalized RF, bipartition counts, shared count, and max RF."""
    A = bipartitions(tree_a, common_leaves)
    B = bipartitions(tree_b, common_leaves)
    shared = A & B
    rf = len(A - B) + len(B - A)
    n = len(common_leaves)
    max_rf = 2 * (n - 3) if n > 3 else 1
    return rf, rf / max_rf, len(A), len(B), len(shared), max_rf


def weighted_rf_branch_score(tree_a, tree_b, common_leaves):
    """
    Branch score / weighted RF distance.

    Uses the union of unrooted splits. Missing split length is treated as 0.
    """
    A = split_lengths(tree_a, common_leaves)
    B = split_lengths(tree_b, common_leaves)
    keys = set(A) | set(B)
    ss = 0.0
    for k in keys:
        ss += (A.get(k, 0.0) - B.get(k, 0.0)) ** 2
    return math.sqrt(ss)


def rooted_clusters(tree, leaves):
    """Return rooted clusters induced by internal nodes."""
    all_leaves = frozenset(leaves)
    n = len(all_leaves)
    clusters = set()

    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            continue
        desc = frozenset(t.name for t in node.get_terminals() if t.name in all_leaves)
        if len(desc) < 2 or len(desc) >= n:
            continue
        clusters.add(desc)

    return clusters


def rooted_cluster_distance(tree_a, tree_b, common_leaves):
    """Rooted cluster distance: symmetric difference of rooted clusters."""
    A = rooted_clusters(tree_a, common_leaves)
    B = rooted_clusters(tree_b, common_leaves)
    shared = A & B
    dist = len(A - B) + len(B - A)
    max_possible = 2 * (len(common_leaves) - 2) if len(common_leaves) > 2 else 1
    return dist, dist / max_possible, len(A), len(B), len(shared), max_possible


def parent_map(tree):
    parents = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parents[child] = clade
    return parents


def patristic_matrix(tree, names):
    """Compute pairwise tip-to-tip path length matrix."""
    names = list(names)
    name_to_tip = {t.name: t for t in tree.get_terminals()}
    parents = parent_map(tree)
    depths = tree.depths()

    def ancestors_to_root(tip):
        chain = [tip]
        while chain[-1] in parents:
            chain.append(parents[chain[-1]])
        return chain

    anc = {n: ancestors_to_root(name_to_tip[n]) for n in names}
    anc_sets = {n: set(anc[n]) for n in names}

    D = np.zeros((len(names), len(names)), dtype=float)

    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            common = anc_sets[a] & anc_sets[b]
            mrca = max(common, key=lambda x: depths.get(x, 0.0))
            d = depths[name_to_tip[a]] + depths[name_to_tip[b]] - 2 * depths[mrca]
            D[i, j] = D[j, i] = d

    return D


def pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x - x.mean()
    y = y - y.mean()
    denom = math.sqrt(float((x * x).sum() * (y * y).sum()))
    return float((x * y).sum() / denom) if denom else float("nan")


def rank_array(x):
    """Average ranks with no scipy dependency."""
    x = np.asarray(x)
    order = np.argsort(x)
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and x[order[j]] == x[order[i]]:
            j += 1
        rank = (i + j - 1) / 2.0 + 1
        ranks[order[i:j]] = rank
        i = j
    return ranks


def mantel_test(D1, D2, n_perm=N_MANTEL_PERMUTATIONS, seed=RANDOM_SEED):
    """Mantel test on upper-triangle distances."""
    rng = np.random.default_rng(seed)
    n = D1.shape[0]
    tri = np.triu_indices(n, k=1)
    x = D1[tri]
    y = D2[tri]

    r_pearson = pearson(x, y)
    r_spearman = pearson(rank_array(x), rank_array(y))

    pearson_extreme = 0
    spearman_extreme = 0
    rx = rank_array(x)

    for _ in range(n_perm):
        perm = rng.permutation(n)
        y_perm = D2[np.ix_(perm, perm)][tri]
        rp = pearson(x, y_perm)
        rs = pearson(rx, rank_array(y_perm))
        if abs(rp) >= abs(r_pearson):
            pearson_extreme += 1
        if abs(rs) >= abs(r_spearman):
            spearman_extreme += 1

    return (
        r_pearson,
        (pearson_extreme + 1) / (n_perm + 1),
        r_spearman,
        (spearman_extreme + 1) / (n_perm + 1),
    )


def quartet_topology_from_dist(D, a, b, c, d):
    """
    Infer quartet topology from pairwise distances.

    For additive tree distances, the smallest of:
      d(a,b)+d(c,d), d(a,c)+d(b,d), d(a,d)+d(b,c)
    gives the split. Ties are treated as unresolved.
    """
    vals = [
        D[a, b] + D[c, d],
        D[a, c] + D[b, d],
        D[a, d] + D[b, c],
    ]

    m = min(vals)
    winners = [i for i, v in enumerate(vals) if abs(v - m) < 1e-10]
    return winners[0] if len(winners) == 1 else -1


def quartet_distance(D1, D2):
    """Exact quartet distance by enumerating all C(n,4) quartets."""
    n = D1.shape[0]
    total = 0
    different = 0
    same = 0
    unresolved = 0

    for a, b, c, d in itertools.combinations(range(n), 4):
        total += 1
        q1 = quartet_topology_from_dist(D1, a, b, c, d)
        q2 = quartet_topology_from_dist(D2, a, b, c, d)

        if q1 == -1 or q2 == -1:
            unresolved += 1
        elif q1 == q2:
            same += 1
        else:
            different += 1

    resolved = same + different
    return {
        "quartets_total": total,
        "quartets_same": same,
        "quartets_different": different,
        "quartets_unresolved_or_tied": unresolved,
        "quartet_distance_normalized_all": different / total if total else float("nan"),
        "quartet_distance_normalized_resolved": different / resolved if resolved else float("nan"),
    }


def mrca_clade(tree, leaf_names):
    name_to_node = {t.name: t for t in tree.get_terminals()}
    targets = [name_to_node[n] for n in leaf_names if n in name_to_node]

    if len(targets) == 0:
        return None

    parents = parent_map(tree)

    def ancestors(t):
        chain = [t]
        while chain[-1] in parents:
            chain.append(parents[chain[-1]])
        return chain

    common = set(ancestors(targets[0]))
    for t in targets[1:]:
        common &= set(ancestors(t))

    if not common:
        return None

    depths = tree.depths()
    return max(common, key=lambda n: depths.get(n, 0.0))


def is_monophyletic(tree, leaf_names):
    """True iff the MRCA of leaf_names contains no other leaves."""
    leaf_names = list(leaf_names)
    if len(leaf_names) < 2:
        return None

    mrca = mrca_clade(tree, leaf_names)
    if mrca is None:
        return None

    mrca_leaves = {t.name for t in mrca.get_terminals()}
    requested = set(leaf_names)
    return len(mrca_leaves - requested) == 0


def monophyly_table(tree16, treeaai, meta, common, group_col, root_method):
    """Build per-group monophyly table for one taxonomic column."""
    id_col = meta.columns[0]
    submeta = meta[meta[id_col].isin(common)].copy()

    rows = []
    for grp, sub in submeta.groupby(group_col):
        gcas = list(sub[id_col])
        if len(gcas) < 2:
            rows.append({
                "rooting_method": root_method,
                "rank": group_col,
                group_col: grp,
                "group": grp,
                "n_genomes": len(gcas),
                "monophyletic_16S": "N/A (n<2)",
                "monophyletic_FastAAI": "N/A (n<2)",
                "same_result": "N/A",
            })
            continue

        m16 = is_monophyletic(tree16, gcas)
        mai = is_monophyletic(treeaai, gcas)

        rows.append({
            "rooting_method": root_method,
            "rank": group_col,
            group_col: grp,
            "group": grp,
            "n_genomes": len(gcas),
            "monophyletic_16S": "yes" if m16 else "no",
            "monophyletic_FastAAI": "yes" if mai else "no",
            "same_result": "yes" if m16 == mai else "no",
        })

    return pd.DataFrame(rows).sort_values(["n_genomes", "group"], ascending=[False, True])


def evaluate_one_rooting(method, meta):
    """Compute all metrics for one root method."""
    print(f"\n=== Rooting method: {method} ===")

    t16s = tree_lib.parse_newick(TREE_16S)
    taai = tree_lib.parse_newick(TREE_AAI)

    t16s, root_note_16 = apply_rooting(t16s, method, meta)
    taai, root_note_aai = apply_rooting(taai, method, meta)

    Path(TREE_OUTDIR).mkdir(parents=True, exist_ok=True)
    with open(f"{TREE_OUTDIR}/16S_{method}.treefile", "w") as f:
        f.write(write_newick(t16s) + "\n")
    with open(f"{TREE_OUTDIR}/FastAAI_{method}.treefile", "w") as f:
        f.write(write_newick(taai) + "\n")

    leaves_16s = leaf_names(t16s)
    leaves_aai = leaf_names(taai)
    common = sorted(leaves_16s & leaves_aai)

    rf, nrf, nA, nB, shared, max_rf = robinson_foulds(t16s, taai, common)
    rcd, nrcd, rcA, rcB, rc_shared, rc_max = rooted_cluster_distance(t16s, taai, common)
    branch_score = weighted_rf_branch_score(t16s, taai, common)

    D16 = patristic_matrix(t16s, common)
    DAAI = patristic_matrix(taai, common)

    # Save matrices only once for the biologically preferred Archaea-rooted version.
    if method == "archaea_outgroup":
        pd.DataFrame(D16, index=common, columns=common).to_csv(
            f"{OUTDIR}/patristic_distance_matrix_16S_archaea_outgroup.tsv", sep="\t"
        )
        pd.DataFrame(DAAI, index=common, columns=common).to_csv(
            f"{OUTDIR}/patristic_distance_matrix_FastAAI_archaea_outgroup.tsv", sep="\t"
        )

    mantel_p, mantel_pp, mantel_s, mantel_sp = mantel_test(D16, DAAI)
    q = quartet_distance(D16, DAAI)

    mono_domain = monophyly_table(t16s, taai, meta, common, "domain", method)
    mono_kingdom = monophyly_table(t16s, taai, meta, common, "kingdom", method)
    mono_phylum = monophyly_table(t16s, taai, meta, common, "phylum", method)

    # Headline group counts only for multi-genome groups.
    multi_kingdom = mono_kingdom[
        (~mono_kingdom["monophyletic_16S"].str.startswith("N/A")) &
        (~mono_kingdom["group"].astype(str).isin(["Unknown", "nan", "NaN", "None", ""]))
    ].copy()

    summary = {
        "rooting_method": method,
        "root_note_16S": root_note_16,
        "root_note_FastAAI": root_note_aai,
        "n_common_leaves": len(common),
        "RF_distance_unrooted": rf,
        "max_possible_RF": max_rf,
        "normalised_RF": round(nrf, 6),
        "bipartitions_16S": nA,
        "bipartitions_FastAAI": nB,
        "shared_bipartitions": shared,
        "rooted_cluster_distance": rcd,
        "max_possible_rooted_cluster_distance": rc_max,
        "normalised_rooted_cluster_distance": round(nrcd, 6),
        "rooted_clusters_16S": rcA,
        "rooted_clusters_FastAAI": rcB,
        "shared_rooted_clusters": rc_shared,
        "weighted_RF_branch_score_distance": round(branch_score, 8),
        "quartets_total": q["quartets_total"],
        "quartets_same": q["quartets_same"],
        "quartets_different": q["quartets_different"],
        "quartets_unresolved_or_tied": q["quartets_unresolved_or_tied"],
        "normalised_quartet_distance_all": round(q["quartet_distance_normalized_all"], 6),
        "normalised_quartet_distance_resolved": round(q["quartet_distance_normalized_resolved"], 6),
        "mantel_pearson_r": round(mantel_p, 6),
        "mantel_pearson_p": round(mantel_pp, 6),
        "mantel_spearman_rho": round(mantel_s, 6),
        "mantel_spearman_p": round(mantel_sp, 6),
        "multi_genome_kingdoms": len(multi_kingdom),
        "kingdoms_monophyletic_16S": int((multi_kingdom["monophyletic_16S"] == "yes").sum()),
        "kingdoms_monophyletic_FastAAI": int((multi_kingdom["monophyletic_FastAAI"] == "yes").sum()),
        "kingdoms_monophyletic_in_both": int(
            ((multi_kingdom["monophyletic_16S"] == "yes") &
             (multi_kingdom["monophyletic_FastAAI"] == "yes")).sum()
        ),
    }

    return summary, mono_domain, mono_kingdom, mono_phylum


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(TREE_OUTDIR, exist_ok=True)
    os.makedirs(f"{OUTDIR}/tables", exist_ok=True)

    sep = "," if METADATA.lower().endswith(".csv") else "\t"
    meta = pd.read_csv(METADATA, sep=sep)
    meta = harmonize_metadata(meta)
    meta.to_csv(f"{OUTDIR}/metadata_harmonized.tsv", sep="\t", index=False)

    all_summaries = []
    mono_domains = []
    mono_kingdoms = []
    mono_phyla = []

    for method in ROOTING_METHODS:
        summary, mdomain, mkingdom, mphylum = evaluate_one_rooting(method, meta)
        all_summaries.append(summary)
        mono_domains.append(mdomain)
        mono_kingdoms.append(mkingdom)
        mono_phyla.append(mphylum)

        # Root-method-specific metrics table.
        pd.DataFrame([summary]).to_csv(
            f"{OUTDIR}/tables/tree_comparison_metrics_{method}.tsv",
            sep="\t",
            index=False,
        )
        mdomain.to_csv(f"{OUTDIR}/tables/per_domain_monophyly_{method}.tsv", sep="\t", index=False)
        mkingdom.to_csv(f"{OUTDIR}/tables/per_kingdom_monophyly_{method}.tsv", sep="\t", index=False)
        mphylum.to_csv(f"{OUTDIR}/tables/per_phylum_monophyly_{method}.tsv", sep="\t", index=False)

    metrics_df = pd.DataFrame(all_summaries)
    metrics_df.to_csv(f"{OUTDIR}/tree_comparison_metrics_by_rooting.tsv", sep="\t", index=False)

    # Keep v2-compatible one-row table using the biologically preferred Archaea-rooted results.
    metrics_df[metrics_df["rooting_method"] == "archaea_outgroup"].to_csv(
        f"{OUTDIR}/tree_comparison_metrics.tsv",
        sep="\t",
        index=False,
    )

    pd.concat(mono_domains, ignore_index=True).to_csv(
        f"{OUTDIR}/per_domain_monophyly_by_rooting.tsv", sep="\t", index=False
    )
    pd.concat(mono_kingdoms, ignore_index=True).to_csv(
        f"{OUTDIR}/per_kingdom_monophyly_by_rooting.tsv", sep="\t", index=False
    )
    pd.concat(mono_phyla, ignore_index=True).to_csv(
        f"{OUTDIR}/per_phylum_monophyly_by_rooting.tsv", sep="\t", index=False
    )

    # v2-compatible copies for Archaea-outgroup tables.
    pd.concat(mono_domains, ignore_index=True).query("rooting_method == 'archaea_outgroup'").to_csv(
        f"{OUTDIR}/per_domain_monophyly.tsv", sep="\t", index=False
    )
    pd.concat(mono_kingdoms, ignore_index=True).query("rooting_method == 'archaea_outgroup'").to_csv(
        f"{OUTDIR}/per_kingdom_monophyly.tsv", sep="\t", index=False
    )
    pd.concat(mono_phyla, ignore_index=True).query("rooting_method == 'archaea_outgroup'").to_csv(
        f"{OUTDIR}/per_phylum_monophyly.tsv", sep="\t", index=False
    )

    print("\n[write] output/tree_comparison_metrics_by_rooting.tsv")
    print(metrics_df.to_string(index=False))


if __name__ == "__main__":
    main()
