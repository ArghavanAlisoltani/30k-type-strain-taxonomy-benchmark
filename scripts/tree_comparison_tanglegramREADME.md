# 16S vs FastAAI tree-comparison bundle

## Tanglegram leaf ordering (16S is the fixed reference)

The 16S (left) tree in both tanglegram scripts is the **reference**: its leaf
order is fixed to the natural post-rooting order — the exact order the circular
branch-length tree uses — and it is **never rotated**. Only the FastAAI (right)
tree is rotated to minimize line crossings. This keeps the 16S display identical
to the circular tree (e.g. the Archaea block stays contiguous and per-kingdom
blocks stay together) and leaves the 16S structure untouched.

Note: some groups (e.g. kingdom *Nanobdellati*) are **not monophyletic** in the
16S tree, so their tips are only *adjacent by rotation*, not a true clade. Fixing
the 16S order to the circular-tree order is what makes that adjacency consistent
between the two figures. The trade-off is more crossing lines on the FastAAI
side, which is expected when the reference tree is held fixed.

## Latest revision

1. **New metadata** `inputs/metadata_with_official_phyla_and_woese_groups.csv`
   is now the source table (replaces `metadata_with_kingdom.tsv`, which is kept
   only for provenance). It provides `official_phylum` and a new
   `carl_woese_historical_group` column, and — importantly — **assigns a kingdom
   to every genome**.
   - *Why the old figures had an "Unknown" kingdom:* the previous
     `metadata_with_kingdom.tsv` had literal `NA` in the kingdom column for three
     genomes (*Dictyoglomus thermophilum*, *Desulfobacter postgatei*,
     *Halobacterium salinarum*). The harmonizer turned `NA -> "Unknown"`, so a
     grey "Unknown" band appeared in the strips. The new CSV fills those in, so
     there is no "Unknown" kingdom anymore. (It also fixes *Methanobacteriati*,
     which now reads as monophyletic in both trees in the step-5 heatmap.)
2. **Fourth metadata strip "W" = Carl Woese historical group.** The step-4
   strip tanglegrams now show D / K / P / **W** strips on each side, where W is
   `carl_woese_historical_group`, colored by `WOESE_COLORS`.
3. **Black tip labels.** Tip-label text is now black (`TIP_LABEL_COLOR`); it is
   no longer tinted by domain. (Connecting lines remain domain-colored, set by
   `LINE_COLOR_BY`.)
4. **User-supplied palettes** for domain / kingdom / phylum are applied verbatim
   in `tanglegram_colors.py`.

## What `tanglegram_colors.py` controls

- `DOMAIN_COLORS`, `KINGDOM_COLORS`, `PHYLUM_COLORS`, `WOESE_COLORS` — explicit
  `name -> hex` maps. Edit a hex to recolor that taxon everywhere.
- `STATUS_COLORS` — step-5 monophyly heatmap (neither / 16S / FastAAI / both).
- `TIP_LABEL_COLOR` — tip-label text color (black).
- `LINE_COLOR_BY` — column that colors the connecting lines (`"domain"` default;
  can be `"kingdom"`, `"phylum"`, or `"woese"`).
- `MISSING_COLOR` — fallback for any value not listed or blank/NA.

## Recoloring workflow

```text
1. Open tanglegram_colors.py
2. Change the hex next to any domain / kingdom / phylum / Woese group (or a STATUS_COLORS entry)
3. Re-run:  python3 step3_...  python3 step4_...  python3 step5_...
```

## Input files

```text
inputs/16S_full.treefile
inputs/all_rnammer_16S.fasta
inputs/fastaai_nj.treefile
inputs/metadata_with_official_phyla_and_woese_groups.csv   <- current metadata
inputs/metadata_with_kingdom.tsv                           <- kept for provenance
```

## Run everything

```bash
bash run_all.sh
```

## Main workflow

### Step 1 — `python3 step1_reduce_16s.py`
Reduces the 110-tip 16S tree to one representative 16S copy per genome.

### Step 2 — `python3 step2_compare_all_rootings.py`
Computes metrics and monophyly under four rooting strategies.

Main outputs:

```text
output/tree_comparison_metrics.tsv
output/tree_comparison_metrics_by_rooting.tsv
output/per_domain_monophyly.tsv
output/per_kingdom_monophyly.tsv
output/per_phylum_monophyly.tsv
output/per_domain_monophyly_by_rooting.tsv
output/per_kingdom_monophyly_by_rooting.tsv
output/per_phylum_monophyly_by_rooting.tsv
```

### Step 3 — `python3 step3_tanglegrams_all_rootings.py`
Standard tanglegrams for all rooting methods. **Line/label colors come from
`tanglegram_colors.py` (`LINE_COLOR_BY`).**

```text
output/tanglegrams/tanglegram_archaea_outgroup.png
output/tanglegrams/tanglegram_unrooted.png
output/tanglegrams/tanglegram_selected_root_GCA_021654395_1.png
output/tanglegrams/tanglegram_midpoint.png
```

v2-compatible copies: `output/tanglegram_archaea_rooted.png` / `.pdf`.

### Step 4 — `python3 step4_tanglegram_metadata_strips_all_rootings.py`
Tanglegrams with Domain/Kingdom/Phylum/**Woese** heatmap strips (D/K/P/W). **All
four strip colors come from the manual maps in `tanglegram_colors.py`.**

```text
output/tanglegrams_metadata_strips/tanglegram_metadata_strips_archaea_outgroup.png
output/tanglegrams_metadata_strips/tanglegram_metadata_strips_unrooted.png
output/tanglegrams_metadata_strips/tanglegram_metadata_strips_selected_root_GCA_021654395_1.png
output/tanglegrams_metadata_strips/tanglegram_metadata_strips_midpoint.png
```

v2-compatible copies: `output/tanglegram_with_metadata_strips_archaea_rooted.png` / `.pdf`.

### Step 5 — `python3 step5_metric_summary_figures.py`
Summary plots comparing metrics across rooting methods. **Heatmap status colors
come from `STATUS_COLORS`; kingdom labels are tinted with `KINGDOM_COLORS`.**

```text
output/figures/rooting_metric_comparison.png
output/figures/mantel_correlation_by_rooting.png
output/figures/monophyly_kingdom_heatmap.png
output/figures/tanglegram_crossings_by_rooting.png
```

## Interpretation guide

- **RF distance** is the standard unrooted split comparison; unchanged by rerooting.
- **Rooted cluster distance** is root-sensitive — use it for root-dependent clade interpretation.
- **Quartet distance** summarizes four-taxon topology differences; less sensitive than RF.
- **Mantel correlation** compares branch-length-weighted pairwise patristic distances.
- **Monophyly tables** should be read with the root method in mind.
- **Archaea-outgroup rooting** is the recommended primary biological interpretation
  because both trees recover a clear Archaea/Bacteria split.

