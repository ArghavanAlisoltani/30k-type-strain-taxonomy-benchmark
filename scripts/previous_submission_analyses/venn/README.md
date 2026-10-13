# 4-way Venn — pairs failing the species-boundary threshold

## What this is
For each of the four methods in `exact_pair_status_table.csv`, we take the pairs
labelled `FAIL_THRESHOLD` (the tool ran successfully but the score fell below
the species-level cutoff) and compute *exact* pair-level intersections
between the four sets.

`FAIL_INITIAL` is excluded on purpose — those pairs had no comparable value
produced (e.g. no shared 16S hit), so they aren't "disagreements", they're
"couldn't measure".

## Sanity check — the four totals reproduce your expected numbers
| Method            | Expected | Found `FAIL_THRESHOLD` |
|-------------------|---------:|-----------------------:|
| 16S               |     349  |               **349**  |
| Mash              |   1,297  |             **1,297**  |
| FastANI           |   2,026  |             **2,026**  |
| FastAAI / Jaccard |   1,271  |             **1,271**  |

Union of all four sets = **2,815** unique pairs.

## Region breakdown (15 non-empty regions)

| Region                                                    | Pairs |
|-----------------------------------------------------------|------:|
| 16S only                                                  |   197 |
| Mash only                                                 |     0 |
| FastANI only                                              |   606 |
| FastAAI/Jaccard only                                      |   561 |
| 16S ∩ Mash                                                |     0 |
| 16S ∩ FastANI                                             |    25 |
| 16S ∩ FastAAI/Jaccard                                     |    28 |
| Mash ∩ FastANI                                            |   683 |
| Mash ∩ FastAAI/Jaccard                                    |     2 |
| FastANI ∩ FastAAI/Jaccard                                 |    97 |
| 16S ∩ Mash ∩ FastANI                                      |    33 |
| 16S ∩ Mash ∩ FastAAI/Jaccard                              |     1 |
| 16S ∩ FastANI ∩ FastAAI/Jaccard                           |     4 |
| Mash ∩ FastANI ∩ FastAAI/Jaccard                          |   517 |
| **All four**                                              |    61 |

Two things jump out biologically:
- **Mash never disagrees alone.** Every Mash-flagged pair is also flagged by
  at least one other tool.
- **Mash & FastANI overlap heavily (683 + 517 + 33 + 61 = 1,294 shared),**
  which is the expected behaviour of two nucleotide-based ANI-like methods.
- **16S dissents most independently** (197 pairs seen only by 16S) — the
  16S signal captures divergence that whole-genome methods miss.

## Files

| File | What it is |
|---|---|
| `venn_fail_threshold.png`  | 4-way Venn (300 dpi raster) |
| `venn_fail_threshold.pdf`  | 4-way Venn (vector) |
| `venn_region_counts.csv`   | Size of each of the 15 non-empty regions |
| `venn_pair_regions.csv`    | Every one of the 2,815 pair_ids labelled with its region |
| `venn_fail_threshold.py`   | Reproducible script (edit + rerun) |

## How to reproduce
```bash
pip install pandas matplotlib venn
python venn_fail_threshold.py
```
The script reads `exact_pair_status_table.csv` (path defined at the top),
verifies the four totals against your expected numbers, then writes the
figure and both tables.

