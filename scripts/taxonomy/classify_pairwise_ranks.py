#!/usr/bin/env python3
"""Assign pairwise scores to exclusive taxonomic categories from a frozen TSV.

Candidate documented implementation of the colleague's pseudocode, not a verified
historical source script. Uses explicit lower-rank disagreement; missing labels
never count as disagreement. Preserves every input row and diagnostic status.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

RANKS = ("species", "genus", "family", "order", "class", "phylum", "kingdom", "domain")
MISSING = {"", "none", "na", "n/a", "nan", "null"}


def norm(x: object) -> str | None:
    s = str(x).strip() if x is not None else ""
    return s if s and s.casefold() not in MISSING else None


def classify(a: dict[str, str], b: dict[str, str]) -> tuple[str, str]:
    """Return deepest proven exclusive category and a diagnostic status.

    A category higher than species requires a known difference at its immediate
    child rank. Unknown lower labels imply unresolved status, not 'different'.
    """
    # Consistent lineage: equal lower-rank assignments cannot have contradictory
    # assigned ancestors. Prevent silent misclassification of inconsistent rows.
    for i, child in enumerate(RANKS):
        child_a, child_b = norm(a.get(child)), norm(b.get(child))
        if child_a is None or child_b is None or child_a != child_b:
            continue
        for ancestor in RANKS[i + 1:]:
            higher_a, higher_b = norm(a.get(ancestor)), norm(b.get(ancestor))
            if higher_a is not None and higher_b is not None and higher_a != higher_b:
                return "", "inconsistent_lineage"
    if all(norm(a.get(r)) == norm(b.get(r)) and norm(a.get(r)) is not None
           for r in ("genus", "species")):
        return "species", "classified"
    for i in range(1, len(RANKS)):
        rank, child = RANKS[i], RANKS[i - 1]
        x, y = norm(a.get(rank)), norm(b.get(rank))
        lower_x, lower_y = norm(a.get(child)), norm(b.get(child))
        if x is None or y is None or x != y:
            continue
        if lower_x is None or lower_y is None:
            return "", "unresolved_lower_rank"
        if lower_x != lower_y:
            return rank, "classified"
    return "", "not_classifiable"


def read_taxonomy(path: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "accession" not in reader.fieldnames:
            raise ValueError("Taxonomy TSV needs an 'accession' header")
        missing = set(RANKS) - set(reader.fieldnames)
        if missing:
            raise ValueError("Taxonomy TSV lacks ranks: " + ", ".join(sorted(missing)))
        for row in reader:
            key = norm(row["accession"])
            if not key or key in result:
                raise ValueError(f"Empty/duplicate assembly accession: {key!r}")
            result[key] = row
    return result


def run(taxonomy: Path, scores: Path, output: Path) -> dict[str, int]:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    if output.resolve() in (taxonomy.resolve(), scores.resolve()):
        raise ValueError("Output must differ from input paths")
    taxa = read_taxonomy(taxonomy)
    stats: dict[str, int] = {}
    # Read/validate input header before creating output.
    with scores.open(encoding="utf-8-sig", newline="") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        needed = {"genome1", "genome2", "score"}
        if not reader.fieldnames or not needed.issubset(reader.fieldnames):
            raise ValueError("Scores TSV needs genome1, genome2, score headers")
        fieldnames = list(reader.fieldnames)
    with scores.open(encoding="utf-8-sig", newline="") as inp, output.open("x", encoding="utf-8", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=[*fieldnames, "exclusive_rank", "status"], delimiter="\t")
        writer.writeheader()
        for row in reader:
            a, b = row["genome1"], row["genome2"]
            if a == b:
                rank, status = "", "self_pair"
            elif a not in taxa or b not in taxa:
                rank, status = "", "unknown_accession"
            elif norm(row["score"]) is None:
                rank, status = "", "missing_numeric_score"
            else:
                try:
                    score = float(row["score"])
                except ValueError:
                    rank, status = "", "invalid_numeric_score"
                else:
                    rank, status = classify(taxa[a], taxa[b]) if math.isfinite(score) else ("", "invalid_numeric_score")
            writer.writerow({**row, "exclusive_rank": rank, "status": status})
            stats[status] = stats.get(status, 0) + 1
    return stats


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--taxonomy", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    for status, n in sorted(run(args.taxonomy, args.scores, args.output).items()):
        print(f"{status}: {n}")
    print("Candidate implementation only. Compare with historical rank-counts before publication.")


if __name__ == "__main__":
    main()
