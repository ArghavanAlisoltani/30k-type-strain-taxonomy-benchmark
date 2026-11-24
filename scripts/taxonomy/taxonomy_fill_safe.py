#!/usr/bin/env python3
"""Propose missing upper-rank assignments in a GenID dictionary.

No network access. Dry-run by default. Existing assigned labels are preserved.
"""
from __future__ import annotations

import argparse
import ast
import csv
from collections import defaultdict
from pathlib import Path
from pprint import pformat

RANKS = ("kingdom", "phylum", "class", "order", "family")
MISSING = {"", "none", "na", "n/a", "nan", "null"}


def clean(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"Expected taxon string or None, got {type(value).__name__}")
    result = value.strip()
    return result if result and result.casefold() not in MISSING else None


def load_genid(path: Path) -> dict[str, dict[str, object]]:
    content = path.read_text(encoding="utf-8-sig").strip()
    if content.startswith("GenID="):
        content = content[len("GenID="):].strip()
    data = ast.literal_eval(content)  # literal values only; never execute input
    if not isinstance(data, dict):
        raise ValueError("GenID must contain a Python dictionary literal")
    for accession, record in data.items():
        if not isinstance(accession, str) or not isinstance(record, dict):
            raise ValueError("Every GenID key must be a string with a dictionary value")
        for rank in (*RANKS, "genus", "species"):
            clean(record.get(rank))
    return data


def lookup_keys(record: dict[str, object]) -> list[tuple[str, tuple[str, ...]]]:
    genus = clean(record.get("genus"))
    species = clean(record.get("species"))
    keys: list[tuple[str, tuple[str, ...]]] = []
    if genus and species:
        keys.append(("species_lookup", (genus.casefold(), species.casefold())))
    if genus:
        keys.append(("genus_lookup", (genus.casefold(),)))
    return keys


def proposed_fills(data: dict[str, dict[str, object]]) -> list[dict[str, str]]:
    # Build lookup labels from input records without cascading proposed values.
    indexes: dict[str, dict[tuple[str, ...], dict[str, set[str]]]] = {
        "species_lookup": defaultdict(lambda: defaultdict(set)),
        "genus_lookup": defaultdict(lambda: defaultdict(set)),
    }
    for record in data.values():
        for kind, key in lookup_keys(record):
            for rank in RANKS:
                value = clean(record.get(rank))
                if value:
                    indexes[kind][key][rank].add(value)

    changes: list[dict[str, str]] = []
    for accession, record in data.items():
        for rank in RANKS:
            if clean(record.get(rank)) is not None:
                continue
            # A disagreement at a more-specific reference level is not
            # overridden with a less-specific proposal.
            for kind, key in lookup_keys(record):
                possible = indexes[kind].get(key, {}).get(rank, set())
                if len(possible) > 1:
                    break  # conflicting references: do not impute
                if len(possible) == 1:
                    taxon = next(iter(possible))  # IMPORTANT: the label, not count
                    changes.append({"accession": accession, "rank": rank,
                                    "previous_value": "", "proposed_value": taxon,
                                    "evidence": kind, "status": "inferred_unverified"})
                    break
    return changes


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True, help="Path to original legacy GenID dictionary")
    ap.add_argument("--output", type=Path, help="New proposed GenID output path")
    ap.add_argument("--audit", type=Path, help="New proposed changes TSV path")
    ap.add_argument("--apply", action="store_true", help="Write output and audit; dry-run otherwise")
    args = ap.parse_args()
    data = load_genid(args.input)
    changes = proposed_fills(data)
    print(f"Input assemblies: {len(data):,}; upper-rank proposals: {len(changes):,}")
    print("No genus/species labels are filled or modified.")
    if not args.apply:
        print("DRY RUN: no files written. Use --apply with --output and --audit.")
        return
    if not args.output or not args.audit:
        ap.error("--apply requires --output and --audit")
    paths = [args.input.resolve(), args.output.resolve(), args.audit.resolve()]
    if len(set(paths)) != len(paths):
        ap.error("Input, output, and audit must be three distinct paths")
    if args.output.exists() or args.audit.exists():
        ap.error("Refusing to overwrite existing output or audit")
    for item in changes:
        data[item["accession"]][item["rank"]] = item["proposed_value"]
    # Explicit new files only; never change the original input.
    with args.output.open("x", encoding="utf-8") as fh:
        fh.write("GenID=" + pformat(data, sort_dicts=False) + "\n")
    with args.audit.open("x", encoding="utf-8", newline="") as fh:
        fields = ("accession", "rank", "previous_value", "proposed_value", "evidence", "status")
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(changes)
    print(f"Wrote proposed outputs: {args.output} and {args.audit}")
    print("IMPORTANT: these proposed fills require validation before scientific use.")


if __name__ == "__main__":
    main()
