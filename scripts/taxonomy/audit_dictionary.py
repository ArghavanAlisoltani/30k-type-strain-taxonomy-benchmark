#!/usr/bin/env python3
"""Audit a legacy GenID dictionary against a frozen accession manifest and taxonomy TSV.

Only use a trusted-size local dictionary file. It is parsed as a Python literal, never executed.
Dry-run by default. --export produces a subset TSV aligned to the manifest without
changing the source dictionary or authoritative published taxonomy.
"""
from __future__ import annotations

import argparse
import ast
import csv
import re
from collections import Counter
from pathlib import Path

RANKS = ('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
ACCESSION = re.compile(r'^GCA_\d{9}\.\d+$')


def load_dictionary(path: Path) -> dict[str, dict[str, str | None]]:
    content = path.read_text(encoding='utf-8-sig').strip()
    if content.startswith('GenID='):
        content = content[len('GenID='):].strip()
    data = ast.literal_eval(content)
    if not isinstance(data, dict):
        raise ValueError('Expected GenID={...} dictionary literal')
    for accession, fields in data.items():
        if not isinstance(accession, str) or not ACCESSION.fullmatch(accession):
            raise ValueError(f'Invalid accession: {accession!r}')
        if not isinstance(fields, dict) or set(fields) != set(RANKS):
            raise ValueError(f'{accession}: expected exactly eight taxonomy ranks')
        if any(value is not None and not isinstance(value, str) for value in fields.values()):
            raise ValueError(f'{accession}: ranks must contain string or None')
    return data


def load_manifest(path: Path) -> list[str]:
    ids = [x.strip() for x in path.read_text(encoding='utf-8').splitlines() if x.strip()]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError('Manifest is empty or contains duplicate accessions')
    if any(not ACCESSION.fullmatch(x) for x in ids):
        raise ValueError('Manifest contains an invalid GCA accession.version')
    return ids


def compare_reference(data: dict[str, dict], ids: list[str], ref: Path) -> tuple[int, list[dict[str, str]]]:
    with ref.open(newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError('Reference TSV has no header')
        key = 'accession' if 'accession' in reader.fieldnames else 'GCA' if 'GCA' in reader.fieldnames else None
        if not key:
            raise ValueError('Reference TSV needs accession or GCA column')
        ranks = [r for r in RANKS if r in reader.fieldnames]
        if not ranks:
            raise ValueError('No matching taxonomy rank columns in reference TSV')
        expected = set(ids)
        mismatches: list[dict[str, str]] = []
        seen = set()
        for row in reader:
            accession = row[key]
            if accession in seen:
                raise ValueError(f'Duplicate reference accession: {accession}')
            seen.add(accession)
            if accession not in expected:
                continue
            for rank in ranks:
                refval = (row[rank] or '').strip()
                orig = data[accession].get(rank)
                dictval = orig.strip() if isinstance(orig, str) else ''
                if dictval != refval:
                    mismatches.append({'accession': accession, 'rank': rank,
                                       'dictionary': dictval, 'reference': refval})
        missing = expected - seen
        if missing:
            raise ValueError(f'Reference TSV is missing {len(missing)} manifest IDs; example {sorted(missing)[0]}')
        return len(ranks), mismatches


def audit(dictionary: Path, manifest: Path, reference: Path | None = None,
          export: Path | None = None) -> dict:
    data = load_dictionary(dictionary)
    ids = load_manifest(manifest)
    selected = set(ids)
    extras, missing = sorted(set(data) - selected), sorted(selected - set(data))
    summary = {'dictionary_entries': len(data), 'manifest_entries': len(ids),
               'extra_dictionary_accessions': extras, 'missing_dictionary_accessions': missing}
    if missing:
        raise ValueError(f'Dictionary lacks {len(missing)} manifest accessions; examples: {missing[:5]}')
    if reference:
        rank_count, mismatches = compare_reference(data, ids, reference)
        summary.update({'reference_ranks_compared': rank_count,
                        'reference_field_mismatches': len(mismatches),
                        'mismatch_examples': mismatches[:8]})
    if export is not None:
        if export.exists() or export.resolve() in {dictionary.resolve(), manifest.resolve()} or (reference is not None and export.resolve() == reference.resolve()):
            raise ValueError('Export path must be a new file distinct from every input')
        with export.open('x', newline='', encoding='utf-8') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(('accession', *RANKS))
            for accession in ids:
                record = data[accession]
                writer.writerow((accession, *(record.get(rank) or '' for rank in RANKS)))
        summary['exported'] = str(export)
    return summary


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dictionary', required=True, type=Path)
    p.add_argument('--manifest', required=True, type=Path)
    p.add_argument('--reference-tsv', type=Path)
    p.add_argument('--export-tsv', type=Path)
    args = p.parse_args()
    try:
        result = audit(args.dictionary, args.manifest, args.reference_tsv, args.export_tsv)
    except (ValueError, SyntaxError, OSError, MemoryError) as exc:
        p.exit(1, f'ERROR: {exc}\n')
    import json
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
