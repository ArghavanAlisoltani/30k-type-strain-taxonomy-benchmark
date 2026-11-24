#!/usr/bin/env python3
"""Convert Mash `triangle` relaxed-PHYLIP lower triangle into three-column pairs."""
from __future__ import annotations

import argparse
import math
from pathlib import Path


def convert(source: Path, destination: Path) -> dict[str, int]:
    if source.resolve() == destination.resolve():
        raise ValueError('Input and output must differ')
    if not source.is_file() or destination.exists():
        raise ValueError('Input must exist and output must not already exist')
    temp = destination.with_name(destination.name + '.partial')
    if temp.exists():
        raise ValueError(f'Remove or inspect previous partial output first: {temp}')
    ids: list[str] = []
    seen_ids: set[str] = set()
    count: int | None = None
    edges = 0
    try:
        with source.open(encoding='utf-8') as inp, temp.open('x', encoding='utf-8') as out:
            for n, line in enumerate(inp, 1):
                if not line.strip():
                    continue
                cells = line.split()
                if count is None:
                    if len(cells) != 1 or not cells[0].isdigit():
                        raise ValueError(f'{source}:{n}: first line must be number of genomes')
                    count = int(cells[0])
                    if count < 1:
                        raise ValueError(f'{source}:{n}: number of genomes must be positive')
                    continue
                index = len(ids)
                if index >= count:
                    raise ValueError(f'{source}:{n}: more matrix rows than declared')
                if len(cells) != index + 1:
                    raise ValueError(f'{source}:{n}: expected name and {index} distances, got {len(cells)} columns')
                genome = cells[0]
                if genome in seen_ids:
                    raise ValueError(f'{source}:{n}: duplicate genome ID {genome}')
                for previous, raw in zip(ids, cells[1:]):
                    try:
                        value = float(raw)
                    except ValueError as exc:
                        raise ValueError(f'{source}:{n}: invalid numeric distance {raw}') from exc
                    if not math.isfinite(value) or not 0 <= value <= 1:
                        raise ValueError(f'{source}:{n}: Mash distance outside [0,1]: {raw}')
                    out.write(f'{genome}\t{previous}\t{raw}\n')
                    edges += 1
                ids.append(genome)
                seen_ids.add(genome)
        if count is None or len(ids) != count:
            raise ValueError(f'{source}: declared {count} genomes but found {len(ids)}')
        expected = count * (count - 1) // 2
        if edges != expected:
            raise ValueError(f'{source}: expected {expected} distances, found {edges}')
        temp.replace(destination)
        return {'genomes': count, 'pairs': edges}
    except BaseException:
        temp.unlink(missing_ok=True)
        raise


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    try:
        stats = convert(a.input, a.output)
    except (ValueError, OSError) as e:
        p.exit(1, f'ERROR: {e}\n')
    print(f"Converted {stats['genomes']} genomes and {stats['pairs']} unique pairs: {a.output}")


if __name__ == '__main__':
    main()
