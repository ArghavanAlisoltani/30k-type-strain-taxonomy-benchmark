#!/usr/bin/env python3
"""Run FastANI v1.33 in bounded query/reference blocks and merge without changing directional counts."""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

EXTENSIONS = ('.fna', '.fa', '.fasta', '.fna.gz', '.fa.gz', '.fasta.gz')


def genomes_from(path: Path, recursive: bool = False) -> list[str]:
    if not path.is_absolute() or not path.exists():
        raise ValueError(f'Genome source must exist and have an absolute path: {path}')
    if path.is_dir():
        iterator = path.rglob('*') if recursive else path.iterdir()
        files = sorted(p.resolve() for p in iterator if p.is_file() and p.name.endswith(EXTENSIONS))
        records = [str(p) for p in files]
    else:
        records = []
        for line_no, raw in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
            value = raw.strip()
            if not value:
                continue
            target = Path(value)
            if not target.is_absolute():
                raise ValueError(f'{path}:{line_no}: genome path must be absolute: {value}')
            records.append(str(target.resolve()))
    if not records:
        raise ValueError(f'No genome files found in {path}')
    if len(records) != len(set(records)):
        raise ValueError(f'Duplicate genome paths in {path}')
    for value in records:
        if '\t' in value or '\n' in value or '\r' in value:
            raise ValueError(f'Path contains tab or newline: {value}')
        p = Path(value)
        if not p.is_file() or not p.name.endswith(EXTENSIONS):
            raise ValueError(f'Invalid or missing genome FASTA: {p}')
    return records


def digest(records: list[str]) -> str:
    return hashlib.sha256(('\n'.join(records) + '\n').encode('utf-8')).hexdigest()


def blocks(records: list[str], size: int) -> list[list[str]]:
    return [records[i:i + size] for i in range(0, len(records), size)]


def jobs_for(nq: int, nr: int, symmetric: bool) -> list[tuple[int, int]]:
    return [(i, j) for i in range(nq) for j in range(i if symmetric else 0, nr)]


def key(i: int, j: int) -> str:
    return f'q{i:04d}_r{j:04d}'


def write_exact(path: Path, content: str) -> None:
    if path.exists():
        if path.read_text(encoding='utf-8') != content:
            raise ValueError(f'Existing file conflicts with run definition: {path}')
    else:
        path.write_text(content, encoding='utf-8')


def record(line: str, location: str, qs: set[str], rs: set[str]) -> tuple[str, str, str, str, str]:
    fields = line.rstrip('\n\r').split('\t')
    if len(fields) != 5:
        raise ValueError(f'{location}: FastANI must have exactly 5 tab-separated columns')
    q, r, ani, shared, total = fields
    if q not in qs or r not in rs:
        raise ValueError(f'{location}: query/reference path is outside the corresponding list')
    try:
        a = float(ani)
        m = int(shared)
        t = int(total)
    except ValueError as exc:
        raise ValueError(f'{location}: invalid ANI or fragment counts') from exc
    if not (math.isfinite(a) and 0 <= a <= 100 and 0 <= m <= t and t > 0):
        raise ValueError(f'{location}: ANI or fragment counts outside valid ranges')
    return q, r, ani, shared, total


def validate_output(path: Path, qs: set[str], rs: set[str]) -> int:
    n = 0
    with path.open(encoding='utf-8') as fh:
        for lineno, line in enumerate(fh, 1):
            record(line, f'{path}:{lineno}', qs, rs)
            n += 1
    return n


def run_one(task: tuple[int, int], args: argparse.Namespace,
            qblocks: list[list[str]], rblocks: list[list[str]], binary: str) -> None:
    i, j = task
    label = key(i, j)
    output = args.out / 'blocks' / (label + '.tsv')
    done = args.out / 'blocks' / (label + '.done')
    if done.exists():
        if not output.exists():
            raise RuntimeError(f'{label}: completion marker exists but output is missing')
        validate_output(output, set(qblocks[i]), set(rblocks[j]))
        return
    temporary = args.out / 'blocks' / (label + '.part')
    log = args.out / 'logs' / (label + '.log')
    if temporary.exists():
        temporary.unlink()  # known, non-authoritative incomplete output from this task
    command = [binary, '--ql', str(args.out / 'lists' / f'q{i:04d}.txt'),
               '--rl', str(args.out / 'lists' / f'r{j:04d}.txt'),
               '--fragLen', str(args.frag_len), '--minFraction', str(args.min_fraction),
               '-t', str(args.threads), '-o', str(temporary)]
    with log.open('w', encoding='utf-8') as fh:
        fh.write(' '.join(command) + '\n')
        fh.flush()
        process = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT, check=False)
    if process.returncode != 0:
        raise RuntimeError(f'{label}: FastANI exit code {process.returncode}; see {log}')
    if not temporary.exists():
        temporary.touch()  # successful FastANI run may have zero passing pairs
    validate_output(temporary, set(qblocks[i]), set(rblocks[j]))
    temporary.replace(output)
    done.write_text('ok\n', encoding='utf-8')


def merge(args: argparse.Namespace, config: dict, qblocks: list[list[str]],
          rblocks: list[list[str]]) -> dict:
    symmetric = config['mode'] == 'symmetric'
    jobs = jobs_for(len(qblocks), len(rblocks), symmetric)
    for i, j in jobs:
        label = key(i, j)
        if not (args.out / 'blocks' / (label + '.done')).exists() or not (args.out / 'blocks' / (label + '.tsv')).exists():
            raise RuntimeError(f'Cannot merge: missing completed block {label}')
    merged = args.out / 'fastani_pairs.tsv'
    summary_path = args.out / 'fastani_summary.json'
    if merged.exists() and summary_path.exists():
        previous = json.loads(summary_path.read_text(encoding='utf-8'))
        if previous.get('mode') != config['mode'] or previous.get('query_genomes') != config['query_genomes'] or previous.get('reference_genomes') != config['reference_genomes']:
            raise ValueError('Existing merged output does not match run configuration')
        return previous
    if merged.exists() or summary_path.exists():
        raise ValueError('Partial merged outputs already exist; inspect before retrying')
    temporary = args.out / 'fastani_pairs.tsv.part'
    if temporary.exists():
        raise ValueError(f'Incomplete merge exists; inspect before retrying: {temporary}')
    raw_count = 0
    retained = 0
    skipped_self = 0
    excluded_mirrors = 0
    with temporary.open('w', encoding='utf-8') as out:
        for i, j in jobs:
            path = args.out / 'blocks' / (key(i, j) + '.tsv')
            qset, rset = set(qblocks[i]), set(rblocks[j])
            if symmetric and i == j:
                selected: dict[tuple[str, str], tuple[str, str, str, str, str]] = {}
                with path.open(encoding='utf-8') as raw:
                    for lineno, line in enumerate(raw, 1):
                        q, r, ani, shared, total = record(line, f'{path}:{lineno}', qset, rset)
                        raw_count += 1
                        if q == r:
                            skipped_self += 1
                            continue
                        pkey = (min(q, r), max(q, r))
                        if pkey in selected:
                            excluded_mirrors += 1
                            # Keep the lexicographic query -> reference direction;
                            # NEVER swap labels without swapping direction-dependent values.
                            if q < r and selected[pkey][0] > selected[pkey][1]:
                                selected[pkey] = (q, r, ani, shared, total)
                        else:
                            selected[pkey] = (q, r, ani, shared, total)
                for value in selected.values():
                    out.write('\t'.join(value) + '\n')
                    retained += 1
            else:
                seen: set[tuple[str, str]] = set()
                with path.open(encoding='utf-8') as raw:
                    for lineno, line in enumerate(raw, 1):
                        q, r, ani, shared, total = record(line, f'{path}:{lineno}', qset, rset)
                        raw_count += 1
                        if q == r:
                            skipped_self += 1
                            continue
                        pkey = (q, r) if not symmetric else (min(q, r), max(q, r))
                        if pkey in seen:
                            raise ValueError(f'{path}:{lineno}: duplicate pair in an off-diagonal/cross block')
                        seen.add(pkey)
                        out.write('\t'.join((q, r, ani, shared, total)) + '\n')
                        retained += 1
    n_q = config['query_genomes']
    n_r = config['reference_genomes']
    # Exclude identical-path self-pairs from the theoretical comparison space.
    q_all = {item for block in qblocks for item in block}
    r_all = {item for block in rblocks for item in block}
    expected = n_q * (n_q - 1) // 2 if symmetric else n_q * n_r - len(q_all & r_all)
    if retained > expected:
        temporary.unlink(missing_ok=True)
        raise ValueError(f'Merged result contains {retained} pairs, exceeding {expected} possible comparisons')
    temporary.replace(merged)
    summary = dict(mode=config['mode'], query_genomes=n_q, reference_genomes=n_r,
                   blocks=len(jobs), expected_pairs=expected, numeric_pairs=retained,
                   no_numeric_output=expected-retained, raw_rows=raw_count,
                   excluded_self_rows=skipped_self, excluded_mirrored_rows=excluded_mirrors,
                   alignment_fraction='matched_fragments/total_query_fragments; directional',
                   format='no header: query_path reference_path ANI_percent matched_fragments total_query_fragments')
    (args.out / 'fastani_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genomes', required=True, type=Path,
                        help='Absolute directory or text file containing absolute FASTA paths')
    parser.add_argument('--references', type=Path,
                        help='Optional absolute directory/list; when given run a cross-list query x reference comparison')
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--fastani', default='fastANI')
    parser.add_argument('--batch-size', type=int, default=500)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--frag-len', type=int, default=3000)
    parser.add_argument('--min-fraction', type=float, default=0.2)
    parser.add_argument('--recursive', action='store_true')
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--execute', action='store_true', help='Run FastANI; default only prints plan')
    action.add_argument('--merge-only', action='store_true', help='Merge existing successful block outputs')
    args = parser.parse_args()
    if not args.out.is_absolute():
        parser.error('--out must be an absolute path')
    if args.batch_size < 1 or args.threads < 1 or args.jobs < 1 or args.frag_len < 1 or not 0 < args.min_fraction <= 1:
        parser.error('Invalid batch-size, threads, jobs, fragment length, or minimum fraction')
    try:
        queries = genomes_from(args.genomes, args.recursive)
        refs = genomes_from(args.references, args.recursive) if args.references else queries
        symmetric = args.references is None
        qblocks = blocks(queries, args.batch_size)
        rblocks = qblocks if symmetric else blocks(refs, args.batch_size)
        config = dict(mode='symmetric' if symmetric else 'cross',
                      query_genomes=len(queries), reference_genomes=len(refs),
                      query_list_sha256=digest(queries), reference_list_sha256=digest(refs),
                      batch_size=args.batch_size, frag_len=args.frag_len,
                      min_fraction=args.min_fraction, threads=args.threads,
                      executable=args.fastani)
        task_count = len(jobs_for(len(qblocks), len(rblocks), symmetric))
        pairs = len(queries)*(len(queries)-1)//2 if symmetric else len(queries)*len(refs)-len(set(queries)&set(refs))
        print(f"Mode={config['mode']}; genomes={len(queries)} x {len(refs)}; jobs={task_count}; potential pairs={pairs:,}; cores at full concurrency={args.jobs*args.threads}")
        if not args.execute and not args.merge_only:
            print('DRY RUN. Use --execute to run; use --merge-only to combine completed blocks.')
            return
        if args.merge_only:
            saved = json.loads((args.out / 'run_config.json').read_text(encoding='utf-8'))
            if saved != config:
                raise ValueError('Run configuration or genome list changed')
        else:
            binary = shutil.which(args.fastani)
            if not binary:
                raise ValueError(f'FastANI executable not found: {args.fastani}')
            if args.out.exists() and not (args.out / 'run_config.json').exists() and any(args.out.iterdir()):
                raise ValueError(f'Output directory is non-empty without a run_config.json: {args.out}')
            for folder in ('lists', 'blocks', 'logs'):
                (args.out / folder).mkdir(parents=True, exist_ok=True)
            saved_config = args.out / 'run_config.json'
            if saved_config.exists():
                if json.loads(saved_config.read_text(encoding='utf-8')) != config:
                    raise ValueError('Run configuration or genome list changed')
            else:
                saved_config.write_text(json.dumps(config, indent=2) + '\n', encoding='utf-8')
            for i, chunk in enumerate(qblocks):
                write_exact(args.out / 'lists' / f'q{i:04d}.txt', '\n'.join(chunk) + '\n')
            for i, chunk in enumerate(rblocks):
                write_exact(args.out / 'lists' / f'r{i:04d}.txt', '\n'.join(chunk) + '\n')
            tasks = jobs_for(len(qblocks), len(rblocks), symmetric)
            if args.jobs == 1:
                for task in tasks:
                    run_one(task, args, qblocks, rblocks, binary)
            else:
                with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
                    results = [pool.submit(run_one, task, args, qblocks, rblocks, binary) for task in tasks]
                    for future in concurrent.futures.as_completed(results):
                        future.result()
        summary = merge(args, config, qblocks, rblocks)
        print(json.dumps(summary, indent=2))
    except (ValueError, RuntimeError, OSError, json.JSONDecodeError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
