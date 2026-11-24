#!/usr/bin/env python3
"""Extract RNAmmer 16S calls, report coordinates and choose a full-length representative."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


def parse_gff(path: Path) -> list[dict]:
    calls = []
    for ln, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) != 9:
            raise ValueError(f"{path}:{ln}: expected 9 GFF columns")
        seqid, source, kind, first, last, raw_score, strand, phase, attribute = cols
        if "16s" not in attribute.casefold() and "ssu" not in attribute.casefold():
            continue
        start, end = int(first), int(last)
        if start < 1 or end < start or strand not in "+-":
            raise ValueError(f"{path}:{ln}: invalid 16S coordinates/strand")
        if raw_score == ".":
            raise ValueError(f"{path}:{ln}: missing RNAmmer score; representative cannot be selected")
        calls.append(dict(seqid=seqid, start=start, end=end, score=float(raw_score), strand=strand))
    return calls


def taxonomy_phyla(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"accession", "phylum"}.issubset(reader.fieldnames):
            raise ValueError("Taxonomy TSV requires accession and phylum columns")
        phyla = {}
        for row in reader:
            if row["accession"] in phyla:
                raise ValueError(f"Duplicate accession in taxonomy: {row['accession']}")
            phyla[row["accession"]] = row["phylum"] or "Unclassified"
    return phyla


def process(genomes: Path, gff_dir: Path, output: Path, minimum: int, maximum: int,
            taxonomy: Path | None = None) -> list[dict]:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    files = sorted(genomes.glob("*.fna"))
    if not files:
        raise ValueError(f"No .fna files in {genomes}")
    for f in files:
        if not (gff_dir / f"{f.stem}.gff").is_file():
            raise FileNotFoundError(f"Missing GFF for {f.stem}")
    phyla = taxonomy_phyla(taxonomy)
    if taxonomy and any(f.stem not in phyla for f in files):
        raise ValueError("An input genome is missing from the taxonomy TSV")
    output.mkdir(parents=True)
    rows, features = [], []
    with (output / "full_length_16s.fna").open("w") as all_fa, (output / "representative_16s.fna").open("w") as best_fa:
        for fasta in files:
            assembly = fasta.stem
            with fasta.open() as handle:
                contigs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            calls = parse_gff(gff_dir / f"{assembly}.gff")
            eligible = []
            per_copy = []
            for i, call in enumerate(calls, 1):
                if call["seqid"] not in contigs:
                    raise ValueError(f"GFF contig {call['seqid']} absent from {fasta}")
                genome = contigs[call["seqid"]].seq
                if call["end"] > len(genome):
                    raise ValueError(f"GFF call beyond contig length: {assembly}, {i}")
                seq = genome[call["start"] - 1:call["end"]]
                if call["strand"] == "-":
                    seq = seq.reverse_complement()
                full = minimum <= len(seq) <= maximum
                feature = dict(accession=assembly, call=i, contig=call["seqid"], start=call["start"],
                               end=call["end"], strand=call["strand"], score=call["score"],
                               length=len(seq), full_length=int(full), selected=0)
                per_copy.append(feature)
                if full:
                    eligible.append((call["score"], call["seqid"], call["start"], i, str(seq)))
                    all_fa.write(f">{assembly}|16S_{i} score={call['score']}\n{seq}\n")
            picked = None
            if eligible:
                picked = sorted(eligible, key=lambda x: (-x[0], x[1], x[2], x[3]))[0]
                best_fa.write(f">{assembly} source_call={picked[3]} score={picked[0]}\n{picked[4]}\n")
                per_copy[picked[3] - 1]["selected"] = 1
            features.extend(per_copy)
            rows.append(dict(accession=assembly, phylum=phyla.get(assembly, ""),
                             all_16s_calls=len(calls), full_length_16s=len(eligible),
                             has_representative=int(bool(eligible)),
                             best_score=(picked[0] if picked else "")))

    tables = [
        ("16s_calls.tsv", features, ["accession", "call", "contig", "start", "end", "strand", "score", "length", "full_length", "selected"]),
        ("16s_recovery_by_genome.tsv", rows, ["accession", "phylum", "all_16s_calls", "full_length_16s", "has_representative", "best_score"]),
    ]
    if taxonomy:
        by_phylum = defaultdict(lambda: dict(genomes=0, with_16s=0, with_full_length_16s=0))
        for row in rows:
            z = by_phylum[row["phylum"]]
            z["genomes"] += 1
            z["with_16s"] += int(row["all_16s_calls"] > 0)
            z["with_full_length_16s"] += row["has_representative"]
        phylum_rows = [dict(phylum=p, **counts) for p, counts in sorted(by_phylum.items())]
        tables.append(("16s_recovery_by_phylum.tsv", phylum_rows,
                       ["phylum", "genomes", "with_16s", "with_full_length_16s"]))
    for filename, data, fieldnames in tables:
        with (output / filename).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(data)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genomes", type=Path, required=True)
    p.add_argument("--gff-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--taxonomy", type=Path, help="Optional frozen TSV for phylum summary")
    p.add_argument("--min-length", type=int, default=1450)
    p.add_argument("--max-length", type=int, default=1700)
    args = p.parse_args()
    if args.min_length < 1 or args.max_length < args.min_length:
        p.error("Invalid length interval")
    rows = process(args.genomes, args.gff_dir, args.out, args.min_length, args.max_length, args.taxonomy)
    print(f"Genomes: {len(rows)}; recovered full-length 16S: {sum(r['has_representative'] for r in rows)}")


if __name__ == "__main__":
    main()
