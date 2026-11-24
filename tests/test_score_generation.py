"""Test file parsing, directional counts, incomplete jobs, and accession filtering."""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def module(filename: str):
    spec = importlib.util.spec_from_file_location(filename.replace('.py', ''), ROOT / {'fastani_blockwise.py': 'scripts/fastani/fastani_blockwise.py', 'mash_triangle_to_pairs.py': 'scripts/mash/mash_triangle_to_pairs.py', 'audit_dictionary.py': 'scripts/taxonomy/audit_dictionary.py'}[filename])
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

fa = module('fastani_blockwise.py')
ma = module('mash_triangle_to_pairs.py')
di = module('audit_dictionary.py')


class FastaniTests(unittest.TestCase):
    def test_dry_run_no_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            g = base / 'GCA_000000001.1.fna'
            g.write_text('>g\nACGT\n')
            result = subprocess.run([sys.executable, str(ROOT / 'scripts/fastani/fastani_blockwise.py'),
                                     '--genomes', str(g.parent), '--out', str(base / 'out')],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('DRY RUN', result.stdout)
            self.assertFalse((base / 'out').exists())

    def test_directional_counts_not_swapped(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            for folder in ('blocks', 'logs', 'lists'):
                (base / folder).mkdir()
            a, b, c = '/data/A.fna', '/data/B.fna', '/data/C.fna'
            (base / 'blocks/q0000_r0000.tsv').write_text(
                f'{b}\t{a}\t90\t8\t12\n{a}\t{b}\t92\t7\t10\n{a}\t{a}\t100\t10\t10\n')
            (base / 'blocks/q0000_r0001.tsv').write_text(
                f'{a}\t{c}\t79\t6\t10\n{b}\t{c}\t85\t7\t10\n')
            (base / 'blocks/q0001_r0001.tsv').write_text('')
            for p in (base / 'blocks').glob('*.tsv'):
                p.with_suffix('.done').write_text('ok\n')
            class Args: out=base
            summary = fa.merge(Args(), {'mode':'symmetric','query_genomes':3,'reference_genomes':3},
                               [[a,b],[c]], [[a,b],[c]])
            rows = (base / 'fastani_pairs.tsv').read_text().splitlines()
            self.assertEqual(len(rows), 3)
            self.assertIn(f'{a}\t{b}\t92\t7\t10', rows)
            self.assertEqual(summary['excluded_mirrored_rows'], 1)
            self.assertEqual(summary['excluded_self_rows'], 1)
            self.assertEqual(summary['no_numeric_output'], 0)

    def test_reject_record_outside_lists(self):
        with self.assertRaises(ValueError):
            fa.record('/other\t/ref\t95\t1\t2', 'test', {'/q'}, {'/ref'})

    def test_reject_bad_fragments(self):
        with self.assertRaises(ValueError):
            fa.record('/q\t/ref\t95\t8\t2', 'test', {'/q'}, {'/ref'})

    def test_deterministic_triangle_jobs(self):
        self.assertEqual(fa.jobs_for(3, 3, True), [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)])
        self.assertEqual(len(fa.jobs_for(3, 2, False)), 6)

    def test_failed_fastani_marks_no_success(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            genomes = base / 'genomes'; genomes.mkdir()
            for i in range(2):
                (genomes / f'GCA_00000000{i+1}.1.fna').write_text('>g\nACGT\n')
            failing = base / 'fail_fastani.sh'
            failing.write_text('#!/usr/bin/env bash\nexit 17\n')
            failing.chmod(0o700)
            dest = base / 'run'
            result = subprocess.run([sys.executable, str(ROOT / 'scripts/fastani/fastani_blockwise.py'),
                                     '--genomes', str(genomes), '--out', str(dest),
                                     '--fastani', str(failing), '--batch-size', '2', '--execute'],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(list((dest / 'blocks').glob('*.done')))
            self.assertTrue(list((dest / 'logs').glob('*.log')))


class MashTests(unittest.TestCase):
    def test_lower_triangle(self):
        with tempfile.TemporaryDirectory() as temp:
            src=Path(temp)/'in.txt'; out=Path(temp)/'out.tsv'
            src.write_text('3\nA\nB 0.1\nC 0.2 0.3\n')
            s=ma.convert(src,out)
            self.assertEqual(s, {'genomes':3,'pairs':3})
            self.assertEqual(out.read_text().splitlines(), ['B\tA\t0.1','C\tA\t0.2','C\tB\t0.3'])

    def test_invalid_triangle_no_output(self):
        with tempfile.TemporaryDirectory() as temp:
            src=Path(temp)/'in.txt'; out=Path(temp)/'out.tsv'
            src.write_text('3\nA\nB 0.1\nC 0.2\n')
            with self.assertRaises(ValueError): ma.convert(src,out)
            self.assertFalse(out.exists())
            self.assertFalse((Path(temp)/'out.tsv.partial').exists())

    def test_reject_mash_dist_five_columns(self):
        with tempfile.TemporaryDirectory() as temp:
            src=Path(temp)/'in.txt';out=Path(temp)/'out.tsv'
            src.write_text('A\tB\t0.1\t0\t99/100\n')
            with self.assertRaises(ValueError): ma.convert(src,out)


class DictionaryTests(unittest.TestCase):
    def test_extra_accessions_not_silently_included(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp); src=base/'dictionary.py'; man=base/'ids.txt'; exported=base/'result.tsv'
            id1='GCA_000000001.1'; id2='GCA_000000002.1'
            row=dict.fromkeys(di.RANKS,'X')
            src.write_text('GenID='+repr({id1:row,id2:row}))
            man.write_text(id1+'\n')
            result=di.audit(src, man, export=exported)
            self.assertEqual(result['extra_dictionary_accessions'],[id2])
            self.assertEqual(len(exported.read_text().splitlines()),2)
            self.assertNotIn(id2,exported.read_text())

    def test_reference_mismatch_reported(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp); src=base/'dictionary.py'; man=base/'ids.txt'; ref=base/'reference.tsv'
            id1='GCA_000000001.1'
            row=dict.fromkeys(di.RANKS,'A')
            src.write_text('GenID='+repr({id1:row}))
            man.write_text(id1+'\n')
            ref.write_text('accession\tdomain\n'+id1+'\tB\n')
            result=di.audit(src,man,reference=ref)
            self.assertEqual(result['reference_field_mismatches'],1)

    def test_code_execution_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp); src=base/'dictionary.py'
            src.write_text("GenID=__import__('os').system('exit 4')")
            with self.assertRaises(ValueError): di.load_dictionary(src)


if __name__ == '__main__':
    unittest.main()
