"""Unit checks of upstream processing functions; not a genome benchmark."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1] / "scripts" / "16s_rRNA"


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


proc = load("process_rnammer")
info = load("conservation")


class RNAmmerProcessingTest(unittest.TestCase):
    def test_gff_score_and_copy_selection(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            (p / "genomes").mkdir(); (p / "gff").mkdir()
            (p / "genomes" / "GCA_000000000.1.fna").write_text(">chr\n" + "ACGT" * 400 + "\n")
            (p / "gff" / "GCA_000000000.1.gff").write_text(
                "chr\tRNAmmer-1.2\trRNA\t1\t1500\t60\t+\t.\t16s_rRNA\n"
                "chr\tRNAmmer-1.2\trRNA\t50\t1549\t90\t-\t.\t16s_rRNA\n")
            res = proc.process(p / "genomes", p / "gff", p / "out", 1450, 1700)
            self.assertEqual(res[0]["full_length_16s"], 2)
            self.assertEqual(res[0]["best_score"], 90)
            self.assertIn("source_call=2", (p / "out" / "representative_16s.fna").read_text())

    def test_missing_gff_fails(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder); (p / "g").mkdir(); (p / "r").mkdir()
            (p / "g" / "a.fna").write_text(">contig\nACTG\n")
            with self.assertRaises(FileNotFoundError):
                proc.process(p / "g", p / "r", p / "out", 2, 6)


class InformationTest(unittest.TestCase):
    def test_conserved_and_variable(self):
        aln = MultipleSeqAlignment([SeqRecord(Seq("AA")), SeqRecord(Seq("AC")), SeqRecord(Seq("AG")), SeqRecord(Seq("AT"))])
        positions, gaps, raw, smoothed = info.compute(aln, window=1)
        self.assertEqual(positions, [1,2])
        self.assertAlmostEqual(raw[0], 2.0)
        self.assertAlmostEqual(raw[1], 0.0)

    def test_normalize_u_is_optional(self):
        aln = MultipleSeqAlignment([SeqRecord(Seq("T")), SeqRecord(Seq("U"))])
        self.assertAlmostEqual(info.compute(aln,window=1)[2][0],1.0)
        self.assertAlmostEqual(info.compute(aln,window=1,normalize_u=True)[2][0],2.0)


if __name__ == "__main__":
    unittest.main()
