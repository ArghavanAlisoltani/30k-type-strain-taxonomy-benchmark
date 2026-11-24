"""Unit tests for candidate scripts. These fixtures are not manuscript benchmark data."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts" / "taxonomy"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


tax = load("taxonomy_fill_safe")
rank = load("classify_pairwise_ranks")


class TaxonomyFillTests(unittest.TestCase):
    def test_species_name_stays_unchanged_and_label_not_count(self):
        d = {
            "GCA_a.1": {"genus": "Alpha", "species": "Alpha a", "family": "F1"},
            "GCA_b.1": {"genus": "Alpha", "species": "Alpha b", "family": "F1"},
            "GCA_c.1": {"genus": "Alpha", "species": "Alpha c", "family": None},
        }
        proposals = tax.proposed_fills(d)
        self.assertEqual([x["proposed_value"] for x in proposals], ["F1"])
        self.assertEqual(proposals[0]["rank"], "family")
        self.assertIsNone(d["GCA_c.1"]["family"])

    def test_conflicting_lookup_is_not_filled(self):
        d = {"a": {"genus": "A", "species": "A a", "family": "X"},
             "b": {"genus": "A", "species": "A b", "family": "Y"},
             "c": {"genus": "A", "species": "A c", "family": ""}}
        self.assertEqual(tax.proposed_fills(d), [])

    def test_parse_safe_literal_and_reject_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "g.txt"
            p.write_text("GenID={'x': {'genus': 'A', 'species': 'A a'}}", encoding="utf-8")
            self.assertEqual(len(tax.load_genid(p)), 1)
            p.write_text("GenID=__import__('os').system('echo BAD')", encoding="utf-8")
            with self.assertRaises((ValueError, SyntaxError)):
                tax.load_genid(p)

    def test_missing_placeholder(self):
        self.assertIsNone(tax.clean("N/A"))
        self.assertIsNone(tax.clean("nan"))


class RankAssignmentTests(unittest.TestCase):
    A = {"species": "A a", "genus": "A", "family": "F", "order": "O", "class": "C", "phylum": "P", "kingdom": "K", "domain": "Bacteria"}

    def test_species_identical(self):
        self.assertEqual(rank.classify(self.A, dict(self.A)), ("species", "classified"))

    def test_same_genus_different_species(self):
        b = {**self.A, "species": "A b"}
        self.assertEqual(rank.classify(self.A, b), ("genus", "classified"))

    def test_same_family_different_genus(self):
        b = {**self.A, "genus": "B", "species": "B a"}
        self.assertEqual(rank.classify(self.A, b), ("family", "classified"))

    def test_missing_species_is_not_genus_negative(self):
        b = {**self.A, "species": ""}
        self.assertEqual(rank.classify(self.A, b), ("", "unresolved_lower_rank"))

    def test_inconsistent_lineage_is_flagged(self):
        b = {**self.A, "family": "F2"}
        self.assertEqual(rank.classify(self.A, b), ("", "inconsistent_lineage"))

    def test_pairwise_tsv_one_row_per_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tf, sf, of = (root / f for f in ("tax.tsv", "scores.tsv", "out.tsv"))
            header = "accession\t" + "\t".join(rank.RANKS) + "\n"
            def line(name, d):
                return name + "\t" + "\t".join(d[x] for x in rank.RANKS) + "\n"
            tf.write_text(header + line("GCA_a.1", self.A) + line("GCA_b.1", {**self.A, "species":"A b"}))
            sf.write_text("genome1\tgenome2\tscore\nGCA_a.1\tGCA_b.1\t98.8\n")
            self.assertEqual(rank.run(tf,sf,of), {"classified":1})
            self.assertIn("\tgenus\tclassified", of.read_text())
            with self.assertRaises(FileExistsError):
                rank.run(tf,sf,of)

if __name__ == "__main__":
    unittest.main()
