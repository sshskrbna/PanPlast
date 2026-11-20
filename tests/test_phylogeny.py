"""Tree input, conservative taxon matching, missing gene states and report output."""
import io
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from panplast.cli import main, parse_args
from panplast.genes import Locus
from panplast.phylogeny import build_overlay, load_tip_map, match_taxa, plot_overlay, read_tree


class TestTreeInput(unittest.TestCase):
    def test_newick_and_iqtree_ml_section_preserve_support(self):
        newick = "((A:0.1,B:0.2)88.5/99:0.3,C:0.4);"
        with tempfile.TemporaryDirectory() as td:
            for ext in ("treefile", "contree", "nwk", "iqtree"):
                path = Path(td) / f"run.{ext}"
                text = ("MAXIMUM LIKELIHOOD TREE\nASCII tree here\nTree in newick format:\n\n" + newick +
                        "\n\nCONSENSUS TREE\nTree in newick format:\n\n(A,B,C);" if ext == "iqtree" else newick)
                path.write_text(text)
                with self.subTest(ext=ext):
                    tree = read_tree(path)
                    self.assertEqual([tip.name for tip in tree.get_terminals()], ["A", "B", "C"])
                    self.assertEqual(tree.root.clades[0].name, "88.5/99")
                    self.assertAlmostEqual(tree.root.clades[0].branch_length, .3)

    def test_invalid_trees_fail_with_clear_error(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bad.treefile"
            for text in ("(A,B)", "(A,A);", "A;", "(A:-1,B:1);", "(A:0.1,B:0.1));", "(A,B);(C,D);"):
                with self.subTest(text=text):
                    path.write_text(text)
                    with self.assertRaises(ValueError):
                        read_tree(path)
            path = Path(td) / "bad.iqtree"
            path.write_text("IQ-TREE report without a tree")
            with self.assertRaisesRegex(ValueError, "treefile"):
                read_tree(path)

    def test_tip_map_requires_valid_unique_entries(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "map.tsv"
            path.write_text("tip\tgenome\nA\tg1\nB\tg2\n")
            self.assertEqual(load_tip_map(path), {"A": "g1", "B": "g2"})
            for text in ("name\tgenome\nA\tg1\n", "tip\tgenome\nA\tg1\nA\tg2\n", "tip\tgenome\nA\t\n"):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    load_tip_map(path)

    def test_tree_options_require_a_tree_and_nonempty_genes(self):
        for argv in (["--tree-genes", "rbcL"], ["--tree-map", "map.tsv"],
                     ["--tree", "x.treefile", "--tree-genes", ","]):
            with self.subTest(argv=argv), patch("sys.stderr", new=io.StringIO()), self.assertRaises(SystemExit):
                parse_args(argv)


class TestTaxonMatching(unittest.TestCase):
    def test_ids_accessions_species_unknown_and_genomes_outside_tree(self):
        matrix = pd.DataFrame([[1, 1, 1, 1]], index=["rbcL"],
                              columns=["NC_000001.1", "NC_000002.1", "g3", "g4"])
        labels = {"NC_000001.1": "Alpha species", "NC_000002.1": "Beta species", "g3": "Gamma species"}
        tree = Phylo.read(io.StringIO("(NC_000001,NC_000002.1_Beta_species,Gamma_species,outgroup);"), "newick")
        taxa = match_taxa(tree, matrix, labels).set_index("tip")
        self.assertEqual(taxa.loc["NC_000001", "genome"], "NC_000001.1")
        self.assertEqual(taxa.loc["NC_000002.1_Beta_species", "genome"], "NC_000002.1")
        self.assertEqual(taxa.loc["Gamma_species", "genome"], "g3")
        self.assertEqual(taxa.loc["outgroup", "status"], "unmatched_tip")
        self.assertEqual(taxa.loc["", "status"], "genome_not_in_tree")
        self.assertEqual(taxa.loc["", "genome"], "g4")

    def test_ambiguous_species_needs_explicit_mapping(self):
        matrix = pd.DataFrame([[1, 2]], index=["rbcL"], columns=["g1", "g2"])
        tree = Phylo.read(io.StringIO("(Same_species,g2);"), "newick")
        labels = {"g1": "Same species", "g2": "Same species"}
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            match_taxa(tree, matrix, labels)
        taxa = match_taxa(tree, matrix, labels, tip_map={"Same_species": "g1"})
        self.assertEqual(list(taxa.genome), ["g1", "g2"])

    def test_duplicate_mappings_and_unknown_ids_are_rejected(self):
        matrix = pd.DataFrame([[1]], index=["rbcL"], columns=["g1"])
        tree = Phylo.read(io.StringIO("(a,b);"), "newick")
        for mapping in ({"a": "g1", "b": "g1"}, {"a": "unknown"}, {}):
            with self.subTest(mapping=mapping), self.assertRaises(ValueError):
                match_taxa(tree, matrix, {}, tip_map=mapping)


class TestGeneOverlay(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "run.treefile"
        self.path.write_text("((g2:0.1,outgroup:0.2)99:0.1,g1:0.3);")
        self.matrix = pd.DataFrame({"g1": [1, 2, 1], "g2": [0, 1, 1]}, index=["infA", "rps19", "rbcL"])
        self.summary = pd.DataFrame([dict(genome=g, IRb_bp=20000, length_bp=150000) for g in ("g1", "g2")])

    def test_unknown_is_not_absent_and_missing_counts_use_matched_genomes(self):
        overlay = build_overlay(self.path, self.matrix, {})
        self.assertEqual(list(overlay.copies.index), ["g2", "outgroup", "g1"])
        self.assertEqual(list(overlay.copies.columns), ["infA"])
        self.assertEqual(overlay.copies.at["g2", "infA"], 0)
        self.assertTrue(pd.isna(overlay.copies.at["outgroup", "infA"]))
        self.assertEqual(list(overlay.absences.genome), ["g2"])
        self.assertEqual(list(overlay.absences.status), ["not_annotated"])

    def test_explicit_genes_can_include_globally_absent_gene_and_flagged_loci(self):
        partial = Locus("infA", "infA", "CDS", 0, 50, 1, [(0, 50)], flags={"partial"})
        overlay = build_overlay(self.path, self.matrix, {}, genes=["infA", "ycf2"], loci_by={"g2": [partial]})
        row = overlay.absences.set_index(["genome", "gene"]).loc[("g2", "infA")]
        self.assertEqual(row.status, "not_counted")
        self.assertEqual(row["flags"], "partial")
        self.assertEqual(overlay.copies.at["g1", "ycf2"], 0)
        with self.assertRaisesRegex(ValueError, "Unknown"):
            build_overlay(self.path, self.matrix, {}, genes=["typo"])

    def test_figure_aligns_unknown_and_absent_states_to_tree_order(self):
        from panplast.plots import plt
        overlay = build_overlay(self.path, self.matrix, {}, genes=["infA", "rps19"])
        with patch("panplast.plots._save") as save:
            plot_overlay(overlay, self.summary, Path(self.tmp.name) / "figure")
            fig = save.call_args.args[0]
            try:
                self.assertEqual(fig.axes[-1].images[0].get_array().tolist(), [[0, 1], [-1, -1], [1, 2]])
                self.assertIn("99", [text.get_text() for text in fig.axes[0].texts])
                self.assertEqual(fig.axes[0].get_xlabel(), "Substitutions / site")
            finally:
                plt.close(fig)

    def test_pagination_and_topology_only_trees(self):
        from panplast.plots import plt
        self.path.write_text("(g1,g2,outgroup);")
        matrix = pd.DataFrame({"g1": [1] * 31, "g2": [0] * 31}, index=[f"gene{i}" for i in range(31)])
        overlay = build_overlay(self.path, matrix, {})
        figures = []
        def capture(fig, path):
            figures.append(fig)
            return {"png": f"{path}.png"}
        try:
            with patch("panplast.plots._save", side_effect=capture):
                result = plot_overlay(overlay, self.summary, Path(self.tmp.name) / "figure")
            self.assertEqual(len(result), 2)
            self.assertEqual(figures[0].axes[-1].images[0].get_array().shape, (3, 30))
            self.assertIn("Topology", figures[0].axes[0].get_xlabel())
        finally:
            for fig in figures:
                plt.close(fig)

    def test_empty_matrix_and_no_absences_still_render(self):
        for matrix in (self.matrix.loc[["rbcL"]], pd.DataFrame(columns=["g1", "g2"])):
            with self.subTest(genes=list(matrix.index)):
                overlay = build_overlay(self.path, matrix, {})
                self.assertTrue(overlay.absences.empty)
                output = plot_overlay(overlay, self.summary, Path(self.tmp.name) / "empty")
                self.assertTrue(Path(output[0]["png"]).exists())


@unittest.skipUnless(shutil.which("blastn"), "blastn not on PATH")
class TestPhylogenyEndToEnd(unittest.TestCase):
    def test_multiple_trees_export_figures_tables_excel_and_report(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            for i, species in enumerate(("Alpha species", "Beta species", "Gamma species"), 1):
                rec = SeqRecord(Seq("ACGT" * 3000), id=f"NC_00000{i}.1",
                                annotations={"molecule_type": "DNA", "organism": species})
                names = [(100, 1162, "psbA")]
                if i != 2:
                    names.append((2000, 3434, "rbcL"))
                if i == 1:
                    names.append((4000, 4050, "infA"))
                for start, end, name in names:
                    for kind in ("gene", "CDS"):
                        rec.features.append(SeqFeature(FeatureLocation(start, end, strand=1), type=kind,
                                                       qualifiers={"gene": [name]}))
                SeqIO.write(rec, root / f"{rec.id}.gb", "genbank")
            newick = "((sample_a:0.1,Beta_species:0.2)88.5/99:0.1,outgroup:0.5);"
            tree = root / "run.contree"
            tree.write_text(newick)
            report_tree = root / "run.iqtree"
            report_tree.write_text("IQ-TREE report\nTree in newick format:\n\n" + newick + "\n\nOTHER RESULTS\n")
            mapping = root / "map.tsv"
            mapping.write_text("tip\tgenome\nsample_a\tNC_000001.1\n")
            out = root / "res"
            main(["--genbank", str(root), "--ir-source", "genbank", "--maps", "0", "--out", str(out),
                  "--tree", str(tree), str(report_tree), "--tree-map", str(mapping),
                  "--tree-genes", "rbcL,infA"])
            for i in (1, 2):
                for ext in ("png", "pdf", "svg"):
                    self.assertTrue((out / "figures" / f"phylogeny_{i}.{ext}").exists())
            taxa = pd.read_csv(out / "tables" / "tree_taxa.tsv", sep="\t")
            self.assertEqual(len(taxa[taxa.status == "matched"]), 4)
            self.assertEqual(len(taxa[taxa.status == "unmatched_tip"]), 2)
            self.assertEqual(len(taxa[taxa.status == "genome_not_in_tree"]), 2)
            copies = pd.read_csv(out / "tables" / "tree_copies.tsv", sep="\t")
            self.assertEqual(len(copies[copies.tip == "outgroup"]), 4)
            self.assertTrue(copies.loc[copies.tip == "outgroup", "copies"].isna().all())
            absences = pd.read_csv(out / "tables" / "tree_absences.tsv", sep="\t")
            self.assertEqual(set(absences.status), {"not_counted", "not_annotated"})
            with pd.ExcelFile(out / "panplast_results.xlsx") as workbook:
                self.assertTrue({"Tree taxa", "Tree absences", "Tree gene copies"} <= set(workbook.sheet_names))
            report = (out / "report.html").read_text()
            self.assertIn('id="phylogeny"', report)
            self.assertIn("not reconstructed evolutionary loss events", report)
            summary = pd.read_csv(out / "tables" / "summary.tsv", sep="\t")
            self.assertEqual(len(summary), 3)  # The tree does not change the pan-genome cohort.


if __name__ == "__main__":
    unittest.main()
