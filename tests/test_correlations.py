"""Numerical checks and end-to-end coverage for genome correlation plots."""
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from panplast.cli import main
from panplast.correlations import correlation_results, fit_regression, scatter_data, scatter_plots
from panplast.genes import Locus


class TestRegression(unittest.TestCase):
    def test_exact_line_and_paired_missing_values(self):
        fit = fit_regression([1, 2, 3, 4, np.nan, 6], [5, 8, 11, 14, 17, np.inf])
        self.assertEqual(fit.summary["n"], 4)
        self.assertEqual(fit.summary["n_excluded"], 2)
        self.assertEqual(fit.summary["status"], "ok")
        self.assertAlmostEqual(fit.summary["slope"], 3)
        self.assertAlmostEqual(fit.summary["intercept"], 2)
        self.assertAlmostEqual(fit.summary["pearson_r"], 1)
        np.testing.assert_allclose(fit.lower, fit.mean)
        np.testing.assert_allclose(fit.upper, fit.mean)

    def test_regression_and_mean_confidence_interval_against_reference(self):
        # For x=[0,1,2,3], y=[0,1,1,3]: slope=.9, intercept=-.1,
        # SSE=.7, Sxx=5; Student t(.975,2)=4.302652729696142.
        fit = fit_regression([0, 1, 2, 3], [0, 1, 1, 3])
        self.assertAlmostEqual(fit.summary["slope"], .9)
        self.assertAlmostEqual(fit.summary["intercept"], -.1)
        self.assertAlmostEqual(fit.summary["pearson_r"], np.sqrt(4.05 / 4.75))
        self.assertAlmostEqual(fit.summary["p_value"], 1 - np.sqrt(4.05 / 4.75))
        expected_half_width = 4.302652729696142 * np.sqrt(.35 * (.25 + 1.5 ** 2 / 5))
        self.assertAlmostEqual(fit.upper[0] - fit.mean[0], expected_half_width)
        self.assertAlmostEqual(fit.mean[0] - fit.lower[0], expected_half_width)

    def test_small_or_constant_samples_have_no_inferential_statistics(self):
        for x, y, status in (([], [], "fewer_than_3_pairs"), ([1, 2], [2, 3], "fewer_than_3_pairs"),
                             ([1, 1, 1], [1, 2, 3], "constant_x"),
                             ([1, 2, 3], [4, 4, 4], "constant_y")):
            with self.subTest(status=status):
                fit = fit_regression(x, y)
                self.assertEqual(fit.summary["status"], status)
                self.assertTrue(np.isnan(fit.summary["pearson_r"]))
                self.assertTrue(np.isnan(fit.summary["p_value"]))
                self.assertEqual(len(fit.x), 0)


class TestScatterData(unittest.TestCase):
    def setUp(self):
        # Genome order differs between summary and matrix: align by ID, not position.
        self.summary = pd.DataFrame(dict(genome=["g2", "g1", "g3"], species=["Beta b", "Alpha a", "Gamma c"],
                                         group=["B", "A", "A"], length_bp=[120000, 100000, 130000],
                                         IRb_bp=[20000, np.nan, 30000]))
        self.matrix = pd.DataFrame({"g1": [2, 1], "g2": [0, 1], "g3": [1, 1]}, index=["ycf2", "rbcL"])

    def test_copy_count_unique_count_and_highlights_are_aligned(self):
        data = scatter_data(self.summary, self.matrix, highlight_genomes=["g2"])
        self.assertEqual(list(data.gene_count), [1, 3, 2])
        self.assertEqual(list(data.highlight_copies), [0, 2, 1])
        self.assertEqual(list(data.marked), [True, False, False])
        unique = scatter_data(self.summary, self.matrix, gene_count="unique")
        self.assertEqual(list(unique.gene_count), [1, 2, 2])
        self.assertEqual(data.loc[0, "genome_kb"], 120)
        self.assertTrue(pd.isna(data.loc[1, "ir_kb"]))

    def test_missing_ir_only_affects_ir_comparisons(self):
        table = correlation_results(scatter_data(self.summary, self.matrix)).set_index("comparison")
        self.assertEqual(table.loc["genes_vs_genome", "n"], 3)
        self.assertEqual(table.loc["genome_vs_ir", "n"], 2)
        self.assertEqual(table.loc["genes_vs_ir", "n_excluded"], 1)

    def test_partial_only_gene_is_not_a_counted_copy(self):
        loc = Locus("infA", "infA", "CDS", 0, 50, 1, [(0, 50)], flags={"partial"})
        data = scatter_data(self.summary, self.matrix, highlight_gene="infA", loci_by={"g1": [loc]})
        self.assertTrue((data.highlight_copies == 0).all())
        disabled = scatter_data(self.summary, self.matrix, highlight_gene=None)
        self.assertTrue(disabled.highlight_copies.isna().all())

    def test_invalid_highlight_options_fail(self):
        for kwargs in (dict(highlight_genomes=["typo"]), dict(highlight_gene="typo"), dict(gene_count="typo")):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                scatter_data(self.summary, self.matrix, **kwargs)

    def test_figures_use_calculated_statistics_and_highlight_layers(self):
        from panplast.plots import plt
        data = scatter_data(self.summary.fillna({"IRb_bp": 10000}), self.matrix, highlight_genomes=["g2"])
        figures = []
        def capture(fig, path):
            figures.append(fig)
            return {"png": str(path) + ".png"}
        try:
            with patch("panplast.plots._save", side_effect=capture):
                result = scatter_plots(data, Path("unused"))
            self.assertEqual(set(result), {"genes_vs_genome", "genome_vs_ir", "genes_vs_ir"})
            ax = figures[0].axes[0]
            self.assertIn("r = -0.65", " ".join(t.get_text() for t in ax.texts))
            self.assertTrue(any(line.get_linestyle() == "--" for line in ax.lines))
            # 3 ordinary markers, a CI band, a red ring and an orange square.
            self.assertEqual(len(ax.collections), 6)
        finally:
            for fig in figures:
                plt.close(fig)

    def test_no_ir_and_constant_counts_still_export(self):
        with tempfile.TemporaryDirectory() as td:
            summary = self.summary.copy()
            summary["IRb_bp"] = np.nan
            matrix = self.matrix * 0
            data = scatter_data(summary, matrix, highlight_gene=None)
            figures = scatter_plots(data, Path(td))
            for figure in figures.values():
                for path in figure.values():
                    self.assertTrue(Path(path).exists())


@unittest.skipUnless(shutil.which("blastn"), "blastn not on PATH")
class TestCorrelationPipeline(unittest.TestCase):
    def test_export_tables_excel_and_report_without_tree(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            for i, length in enumerate((10000, 12000, 14000), 1):
                rec = SeqRecord(Seq("ACGT" * (length // 4)), id=f"g{i}",
                                annotations={"molecule_type": "DNA", "organism": f"Species {i}"})
                for start, end, name in [(100, 1162, "psbA")] + ([(2000, 3434, "rbcL")] if i > 1 else []):
                    for kind in ("gene", "CDS"):
                        rec.features.append(SeqFeature(FeatureLocation(start, end, strand=1), type=kind,
                                                       qualifiers={"gene": [name]}))
                SeqIO.write(rec, root / f"g{i}.gb", "genbank")
            meta = root / "meta.tsv"
            meta.write_text("accession\tspecies\tgroup\ng1\tSpecies 1\tA\ng2\tSpecies 2\tB\ng3\tSpecies 3\tB\n")
            out = root / "out"
            main(["--genbank", str(root), "--meta", str(meta), "--ir-source", "genbank", "--maps", "0",
                  "--highlight-gene", "rbcL", "--highlight-genomes", "g1", "--gene-count", "unique", "--out", str(out)])
            for name in ("genes_vs_genome", "genome_vs_ir", "genes_vs_ir"):
                for ext in ("png", "pdf", "svg"):
                    self.assertTrue((out / "figures" / f"{name}.{ext}").exists())
            table = pd.read_csv(out / "tables/correlations.tsv", sep="\t").set_index("comparison")
            self.assertEqual(table.loc["genes_vs_genome", "n"], 3)
            self.assertEqual(table.loc["genome_vs_ir", "n"], 0)
            points = pd.read_csv(out / "tables/scatter_points.tsv", sep="\t").set_index("genome")
            self.assertEqual(points.loc["g1", "highlight_copies"], 0)
            self.assertEqual(points.loc["g2", "gene_count"], 2)
            with pd.ExcelFile(out / "panplast_results.xlsx") as book:
                self.assertTrue({"Scatterplot data", "Correlations"} <= set(book.sheet_names))
            html = (out / "report.html").read_text()
            self.assertTrue('id="correlations"' in html)
            self.assertLess(html.index('<section id="correlations">'), html.index('<section id="cats">'))
            self.assertTrue('href="correlations.html"' in html)
            standalone = (out / "correlations.html").read_text()
            self.assertEqual(standalone.count('src="data:image/png;base64,'), 3)
            self.assertTrue('href="report.html"' in standalone)
            self.assertTrue('id="phylogeny"' not in standalone)
            self.assertTrue("without phylogenetic correction" in html)
            self.assertFalse('id="phylogeny"' in html)


if __name__ == "__main__":
    unittest.main()
