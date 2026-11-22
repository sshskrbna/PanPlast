"""Shared species/group colouring, with and without metadata."""
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
from matplotlib.colors import to_rgba

from panplast.cli import main, parse_args
from panplast.correlations import scatter_data
from panplast.phylogeny import build_overlay, plot_overlay
from panplast.styling import taxon_styles


class TestTaxonStyles(unittest.TestCase):
    def test_same_species_shares_color_across_different_groups(self):
        styles = taxon_styles(["g1", "g2", "g3"], {"g1": "Alpha species", "g2": "Alpha_species", "g3": "Beta species"},
                              {"g1": "A", "g2": "B", "g3": "A"}, color_by="species").set_index("genome")
        self.assertEqual(styles.at["g1", "color"], styles.at["g2", "color"])
        self.assertNotEqual(styles.at["g1", "color"], styles.at["g3", "color"])

    def test_without_metadata_uses_species_or_ids_and_can_hide_names(self):
        styles = taxon_styles(["g1", "g2"], {"g1": "Alpha species"}, color_by="species",
                              label_by="accession").set_index("genome")
        self.assertEqual(styles.at["g1", "species"], "Alpha species")
        self.assertEqual(styles.at["g2", "species"], "g2")
        self.assertEqual(styles.at["g1", "display_label"], "g1")
        self.assertNotEqual(styles.at["g1", "color"], styles.at["g2", "color"])

    def test_custom_colors_are_validated_for_the_selected_category(self):
        labels = {"g1": "Alpha species", "g2": "Beta species"}
        groups = {"g1": "A", "g2": "A"}
        meta = {"g1": {"color": "red"}, "g2": {"color": "blue"}}
        species = taxon_styles(labels, labels, groups, meta, color_by="species").set_index("genome")
        self.assertEqual(species.at["g1", "color"], "#ff0000ff")
        with self.assertRaisesRegex(ValueError, "Conflicting"):
            taxon_styles(labels, labels, groups, meta, color_by="group")
        uniform = taxon_styles(labels, labels, groups, meta, color_by="none")
        self.assertEqual(uniform.color.nunique(), 1)

    def test_cli_metadata_is_optional(self):
        args = parse_args(["--color-by", "species", "--label-by", "accession"])
        self.assertIsNone(args.meta)
        self.assertEqual(args.color_by, "species")
        self.assertEqual(args.label_by, "accession")

    def test_tree_subset_keeps_scatter_colors_and_matching_ignores_label_mode(self):
        from panplast.plots import plt
        labels = {"g1": "Alpha species", "g2": "Beta species", "g3": "Gamma species"}
        styles = taxon_styles(labels, labels, color_by="species", label_by="accession")
        matrix = pd.DataFrame({"g1": [1, 1], "g2": [0, 1], "g3": [1, 1]}, index=["infA", "rbcL"])
        summary = pd.DataFrame([dict(genome=g, species=labels[g], length_bp=150000, IRb_bp=20000) for g in labels])
        points = scatter_data(summary, matrix, styles=styles).set_index("genome")
        with tempfile.TemporaryDirectory() as td:
            tree = Path(td) / "test.treefile"
            tree.write_text("(Beta_species:0.1,Gamma_species:0.2,outgroup:0.3);")
            overlay = build_overlay(tree, matrix, labels, styles=styles)
            taxa = overlay.taxa[overlay.taxa.status == "matched"].set_index("genome")
            self.assertEqual(taxa.at["g2", "color"], points.at["g2", "color"])
            self.assertEqual(taxa.at["g3", "color"], points.at["g3", "color"])
            self.assertEqual(taxa.at["g2", "display_label"], "g2")
            with patch("panplast.plots._save") as save:
                plot_overlay(overlay, summary, Path(td) / "plot")
                fig = save.call_args.args[0]
                try:
                    np.testing.assert_allclose(fig.axes[0].collections[0].get_facecolors()[0],
                                               to_rgba(points.at["g2", "color"]))
                    shown = [text.get_text() for text in fig.axes[1].texts]
                    self.assertEqual(shown, ["g2", "g3", "outgroup [no data]"])
                finally:
                    plt.close(fig)


@unittest.skipUnless(shutil.which("blastn"), "blastn not on PATH")
class TestMetadataColorPipeline(unittest.TestCase):
    def test_species_colors_with_metadata_and_without_it(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            labels = {f"g{i}": name for i, name in enumerate(("Alpha species", "Beta species", "Gamma species"), 1)}
            for genome, species in labels.items():
                rec = SeqRecord(Seq("ACGT" * 3000), id=genome,
                                annotations={"molecule_type": "DNA", "organism": species})
                for kind in ("gene", "CDS"):
                    rec.features.append(SeqFeature(FeatureLocation(100, 1534, strand=1), type=kind,
                                                   qualifiers={"gene": ["rbcL"]}))
                SeqIO.write(rec, root / f"{genome}.gb", "genbank")
            tree = root / "tree.nwk"
            tree.write_text("(Alpha_species:0.1,Beta_species:0.2,Gamma_species:0.3);")
            meta = root / "meta.tsv"
            meta.write_text("accession\tspecies\tgroup\tcolor\ng1\tAlpha species\tOne group\t#ff0000\n"
                            "g2\tBeta species\tOne group\t#0000ff\ng3\tGamma species\tOne group\t#00a040\n")
            for use_meta in (True, False):
                with self.subTest(use_meta=use_meta):
                    out = root / ("with_meta" if use_meta else "without_meta")
                    argv = ["--genbank", str(root), "--tree", str(tree), "--ir-source", "genbank", "--maps", "0",
                            "--color-by", "species", "--label-by", "accession", "--out", str(out)]
                    if use_meta:
                        argv += ["--meta", str(meta)]
                    main(argv)
                    style = pd.read_csv(out / "tables/taxon_styles.tsv", sep="\t").set_index("genome")
                    points = pd.read_csv(out / "tables/scatter_points.tsv", sep="\t").set_index("genome")
                    taxa = pd.read_csv(out / "tables/tree_taxa.tsv", sep="\t").set_index("genome")
                    self.assertEqual(style.color.nunique(), 3)
                    for genome in labels:
                        self.assertEqual(points.at[genome, "color"], taxa.at[genome, "color"])
                        self.assertEqual(taxa.at[genome, "display_label"], genome)
                    if use_meta:
                        self.assertEqual(style.at["g1", "color"], "#ff0000ff")
                    with pd.ExcelFile(out / "panplast_results.xlsx") as book:
                        self.assertIn("Taxon colours", book.sheet_names)


if __name__ == "__main__":
    unittest.main()
