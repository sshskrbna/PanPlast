"""Tests on synthetic plastomes. Run with `python -m pytest` or `python -m unittest`."""
import random
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from panplast import analysis
from panplast.genes import dedupe, extract_loci, genes_from_record, load_synonyms, read_genbank, resolve_name
from panplast.ir import detect_quadripartite

HAVE_BLAST = shutil.which("blastn") is not None
SYN = load_synonyms(None)


def random_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def synthetic_plastome(seed=1, lsc=9000, ir=3000, ssc=2500):
    """LSC + IRb + SSC + IRa with IRa = reverse complement of IRb."""
    L, I, S = random_seq(lsc, seed), random_seq(ir, seed + 100), random_seq(ssc, seed + 200)
    return L + I + S + str(Seq(I).reverse_complement()), (lsc, ir, ssc)


def gene(start, end, name, ftype="CDS", strand=1, **qual):
    loc = FeatureLocation(start, end, strand=strand)
    q = {"gene": [name], **{k: [v] for k, v in qual.items()}}
    return [SeqFeature(loc, type="gene", qualifiers=dict(q)), SeqFeature(loc, type=ftype, qualifiers={"gene": [name]})]


class TestNames(unittest.TestCase):
    def test_synonyms(self):
        self.assertEqual(resolve_name("pafI", "CDS", SYN)[0], "ycf3")
        self.assertEqual(resolve_name("pbf1", "CDS", SYN)[0], "psbN")
        self.assertEqual(resolve_name("trnfM", "tRNA", SYN)[0], "trnfM-CAU")
        self.assertEqual(resolve_name("trnH-gug", "tRNA", SYN)[0], "trnH-GUG")

    def test_fragment(self):
        name, flags, _ = resolve_name("rps12-fragment", "CDS", SYN)
        self.assertEqual(name, "rps12")
        self.assertIn("fragment", flags)

    def test_cds_candidates_drop_trna(self):
        self.assertEqual(resolve_name("trnI-GAU | ycf68", "CDS", SYN)[0], "ycf68")

    def test_aragorn_miscall_loses_to_blat(self):
        # GeSeq: ARAGORN calls intron-containing trnG-UCC "trnS-CGA"; BLAT names it trnG-UCC
        seq, _ = synthetic_plastome()
        rec = SeqRecord(Seq(seq), id="X")
        rec.features = (gene(100, 172, "trnS-CGA", "tRNA", annotator="ARAGORN") +
                        gene(101, 172, "trnG-UCC", "tRNA", annotator="blatN",
                             info="blatN_hit trnG-UCC_NC_050409.1, position 1 - 71, psl score 70.5, "
                                  "coverage 100.00%, match 100.00%"))
        loci = dedupe(extract_loci(rec, SYN))
        self.assertEqual([x.name for x in loci], ["trnG-UCC"])

    def test_trnM_vs_trnfM(self):
        seq, _ = synthetic_plastome()
        rec = SeqRecord(Seq(seq), id="X")
        rec.features = gene(200, 274, "trnM-CAU | trnfM-CAU", "tRNA", annotator="blatN, tRNAscan-SE_v2; merged",
                            info="annotated by tRNAscan-SE; blatN_hit trnfM-CAU_NC_1.1, position 1 - 74, "
                                 "psl score 100.0, coverage 100.00%, match 100.00%; merged")
        self.assertEqual(genes_from_record(rec, SYN)[0].name, "trnfM-CAU")

    def test_remnant_cds_is_partial(self):
        seq, _ = synthetic_plastome()
        rec = SeqRecord(Seq(seq), id="X")
        rec.features = gene(300, 470, "accD")      # 170 bp of a ~1500 bp gene
        loc = genes_from_record(rec, SYN)[0]
        self.assertIn("partial", loc.flags)
        self.assertFalse(loc.counted)

    def test_matK_inside_trnK_not_merged(self):
        seq, _ = synthetic_plastome()
        rec = SeqRecord(Seq(seq), id="X")
        rec.features = gene(1000, 3600, "trnK-UUU", "tRNA") + gene(1300, 2830, "matK")
        self.assertEqual(sorted(x.name for x in genes_from_record(rec, SYN)), ["matK", "trnK-UUU"])


class TestGenBankRepair(unittest.TestCase):
    def test_zero_end_location(self):
        seq, _ = synthetic_plastome()
        text = (f"LOCUS       TEST1                {len(seq)} bp    DNA     circular     01-JAN-2026\n"
                "FEATURES             Location/Qualifiers\n"
                f"     intron          {len(seq) + 1}..0\n"
                "                     /gene=\"rps12\"\n"
                "     gene            100..400\n"
                "                     /gene=\"psbA\"\n"
                "ORIGIN\n")
        body = "".join(f"{i + 1:>9} {seq[i:i + 60]}\n" for i in range(0, len(seq), 60))
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "t.gb"
            p.write_text(text + body + "//\n")
            rec = read_genbank(p)
        self.assertEqual(len(rec.seq), len(seq))


@unittest.skipUnless(HAVE_BLAST, "blastn not on PATH")
class TestInvertedRepeat(unittest.TestCase):
    def test_detect(self):
        seq, (lsc, ir, ssc) = synthetic_plastome()
        q = detect_quadripartite(seq, min_ir=1000)
        self.assertIsNotNone(q)
        self.assertAlmostEqual(q.region_len("IRb"), ir, delta=5)
        self.assertAlmostEqual(q.region_len("LSC"), lsc, delta=10)
        self.assertAlmostEqual(q.region_len("SSC"), ssc, delta=10)

    def test_ir_across_origin(self):
        seq, (lsc, ir, ssc) = synthetic_plastome()
        rot = seq[-700:] + seq[:-700]          # IRa now wraps the origin
        q = detect_quadripartite(rot, min_ir=1000)
        self.assertAlmostEqual(q.region_len("IRa"), ir, delta=5)


class TestCategories(unittest.TestCase):
    def test_thresholds(self):
        m = pd.DataFrame({f"g{i}": [1, 1 if i else 0, 1 if i < 5 else 0, 1 if i == 0 else 0] for i in range(20)},
                         index=["all", "minus_one", "quarter", "one"])
        c = analysis.categorize(m, 1.0, 0.95, 0.15).set_index("gene")["category"].astype(str)
        self.assertEqual(c["all"], "core")
        self.assertEqual(c["minus_one"], "soft-core")   # 19/20 = 95 %
        self.assertEqual(c["quarter"], "shell")
        self.assertEqual(c["one"], "cloud")


@unittest.skipUnless(HAVE_BLAST, "blastn not on PATH")
class TestEndToEnd(unittest.TestCase):
    def test_pipeline(self):
        from panplast.cli import main
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            for k in range(3):
                seq, (lsc, ir, ssc) = synthetic_plastome(seed=k + 1, ir=6000)
                rec = SeqRecord(Seq(seq), id=f"NC_00000{k}.1", annotations={"molecule_type": "DNA", "topology": "circular"})
                feats = gene(500, 1562, "psbA", strand=-1) + gene(3000, 4434, "rbcL")
                feats += gene(lsc + 200, lsc + 479, "rps19") + gene(lsc + ir + 2 * ssc // 3, lsc + ir + ssc + 200, "ndhF")
                if k:
                    feats += gene(6000, 6234, "infA")
                rec.features = feats
                SeqIO.write(rec, td / f"NC_00000{k}.1_Testus_{k}.gb", "genbank")
            out = td / "res"
            main(["--genbank", str(td), "--out", str(out), "--maps", "1", "--soft-core", "0.6"])
            cats = pd.read_csv(out / "tables" / "categories.tsv", sep="\t").set_index("gene")
            self.assertEqual(cats.at["rbcL", "category"], "core")
            self.assertEqual(cats.at["infA", "category"], "soft-core")
            self.assertTrue((out / "report.html").exists())
            self.assertTrue((out / "panplast_results.xlsx").exists())
            self.assertTrue((out / "figures" / "ir_junctions.pdf").exists())
            j = pd.read_csv(out / "tables" / "junction_genes.tsv", sep="\t")
            self.assertIn("ndhF", set(j.gene))


if __name__ == "__main__":
    unittest.main()
