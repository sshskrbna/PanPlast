"""Cross-genome reconciliation of tRNA names by sequence homology.

ARAGORN reads the wrong anticodon of some intron-containing tRNAs
(trnG-UCC -> "trnS-CGA", trnI-GAU -> "trnT-CGU"). Where BLAT gave no better
name, the locus is compared with the well-named tRNAs of the other genomes
and takes the name of its best homolog.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from Bio.Seq import Seq

from .genes import PLASTID_TRNAS, Locus
from .text import t


def _locus_seq(seq: str, loc: Locus) -> str:
    s = "".join(seq[a:b] for a, b in loc.parts)
    return str(Seq(s).reverse_complement()) if loc.strand == -1 else s


def reconcile_trna_names(loci_by_genome: dict[str, list[Locus]], seqs: dict[str, str],
                         min_cov=0.9, log=print) -> list[dict]:
    queries, lib = [], []
    for g, loci in loci_by_genome.items():
        for i, x in enumerate(loci):
            if x.gtype != "tRNA" or not x.counted:
                continue
            if x.name in PLASTID_TRNAS:
                lib.append((f"{g}|{i}|{x.name}", _locus_seq(seqs[g], x)))
            else:
                queries.append((g, i, x))
    if not queries or not lib:
        return []
    changes = []
    with tempfile.TemporaryDirectory() as td:
        qf, sf = Path(td) / "q.fa", Path(td) / "s.fa"
        qf.write_text("".join(f">{g}|{i}\n{_locus_seq(seqs[g], x)}\n" for g, i, x in queries))
        sf.write_text("".join(f">{k}\n{s}\n" for k, s in lib))
        out = subprocess.run(["blastn", "-query", str(qf), "-subject", str(sf), "-task", "blastn", "-word_size", "7", "-evalue", "1e-5",
                              "-outfmt", "6 qseqid sseqid pident nident qlen slen bitscore"],
                             capture_output=True, text=True, check=True).stdout
    best = {}
    for line in out.splitlines():
        q, s, pid, nid, ql, sl, bits = line.split("\t")
        # identical bases over the homolog's length: tolerant to the mis-placed
        # intron borders that come with a wrong ARAGORN call
        cov = int(nid) / int(sl)
        if cov < min_cov or int(nid) / int(ql) < 0.7:
            continue
        if q not in best or float(bits) > best[q][1]:
            best[q] = (s.split("|")[2], float(bits), float(pid), cov)
    for g, i, x in queries:
        hit = best.get(f"{g}|{i}")
        if hit and hit[0] != x.name:
            changes.append(dict(genome=g, start=x.start + 1, old=x.name, new=hit[0],
                                identity=round(hit[2], 1), coverage=round(100 * hit[3], 1)))
            x.candidates = list(dict.fromkeys([hit[0], x.name] + x.candidates))
            x.name = hit[0]
            x.flags.add("renamed_by_homology")
        elif not hit:
            x.flags.add("nonstandard_tRNA")
    if changes:
        log(t("log_renamed", n=len(changes)))
    return changes
