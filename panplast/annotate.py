"""Annotation of raw FASTA plastomes.

Two back-ends:
  * PGA (Qu et al. 2019) if `--pga /path/to/PGA.pl` is given;
  * built-in reference transfer by blastn (default): every gene of the
    reference GenBank(s) is searched in the new genome and placed where it
    hits. It reproduces gene presence/copy number and positions (for the
    boundary plot) well, but exon/intron borders are not refined, so for
    submission-quality annotation use GeSeq or PGA.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from Bio.SeqRecord import SeqRecord

from .genes import genes_from_record, load_synonyms, read_genbank
from .text import t


def build_reference_library(ref_files: list[Path], syn) -> list[tuple[str, str, str, str]]:
    """(query_id, gene_name, gene_type, sequence) for every counted reference locus."""
    lib = []
    for rf in ref_files:
        rec = read_genbank(rf)
        seq = rec.seq
        for i, loc in enumerate(genes_from_record(rec, syn)):
            if not loc.counted:
                continue
            spans = loc.parts if loc.flags & {"trans_spliced", "wraps_origin"} else [(loc.start, loc.end)]
            for j, (s, e) in enumerate(spans):
                sub = seq[s:e]
                if loc.strand == -1:
                    sub = sub.reverse_complement()
                if len(sub) < 30:
                    continue
                part = f"p{j + 1}" if len(spans) > 1 else ""
                lib.append((f"{rec.id}|{i}{part}|{loc.name}", loc.name, loc.gtype, str(sub)))
    return lib


def transfer_annotate(fasta_rec: SeqRecord, lib, min_cov=0.6, min_ident=75.0) -> SeqRecord:
    L = len(fasta_rec.seq)
    target = str(fasta_rec.seq) + str(fasta_rec.seq)[:20000]   # circular wrap
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.fa"
        t = Path(td) / "t.fa"
        q.write_text("".join(f">{qid}\n{s}\n" for qid, _, _, s in lib))
        t.write_text(f">t\n{target}\n")
        cmd = ["blastn", "-query", str(q), "-subject", str(t), "-evalue", "1e-5", "-word_size", "11",
               "-outfmt", "6 qseqid qstart qend qlen sstart send pident bitscore"]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    meta = {qid: (name, gtype) for qid, name, gtype, _ in lib}
    hits = []
    for line in out.splitlines():
        qid, qs, qe, ql, ss, se, pid, bits = line.split("\t")
        qs, qe, ql, ss, se = int(qs), int(qe), int(ql), int(ss), int(se)
        cov = (qe - qs + 1) / ql
        if cov < min_cov or float(pid) < min_ident:
            continue
        strand = 1 if ss <= se else -1
        s0, s1 = (ss, se) if strand == 1 else (se, ss)
        # extend to full reference length where the alignment is partial
        lpad, rpad = qs - 1, ql - qe
        if strand == 1:
            s0, s1 = s0 - lpad, s1 + rpad
        else:
            s0, s1 = s0 - rpad, s1 + lpad
        s0 = max(1, s0)
        if s0 > L:   # hit fully in the wrap copy -> duplicate of an earlier one
            continue
        name, gtype = meta[qid]
        part = "p" in qid.split("|")[1]
        hits.append(dict(name=name, gtype=gtype, s=s0 - 1, e=s1, strand=strand, pid=float(pid),
                         cov=cov, bits=float(bits), ref=qid.split("|")[0], part=part))
    # best hit per locus: greedy by bitscore, one per overlapping region per gene name
    hits.sort(key=lambda h: -h["bits"])
    placed = []
    for h in hits:
        clash = False
        for p in placed:
            ov = min(h["e"], p["e"]) - max(h["s"], p["s"])
            if ov > 0.5 * min(h["e"] - h["s"], p["e"] - p["s"]):
                if p["name"] == h["name"] or (p["gtype"] == h["gtype"] and ov > 0.8 * (h["e"] - h["s"])):
                    clash = True
                    break
        if not clash:
            placed.append(h)
    feats = []
    by_name = {}
    for h in placed:
        if h["part"]:   # trans-spliced parts are joined per gene later
            by_name.setdefault(h["name"], []).append(h)
            continue
        feats.extend(_feature(h, L))
    for name, parts in by_name.items():
        locs = [FeatureLocation(p["s"] % L, min(p["e"], L), strand=p["strand"]) for p in parts]
        loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)
        note = f"transferred by blastn from {parts[0]['ref']}; trans-spliced parts"
        feats.append(SeqFeature(loc, type="gene", qualifiers={"gene": [name], "annotator": ["panplast-blast"], "note": [note]}))
        feats.append(SeqFeature(loc, type=parts[0]["gtype"], qualifiers={"gene": [name]}))
    rec = SeqRecord(Seq(str(fasta_rec.seq)), id=fasta_rec.id, name=fasta_rec.id[:16],
                    description=fasta_rec.description,
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    rec.features = sorted(feats, key=lambda f: int(f.location.start))
    return rec


def _feature(h, L):
    s, e = h["s"], h["e"]
    if e > L:   # wraps the origin
        loc = CompoundLocation([FeatureLocation(s, L, strand=h["strand"]),
                                FeatureLocation(0, e - L, strand=h["strand"])])
    else:
        loc = FeatureLocation(s, e, strand=h["strand"])
    note = f"transferred by blastn from {h['ref']}, identity {h['pid']:.1f}%, coverage {100 * h['cov']:.0f}%"
    return [SeqFeature(loc, type="gene", qualifiers={"gene": [h["name"]], "annotator": ["panplast-blast"], "note": [note]}),
            SeqFeature(loc, type=h["gtype"], qualifiers={"gene": [h["name"]]})]


def annotate_fastas(fastas: list[Path], refs: list[Path], outdir: Path, syn_path=None,
                    pga: str | None = None, log=print) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    if pga:
        return _run_pga(fastas, refs, outdir, pga, log)
    syn = load_synonyms(syn_path)
    lib = build_reference_library(refs, syn)
    log(t("log_lib", n=len(lib), r=len(refs)))
    out = []
    for fa in fastas:
        rec = SeqIO.read(fa, "fasta")
        ann = transfer_annotate(rec, lib)
        flat = ann.features
        path = outdir / f"{fa.stem}.gb"
        SeqIO.write(ann, path, "genbank")
        out.append(path)
        log(t("log_transferred", f=fa.name, n=sum(f.type == "gene" for f in flat)))
    return out


def _run_pga(fastas, refs, outdir, pga, log):
    """Run PGA.pl (needs perl + BLAST+). Not exercised in tests of this pipeline."""
    with tempfile.TemporaryDirectory() as td:
        rdir, tdir = Path(td) / "ref", Path(td) / "target"
        rdir.mkdir()
        tdir.mkdir()
        for r in refs:
            (rdir / r.name).write_bytes(r.read_bytes())
        for f in fastas:
            (tdir / f"{f.stem}.fasta").write_bytes(f.read_bytes())
        cmd = ["perl", pga, "-r", str(rdir), "-t", str(tdir), "-o", str(outdir)]
        log("  " + " ".join(cmd))
        subprocess.run(cmd, check=True)
    return sorted(outdir.glob("*.gb"))
