"""Inverted-repeat detection by self-BLAST and the four LSC/IR/SSC junctions."""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

REGIONS = ("LSC", "IRb", "SSC", "IRa")
JUNCTIONS = ("JLB", "JSB", "JSA", "JLA")   # LSC/IRb, IRb/SSC, SSC/IRa, IRa/LSC
JUNCTION_SIDES = {"JLB": ("LSC", "IRb"), "JSB": ("IRb", "SSC"),
                  "JSA": ("SSC", "IRa"), "JLA": ("IRa", "LSC")}


@dataclass
class Quadripartite:
    length: int
    regions: dict[str, tuple[int, int]]   # 0-based half-open, may wrap (start > end)
    junctions: dict[str, int]             # 0-based index of the first base right of the junction
    ir_identity: float = 100.0
    source: str = "blast"

    def region_len(self, name: str) -> int:
        s, e = self.regions[name]
        return e - s if e >= s else self.length - s + e

    def region_of(self, pos: int) -> str:
        for name in REGIONS:
            s, e = self.regions[name]
            if (s <= pos < e) if s <= e else (pos >= s or pos < e):
                return name
        return "?"


def _blast_inverted(seq: str, workdir: Path) -> list[tuple[int, int, int, int, float]]:
    fa = workdir / "self.fa"
    fa.write_text(">s\n" + seq + "\n")
    cmd = ["blastn", "-query", str(fa), "-subject", str(fa), "-strand", "minus",
           "-evalue", "1e-50", "-dust", "no", "-outfmt", "6 qstart qend sstart send pident length"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    hits = []
    for line in out.splitlines():
        qs, qe, ss, se, pid, ln = line.split("\t")
        qs, qe, ss, se = int(qs), int(qe), int(ss), int(se)
        if int(ln) < 300:
            continue
        # minus strand: ss > se. Keep only one orientation of the symmetric pair.
        if qe < se:
            hits.append((qs, qe, se, ss, float(pid)))   # q-range, s-range (ascending)
    return hits


def _chain(hits, max_gap=1500):
    """Chain collinear inverted HSPs (q increases while s decreases)."""
    hits = sorted(hits)
    best, best_len = None, 0
    for i, h in enumerate(hits):
        chain = [h]
        for h2 in hits[i + 1:]:
            last = chain[-1]
            if 0 <= h2[0] - last[1] <= max_gap and 0 <= last[2] - h2[3] <= max_gap:
                chain.append(h2)
        qlen = chain[-1][1] - chain[0][0] + 1
        if qlen > best_len:
            best, best_len = chain, qlen
    if not best:
        return None
    q0, q1 = best[0][0], best[-1][1]
    s0, s1 = best[-1][2], best[0][3]
    pid = sum(h[4] * (h[1] - h[0]) for h in best) / sum(h[1] - h[0] for h in best)
    return q0 - 1, q1, s0 - 1, s1, pid   # 0-based half-open


def detect_quadripartite(seq: str, min_ir: int = 5000) -> Quadripartite | None:
    """Find the IR pair; tries two rotations so an IR spanning the origin is found too."""
    L = len(seq)
    best = None
    with tempfile.TemporaryDirectory() as td:
        for shift in (0, L // 3):
            rot = seq[shift:] + seq[:shift]
            ch = _chain(_blast_inverted(rot, Path(td)))
            if ch is None:
                continue
            a, b, c, d, pid = ch
            if b - a < min_ir or b > c:
                continue
            if best is None or (b - a) > (best[1] - best[0]):
                best = (a, b, c, d, pid, shift)
    if best is None:
        return None
    a, b, c, d, pid, shift = best
    # single-copy regions in rotated coordinates: S1 = [b, c), S2 = [d, a) wrapping
    s1 = c - b
    s2 = L - d + a
    if s1 >= s2:   # S1 is LSC -> after LSC comes IR2 (IRb), then S2 (SSC), then IR1 (IRa)
        reg = {"LSC": (b, c), "IRb": (c, d), "SSC": (d, a), "IRa": (a, b)}
    else:          # S2 is LSC -> IR1 is IRb, S1 is SSC, IR2 is IRa
        reg = {"LSC": (d, a), "IRb": (a, b), "SSC": (b, c), "IRa": (c, d)}
    reg = {k: ((s + shift) % L, (e + shift) % L) for k, (s, e) in reg.items()}
    for k, (s, e) in reg.items():   # avoid 0-length looking like full wrap
        if e == 0:
            reg[k] = (s, L)
    junc = {"JLB": reg["IRb"][0], "JSB": reg["SSC"][0], "JSA": reg["IRa"][0], "JLA": reg["LSC"][0]}
    return Quadripartite(L, reg, junc, round(pid, 2), "blast")


def from_annotation(regions: dict[str, tuple[int, int]], L: int) -> Quadripartite | None:
    if not all(r in regions for r in REGIONS):
        return None
    reg = {k: (s, e) for k, (s, e) in regions.items() if k in REGIONS}
    junc = {"JLB": reg["IRb"][0], "JSB": reg["SSC"][0], "JSA": reg["IRa"][0], "JLA": reg["LSC"][0] % L}
    return Quadripartite(L, reg, junc, float("nan"), "genbank")


def gc(seq: str) -> float:
    s = seq.upper()
    acgt = sum(s.count(x) for x in "ACGT")
    return 100.0 * (s.count("G") + s.count("C")) / acgt if acgt else float("nan")


def region_seq(seq: str, span: tuple[int, int]) -> str:
    s, e = span
    return seq[s:e] if e >= s else seq[s:] + seq[:e]
