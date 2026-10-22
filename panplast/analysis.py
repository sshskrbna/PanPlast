"""Per-genome summaries, junction genes and core / soft-core statistics."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .genes import Locus, functional_group
from .text import t
from .ir import JUNCTIONS, JUNCTION_SIDES, REGIONS, Quadripartite, gc, region_seq

CATEGORY_ORDER = ["core", "soft-core", "shell", "cloud"]


def _rel(pos: int, j: int, L: int) -> int:
    """Signed circular offset of pos relative to junction j."""
    d = (pos - j) % L
    return d - L if d > L // 2 else d


def junction_genes(loci: list[Locus], q: Quadripartite, window: int = 3000):
    """Genes crossing or flanking each junction (IRscope-like)."""
    rows = []
    L = q.length
    for jn in JUNCTIONS:
        j = q.junctions[jn]
        left_reg, right_reg = JUNCTION_SIDES[jn]
        near = []
        for loc in loci:
            if loc.flags & {"unresolved", "excluded"}:
                continue
            spans = loc.parts if loc.flags & {"trans_spliced", "wraps_origin"} else [(loc.start, loc.end)]
            for s, e in spans:
                a = _rel(s, j, L)
                b = a + (e - s)
                if b < -window or a > window:
                    continue
                if a < 0 < b:
                    rel, dist = t("crosses"), 0
                elif b <= 0:
                    rel, dist = t("in_region", r=left_reg), -b
                else:
                    rel, dist = t("in_region", r=right_reg), a
                near.append(dict(junction=jn, gene=loc.name, fragment=bool(loc.flags & {"fragment", "partial", "pseudo"}),
                                 strand="+" if loc.strand == 1 else "−", rel_start=a, rel_end=b,
                                 relation=rel, distance_bp=dist,
                                 bp_left=max(0, min(b, 0) - a) if a < 0 else 0,
                                 bp_right=max(0, b - max(a, 0)) if b > 0 else 0,
                                 gene_start=s + 1, gene_end=e))
        # keep crossing genes + the nearest 2 on each side
        cross = [n for n in near if n["relation"] == t("crosses")]
        left = sorted([n for n in near if n["rel_end"] <= 0], key=lambda n: n["distance_bp"])[:2]
        right = sorted([n for n in near if n["rel_start"] >= 0], key=lambda n: n["distance_bp"])[:2]
        rows += cross + left + right
    return rows


def genome_summary(gid, label, seq, loci, q: Quadripartite | None, q_gb: Quadripartite | None, src):
    counted = [x for x in loci if x.counted]
    uniq = {x.name for x in counted}
    by_type = {t: {x.name for x in counted if x.gtype == t} for t in ("CDS", "tRNA", "rRNA")}
    dup = set()
    if q:
        for x in counted:
            if q.region_of(x.start) in ("IRa", "IRb") and sum(y.name == x.name for y in counted) >= 2:
                dup.add(x.name)
    n_bp = seq.upper().count("N")
    row = dict(genome=gid, species=label, source=src, length_bp=len(seq), GC_total=round(gc(seq), 2),
               N_bp=n_bp, N_percent=round(100 * n_bp / max(1, len(seq)), 2))
    for r in REGIONS:
        row[f"{r}_bp"] = q.region_len(r) if q else np.nan
    for r in REGIONS:
        row[f"GC_{r}"] = round(gc(region_seq(seq, q.regions[r])), 2) if q else np.nan
    row["IR_identity_%"] = q.ir_identity if q else np.nan
    if q_gb:
        row["IR_bp_GenBank"] = q_gb.region_len("IRb")
        row["ΔIR_vs_GenBank_bp"] = (q.region_len("IRb") - q_gb.region_len("IRb")) if q else np.nan
    row.update({"genes_unique": len(uniq), "PCG_unique": len(by_type["CDS"]),
                "tRNA_unique": len(by_type["tRNA"]), "rRNA_unique": len(by_type["rRNA"]),
                "genes_total_with_IR_copies": len(counted), "duplicated_in_IR": len(dup),
                "fragments": sum("fragment" in x.flags for x in loci),
                "ambiguous_names": sum(bool(x.flags & {"ambiguous", "name_conflict"}) for x in loci)})
    return row


def presence_matrix(loci_by_genome: dict[str, list[Locus]]) -> pd.DataFrame:
    genomes = list(loci_by_genome)
    counts = {}
    for g, loci in loci_by_genome.items():
        last_end = {}
        for x in sorted((x for x in loci if x.counted), key=lambda x: x.start):
            key = (x.name, x.strand)
            if key in last_end and x.start - last_end[key] < 3000:
                last_end[key] = max(last_end[key], x.end)
                continue   # piece of the same copy (split exons etc.)
            last_end[key] = x.end
            counts.setdefault(x.name, {}).setdefault(g, 0)
            counts[x.name][g] += 1
    m = pd.DataFrame(counts).T.reindex(columns=genomes).fillna(0).astype(int)
    return m


def categorize(m: pd.DataFrame, core=1.0, soft=0.95, shell=0.15) -> pd.DataFrame:
    n = m.shape[1]
    present = (m > 0).sum(axis=1)
    frac = present / n

    def cat(f):
        if f >= core - 1e-9:
            return "core"
        if f >= soft:
            return "soft-core"
        if f >= shell:
            return "shell"
        return "cloud"

    from .genes import gene_type_from_name
    df = pd.DataFrame({
        "gene": m.index,
        "type": [gene_type_from_name(g) for g in m.index],
        "functional_group": [functional_group(g) for g in m.index],
        "n_genomes": present.values,
        "percent": (100 * frac).round(1).values,
        "category": [cat(f) for f in frac],
        "copies_min": m.where(m > 0).min(axis=1).fillna(0).astype(int).values,
        "copies_max": m.max(axis=1).values,
        "absent_in": [", ".join(m.columns[m.loc[g] == 0]) for g in m.index],
    })
    df["category"] = pd.Categorical(df["category"], CATEGORY_ORDER, ordered=True)
    return df.sort_values(["category", "type", "gene"]).reset_index(drop=True)


def rarefaction(m: pd.DataFrame, perms=200, seed=1):
    pres = (m.values > 0)
    n = pres.shape[1]
    rng = np.random.default_rng(seed)
    core = np.zeros((perms, n))
    pan = np.zeros((perms, n))
    for p in range(perms):
        order = rng.permutation(n)
        acc_core = np.ones(pres.shape[0], bool)
        acc_pan = np.zeros(pres.shape[0], bool)
        for k, gi in enumerate(order):
            acc_core &= pres[:, gi]
            acc_pan |= pres[:, gi]
            core[p, k] = acc_core.sum()
            pan[p, k] = acc_pan.sum()
    return core, pan


def gene_lengths(loci_by_genome) -> pd.DataFrame:
    rows = []
    for g, loci in loci_by_genome.items():
        for x in loci:
            rows.append(dict(genome=g, gene=x.name, raw_name=x.raw, type=x.gtype,
                             start=x.start + 1, end=x.end, strand="+" if x.strand == 1 else "−",
                             length=x.length, counted=x.counted, flags=",".join(sorted(x.flags)),
                             candidates=" | ".join(x.candidates), annotator=x.annotator))
    df = pd.DataFrame(rows)
    med = df[df.counted].groupby("gene")["length"].median()
    df["median_len_all_genomes"] = df["gene"].map(med)
    short = df.counted & (df["length"] < 0.6 * df["median_len_all_genomes"]) & (df["type"] == "CDS")
    df.loc[short, "flags"] = df.loc[short, "flags"].str.cat(["short_vs_median"] * short.sum(), sep=",").str.strip(",")
    return df


def orientation_issues(loci: list[Locus], q: Quadripartite) -> list[str]:
    """Usual angiosperm layout: psbA/trnH sit at JLA (IRa/LSC), ndhF at the IRb/SSC side.
    The opposite means the LSC or SSC was deposited reverse-complemented."""
    out = []
    L = q.length

    def dist(name, jn):
        ds = [min(abs(_rel(x.start, q.junctions[jn], L)), abs(_rel(x.end, q.junctions[jn], L)))
              for x in loci if x.name == name and x.counted]
        return min(ds) if ds else None

    a, b = dist("psbA", "JLA"), dist("psbA", "JLB")
    if a is not None and b is not None and b < a:
        out.append(t("qc_lsc_rev"))
    a, b = dist("ndhF", "JSB"), dist("ndhF", "JSA")
    if a is not None and b is not None and b < a:
        out.append(t("qc_ssc_rev"))
    return out
