"""Gene extraction from GenBank records and name normalisation.

GeSeq / PGA / NCBI name the same plastid gene in different ways
("pafI" vs "ycf3", "trnM-CAU | trnfM-CAU", "rps12-fragment", ...).
Everything here turns a GenBank record into a clean list of loci with
one canonical name each, plus QC flags explaining every decision.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

from .text import t

# Default synonyms: GeSeq/Chloë names -> traditional plastome names.
DEFAULT_SYNONYMS = {
    "pafI": "ycf3", "pafII": "ycf4", "pbf1": "psbN", "clpP1": "clpP",
    "lhbA": "psbZ", "ycf9": "psbZ", "ycf10": "cemA", "ycf5": "ccsA",
    "trnfM": "trnfM-CAU", "trnK": "trnK-UUU", "trnH": "trnH-GUG",
    "trnW": "trnW-CCA", "trnF": "trnF-GAA", "trnC": "trnC-GCA",
    "trnD": "trnD-GUC", "trnN": "trnN-GUU", "trnQ": "trnQ-UUG",
    "trnY": "trnY-GUA", "trnA": "trnA-UGC", "trnE": "trnE-UUC",
}

PREDICTORS = ("ARAGORN", "tRNAscan", "Chloe", "Chloë", "HMMER")
TRN_RE = re.compile(r"^trn(fM|[A-Z])(?:-?([ACGUT]{3}))?$", re.I)
RRN_RE = re.compile(r"^(?:rrn)?(4\.5|5|16|23)S?(?:\s*rRNA)?$", re.I)
HIT_RE = re.compile(r"blat([NX])_hit\s+(\S+?),\s*position[^,]*,\s*psl score\s+([\d.]+),\s*coverage\s+([\d.]+)%")

# tRNA genes known from land-plant plastomes. A tRNA name outside this list
# (e.g. ARAGORN calling intron-containing trnG-UCC "trnS-CGA") loses to a
# listed candidate supported by BLAT.
PLASTID_TRNAS = {
    "trnA-UGC", "trnC-GCA", "trnD-GUC", "trnE-UUC", "trnF-GAA", "trnfM-CAU", "trnG-GCC",
    "trnG-UCC", "trnH-GUG", "trnI-CAU", "trnI-GAU", "trnK-UUU", "trnL-CAA", "trnL-UAA",
    "trnL-UAG", "trnM-CAU", "trnN-GUU", "trnP-GGG", "trnP-UGG", "trnQ-UUG", "trnR-ACG",
    "trnR-CCG", "trnR-UCU", "trnS-GCU", "trnS-GGA", "trnS-UGA", "trnT-GGU", "trnT-UGU",
    "trnV-GAC", "trnV-UAC", "trnW-CCA", "trnY-GUA", "trnsU-UCA",
}

# Typical full CDS lengths (bp, Nicotiana/Oryza-like). A CDS shorter than
# --min-cds-frac of this is a remnant (e.g. accD/ycf1/ycf2 in grasses).
EXPECTED_CDS_LEN = {
    "accD": 1500, "atpA": 1524, "atpB": 1497, "atpE": 402, "atpF": 555, "atpH": 246, "atpI": 744,
    "ccsA": 942, "cemA": 690, "clpP": 591, "infA": 234, "matK": 1530, "ndhA": 1092, "ndhB": 1533,
    "ndhC": 363, "ndhD": 1503, "ndhE": 306, "ndhF": 2223, "ndhG": 531, "ndhH": 1182, "ndhI": 504,
    "ndhJ": 477, "ndhK": 678, "petA": 963, "petB": 648, "petD": 483, "petG": 114, "petL": 96,
    "petN": 90, "psaA": 2253, "psaB": 2205, "psaC": 246, "psaI": 111, "psaJ": 135, "psbA": 1062,
    "psbB": 1527, "psbC": 1422, "psbD": 1062, "psbE": 252, "psbF": 120, "psbH": 222, "psbI": 111,
    "psbJ": 123, "psbK": 186, "psbL": 117, "psbM": 105, "psbN": 132, "psbT": 105, "psbZ": 189,
    "rbcL": 1434, "rpl14": 369, "rpl16": 408, "rpl2": 822, "rpl20": 387, "rpl22": 400, "rpl23": 282,
    "rpl32": 168, "rpl33": 201, "rpl36": 114, "rpoA": 1014, "rpoB": 3213, "rpoC1": 2046,
    "rpoC2": 4176, "rps11": 417, "rps12": 372, "rps14": 303, "rps15": 264, "rps16": 267,
    "rps18": 306, "rps19": 279, "rps2": 711, "rps3": 657, "rps4": 606, "rps7": 468, "rps8": 405,
    "ycf1": 5700, "ycf2": 6840, "ycf3": 507, "ycf4": 555,
}

FUNCTIONAL_GROUPS = [
    # (text key, test)
    ("fg_psa", lambda g: g.startswith("psa")),
    ("fg_psb", lambda g: g.startswith("psb")),
    ("fg_pet", lambda g: g.startswith("pet")),
    ("fg_atp", lambda g: g.startswith("atp")),
    ("fg_ndh", lambda g: g.startswith("ndh")),
    ("fg_rbc", lambda g: g.startswith("rbc")),
    ("fg_rpo", lambda g: g.startswith("rpo")),
    ("fg_rps", lambda g: g.startswith("rps")),
    ("fg_rpl", lambda g: g.startswith("rpl")),
    ("fg_trn", lambda g: g.startswith("trn")),
    ("fg_rrn", lambda g: g.startswith("rrn")),
    ("fg_ycf34", lambda g: g in ("ycf3", "ycf4")),
    ("fg_misc", lambda g: g in ("matK", "clpP", "accD", "ccsA", "cemA", "infA")),
    ("fg_ycf", lambda g: g.startswith(("ycf", "orf"))),
]


def functional_group_key(name: str) -> str:
    for key, test in FUNCTIONAL_GROUPS:
        if test(name):
            return key
    return "fg_other"


def functional_group(name: str) -> str:
    return t(functional_group_key(name))


def gene_type_from_name(name: str) -> str:
    if name.startswith("trn"):
        return "tRNA"
    if name.startswith("rrn"):
        return "rRNA"
    return "CDS"


def load_synonyms(path: str | None) -> dict[str, str]:
    syn = dict(DEFAULT_SYNONYMS)
    if path:
        with open(path) as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if not row or row[0].startswith("#") or len(row) < 2:
                    continue
                syn[row[0].strip()] = row[1].strip()
    return syn


@dataclass
class Locus:
    name: str                      # canonical name
    raw: str                       # name as written in the GenBank
    gtype: str                     # CDS / tRNA / rRNA
    start: int                     # 0-based, inclusive (span of all parts)
    end: int                       # 0-based, exclusive
    strand: int
    parts: list[tuple[int, int]]   # exon parts, 0-based half-open
    annotator: str = ""
    flags: set[str] = field(default_factory=set)
    candidates: list[str] = field(default_factory=list)
    evidence: dict[str, float] = field(default_factory=dict)   # BLAT psl score per name
    first: str = ""                # first name in the raw string (GeSeq: predictor's call)

    @property
    def length(self) -> int:
        return sum(e - s for s, e in self.parts)

    @property
    def counted(self) -> bool:
        """Does this locus count toward gene presence/absence?"""
        return not (self.flags & {"fragment", "unresolved", "pseudo", "partial", "excluded"})


def _norm_one(name: str, syn: dict[str, str]) -> tuple[str, bool]:
    """Normalise one candidate name. Returns (name, is_fragment)."""
    s = name.strip().replace(" ", "")
    frag = False
    for suffix in ("-fragment", "_fragment", "-fr"):
        if s.endswith(suffix):
            s, frag = s[: -len(suffix)], True
    s = syn.get(s, s)
    m = TRN_RE.match(s)
    if m:
        aa, codon = m.group(1), m.group(2)
        aa = "fM" if aa.lower() == "fm" else aa.upper()
        s = f"trn{aa}" + (f"-{codon.upper().replace('T', 'U')}" if codon else "")
        s = syn.get(s, s)
    m = RRN_RE.match(s)
    if m and (s.lower().startswith("rrn") or "rrna" in s.lower()):
        s = f"rrn{m.group(1)}"
    return s, frag


def _codon(name: str) -> str | None:
    m = TRN_RE.match(name)
    return m.group(2) if m else None


def _parse_evidence(info: str, syn) -> tuple[dict[str, float], dict[str, float]]:
    """BLAT hits in GeSeq /info -> ({name: psl score}, {name: coverage %})."""
    psl, cov = {}, {}
    for m in HIT_RE.finditer(" ".join(info.split())):
        name, frag = _norm_one(m.group(2).split("_")[0], syn)
        if frag:
            continue
        psl[name] = max(psl.get(name, 0.0), float(m.group(3)))
        cov[name] = max(cov.get(name, 0.0), float(m.group(4)))
    return psl, cov


def resolve_name(raw: str, gtype: str, syn: dict[str, str]) -> tuple[str, set[str], list[str]]:
    """Split a GeSeq-style 'A | B | C' string into normalised candidates."""
    flags: set[str] = set()
    normed = [_norm_one(c, syn) for c in raw.split("|") if c.strip()]
    full = [n for n, fr in normed if not fr]
    if not full:
        return normed[0][0], {"fragment"}, list(dict.fromkeys(n for n, _ in normed))
    cands = list(dict.fromkeys(full))
    if gtype == "tRNA":
        cands = [c for c in cands if c.startswith("trn")] or cands
    elif gtype == "CDS":
        cands = [c for c in cands if not c.startswith(("trn", "rrn"))] or cands
    if len(cands) > 1:
        flags.add("ambiguous")
    return cands[0], flags, cands


def score_name(name: str, gtype: str, members: list["Locus"]) -> float:
    s = 0.0
    if gtype == "tRNA":
        s += 3.0 if name in PLASTID_TRNAS else 0.0
        s += 0.5 if _codon(name) else 0.0
    s += 2.0 * max((m.evidence.get(name, 0.0) for m in members), default=0.0) / 100.0
    if any(m.first == name and any(p.lower() in m.annotator.lower() for p in PREDICTORS) for m in members):
        s += 1.0
    if any(m.name == name and len(m.candidates) == 1 for m in members):
        s += 0.5
    return s


def _parts(feature) -> list[tuple[int, int]]:
    return sorted((int(p.start), int(p.end)) for p in feature.location.parts)


def _qual(feature, key: str) -> str:
    return feature.qualifiers.get(key, [""])[0]


def extract_loci(record, syn: dict[str, str]) -> list[Locus]:
    """All gene loci of a record, before de-duplication."""
    typed = []
    genes = []
    for f in record.features:
        if f.location is None:
            continue
        if f.type in ("CDS", "tRNA", "rRNA"):
            typed.append(f)
        elif f.type == "gene":
            genes.append(f)

    def make(f, gtype, info=""):
        raw = _qual(f, "gene") or _qual(f, "locus_tag") or _qual(f, "product")
        if not raw:
            return None
        if gtype is None:
            gtype = gene_type_from_name(_norm_one(raw.split("|")[0], syn)[0])
        name, flags, cands = resolve_name(raw, gtype, syn)
        if "pseudo" in f.qualifiers or "pseudogene" in f.qualifiers:
            flags.add("pseudo")
        parts = _parts(f)
        strand = f.location.strand or 1
        psl, _ = _parse_evidence(info or _qual(f, "info"), syn)
        return Locus(name, raw, gtype, min(p[0] for p in parts), max(p[1] for p in parts),
                     strand, parts, _qual(f, "annotator"), flags, cands, psl,
                     _norm_one(raw.split("|")[0], syn)[0])

    loci = []
    if genes:
        for g in genes:
            raw = _qual(g, "gene")
            gs, ge = int(g.location.start), int(g.location.end)
            # typed feature of the same gene; prefer identical span, then largest overlap
            match, best = None, -1
            for t in typed:
                ts, te = int(t.location.start), int(t.location.end)
                if _qual(t, "gene") == raw and ts < ge and te > gs:
                    score = (ts == gs and te == ge) * 10**9 + min(te, ge) - max(ts, gs)
                    if score > best:
                        match, best = t, score
            loc = make(match or g, match.type if match else None, _qual(g, "info"))
            if loc is None:
                continue
            if match is not None and "pseudo" in g.qualifiers:
                loc.flags.add("pseudo")
            if not loc.annotator:
                loc.annotator = _qual(g, "annotator")
            loci.append(loc)
        # typed features that have no gene feature (some NCBI records)
        for t in typed:
            raw = _qual(t, "gene")
            if raw and not any(_qual(g, "gene") == raw and int(g.location.start) < int(t.location.end)
                               and int(g.location.end) > int(t.location.start) for g in genes):
                loc = make(t, t.type)
                if loc:
                    loci.append(loc)
    else:
        for t in typed:
            loc = make(t, t.type)
            if loc:
                loci.append(loc)
    # trans-spliced genes (rps12) have parts far apart: keep the real parts,
    # but the 'span' of such a locus is meaningless -> mark it.
    L = len(record.seq)
    for loc in loci:
        if loc.end - loc.start > 20000:
            touches = any(s == 0 for s, _ in loc.parts) and any(e == L for _, e in loc.parts)
            loc.flags.add("wraps_origin" if touches and loc.length < 20000 else "trans_spliced")
    return loci


def _overlap(a: Locus, b: Locus) -> int:
    o = 0
    for s1, e1 in a.parts:
        for s2, e2 in b.parts:
            o += max(0, min(e1, e2) - max(s1, s2))
    return o


def dedupe(loci: list[Locus]) -> list[Locus]:
    """Merge loci that describe the same feature (same type, mutual overlap >= 50 %)
    and give the merged locus the best-supported name."""
    loci = sorted(loci, key=lambda x: (x.start, x.end))
    clusters: list[list[Locus]] = []
    def same(x: Locus, y: Locus) -> bool:
        if x.gtype != y.gtype or x.strand != y.strand:
            return False
        o = _overlap(x, y)
        short = min(x.length, y.length)
        if o >= 0.5 * x.length and o >= 0.5 * y.length:
            return True
        # GeSeq sometimes gives a tRNA both as exons and as one whole-span feature
        if o >= 0.9 * short and set(x.candidates) & set(y.candidates):
            return True
        # two RNA genes never overlap: same-strand overlapping tRNA/rRNA calls are one gene
        return x.gtype in ("tRNA", "rRNA") and o >= 0.5 * short

    for loc in loci:
        for cl in clusters:
            if any(same(m, loc) for m in cl):
                cl.append(loc)
                break
        else:
            clusters.append([loc])
    out = []
    for cl in clusters:
        if all("fragment" in x.flags for x in cl):
            out.append(max(cl, key=lambda x: x.length))
            continue
        good = [x for x in cl if "fragment" not in x.flags]
        gtype = good[0].gtype
        names = []
        for x in good:
            names += x.candidates + [n for n in x.evidence if gene_type_from_name(n) == gtype]
        names = list(dict.fromkeys(names))
        scored = sorted(((score_name(n, gtype, good), n) for n in names), reverse=True)
        chosen = scored[0][1]
        holders = [x for x in good if chosen in x.candidates] or good
        best = max(holders, key=lambda x: (len(x.parts), any(p.lower() in x.annotator.lower() for p in PREDICTORS), x.length))
        flags = set().union(*(x.flags for x in good)) - {"ambiguous"}
        if not all("pseudo" in x.flags for x in good):
            flags.discard("pseudo")
        if len(names) > 1:
            flags.add("renamed" if chosen != best.first else "ambiguous")
        if gtype == "tRNA" and not _codon(chosen):
            flags.add("unresolved")
        else:
            flags.discard("unresolved")
        best.name, best.flags = chosen, flags
        best.candidates = [n for _, n in scored]
        best.annotator = "; ".join(dict.fromkeys(x.annotator for x in good if x.annotator))
        out.append(best)
    return out


def genes_from_record(record, syn: dict[str, str], exclude: set[str] = frozenset(),
                      include_pseudo: bool = False, min_cds_frac: float = 0.5) -> list[Locus]:
    loci = dedupe(extract_loci(record, syn))
    for loc in loci:
        exp = EXPECTED_CDS_LEN.get(loc.name)
        if loc.gtype == "CDS" and exp and loc.length < min_cds_frac * exp:
            loc.flags.add("partial")
        if loc.name in exclude:
            loc.flags.add("excluded")
        if include_pseudo:
            loc.flags.discard("pseudo")
    return sorted(loci, key=lambda x: x.start)


def regions_from_record(record) -> dict[str, tuple[int, int]]:
    """LSC/IRb/SSC/IRa as annotated in the GenBank (GeSeq/OGDRAW, NCBI notes)."""
    out = {}
    for f in record.features:
        if f.type not in ("misc_feature", "repeat_region"):
            continue
        txt = " ".join(f.qualifiers.get("note", []) + f.qualifiers.get("rpt_type", [])).lower()
        key = None
        if "large single copy" in txt or re.search(r"\blsc\b", txt):
            key = "LSC"
        elif "small single copy" in txt or re.search(r"\bssc\b", txt):
            key = "SSC"
        elif re.search(r"(inverted repeat|\bir)\s*-?\s*b\b|irb", txt):
            key = "IRb"
        elif re.search(r"(inverted repeat|\bir)\s*-?\s*a\b|ira", txt):
            key = "IRa"
        if key and key not in out:
            parts = list(f.location.parts)
            out[key] = (int(parts[0].start), int(parts[-1].end))
    return out


_LOC_LINE = re.compile(r"^(     \S+\s+)(complement\()?(\d+)\.\.(\d+)(\)?)\s*$")


def read_genbank(path):
    """SeqIO.read with repairs for broken origin-wrapping locations (GeSeq writes '136376..0')."""
    import io
    import warnings
    from Bio import BiopythonParserWarning, SeqIO
    lines = Path(path).read_text().splitlines()
    m = re.match(r"LOCUS\s+\S+\s+(\d+)\s+bp", lines[0]) if lines else None
    L = int(m.group(1)) if m else None
    fixed = []
    for line in lines:
        mm = _LOC_LINE.match(line) if L else None
        if mm:
            pre, comp, a, b, close = mm.groups()
            a, b = int(a), int(b)
            if b == 0 and a > L:          # empty feature past the end -> drop the location line's feature
                line = f"{pre}{comp or ''}{L}..{L}{close}"
            elif b < a:                   # wraps the origin
                inner = f"join({a}..{L},1..{b})" if b > 0 else f"{a}..{L}"
                line = f"{pre}complement({inner})" if comp else f"{pre}{inner}"
        fixed.append(line)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonParserWarning)
        return SeqIO.read(io.StringIO("\n".join(fixed) + "\n"), "genbank")


ACC_RE = re.compile(r"([A-Z]{1,2}_?\d{5,9}\.\d+)")


def accession_of(text: str) -> str | None:
    m = ACC_RE.search(text)
    return m.group(1) if m else None


def find_files(paths: list[str], exts: tuple[str, ...]) -> list[Path]:
    out = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            out += sorted(x for x in p.iterdir() if x.suffix.lower() in exts and not x.name.startswith("."))
        elif p.exists():
            out.append(p)
    return out
