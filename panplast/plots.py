"""Figures: IR junctions (IRscope-like), gene matrix, core/soft-core, region lengths, circular maps."""
from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, Rectangle, Wedge, Polygon

from .analysis import CATEGORY_ORDER
from .genes import FUNCTIONAL_GROUPS, functional_group_key
from .text import t
from .ir import JUNCTIONS, JUNCTION_SIDES, REGIONS

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": "#8a8f98",
    "axes.labelcolor": "#2b2f36", "xtick.color": "#5b616b", "ytick.color": "#5b616b",
    "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

INK, MUTED, GRID = "#22262d", "#6b717c", "#e3e6ea"
REGION_COLORS = {"LSC": "#4f7cac", "IRb": "#d08c3c", "IRa": "#d08c3c", "SSC": "#5a9e6f"}
CATEGORY_COLORS = {"core": "#1f5f75", "soft-core": "#5fa8bf", "shell": "#d9a441", "cloud": "#b9bec6"}
TYPE_COLORS = {"CDS": "#3d6fa3", "tRNA": "#c9733a", "rRNA": "#7b5ea7"}
GROUP_COLORS = dict(zip([g for g, _ in FUNCTIONAL_GROUPS] + ["fg_other"], [
    "#2e6b3a", "#5aa05a", "#8bbf4d", "#b5c832", "#e0c341", "#3e8f8a",
    "#c75f4e", "#e3956b", "#b04a6b", "#6f5aa8", "#9c7fc9", "#7aa6c9", "#8d6e4f", "#a9adb4", "#cfd2d6"]))


def _save(fig, path_base, formats=("png", "pdf", "svg")):
    out = {}
    for ext in formats:
        p = f"{path_base}.{ext}"
        fig.savefig(p, dpi=220 if ext == "png" else None, bbox_inches="tight", facecolor="white")
        out[ext] = p
    plt.close(fig)
    return out


def _italic_label(species: str, gid: str) -> str:
    return f"{species}  ({gid})" if species and species != gid else gid


# --------------------------------------------------------------------------- junctions
def junction_plot(order, labels, quads, jrows, path_base, window=1200):
    n = len(order)
    row_h = 0.62
    fig_h = 0.9 + n * row_h
    fig = plt.figure(figsize=(17, fig_h))
    label_w = 0.20
    gap = 0.008
    pw = (1 - label_w - 4 * gap) / 4
    top = 1 - 0.6 / fig_h
    bottom = 0.25 / fig_h
    axes = []
    for k, jn in enumerate(JUNCTIONS):
        ax = fig.add_axes([label_w + k * (pw + gap), bottom, pw, top - bottom])
        ax.set_xlim(-window, window)
        ax.set_ylim(n - 0.5, -0.5)
        ax.axis("off")
        l, r = JUNCTION_SIDES[jn]
        ax.text(0, -0.85, f"{jn}  ({l}/{r})", ha="center", va="bottom", fontsize=10, fontweight="bold", color=INK)
        ax.axvline(0, color="#30343b", lw=0.8, ls=(0, (3, 2)), zorder=1)
        axes.append(ax)
    lab_ax = fig.add_axes([0.005, bottom, label_w - 0.01, top - bottom])
    lab_ax.set_xlim(0, 1)
    lab_ax.set_ylim(n - 0.5, -0.5)
    lab_ax.axis("off")

    by_g = {}
    for r in jrows:
        by_g.setdefault((r["genome"], r["junction"]), []).append(r)

    band_h = 0.16
    for i, g in enumerate(order):
        q = quads.get(g)
        sp = labels.get(g, g)
        lab_ax.text(0.0, i - 0.12, sp, ha="left", va="center", fontsize=8, style="italic", color=INK)
        if q:
            txt = (f"{g} · {q.length:,} bp · LSC {q.region_len('LSC'):,} · "
                   f"IR {q.region_len('IRb'):,} · SSC {q.region_len('SSC'):,}").replace(",", " ")
        else:
            txt = f"{g} · {t('ir_none')}"
        lab_ax.text(0.0, i + 0.17, txt, ha="left", va="center", fontsize=6.2, color=MUTED)
        if i % 2 == 0:
            for ax in axes:
                ax.add_patch(Rectangle((-window, i - 0.5), 2 * window, 1, color="#f4f6f8", zorder=0, lw=0))
        if not q:
            continue
        for ax, jn in zip(axes, JUNCTIONS):
            l, r = JUNCTION_SIDES[jn]
            ax.add_patch(Rectangle((-window, i - band_h / 2), window, band_h, color=REGION_COLORS[l], lw=0, zorder=2))
            ax.add_patch(Rectangle((0, i - band_h / 2), window, band_h, color=REGION_COLORS[r], lw=0, zorder=2))
            ax.text(-window + 20, i, l, ha="left", va="center", fontsize=5.5, color="white", zorder=3, fontweight="bold")
            ax.text(window - 20, i, r, ha="right", va="center", fontsize=5.5, color="white", zorder=3, fontweight="bold")
            genes = sorted(by_g.get((g, jn), []), key=lambda x: x["rel_start"])
            occupied = {"+": [[], [], []], "−": [[], [], []]}   # label x-intervals per level
            char_w = 0.031 * window
            for gr in genes:
                a, b = max(gr["rel_start"], -window), min(gr["rel_end"], window)
                if b <= -window or a >= window:
                    continue
                up = gr["strand"] == "+"
                y0 = i - band_h / 2 - 0.15 if up else i + band_h / 2
                color = GROUP_COLORS.get(functional_group_key(gr["gene"]), "#999")
                ax.add_patch(Rectangle((a, y0), max(b - a, 6), 0.15, facecolor=color,
                                       edgecolor="white" if not gr["fragment"] else color, lw=0.4,
                                       hatch="////" if gr["fragment"] else None,
                                       alpha=0.55 if gr["fragment"] else 1, zorder=4))
                cx = (a + b) / 2
                name = gr["gene"] + ("*" if gr["fragment"] else "")
                if gr["relation"] == t("crosses"):
                    lab = f"{name} {gr['bp_left']}|{gr['bp_right']}"
                else:
                    lab = f"{name} {gr['distance_bp']} bp"
                half = len(lab) * char_w / 2
                cx = min(max(cx, -window + half), window - half)
                lv = 0
                for lv in range(3):
                    if all(cx + half < o0 or cx - half > o1 for o0, o1 in occupied[gr["strand"]][lv]):
                        break
                occupied[gr["strand"]][lv].append((cx - half - 0.01 * window, cx + half + 0.01 * window))
                dy = 0.085 + 0.095 * lv
                ly = (y0 - dy + 0.03) if up else (y0 + 0.15 + dy - 0.03)
                ax.text(cx, ly, lab, ha="center", va="center", fontsize=5.4, style="italic",
                        color=INK, zorder=6, clip_on=False)
    handles = [Patch(color=REGION_COLORS[r], label=r) for r in ("LSC", "IRb", "SSC")]
    handles += [Patch(facecolor="#999", alpha=0.55, hatch="////", label=t("f_fragment"))]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.05 / fig_h), ncol=4,
               frameon=False, fontsize=7.5)
    fig.text(0.995, 1 - 0.12 / fig_h, t("f_both_strands", w=window),
             ha="right", va="top", fontsize=7, color=MUTED)
    return _save(fig, path_base)


# --------------------------------------------------------------------------- regions
def region_lengths_plot(order, labels, summary, path_base):
    df = summary.set_index("genome").loc[order]
    n = len(order)
    fig, ax = plt.subplots(figsize=(10, 0.9 + 0.22 * n))
    left = np.zeros(n)
    for r in REGIONS:
        vals = df[f"{r}_bp"].fillna(0).values
        ax.barh(range(n), vals, left=left, color=REGION_COLORS[r], edgecolor="white", lw=0.8,
                height=0.72, label=r if r != "IRa" else None)
        left += vals
    ax.set_yticks(range(n), [labels.get(g, g) for g in order], fontsize=7, style="italic")
    ax.invert_yaxis()
    ax.set_xlabel(t("f_length_bp"))
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    lo = np.nanmin(df["length_bp"]) * 0.0
    ax.set_xlim(lo, np.nanmax(df["length_bp"]) * 1.01)
    ax.legend(ncol=3, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=8)
    ax.set_title(t("f_regions_title"), loc="left", fontsize=10, pad=18, color=INK)
    return _save(fig, path_base)


# --------------------------------------------------------------------------- gene matrix
def gene_heatmap(matrix, cats, order, labels, path_base, only_variable=False):
    cat_of = dict(zip(cats.gene, cats.category.astype(str)))
    genes = [g for g in cats.gene if (not only_variable or cat_of[g] != "core")]
    if not genes:
        return None
    m = matrix.loc[genes, order].clip(upper=3).values
    nG, nS = m.shape
    fig_h = max(3.4, 1.9 + 0.13 * nG)
    fig_w = max(7.0, 3.0 + 0.16 * nS)
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([2.0 / fig_w, 1.4 / fig_h, (fig_w - 2.6) / fig_w, (fig_h - 2.4) / fig_h])
    cmap = ListedColormap(["#f1f2f4", "#8fbcd4", "#2f6f95", "#163c57"])
    ax.imshow(m, aspect="auto", cmap=cmap, norm=BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5], 4), interpolation="nearest")
    ax.set_xticks(range(nS), [labels.get(g, g) for g in order], rotation=90, fontsize=6, style="italic")
    ax.set_yticks(range(nG), genes, fontsize=6, style="italic")
    ax.set_xticks(np.arange(-.5, nS), minor=True)
    ax.set_yticks(np.arange(-.5, nG), minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    # category band
    bx = fig.add_axes([ax.get_position().x0 - 0.75 / fig_w, ax.get_position().y0, 0.12 / fig_w, ax.get_position().height])
    cvals = [CATEGORY_ORDER.index(cat_of[g]) for g in genes]
    bx.imshow(np.array(cvals)[:, None], aspect="auto", interpolation="nearest",
              cmap=ListedColormap([CATEGORY_COLORS[c] for c in CATEGORY_ORDER]), vmin=-.5, vmax=3.5)
    bx.axis("off")
    h1 = [Patch(color=c, label=l) for c, l in zip(cmap.colors, [t("f_absent"), t("f_1copy"), t("f_2copy"), t("f_3copy")])]
    h2 = [Patch(color=CATEGORY_COLORS[c], label=c) for c in CATEGORY_ORDER]
    l1 = fig.legend(handles=h1, loc="upper left", ncol=4, frameon=False, fontsize=7, title=t("f_copies"),
                    title_fontsize=7, alignment="left", bbox_to_anchor=(0.01, 1 - 0.02 / fig_h))
    fig.legend(handles=h2, loc="upper left", ncol=4, frameon=False, fontsize=7, title=t("f_cat_band"),
               title_fontsize=7, alignment="left", bbox_to_anchor=(0.45, 1 - 0.02 / fig_h))
    fig.add_artist(l1)
    title = t("f_var_title") if only_variable else t("f_all_title")
    fig.text(0.01, 1 - 0.62 / fig_h, title, fontsize=10, color=INK, va="top")
    return _save(fig, path_base)


# --------------------------------------------------------------------------- pan-genome summary
def pangenome_plot(cats, matrix, rare_core, rare_pan, path_base, thresholds):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), gridspec_kw=dict(width_ratios=[1.1, 1.2, 1.2], wspace=0.32))
    ax = axs[0]
    tab = cats.groupby(["category", "type"], observed=False).size().unstack(fill_value=0).reindex(CATEGORY_ORDER).fillna(0)
    left = np.zeros(len(CATEGORY_ORDER))
    for tp in ("CDS", "tRNA", "rRNA"):
        vals = tab[tp].values if tp in tab else np.zeros(4)
        ax.barh(range(4), vals, left=left, color=TYPE_COLORS[tp], edgecolor="white", lw=1, height=0.62, label=tp)
        left += vals
    core_t, soft_t, shell_t = thresholds
    ylabels = [f"core\n(100 %)" if core_t >= 1 else f"core\n(≥{core_t:.0%})", f"soft core\n({soft_t:.0%}–{core_t:.0%})",
               f"shell\n({shell_t:.0%}–{soft_t:.0%})", f"cloud\n(<{shell_t:.0%})"]
    ax.set_yticks(range(4), ylabels, fontsize=8)
    ax.invert_yaxis()
    for k, v in enumerate(left):
        ax.text(v + max(left) * 0.015, k, f"{int(v)}", va="center", fontsize=9, color=INK, fontweight="bold")
    ax.set_xlim(0, max(left) * 1.15 + 1)
    ax.set_xlabel(t("f_ngenes"))
    ax.legend(frameon=False, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=8)
    ax.set_title(t("f_cats_title"), loc="left", pad=20, color=INK)

    ax = axs[1]
    n = matrix.shape[1]
    counts = cats.n_genomes.value_counts().reindex(range(1, n + 1), fill_value=0)
    colors = []
    for k in counts.index:
        f = k / n
        colors.append(CATEGORY_COLORS["core" if f >= core_t - 1e-9 else "soft-core" if f >= soft_t else
                                      "shell" if f >= shell_t else "cloud"])
    ax.bar(counts.index, counts.values, color=colors, width=0.8, edgecolor="white", lw=0.5)
    ax.set_xlabel(t("f_in_k"))
    ax.set_ylabel(t("f_ngenes"))
    ax.set_xlim(0.3, n + 0.7)
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(t("f_spectrum"), loc="left", pad=20, color=INK)

    ax = axs[2]
    x = np.arange(1, n + 1)
    for arr, c, lab in ((rare_pan, "#d08c3c", t("f_pan")), (rare_core, "#1f5f75", "core")):
        mu = arr.mean(0)
        lo, hi = np.percentile(arr, [2.5, 97.5], axis=0)
        ax.fill_between(x, lo, hi, color=c, alpha=0.18, lw=0)
        ax.plot(x, mu, color=c, lw=2, label=lab)
        ax.text(n + 0.6, mu[-1], f"{mu[-1]:.0f}", color=c, va="center", fontsize=8, fontweight="bold")
    ax.set_xlabel(t("f_ngenomes"))
    ax.set_ylabel(t("f_ngenes"))
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=8)
    ax.set_title(t("f_rare"), loc="left", pad=20, color=INK)
    return _save(fig, path_base)


# --------------------------------------------------------------------------- circular map
def circular_map(gid, label, seq, loci, q, path_base, formats=("png", "pdf")):
    L = len(seq)
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-1.32, 1.32)
    ax.set_ylim(-1.32, 1.32)
    ax.set_aspect("equal")
    ax.axis("off")

    def ang(pos):   # degrees, clockwise from top
        return 90 - 360 * pos / L

    r0 = 0.82
    ax.add_patch(Wedge((0, 0), r0 + 0.003, 0, 360, width=0.006, color="#30343b"))
    for x in loci:
        if x.flags & {"unresolved", "excluded"}:
            continue
        col = GROUP_COLORS.get(functional_group_key(x.name), "#bbb")
        out = x.strand == 1
        r_out = r0 + 0.075 if out else r0
        for s, e in x.parts:
            a1, a2 = ang(e), ang(s)
            if a2 - a1 < 0.25:
                a1 = a2 - 0.25
            ax.add_patch(Wedge((0, 0), r_out, a1, a2, width=0.075, facecolor=col,
                               edgecolor="none", alpha=0.45 if not x.counted else 1))
    # labels
    seen = []
    for x in loci:
        if not x.counted:
            continue
        mid = (x.parts[0][0] + x.parts[0][1]) / 2 if "trans_spliced" in x.flags else (x.start + x.end) / 2
        a = math.radians(ang(mid))
        out = x.strand == 1
        rr = r0 + 0.095 if out else r0 - 0.095
        deg = math.degrees(a) % 360
        flip = 90 < deg < 270
        rot = deg + 180 if flip else deg
        ha = ("right" if flip else "left") if out else ("left" if flip else "right")
        ax.text(rr * math.cos(a), rr * math.sin(a), x.name, rotation=rot, rotation_mode="anchor",
                ha=ha, va="center", fontsize=4.6, style="italic", color=INK)
    # GC content ring
    win = max(200, L // 600)
    gcs = []
    for i in range(0, L, win):
        w = seq[i:i + win].upper()
        acgt = sum(w.count(c) for c in "ACGT") or 1
        gcs.append((w.count("G") + w.count("C")) / acgt)
    gcs = np.array(gcs)
    rg = 0.52
    ax.add_patch(Wedge((0, 0), rg + 0.13, 0, 360, width=0.13, color="#f2f3f5"))
    for k, v in enumerate(gcs):
        a1, a2 = ang(min((k + 1) * win, L)), ang(k * win)
        ax.add_patch(Wedge((0, 0), rg + 0.13 * v / max(gcs.max(), 1e-9), a1, a2, width=0.13 * v / max(gcs.max(), 1e-9),
                           color="#7d8590", lw=0))
    # quadripartite ring
    if q:
        rq = 0.70
        for reg in REGIONS:
            s, e = q.regions[reg]
            spans = [(s, e)] if e >= s else [(s, L), (0, e)]
            for a, b in spans:
                ax.add_patch(Wedge((0, 0), rq, ang(b), ang(a), width=0.035, color=REGION_COLORS[reg]))
            midp = (s + q.region_len(reg) / 2) % L
            th = math.radians(ang(midp))
            ax.text((rq - 0.07) * math.cos(th), (rq - 0.07) * math.sin(th), reg, ha="center", va="center",
                    fontsize=8, fontweight="bold", color=REGION_COLORS[reg])
    ax.text(0, 0.07, label, ha="center", va="center", fontsize=11, style="italic", color=INK)
    ax.text(0, -0.02, gid, ha="center", va="center", fontsize=8, color=MUTED)
    gc_total = 100 * (seq.upper().count("G") + seq.upper().count("C")) / max(1, sum(seq.upper().count(c) for c in "ACGT"))
    ax.text(0, -0.1, f"{L:,} {t('f_bp')} · GC {gc_total:.1f} %".replace(",", " "), ha="center", va="center", fontsize=8, color=MUTED)
    used = {functional_group_key(x.name) for x in loci if x.counted}
    handles = [Patch(color=GROUP_COLORS[k], label=t(k)) for k in GROUP_COLORS if k in used]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.01, -0.01), fontsize=6, frameon=False, ncol=1)
    return _save(fig, path_base, formats)
