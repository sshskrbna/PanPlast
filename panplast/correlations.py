"""Descriptive genome-size / gene-content scatterplots and ordinary least squares."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from .genes import EXPECTED_CDS_LEN, PLASTID_TRNAS
from .styling import category_colors


COMPARISONS = (
    ("genes_vs_genome", "genome_kb", "gene_count", "Counted genes vs. plastome size"),
    ("genome_vs_ir", "ir_kb", "genome_kb", "Plastome size vs. IR length"),
    ("genes_vs_ir", "ir_kb", "gene_count", "Counted genes vs. IR length"),
)

def _group_colors(data):
    if "color_label" in data and "color" in data:
        return category_colors(data, "color_label", "color")
    return category_colors(data)


@dataclass
class Regression:
    summary: dict
    x: np.ndarray
    mean: np.ndarray
    lower: np.ndarray
    upper: np.ndarray


def scatter_data(summary, matrix, gene_count="total", highlight_gene=None, highlight_genomes=(), loci_by=None,
                 colors_by_group=None, styles=None):
    """Use counted copies from the same matrix as the main pan-genome analysis."""
    if gene_count not in {"total", "unique"}:
        raise ValueError("gene_count must be 'total' or 'unique'")
    marked = set(highlight_genomes)
    unknown = marked - set(summary.genome)
    if unknown:
        raise ValueError("Unknown --highlight-genomes IDs: " + ", ".join(sorted(unknown)))
    if highlight_gene:
        known = set(matrix.index) | set(EXPECTED_CDS_LEN) | PLASTID_TRNAS
        known |= {loc.name for loci in (loci_by or {}).values() for loc in loci}
        if highlight_gene not in known:
            raise ValueError(f"Unknown --highlight-gene name: {highlight_gene}")
    data = summary[["genome", "species"]].copy()
    data["group"] = summary["group"].fillna("").replace("", "Ungrouped") if "group" in summary else "Ungrouped"
    from matplotlib.colors import to_hex
    data["group_color"] = data.group.map(colors_by_group or {}).fillna("")
    colors = _group_colors(data)
    data["group_color"] = data.group.map({group: to_hex(color, keep_alpha=True) for group, color in colors.items()})
    if styles is not None:
        style_index = styles.set_index("genome")
        for column in ("color_by", "color_label", "color", "label_by", "display_label"):
            data[column] = data.genome.map(style_index[column])
        data["group_color"] = data["color"]  # Kept for compatibility with older exports.
    data["genome_kb"] = summary["length_bp"] / 1000
    data["ir_kb"] = summary["IRb_bp"] / 1000
    counts = matrix.sum(axis=0) if gene_count == "total" else (matrix > 0).sum(axis=0)
    data["gene_count"] = data.genome.map(counts)
    data["gene_count_mode"] = gene_count
    data["highlight_gene"] = highlight_gene or ""
    data["highlight_copies"] = (data.genome.map(matrix.loc[highlight_gene]) if highlight_gene in matrix.index
                                 else 0) if highlight_gene else np.nan
    data["marked"] = data.genome.isin(marked)
    return data


def fit_regression(x, y):
    """Paired finite observations, Pearson r/p and a pointwise 95% mean-response CI.

    Fewer than three observations or a constant variable have no reported
    inferential statistics. Missing IR lengths are never replaced by zero.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    n = len(x)
    result = dict(n=n, n_excluded=int(len(finite) - n), pearson_r=np.nan, p_value=np.nan,
                  slope=np.nan, intercept=np.nan, confidence_level=.95, status="ok")
    empty = np.array([], dtype=float)
    if n < 3:
        result["status"] = "fewer_than_3_pairs"
    elif np.ptp(x) == 0:
        result["status"] = "constant_x"
    elif np.ptp(y) == 0:
        result["status"] = "constant_y"
    if result["status"] != "ok":
        return Regression(result, empty, empty, empty, empty)
    fit = stats.linregress(x, y)
    pearson = stats.pearsonr(x, y)
    result.update(pearson_r=float(pearson.statistic), p_value=float(pearson.pvalue),
                  slope=float(fit.slope), intercept=float(fit.intercept))
    grid = np.linspace(x.min(), x.max(), 200)
    mean = fit.intercept + fit.slope * grid
    residual = y - (fit.intercept + fit.slope * x)
    variance = np.sum(residual ** 2) / (n - 2)
    sxx = np.sum((x - x.mean()) ** 2)
    se = np.sqrt(variance * (1 / n + (grid - x.mean()) ** 2 / sxx))
    half_width = stats.t.ppf(.975, df=n - 2) * se
    return Regression(result, grid, mean, mean - half_width, mean + half_width)


def correlation_results(data):
    rows = []
    for name, x, y, _ in COMPARISONS:
        fit = fit_regression(data[x], data[y])
        rows.append(dict(comparison=name, x=x, y=y, **fit.summary))
    return pd.DataFrame(rows)


def _point_labels(data):
    labels = []
    for row in data.itertuples():
        words = str(getattr(row, "display_label", None) or row.species or row.genome).split()
        label = f"{words[0][0]}. {' '.join(words[1:])}" if len(words) >= 2 else " ".join(words)
        labels.append(label)
    repeated = Counter(labels)
    return [f"{label} ({genome})" if repeated[label] > 1 else label for label, genome in zip(labels, data.genome)]


def _place_labels(ax, data, x, y, colors, reserved):
    """Choose nearby text positions that minimise label / marker overlaps."""
    from matplotlib.transforms import Bbox

    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    positions = ax.transData.transform(data[[x, y]].to_numpy(dtype=float))
    point_boxes = [Bbox.from_bounds(px - 5, py - 5, 10, 10) for px, py in positions]
    occupied = list(reserved)
    offsets = [(dx, dy) for dy in (11, -12, 23, -24, 36, -37, 50, -51) for dx in (0, 16, -16, 34, -34)]

    def area(a, b):
        return max(0, min(a.x1, b.x1) - max(a.x0, b.x0)) * max(0, min(a.y1, b.y1) - max(a.y0, b.y0))

    for row, label in zip(data.itertuples(), _point_labels(data)):
        color = colors[getattr(row, "color_label", row.group)]
        text = ax.annotate(label, (getattr(row, x), getattr(row, y)), xytext=(0, 11),
                           textcoords="offset points", ha="center", va="center", fontsize=7,
                           style="italic", color=color, zorder=5)
        best = None
        for offset in offsets:
            text.set_position(offset)
            box = text.get_window_extent(renderer).expanded(1.05, 1.15)
            overflow = max(0, box.width * box.height - area(box, ax.bbox))
            cost = (100 * overflow + 10 * sum(area(box, other) for other in occupied)
                    + sum(area(box, other) for other in point_boxes) + .03 * np.hypot(*offset))
            if best is None or cost < best[0]:
                best = cost, offset, box
        text.set_position(best[1])
        occupied.append(best[2])
        if np.hypot(*best[1]) > 16:
            ax.annotate("", (getattr(row, x), getattr(row, y)), xytext=best[1], textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", color="#9ba1aa", lw=.45, shrinkA=3, shrinkB=4), zorder=2)


def scatter_plots(data, outdir):
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    from .plots import INK, MUTED, _save, plt

    groups = sorted(data["color_label" if "color_label" in data else "group"].unique())
    colors = _group_colors(data)
    color_by = str(data.color_by.iloc[0]) if "color_by" in data and len(data) else "group"
    color_title = "Species" if color_by == "species" else "Group"
    gene = str(data.highlight_gene.iloc[0]) if len(data) else ""
    mode = str(data.gene_count_mode.iloc[0]) if len(data) else "total"
    ylabel = "Counted gene copies (including IR copies)" if mode == "total" else "Unique counted genes"
    axes_labels = {"genome_kb": "Total plastome size (kb)", "ir_kb": "IRb length (kb)", "gene_count": ylabel}
    figures = {}
    for name, x, y, title in COMPARISONS:
        finite = np.isfinite(data[x].astype(float)) & np.isfinite(data[y].astype(float))
        points = data.loc[finite]
        fit = fit_regression(data[x], data[y])
        height = max(8.2, 3.0 + .25 * len(groups) + (1 if gene else 0))
        fig, ax = plt.subplots(figsize=(13, height))
        fig.subplots_adjust(left=.08, right=.73, bottom=.15, top=.86)
        ax.set_axisbelow(True)
        ax.grid(color="#f0f0f0", lw=.65)
        if fit.summary["status"] == "ok":
            ax.fill_between(fit.x, fit.lower, fit.upper, color="#bdbdbd", alpha=.20, lw=0, zorder=1)
            ax.plot(fit.x, fit.mean, color="#666666", lw=1.0, ls=(0, (4, 4)), zorder=2)
        for row in points.itertuples():
            absent = bool(gene) and row.highlight_copies == 0
            ax.scatter(getattr(row, x), getattr(row, y), marker="^" if absent else "o", s=46,
                       color=colors[getattr(row, "color_label", row.group)], edgecolors="white", linewidths=.35, zorder=3)
            if absent:
                ax.scatter(getattr(row, x), getattr(row, y), marker="o", s=110, facecolors="none",
                           edgecolors="#c0392b", linewidths=1.7, zorder=4)
            if row.marked:
                ax.scatter(getattr(row, x), getattr(row, y), marker="s", s=150, facecolors="none",
                           edgecolors="#e67e22", linewidths=1.7, zorder=4)
        ax.margins(x=.10, y=.18)
        ax.set_xlabel(axes_labels[x], fontsize=10, fontweight="bold")
        ax.set_ylabel(axes_labels[y], fontsize=10, fontweight="bold")
        if y == "gene_count":
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        s = fit.summary
        if s["status"] == "ok":
            annotation = f"r = {s['pearson_r']:.2f}\np = {s['p_value']:.3g}\nn = {s['n']}"
        else:
            reason = "fewer than 3 valid pairs" if s["n"] < 3 else "constant variable"
            annotation = f"n = {s['n']}\nr, p: undefined\n{reason}"
        stat_text = ax.text(.025, .97, annotation, transform=ax.transAxes, va="top", fontsize=9,
                            fontweight="bold", color=INK,
                            bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=3))
        if points.empty:
            ax.text(.5, .5, "No paired data available", transform=ax.transAxes, ha="center", color=MUTED)
        group_legend = ax.legend(handles=[Line2D([], [], marker="o", ls="", color=colors[g], label=g) for g in groups],
                                 title=color_title, loc="upper left", bbox_to_anchor=(1.025, 1),
                                 frameon=False, fontsize=8, title_fontsize=9)
        ax.add_artist(group_legend)
        fig.canvas.draw()
        bottom = ax.transAxes.inverted().transform(group_legend.get_window_extent().get_points())[0, 1] - .04
        special = []
        if gene:
            special += [Line2D([], [], marker="o", ls="", color=INK, label=f"{gene}: counted copy"),
                        Line2D([], [], marker="^", ls="", color=INK, label=f"{gene}: no counted copy")]
        elif data.marked.any():
            special.append(Line2D([], [], marker="s", ls="", markerfacecolor="none", markeredgecolor="#e6862b",
                                  markeredgewidth=1.5, markersize=9, label="Selected genomes"))
        if special:
            ax.legend(handles=special, loc="upper left", bbox_to_anchor=(1.025, bottom), frameon=False,
                      title=f"{gene} status" if gene else "Selection", fontsize=8, title_fontsize=9)
        _place_labels(ax, points, x, y, colors, [stat_text.get_window_extent(fig.canvas.get_renderer())])
        fig.text(.08, .95, title, fontsize=13, fontweight="bold", color=INK, va="top")
        colour_description = "Uniform taxon colour" if color_by == "none" else f"Colour: {color_title.lower()}"
        subtitle = f"{colour_description}; red rings: no counted {gene} copy" if gene else colour_description
        if data.marked.any():
            subtitle += "; orange squares: selected genomes"
        fig.text(.08, .91, subtitle, fontsize=9, color=MUTED, va="top")
        footnote = ("Dashed: OLS regression; shading: pointwise 95% confidence interval for the mean response.\n"
                    "Pearson p is nominal, with no phylogenetic correction.")
        if gene:
            footnote += " Gene status follows the current annotation filters."
        fig.text(.08, .045, footnote,
                 fontsize=7.5, color=MUTED)
        figures[name] = _save(fig, outdir / name)
    return figures
