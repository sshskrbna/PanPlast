"""IQ-TREE input and gene-content overlays; no ancestral loss inference.

Unmatched tree tips retain unknown states, rather than being treated as genomes
with no genes. Taxon matching is deliberately conservative and one-to-one.
"""
from __future__ import annotations

import csv
import io
import math
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import Phylo
from Bio.Phylo.NewickIO import NewickError

from .genes import EXPECTED_CDS_LEN, PLASTID_TRNAS, accession_of


@dataclass
class TreeOverlay:
    tree: object
    source: str
    taxa: pd.DataFrame
    copies: pd.DataFrame          # tips x genes; NaN means unmatched tip
    absences: pd.DataFrame


def read_tree(path):
    """Read one Newick tree, or the ML tree embedded in an IQ-TREE report."""
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig").strip()
    if path.suffix.lower() == ".iqtree":
        marker = re.search(r"Tree in newick format:\s*", text, re.I)
        if not marker:
            raise ValueError(f"{path}: no embedded Newick tree; supply the .treefile or .contree instead")
        # IQ-TREE prints the ML tree first; later sections may contain a consensus.
        text = text[marker.end():].split(";", 1)[0] + ";"
    if not text.endswith(";"):
        raise ValueError(f"{path}: Newick tree must end with ';'")
    try:
        tree = Phylo.read(io.StringIO(text), "newick")
    except (ValueError, IndexError, NewickError) as exc:
        raise ValueError(f"{path}: cannot read a single Newick tree: {exc}") from exc
    tips = tree.get_terminals()
    names = [tip.name for tip in tips]
    if len(tips) < 2 or any(not name for name in names):
        raise ValueError(f"{path}: a tree needs at least two named tips")
    if len(set(names)) != len(names):
        raise ValueError(f"{path}: duplicate tip names; each tip must have a unique name")
    for node in tree.find_clades():
        if node.branch_length is not None and (not math.isfinite(node.branch_length) or node.branch_length < 0):
            raise ValueError(f"{path}: branch lengths must be finite and non-negative")
    return tree


def load_tip_map(path):
    """Optional TSV: tip (exact tree label), genome (exact PanPlast genome ID)."""
    if path is None:
        return {}
    result = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not {"tip", "genome"} <= set(reader.fieldnames or []):
            raise ValueError("--tree-map needs TSV columns 'tip' and 'genome'")
        for row in reader:
            tip, genome = ((row.get(k) or "").strip() for k in ("tip", "genome"))
            if not tip or not genome or tip in result:
                raise ValueError("--tree-map contains an empty value or a duplicate tip")
            result[tip] = genome
    return result


def _species_key(value):
    return " ".join(value.replace("_", " ").split()).casefold()


def match_taxa(tree, matrix, labels, groups=None, tip_map=None):
    """Match IDs, embedded accessions, then unambiguous species labels."""
    groups, tip_map = groups or {}, tip_map or {}
    genomes = list(matrix.columns)
    unknown = set(tip_map.values()) - set(genomes)
    if unknown:
        raise ValueError("--tree-map references unknown genome IDs: " + ", ".join(sorted(unknown)))
    aliases = {}
    for genome in genomes:
        for alias in {genome, re.sub(r"\.\d+$", "", genome)}:
            aliases.setdefault(alias, set()).add(genome)
    species = {}
    for genome in genomes:
        label = labels.get(genome, "")
        if label.strip():
            species.setdefault(_species_key(label), set()).add(genome)
    rows, used = [], set()
    for tip in tree.get_terminals():
        name = tip.name
        if name in tip_map:
            candidates, method = {tip_map[name]}, "explicit_map"
        elif name in matrix.columns:
            candidates, method = {name}, "genome_id"
        else:
            acc = accession_of(name)
            candidates = aliases.get(acc or name, set())
            method = "accession"
            if not candidates:
                candidates = species.get(_species_key(name), set())
                method = "species"
        if len(candidates) > 1:
            raise ValueError(f"Ambiguous tree tip {name!r}; use --tree-map with an exact genome ID")
        genome = next(iter(candidates), None)
        if genome in used:
            raise ValueError(f"Multiple tree tips map to genome {genome!r}; use one tip per genome")
        if genome is not None:
            used.add(genome)
        rows.append(dict(tip=name, genome=genome or "", species=labels.get(genome, ""),
                         group=groups.get(genome, ""), status="matched" if genome is not None else "unmatched_tip",
                         matched_by=method if genome is not None else ""))
    if not used:
        raise ValueError("No tree tips match the input genomes. Use --tree-map (tip, genome) to link their IDs")
    for genome in genomes:
        if genome not in used:
            rows.append(dict(tip="", genome=genome, species=labels.get(genome, ""), group=groups.get(genome, ""),
                             status="genome_not_in_tree", matched_by=""))
    return pd.DataFrame(rows)


def build_overlay(path, matrix, labels, groups=None, tip_map=None, genes=None, loci_by=None):
    tree = read_tree(path)
    taxa = match_taxa(tree, matrix, labels, groups, tip_map)
    matched = taxa[taxa.status == "matched"]
    if genes is None:
        # Select observed pan-genes missing in at least one matched genome.
        genes = list(matrix.index[(matrix[matched.genome] == 0).any(axis=1)])
    else:
        genes = list(dict.fromkeys(genes))
        known = set(matrix.index) | set(EXPECTED_CDS_LEN) | PLASTID_TRNAS
        known |= {loc.name for loci in (loci_by or {}).values() for loc in loci}
        unknown = set(genes) - known
        if unknown:
            raise ValueError("Unknown --tree-genes names: " + ", ".join(sorted(unknown)))
    counts = matrix.reindex(genes, fill_value=0)
    copies = pd.DataFrame(float("nan"), index=[tip.name for tip in tree.get_terminals()], columns=genes)
    absences = []
    for row in matched.itertuples():
        copies.loc[row.tip] = counts[row.genome].to_numpy()
        for gene in genes:
            if counts.at[gene, row.genome] != 0:
                continue
            loci = [loc for loc in (loci_by or {}).get(row.genome, []) if loc.name == gene]
            flags = sorted({flag for loc in loci for flag in loc.flags})
            absences.append(dict(tip=row.tip, genome=row.genome, species=row.species, gene=gene,
                                 status="not_counted" if loci else "not_annotated",
                                 flags=",".join(flags)))
    return TreeOverlay(tree, str(path), taxa, copies,
                       pd.DataFrame(absences, columns=["tip", "genome", "species", "gene", "status", "flags"]))


def plot_overlay(overlay, summary, path_base, genes_per_page=30):
    """Draw aligned phylogram, labels, genome lengths and gene states in pages."""
    import numpy as np
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    from .plots import INK, MUTED, _save, plt

    if genes_per_page < 1:
        raise ValueError("genes_per_page must be positive")
    tree = overlay.tree
    tips = tree.get_terminals()
    by_tip = overlay.taxa[overlay.taxa.tip != ""].set_index("tip")
    n = len(tips)
    y = {tip: i for i, tip in enumerate(tips)}
    for node in tree.get_nonterminals(order="postorder"):
        y[node] = (y[node.clades[0]] + y[node.clades[-1]]) / 2
    branches = [node for node in tree.find_clades() if node is not tree.root]
    use_lengths = all(node.branch_length is not None for node in branches) and any(
        node.branch_length > 0 for node in branches)
    x = {tree.root: 0.0}
    for node in tree.find_clades(order="preorder"):
        for child in node.clades:
            x[child] = x[node] + (child.branch_length if use_lengths else 1)
    xmax = max(x.values()) or 1
    groups = sorted(set(overlay.taxa.loc[overlay.taxa.status == "matched", "group"]) - {""})
    palette = plt.get_cmap("tab20")
    group_colors = {group: palette(i % 20) for i, group in enumerate(groups)}
    stats = summary.set_index("genome")
    tip_labels = []
    for tip in tips:
        row = by_tip.loc[tip.name]
        label = f"{row.species or row.genome} ({row.genome})" if row.genome else f"{tip.name} [no data]"
        if row.species == row.genome or not row.species:
            label = row.genome or f"{tip.name} [no data]"
        tip_labels.append(label)
    label_width = max(2.6, max(map(len, tip_labels)) * 0.052 + 0.6)
    genes = list(overlay.copies.columns)
    pages = [genes[i:i + genes_per_page] for i in range(0, len(genes), genes_per_page)] or [[]]
    results = []
    legend_rows = math.ceil(len(groups) / 5)
    for page_idx, page in enumerate(pages, 1):
        widths = [3.0, label_width, 1.05, 1.05, max(2.0, 0.29 * len(page))]
        height = max(4.8, n * 0.25 + 2.1 + legend_rows * 0.22)
        fig = plt.figure(figsize=(sum(widths) + 0.8, height))
        bottom = (1.05 + legend_rows * 0.22) / height
        grid = fig.add_gridspec(1, 5, width_ratios=widths, left=0.035, right=0.985,
                               top=1 - 0.85 / height, bottom=bottom, wspace=0.035)
        axes = [fig.add_subplot(grid[0, i]) for i in range(5)]
        for ax in axes:
            ax.set_ylim(n - 0.5, -0.5)
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.tick_params(left=False, labelleft=False)
        ax_tree, ax_labels, ax_ir, ax_size, ax_genes = axes
        for node in tree.find_clades():
            if node.clades:
                ax_tree.plot([x[node], x[node]], [y[node.clades[0]], y[node.clades[-1]]], color=INK, lw=0.65)
            for child in node.clades:
                ax_tree.plot([x[node], x[child]], [y[child], y[child]], color=INK, lw=0.65)
            support = node.name or (f"{node.confidence:g}" if node.confidence is not None else "")
            if node is not tree.root and node.clades and support:
                ax_tree.annotate(support, (x[node], y[node]), xytext=(-2, 3), textcoords="offset points",
                                 fontsize=5, ha="right", color=MUTED)
        for i, tip in enumerate(tips):
            row = by_tip.loc[tip.name]
            color = group_colors.get(row.group, INK) if row.genome else "#aeb4bc"
            # Pastel group colours are suitable for markers but too faint for text.
            label_color = tuple(c * 0.65 for c in color[:3]) if row.genome and row.group in group_colors else INK
            if not row.genome:
                label_color = "#737b86"
            ax_tree.plot([x[tip], xmax * 1.025], [i, i], color="#d8dce2", lw=0.5, ls=":")
            ax_tree.scatter([x[tip]], [i], s=9, color=color, zorder=3)
            ax_labels.text(0, i, tip_labels[i], color=label_color, va="center", fontsize=7, style="italic")
        ax_tree.set_xlim(-0.02 * xmax, xmax * 1.04)
        ax_tree.set_xlabel("Substitutions / site" if use_lengths else "Topology (branch lengths unavailable)", fontsize=7)
        ax_tree.tick_params(axis="x", labelsize=6)
        ax_tree.locator_params(axis="x", nbins=4)
        ax_labels.set_xlim(0, 1)
        ax_labels.axis("off")
        for ax, column, title, color in ((ax_ir, "IRb_bp", "IRb (kb)", "#d08c3c"),
                                         (ax_size, "length_bp", "Genome (kb)", "#5a9e6f")):
            values = []
            for i, tip in enumerate(tips):
                genome = by_tip.loc[tip.name, "genome"]
                value = stats.at[genome, column] / 1000 if genome else float("nan")
                values.append(value)
                if pd.isna(value):
                    ax.text(0, i, "—", va="center", fontsize=7, color=MUTED)
                else:
                    ax.barh(i, value, height=0.6, color=color)
            finite = [v for v in values if pd.notna(v)]
            ax.set_xlim(0, max(finite, default=1) * 1.08 or 1)
            ax.set_title(title, fontsize=8, pad=12)
            ax.tick_params(axis="x", labelsize=6)
            ax.locator_params(axis="x", nbins=3)
        if page:
            values = overlay.copies[page].to_numpy().clip(max=2)
            # -1: no matching genome; 0: no counted copy; 1 or 2: counted copies.
            values = np.nan_to_num(values, nan=-1)
            cmap = ListedColormap(["#e2e5ea", "#c4514d", "#8fbcd4", "#2f6f95"])
            ax_genes.imshow(values, aspect="auto", interpolation="nearest", cmap=cmap,
                            norm=BoundaryNorm([-1.5, -.5, .5, 1.5, 2.5], 4))
            ax_genes.set_xticks(range(len(page)), page, rotation=90, fontsize=7, style="italic")
            ax_genes.set_xticks(np.arange(-.5, len(page)), minor=True)
            ax_genes.set_yticks(np.arange(-.5, n), minor=True)
            ax_genes.grid(which="minor", color="white", lw=0.5)
            ax_genes.tick_params(which="both", length=0)
        else:
            ax_genes.text(.5, .5, "No missing genes\nin matched genomes", transform=ax_genes.transAxes,
                          ha="center", va="center", fontsize=8, color=MUTED)
            ax_genes.axis("off")
        ax_genes.set_title("Gene copies" + (f" · {page_idx}/{len(pages)}" if len(pages) > 1 else ""),
                           fontsize=8, pad=12)
        fig.text(0.035, 1 - 0.22 / height, f"Gene absence on IQ-TREE · {Path(overlay.source).name}",
                 fontsize=11, fontweight="bold", va="top", color=INK)
        fig.text(0.035, 1 - 0.5 / height,
                 "Observed gene content at tips; red = no counted copy. Branch labels retain input support values.",
                 fontsize=7, va="top", color=MUTED)
        legend = [Patch(color=c, label=label) for c, label in
                  (("#c4514d", "Absent / not counted"), ("#8fbcd4", "1 copy"),
                   ("#2f6f95", "≥2 copies"), ("#e2e5ea", "No matching genome"))]
        fig.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.03, (0.34 + legend_rows * .22) / height),
                   ncol=4, frameon=False, fontsize=7)
        if groups:
            fig.legend(handles=[Patch(color=group_colors[g], label=g) for g in groups],
                       loc="lower left", bbox_to_anchor=(0.03, 0.3 / height), ncol=5, frameon=False, fontsize=6)
        fig.text(0.035, 0.13 / height,
                 "Absence can reflect incomplete annotation or assembly; ancestral loss events are not inferred.",
                 fontsize=6.5, color=MUTED)
        base = f"{path_base}_{page_idx}" if len(pages) > 1 else path_base
        results.append(_save(fig, base))
    return results
