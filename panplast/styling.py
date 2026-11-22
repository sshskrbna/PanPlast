"""One taxon colour key shared by figures and phylogenies."""
from __future__ import annotations

import pandas as pd
from matplotlib.colors import hsv_to_rgb, to_hex, to_rgba


BRIGHT_PALETTE = (
    "#4363d8", "#c7003e", "#009dcc", "#f58231", "#800000",
    "#0277bd", "#00a040", "#1a237e", "#00796b", "#6a1fa8",
    "#bf5000", "#33691e", "#911eb4", "#5d4037", "#a08000",
    "#f032e6", "#5d8a00", "#e6194b", "#37474f", "#3cb44b",
)


def category_colors(data, label_column="group", color_column="group_color"):
    """Resolve optional user colours, rejecting conflicting category assignments."""
    colors = {}
    for i, label in enumerate(sorted(data[label_column].unique())):
        supplied = set()
        if color_column in data:
            for value in data.loc[data[label_column] == label, color_column].dropna():
                if str(value).strip():
                    try:
                        supplied.add(to_rgba(str(value).strip()))
                    except ValueError as exc:
                        raise ValueError(f"Invalid metadata color for {label!r}: {value!r}") from exc
        if len(supplied) > 1:
            raise ValueError(f"Conflicting metadata colors for {label!r}; use one color per category")
        if supplied:
            colors[label] = supplied.pop()
        elif str(label).lower() in {"unknown", "ungrouped", "unassigned", "—", "all genomes"}:
            colors[label] = to_rgba("#37474f")
        elif i < len(BRIGHT_PALETTE):
            colors[label] = to_rgba(BRIGHT_PALETTE[i])
        else:
            colors[label] = to_rgba(hsv_to_rgb(((i * .61803398875) % 1, .85, .75)))
    return colors


def taxon_styles(genomes, labels, groups=None, metadata=None, color_by="group", label_by="species"):
    """Build the full-dataset colour key once so subsets of a tree cannot shift it."""
    if color_by not in {"species", "group", "none"}:
        raise ValueError("color_by must be species, group or none")
    if label_by not in {"species", "accession"}:
        raise ValueError("label_by must be species or accession")
    groups, metadata = groups or {}, metadata or {}
    rows, canonical_species = [], {}
    for genome in genomes:
        species = str(labels.get(genome) or genome).strip() or genome
        # Underscores and spaces are interchangeable when identifying a species.
        species_key = " ".join(species.replace("_", " ").split()).casefold()
        canonical_species.setdefault(species_key, species)
        group = groups.get(genome) or "Ungrouped"
        category = canonical_species[species_key] if color_by == "species" else group
        if color_by == "none":
            category = "All genomes"
        rows.append(dict(genome=genome, species=species, group=group, color_by=color_by,
                         color_label=category, color=metadata.get(genome, {}).get("color", "") if color_by != "none" else "",
                         label_by=label_by, display_label=species if label_by == "species" else genome))
    frame = pd.DataFrame(rows, columns=["genome", "species", "group", "color_by", "color_label", "color",
                                        "label_by", "display_label"])
    palette = category_colors(frame, "color_label", "color")
    frame["color"] = frame.color_label.map({key: to_hex(value, keep_alpha=True) for key, value in palette.items()})
    return frame
