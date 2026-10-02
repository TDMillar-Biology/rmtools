#!/usr/bin/env python3
"""
Plot RepeatMasker annotations along a contig.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.patches import Patch
from matplotlib import ticker as mticker
from .universal import parse_region, load_sizes


def load_data(path: Path, contig = None):
    df = pd.read_csv(path, sep="\t", keep_default_na=False, dtype={"chrom": str})
    if contig:
        return df[df["chrom"] == contig].sort_values("start")
    else:
        return df.sort_values("start")


def choose_taxonomy(df, level):
    if level == "class":
        return df["repeat_class"]
    elif level == "family":
        return df["repeat_class"] + "/" + df["repeat_family"]
    elif level == "name":
        return df["repeat_name"]
    else:
        raise ValueError(level)


def plot_raw_intervals(df, taxonomy_col, ax, color_map):
    for i, row in df.iterrows():
        ax.broken_barh(
            [(row.start, row.end - row.start)],
            (0, 1),
            facecolors=color_map[taxonomy_col.loc[i]]
        )

    ax.set_ylim(0, 1)
    ax.set_yticks([])



def plot_binned(df_bins, ax, color_map):
    '''
    Taxon below refers to a named repeat class, family, unit in the taxonomy
    Not sure if that's the correct way to name that entity
    '''
    taxa = df_bins["taxonomy"].unique()
    taxon_order = [t for t in taxa if t != "Unannotated"] + ["Unannotated"]

    bottoms = defaultdict(int)

    for taxon in taxon_order:
        sub = df_bins[df_bins["taxonomy"] == taxon]
        ax.bar(
            sub["bin_start"],
            sub["coverage"],
            width=sub["bin_end"] - sub["bin_start"],
            bottom=[bottoms[b] for b in sub["bin_start"]],
            color=color_map[taxon],
            align="edge"
        )
        for b, h in zip(sub["bin_start"], sub["coverage"]):
            bottoms[b] += h


def add_legend(ax, color_map, title="Repeat class"):
    """
    Add a legend based on taxonomy → color mapping.
    """
    handles = [
        Patch(facecolor=color, label=label)
        for label, color in color_map.items()
    ]

    ax.legend(
        handles=handles,
        title=title,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False
    )

def make_color_map(categories, cmap=plt.cm.tab20):
    """
    Assign a consistent color to each taxonomy category.
    """
    categories = sorted(set(categories) - {"Unannotated"})
    color_map = {cat: cmap(i % cmap.N) for i, cat in enumerate(categories)}
    color_map["Unannotated"] = '#E6E6E6' ## just off white
    return color_map


def merge_intervals(intervals):
    """
    Merge overlapping intervals.
    intervals: list of (start, end)
    Returns list of merged (start, end)
    """
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [intervals[0]]

    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged


def clip_intervals_to_bin(df, taxonomy_col, bin_start, bin_end):
    """
    Return list of (start, end, taxonomy) clipped to bin boundaries.
    Handles the case of an annotation whose boundaries cannot be contained in a single bin
    """
    clipped = []

    for _, row in df.iterrows():
        s = max(row.start, bin_start)
        e = min(row.end, bin_end)
        if s < e:
            clipped.append((s, e, taxonomy_col.loc[_]))

    return clipped


def _bin_intervals(df, taxonomy_col, bin_size, mode, start=None, end=None):
    """Use genome-anchored bins clipped to a half-open plotting extent."""
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")
    start = 0 if start is None else start
    if end is None:
        end = int(df["end"].max()) if not df.empty else start
    if start < 0 or end < start:
        raise ValueError("Invalid binning extent")
    records = []
    for anchor in range((start // bin_size) * bin_size, end, bin_size):
        bin_start, bin_end = max(anchor, start), min(anchor + bin_size, end)
        width = bin_end - bin_start
        window = df[(df.start < bin_end) & (df.end > bin_start)]
        clipped = clip_intervals_to_bin(window, taxonomy_col, bin_start, bin_end)
        repeat_bp = sum(e - s for s, e in merge_intervals([(s, e) for s, e, _ in clipped]))
        class_bp = {}
        for taxon in sorted({c for _, _, c in clipped}):
            intervals = [(s, e) for s, e, c in clipped if c == taxon]
            class_bp[taxon] = (sum(e - s for s, e in intervals) if mode == "sum"
                              else sum(e - s for s, e in merge_intervals(intervals)))
        if mode == "dominant" and class_bp:
            class_bp = {max(class_bp, key=class_bp.get): repeat_bp}
        elif mode == "composition" and class_bp:
            total = sum(class_bp.values())
            class_bp = {c: bp / total * repeat_bp for c, bp in class_bp.items()}
        unannotated = max(width - (sum(class_bp.values()) if mode == "sum" else repeat_bp), 0)
        records.append(dict(bin_start=bin_start, bin_end=bin_end,
                            taxonomy="Unannotated", coverage=unannotated))
        records.extend(dict(bin_start=bin_start, bin_end=bin_end, taxonomy=c, coverage=bp)
                       for c, bp in class_bp.items())
    return pd.DataFrame(records, columns=["bin_start", "bin_end", "taxonomy", "coverage"])


def bin_intervals(df, taxonomy_col, bin_size, start=None, end=None):
    """Sum annotation lengths; overlapping annotations may exceed bin width."""
    return _bin_intervals(df, taxonomy_col, bin_size, "sum", start, end)


def bin_intervals_dominant(df, taxonomy_col, bin_size, start=None, end=None):
    """Assign the union of repeat-covered bases to the dominant class."""
    return _bin_intervals(df, taxonomy_col, bin_size, "dominant", start, end)


def bin_intervals_repeat_composition(df, taxonomy_col, bin_size, start=None, end=None):
    """Scale per-class union lengths proportionally to total repeat union length.

    Contributions are a proportional summary, not exclusive per-base assignments.
    """
    return _bin_intervals(df, taxonomy_col, bin_size, "composition", start, end)


def run_from_cli(args):
    contig, r_start, r_end = parse_region(args.region)

    df = load_data(Path(args.rm), contig)

    if r_start is not None: # subset to coords specified by user
        df = df[
            (df["end"] > r_start) &
            (df["start"] < r_end)
        ].copy()

        # Rebase to local coordinates (if you want relative coordinate system)
        #df["start"] -= r_start
        #df["end"] -= r_start
    sizes = load_sizes(getattr(args, "sizes", None))
    plot_start = 0 if r_start is None else r_start
    plot_end = r_end if r_end is not None else sizes.get(contig)
    if plot_end is None and not df.empty:
        plot_end = int(df.end.max())
    taxonomy_col = choose_taxonomy(df, args.taxonomy)

    categories = taxonomy_col.unique()
    color_map = make_color_map(categories)

    fig, ax = plt.subplots(figsize=(12, 2))

    if args.bin_size is None:
        plot_raw_intervals(df, taxonomy_col, ax, color_map)
    else:
        binned = bin_intervals_dominant(df, taxonomy_col, args.bin_size, start=plot_start, end=plot_end)
        plot_binned(binned, ax, color_map)

    if plot_end is not None and plot_end > plot_start:
        ax.set_xlim(plot_start, plot_end)
    if df.empty:
        ax.text(0.5, 0.5, "No repeat annotations", transform=ax.transAxes, ha="center")
    add_legend(ax, color_map, title=f"Repeat {args.taxonomy}")

    ## Format axis labels
    # x axis
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, pos: f"{x/1e6:.1f}"))
    ax.xaxis.set_major_locator(mticker.MultipleLocator(5e6))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(1e6))
    ax.set_xlabel("Genomic position (Mb)")

    # y axis
    ax.set_ylabel("Repeat coverage (bp)")

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
