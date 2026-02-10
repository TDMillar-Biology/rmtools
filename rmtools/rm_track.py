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
from .universal import parse_region
import pdb

def load_data(path: Path, contig: str):
    df = pd.read_csv(path, sep="\t")
    print(f'Loaded data from {path}. Dataframe shape: {df.shape}')
    return df[df["chrom"] == contig].sort_values("start")


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



def bin_intervals(df, taxonomy_col, bin_size):
    max_pos = df["end"].max()
    bins = range(0, max_pos + bin_size, bin_size)
    records = []

    for bin_start in bins:
        bin_end = bin_start + bin_size
        window = df[(df.start < bin_end) & (df.end > bin_start)]

        total_covered = 0

        for taxon, sub in window.groupby(taxonomy_col):
            covered = (
                sub[["start", "end"]]
                .apply(
                    lambda x: min(x.end, bin_end) - max(x.start, bin_start),
                    axis=1
                )
                .sum()
            )

            total_covered += covered

            records.append({
                "bin_start": bin_start,
                "bin_end": bin_end,
                "taxonomy": taxon,
                "coverage": covered
            })

        # add unannotated explicitly
        records.append({
            "bin_start": bin_start,
            "bin_end": bin_end,
            "taxonomy": "Unannotated",
            "coverage": max(bin_size - total_covered, 0)
        })

    return pd.DataFrame(records)



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
    categories = sorted(categories)
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


def bin_intervals_dominant(df, taxonomy_col, bin_size):
    """
    Bin intervals by dominant repeat class with exclusive base accounting.
    FREQUENTLY VIOLATED ASSUMPTION -- rm output doesnt overlap -- use with caution
    For each bin:
      - compute union of all repeat intervals
      - define unannotated = bin_size - union_size
      - assign repeat portion to dominant repeat class
    """
    max_pos = df["end"].max()
    bins = range(0, max_pos + bin_size, bin_size)
    records = []

    for bin_start in bins:
        bin_end = bin_start + bin_size

        window = df[(df.start < bin_end) & (df.end > bin_start)]
        if window.empty:
            records.append({
                "bin_start": bin_start,
                "bin_end": bin_end,
                "taxonomy": "Unannotated",
                "coverage": bin_size
            })
            continue

        # Clip intervals to bin
        clipped = clip_intervals_to_bin(window, taxonomy_col, bin_start, bin_end)

        # ---- union of all repeat intervals to calculate unannotated bp ----
        all_intervals = [(s, e) for s, e, _ in clipped]
        union = merge_intervals(all_intervals) ##
        
        repeat_bp = sum(e - s for s, e in union)
        unannotated_bp = max(bin_size - repeat_bp, 0)


        # ---- compute unique bp per class (for dominance only) ----
        class_bp = defaultdict(int)
        for cls in set(c for _, _, c in clipped):
            cls_intervals = [(s, e) for s, e, c in clipped if c == cls]
            cls_union = merge_intervals(cls_intervals)
            class_bp[cls] = sum(e - s for s, e in cls_union)

        # ---- emit records ----
        records.append({
            "bin_start": bin_start,
            "bin_end": bin_end,
            "taxonomy": "Unannotated",
            "coverage": unannotated_bp
        })

        if repeat_bp > 0 and class_bp:
            dominant = max(class_bp, key=class_bp.get)
            records.append({
                "bin_start": bin_start,
                "bin_end": bin_end,
                "taxonomy": dominant,
                "coverage": repeat_bp
            })
        #pdb.set_trace()

    return pd.DataFrame(records)

def bin_intervals_repeat_composition(df, taxonomy_col, bin_size):
    """
    Bin intervals and represent repeat composition proportionally within
    the repeat-covered portion of each bin.

    For each bin:
      1. Compute union of all repeat intervals → repeat_bp
      2. Define unannotated = bin_size - repeat_bp
      3. Compute per-class annotated bp (overlaps allowed)
      4. Project class annotation space onto repeat space proportionally

    NOTE:
    - RepeatMasker annotations may overlap.
    - Per-class bp are computed independently.
    - Class contributions are normalized within annotation space and
      scaled to repeat space.
    """
    max_pos = df["end"].max()
    bins = range(0, max_pos + bin_size, bin_size)
    records = []

    for bin_start in bins:
        bin_end = bin_start + bin_size

        window = df[(df.start < bin_end) & (df.end > bin_start)]
        if window.empty:
            records.append({
                "bin_start": bin_start,
                "bin_end": bin_end,
                "taxonomy": "Unannotated",
                "coverage": bin_size
            })
            continue

        # ---- clip intervals to bin ----
        clipped = clip_intervals_to_bin(
            window,
            taxonomy_col,
            bin_start,
            bin_end
        )
        # clipped: List[(start, end, class)]

        # ---- union of all repeat intervals ----
        all_intervals = [(s, e) for s, e, _ in clipped]
        union = merge_intervals(all_intervals)

        repeat_bp = sum(e - s for s, e in union)
        unannotated_bp = max(bin_size - repeat_bp, 0)

        # ---- compute per-class annotated bp (annotation space) ----
        class_bp = defaultdict(int)
        for cls in set(c for _, _, c in clipped):
            cls_intervals = [(s, e) for s, e, c in clipped if c == cls]
            cls_union = merge_intervals(cls_intervals)
            class_bp[cls] = sum(e - s for s, e in cls_union)

        annotation_bp = sum(class_bp.values())

        # ---- emit unannotated ----
        records.append({
            "bin_start": bin_start,
            "bin_end": bin_end,
            "taxonomy": "Unannotated",
            "coverage": unannotated_bp
        })

        # ---- project annotation space → repeat space ----
        if repeat_bp > 0 and annotation_bp > 0:
            for cls, cls_bp in class_bp.items():
                scaled_bp = (cls_bp / annotation_bp) * repeat_bp
                if scaled_bp > 0:
                    records.append({
                        "bin_start": bin_start,
                        "bin_end": bin_end,
                        "taxonomy": cls,
                        "coverage": scaled_bp
                    })

    return pd.DataFrame(records)



def run_from_cli(args):
    contig, r_start, r_end = parse_region(args.region)

    df = load_data(Path(args.input), contig)
    print(f'Contig Specific Dataframe shape: {df.shape}')

    if r_start is not None: # subset to coords specified by user
        df = df[
            (df["end"] > r_start) &
            (df["start"] < r_end)
        ].copy()

        # Rebase to local coordinates (if you want relative coordinate system)
        #df["start"] -= r_start
        #df["end"] -= r_start
    taxonomy_col = choose_taxonomy(df, args.taxonomy)

    categories = taxonomy_col.unique()
    color_map = make_color_map(categories)

    fig, ax = plt.subplots(figsize=(12, 2))

    if args.bin_size is None:
        plot_raw_intervals(df, taxonomy_col, ax, color_map)
    else:
        binned = bin_intervals_dominant(df, taxonomy_col, args.bin_size)
        plot_binned(binned, ax, color_map)

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