#!/usr/bin/env python3
"""
AGP golden-path plotting utilities.

Plots each W component on its own horizontal layer so that
scaffold structure and breakpoints are visually explicit.
"""

from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt
from .universal import parse_region

# ------------------------------------------------------------
# AGP loading
# ------------------------------------------------------------

def load_agp(path: Path):
    """
    Load an AGP file into a DataFrame.

    Assumes standard AGP 2.0 columns.
    """
    cols = [
        "object",
        "obj_beg",
        "obj_end",
        "part_number",
        "component_type",
        "component_id",
        "comp_beg",
        "comp_end",
        "orientation",
    ]

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=cols,
        dtype={
            "object": str,
            "obj_beg": int,
            "obj_end": int,
            "component_type": str,
            "component_id": str,
        }
    )

    return df


def subset_agp(df, contig, start=None, end=None):
    """
    Subset AGP dataframe to a contig and optional region (overlap-based).
    """
    sub = df[df["object"] == contig].copy()

    if start is not None:
        sub = sub[(sub.obj_end > start) & (sub.obj_beg < end)].copy()

    return sub



# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

def plot_agp_layers(
    agp_df,
    ax,
    region_start=None,
    region_end=None,
    bar_height=0.8,
    color="black",
    rebase=True,
):
    """
    Plot AGP W-components as layered bars, clipped to a region.

    Each W component gets its own horizontal lane.
    Only the portion overlapping the region is plotted.
    """
    w_df = agp_df[agp_df["component_type"] == "W"].copy()

    if w_df.empty:
        ax.text(
            0.5, 0.5,
            "No AGP W components",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )
        ax.set_yticks([])
        return

    # Preserve AGP order
    w_df = w_df.reset_index(drop=True)
    w_df["layer"] = range(len(w_df))

    for _, row in w_df.iterrows():
        # Clip to region if provided
        seg_start = row.obj_beg
        seg_end = row.obj_end

        if region_start is not None:
            seg_start = max(seg_start, region_start)
            seg_end = min(seg_end, region_end)

        if seg_start >= seg_end:
            continue  # nothing to plot

        # Optional rebasing so region_start -> 0
        if rebase and region_start is not None:
            seg_start -= region_start

        width = seg_end - seg_start
        y = row.layer

        ax.broken_barh(
            [(seg_start, width)],
            (y, bar_height),
            facecolors=color,
            edgecolors="none",
        )

    ax.set_ylim(0, len(w_df))
    ax.set_yticks([])
    ax.set_ylabel("AGP\ncomponents")


def run_from_cli(args):
    agp_path = Path(args.agp)

    contig, start, end = parse_region(args.region)

    agp_df = load_agp(agp_path)

    agp_sub = subset_agp(
        agp_df,
        contig=contig,
        start=start,
        end=end,
    )

    fig, ax = plt.subplots(figsize=(12, 2))

    plot_agp_layers(
        agp_sub,
        ax,
        region_start=start,
        region_end=end,
        rebase=True,
    )

    ax.set_xlabel("Genomic position (bp)")

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)


