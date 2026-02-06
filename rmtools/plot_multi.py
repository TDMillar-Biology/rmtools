#!/usr/bin/env python3
"""
Plot RepeatMasker annotations for multiple assemblies/contigs
as stacked, left-aligned tracks.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker

from .rm_track import (
    load_data,
    choose_taxonomy,
    bin_intervals_dominant,
    plot_binned,
    make_color_map,
    add_legend,
)

from .universal import parse_region

# ----------------------------
# Control file handling
# ----------------------------

def load_control_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, delim_whitespace=True)
    print(df)
    required = {"path", "contig", "label"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Control file missing columns: {missing}")
    return df


# ----------------------------
# Main plotting routine
# ----------------------------

def plot_multi(
    control_df: pd.DataFrame,
    taxonomy: str,
    bin_size: int,
    figsize=(12, None),
):
    n_tracks = len(control_df)
    height = figsize[1] or 2.2 * n_tracks

    fig, axes = plt.subplots(
        nrows=n_tracks,
        figsize=(figsize[0], height),
        sharex=True,
    )

    if n_tracks == 1:
        axes = [axes]

    for ax, row in zip(axes, control_df.itertuples()):
        # contig may include region: contig[:start-end]
        contig, r_start, r_end = parse_region(row.contig)

        df = load_data(Path(row.path), contig)

        # optional region subsetting
        if r_start is not None:
            df = df[
                (df["end"] > r_start) &
                (df["start"] < r_end)
            ].copy()

        if df.empty:
            ax.text(
                0.5, 0.5,
                "No data",
                transform=ax.transAxes,
                ha="center",
                va="center"
            )
            ax.set_ylabel(row.label, rotation=0, ha="right", va="center")
            continue

        # --------------------------------------------------
        # Left-align: rebase to local coordinates
        # --------------------------------------------------
        left = df["start"].min()
        df["start"] -= left
        df["end"] -= left

        taxonomy_col = choose_taxonomy(df, taxonomy)
        categories = taxonomy_col.unique()
        color_map = make_color_map(categories)

        binned = bin_intervals_dominant(df, taxonomy_col, bin_size)
        plot_binned(binned, ax, color_map)

        

        # Track-style label
        ax.set_ylabel(row.label, rotation=0, ha="right", va="center")
        ax.yaxis.set_label_coords(-0.10, 0.5)

        # Label each x axis
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, pos: f"{x/1e6:.1f}")
        )
        ax.xaxis.set_major_locator(mticker.MultipleLocator(5e6))
        ax.xaxis.set_minor_locator(mticker.MultipleLocator(1e6))

    # ----------------------------
    # Shared x-axis formatting
    # ----------------------------
    ax = axes[-1]
    ax.set_xlabel("Relative position along contig (Mb)")


    # ----------------------------
    # Shared Legend
    # ----------------------------
    ax = axes[0]
    add_legend(ax, color_map, title=f"Repeat Class")

    return fig

def run_from_cli(args):
    control_df = load_control_file(Path(args.control))
    fig = plot_multi(
        control_df,
        taxonomy=args.taxonomy,
        bin_size=args.bin_size,
    )
    fig.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")

