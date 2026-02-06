#!/usr/bin/env python3
"""
Multi-track diagnostic panel plotting.

A panel is a vertically stacked set of tracks (e.g. RM, depth, AGP)
sharing a common x-axis for a single genomic region.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker

from .universal import parse_region
from .rm_track import (
    load_data as load_rm,
    choose_taxonomy,
    bin_intervals,
    plot_binned,
    make_color_map,
)
from .depth_track import load_depth, subset_depth, plot_depth
from .agp_track import load_agp, subset_agp, plot_agp_layers


# ------------------------------------------------------------
# Panel orchestration
# ------------------------------------------------------------

def plot_panel(
    region,
    *,
    rm_path=None,
    depth_path=None,
    agp_path=None,
    rm_taxonomy="class",
    rm_bin_size=50_000,
    depth_bin_size=10_000,
    fig=None,
    gs=None,
):
    """
    Plot a multi-track diagnostic panel for a genomic region.

    Parameters
    ----------
    region : str
        CHROM or CHROM:start-end
    rm_path : Path or None
        RepeatMasker TSV
    depth_path : Path or None
        samtools depth TSV
    agp_path : Path or None
        AGP file
    fig : matplotlib Figure or None
        If provided, plot into this figure
    gs : matplotlib GridSpec or None
        GridSpec slot to draw into

    Returns
    -------
    axes : list[matplotlib Axes]
        Axes objects in top-to-bottom order
    """
    contig, start, end = parse_region(region)

    tracks = []
    heights = []

    if rm_path:
        tracks.append("rm")
        heights.append(3)

    if depth_path:
        tracks.append("depth")
        heights.append(2)

    if agp_path:
        tracks.append("agp")
        heights.append(1)

    if not tracks:
        raise ValueError("At least one of rm_path, depth_path, or agp_path must be provided")

    # Create figure / gridspec if needed
    if fig is None:
        fig = plt.figure(figsize=(12, sum(heights)))

    if gs is None:
        gs = fig.add_gridspec(
            nrows=len(tracks),
            ncols=1,
            height_ratios=heights,
            hspace=0.05,
        )

    axes = []
    ax_map = {}

    for i, track in enumerate(tracks):
        ax = fig.add_subplot(gs[i, 0], sharex=axes[0] if axes else None)
        axes.append(ax)
        ax_map[track] = ax

    # --------------------------------------------------------
    # RepeatMasker track
    # --------------------------------------------------------
    if rm_path:
        df = load_rm(Path(rm_path), contig)

        if start is not None:
            df = df[(df.end > start) & (df.start < end)].copy()

        taxonomy_col = choose_taxonomy(df, rm_taxonomy)
        categories = taxonomy_col.unique()
        color_map = make_color_map(categories)

        binned = bin_intervals(df, taxonomy_col, rm_bin_size)
        plot_binned(binned, ax_map["rm"], color_map)

        ax_map["rm"].set_ylabel("Repeats")

    # --------------------------------------------------------
    # Depth track
    # --------------------------------------------------------
    if depth_path:
        depth = load_depth(Path(depth_path))
        depth_sub = subset_depth(depth, contig, start, end)

        plot_depth(
            depth_sub,
            ax_map["depth"],
            bin_size=depth_bin_size,
            region_start=start,
        )

    # --------------------------------------------------------
    # AGP track
    # --------------------------------------------------------
    if agp_path:
        agp = load_agp(Path(agp_path))
        agp_sub = subset_agp(agp, contig, start, end)

        plot_agp_layers(
            agp_sub,
            ax_map["agp"],
            region_start=start,
            region_end=end,
            rebase=True,
        )

    # --------------------------------------------------------
    # Axis formatting
    # --------------------------------------------------------
    axes[-1].set_xlabel("Genomic position (Mb)")
    axes[-1].xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{x / 1e6:.1f}")
    )
    axes[-1].xaxis.set_major_locator(mticker.MultipleLocator(5e6))
    axes[-1].xaxis.set_minor_locator(mticker.MultipleLocator(1e6))

    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)

    return axes

def run_from_cli(args):
    plot_panel(
        region=args.region,
        rm_path=args.rm,
        depth_path=args.depth,
        agp_path=args.agp,
        rm_taxonomy=args.taxonomy,
        rm_bin_size=args.rm_bin,
        depth_bin_size=args.depth_bin,
    )

    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close()
