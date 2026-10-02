#!/usr/bin/env python3
"""
Multi-track diagnostic panel plotting.

A panel is a vertically stacked set of tracks (e.g. RM, depth, AGP)
sharing a common x-axis for a single genomic region.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker

from .universal import parse_region, load_sizes
from .rm_track import (
    load_data as load_rm,
    choose_taxonomy,
    bin_intervals_repeat_composition,
    plot_binned,
    make_color_map,
    add_legend,
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
    rm_color_map=None,
    contig_length=None,
):
    """
    Plot a multi-track diagnostic panel for a genomic region.

    Parameters
    ----------
    region : str
        CHROM or CHROM:start-end. Tracks retain absolute genomic positions;
        an explicit region fixes the shared x-axis to start/end.
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

    rm_df = load_rm(Path(rm_path), contig) if rm_path else None
    depth_df = load_depth(Path(depth_path)) if depth_path else None
    agp_df = load_agp(Path(agp_path)) if agp_path else None
    extent_end = end if end is not None else contig_length
    if extent_end is None:
        extents = []
        if rm_df is not None and not rm_df.empty:
            extents.append(int(rm_df.end.max()))
        if depth_df is not None:
            sub = subset_depth(depth_df, contig)
            if not sub.empty:
                extents.append(int(sub.pos.max()) + 1)
        if agp_df is not None:
            sub = subset_agp(agp_df, contig)
            if not sub.empty:
                extents.append(int(sub.obj_end.max()))
        extent_end = max(extents) if extents else None

    # --------------------------------------------------------
    # RepeatMasker track
    # --------------------------------------------------------
    if rm_path:
        df = rm_df

        if start is not None:
            df = df[(df.end > start) & (df.start < end)].copy()

        taxonomy_col = choose_taxonomy(df, rm_taxonomy)
        categories = taxonomy_col.unique()
        color_map = rm_color_map if rm_color_map is not None else make_color_map(categories)

        binned = bin_intervals_repeat_composition(df, taxonomy_col, rm_bin_size, start=start, end=extent_end)
        plot_binned(binned, ax_map["rm"], color_map)
        if df.empty:
            ax_map["rm"].text(0.5, 0.5, "No repeat annotations",
                              transform=ax_map["rm"].transAxes, ha="center")

        ax_map["rm"].set_ylabel("Repeats")

        add_legend(ax_map["rm"], color_map, title=f"Repeat {rm_taxonomy}")

    # --------------------------------------------------------
    # Depth track
    # --------------------------------------------------------
    if depth_path:
        depth = depth_df
        depth_sub = subset_depth(depth, contig, start, end)

        plot_depth(
            depth_sub,
            ax_map["depth"],
            bin_size=depth_bin_size,
            region_start=start,
            rebase=False,
        )

    # --------------------------------------------------------
    # AGP track
    # --------------------------------------------------------
    if agp_path:
        agp = agp_df
        agp_sub = subset_agp(agp, contig, start, end)

        plot_agp_layers(
            agp_sub,
            ax_map["agp"],
            region_start=start,
            region_end=end,
            rebase=False,
        )

    # --------------------------------------------------------
    # Axis formatting
    # --------------------------------------------------------
    if extent_end is not None:
        axes[0].set_xlim(0 if start is None else start, extent_end)

    axes[0].set_title(
        contig,
        loc="left",
        fontsize=12,
        fontweight="bold"
        )
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
        contig_length=load_sizes(getattr(args, "sizes", None)).get(parse_region(args.region)[0]),
    )


    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close()
