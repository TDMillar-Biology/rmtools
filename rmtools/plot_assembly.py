#!/usr/bin/env python3

from pathlib import Path
import matplotlib.pyplot as plt

from .plot_panel import plot_panel
from .rm_track import load_data, choose_taxonomy, make_color_map
from .universal import load_sizes


def run_from_cli(args):

    contigs = args.main

    # determine number of tracks
    tracks = []
    if args.rm:
        tracks.append("rm")
    if args.depth:
        tracks.append("depth")
    if args.agp:
        tracks.append("agp")

    tracks_per_contig = len(tracks)

    if not tracks:
        raise ValueError("At least one of --rm, --depth, or --agp must be provided")
    sizes = load_sizes(getattr(args, "sizes", None))
    color_map = None
    if args.rm:
        df = load_data(Path(args.rm))
        df = df[df.chrom.isin(contigs)]
        color_map = make_color_map(choose_taxonomy(df, args.taxonomy).unique())

    fig = plt.figure(figsize=(12, len(contigs) * 4))

    fig.suptitle(
        f"Assembly Overview {args.out}",
        fontsize=16,
        fontweight="bold",
        y=0.995
    )

    gs_outer = fig.add_gridspec(
        nrows=len(contigs),
        ncols=1,
        hspace=0.2
    )

    for i, contig in enumerate(contigs):

        # this row becomes a panel
        panel_spec = gs_outer[i]

        panel_gs = panel_spec.subgridspec(
            tracks_per_contig, 1,
            height_ratios=[{"rm": 3, "depth": 2, "agp": 1}[t] for t in tracks]
        )

        plot_panel(
            region=contig,
            rm_path=args.rm,
            depth_path=args.depth,
            agp_path=args.agp,
            rm_taxonomy=args.taxonomy,
            rm_bin_size=args.rm_bin,
            depth_bin_size=args.depth_bin,
            fig=fig,
            gs=panel_gs,
            rm_color_map=color_map,
            contig_length=sizes.get(contig)
        )

    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close()