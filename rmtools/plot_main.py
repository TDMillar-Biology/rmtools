#!/usr/bin/env python3
"""
Plot RepeatMasker annotations for all contigs in an asm
as stacked, left-aligned tracks.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker

from .rm_track import (
    load_data,
    choose_taxonomy,
    bin_intervals_repeat_composition,
    plot_binned,
    make_color_map,
    add_legend,
    plot_raw_intervals
)

from .universal import load_sizes

def run_from_cli(args):
    df = load_data(Path(args.rm))
    taxonomy_col = choose_taxonomy(df, args.taxonomy)

    categories = taxonomy_col.unique()
    color_map = make_color_map(categories)

    number_plots = len(args.main)

    fig, axes = plt.subplots(number_plots, 1, sharex=True)

    if number_plots == 1:
        axes = [axes]

    sizes = load_sizes(getattr(args, "sizes", None))
    selected = df[df["chrom"].isin(args.main)]
    max_pos = max([sizes.get(c, 0) for c in args.main] +
                  ([int(selected.end.max())] if not selected.empty else [0]))

    for i, contig in enumerate(args.main):
        sub = df[df["chrom"] == contig]
        ax = axes[i]

        if args.bin_size is None:
            plot_raw_intervals(sub, taxonomy_col, ax, color_map)
        else:
            binned = bin_intervals_repeat_composition(
                sub, taxonomy_col, args.bin_size, end=sizes.get(contig)
            )
            plot_binned(binned, ax, color_map)

        # Format axis
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, pos: f"{x/1e6:.1f}")
        )
        ax.xaxis.set_major_locator(mticker.MultipleLocator(5e6))
        ax.xaxis.set_minor_locator(mticker.MultipleLocator(1e6))

        ax.set_title(contig, loc="left", fontsize=10, fontweight="bold")

        if max_pos:
            ax.set_xlim(0, max_pos)

    axes[-1].set_xlabel("Genomic position (Mb)")

    add_legend(fig, color_map, title=f"Repeat {args.taxonomy}")

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
