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
    bin_intervals_dominant,
    plot_binned,
    make_color_map,
    add_legend,
    plot_raw_intervals
)

from .universal import parse_region

def run_from_cli(args):
    
    df = load_data(Path(args.rm))
    taxonomy_col = choose_taxonomy(df, args.taxonomy)

    categories = taxonomy_col.unique()
    color_map = make_color_map(categories)
    number_plots = len(args.main)
    fig, axes = plt.subplots(1,number_plots)
    for i, contig in enumerate(args.main):
        df = df[df["chrom"] == contig]
        ax = axes[i]
        if args.bin_size is None:
            plot_raw_intervals(df, taxonomy_col, ax, color_map)
        else:
            binned = bin_intervals_dominant(df, taxonomy_col, args.bin_size)
            plot_binned(binned, ax, color_map)


        ## Format axis labels
        # x axis
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, pos: f"{x/1e6:.1f}"))
        ax.xaxis.set_major_locator(mticker.MultipleLocator(5e6))
        ax.xaxis.set_minor_locator(mticker.MultipleLocator(1e6))
        ax.set_xlabel("Genomic position (Mb)")

            # y axis
        ax.set_ylabel("Repeat coverage (bp)")

    add_legend(ax, color_map, title=f"Repeat {args.taxonomy}")


    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")