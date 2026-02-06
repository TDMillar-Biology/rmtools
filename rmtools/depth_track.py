#!/usr/bin/env python3
"""
Read depth plotting utilities.

Designed for samtools depth -a output.
Provides region-aware, optionally binned depth tracks
for assembly diagnostics.
"""

from pathlib import Path
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from .universal import parse_region

# ------------------------------------------------------------
# Loading
# ------------------------------------------------------------

def load_depth(path: Path):
    """
    Load samtools depth output.

    Expected columns:
      chrom  position  depth
    """
    df = pd.read_csv(
        path,
        sep="\t",
        names=["chrom", "pos", "depth"],
        dtype={
            "chrom": str,
            "pos": int,
            "depth": int,
        }
    )
    return df


def subset_depth(df, contig, start=None, end=None):
    """
    Subset depth dataframe to contig and optional region.
    """
    sub = df[df["chrom"] == contig].copy()

    if start is not None:
        sub = sub[(sub.pos >= start) & (sub.pos <= end)].copy()

    return sub


# ------------------------------------------------------------
# Binning / smoothing
# ------------------------------------------------------------

def bin_depth(df, bin_size, statistic="mean"):
    """
    Bin depth values.

    statistic: mean | median | sum
    """
    df = df.copy()
    df["bin"] = df.pos // bin_size

    if statistic == "mean":
        agg = df.groupby("bin")["depth"].mean()
    elif statistic == "median":
        agg = df.groupby("bin")["depth"].median()
    elif statistic == "sum":
        agg = df.groupby("bin")["depth"].sum()
    else:
        raise ValueError(f"Unknown statistic: {statistic}")

    out = agg.reset_index()
    out["x"] = out["bin"] * bin_size
    return out


# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

def plot_depth(
    depth_df,
    ax,
    bin_size=None,
    statistic="mean",
    region_start=None,
    rebase=True,
    color="black",
    linewidth=1.0,
):
    """
    Plot read depth.

    If bin_size is provided, depth is binned before plotting.
    """
    df = depth_df.copy()

    if bin_size is not None:
        df = bin_depth(df, bin_size, statistic=statistic)
        x = df["x"]
        y = df["depth"]
    else:
        x = df["pos"]
        y = df["depth"]

    # Optional rebasing
    if rebase and region_start is not None:
        x = x - region_start

    ax.plot(x, y, color=color, linewidth=linewidth)
    ax.set_ylabel("Depth")

def run_from_cli(args):
    contig, start, end = parse_region(args.region)

    depth = load_depth(args.depth)
    depth_sub = subset_depth(depth, contig, start, end)

    fig, ax = plt.subplots(figsize=(12, 2))

    plot_depth(
        depth_sub,
        ax,
        bin_size=args.bin_size,
        region_start=start,
    )
    ax.set_xlabel("Genomic position (bp)")

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)