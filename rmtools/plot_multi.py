"""Plot repeat tracks with one shared taxonomy color map."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from .rm_track import (load_data, choose_taxonomy, bin_intervals_repeat_composition,
                       plot_binned, make_color_map, add_legend)
from .universal import parse_region, load_sizes


def load_control_file(path: Path):
    df = pd.read_csv(path, sep=r"\s+", keep_default_na=False)
    missing = {"path", "contig", "label"} - set(df.columns)
    if missing:
        raise ValueError(f"Control file missing columns: {missing}")
    if df.empty:
        raise ValueError("Control file has no tracks")
    return df


def plot_multi(control_df, taxonomy, bin_size, figsize=(12, None)):
    if control_df.empty:
        raise ValueError("Control file has no tracks")
    prepared = []
    categories = set()
    for row in control_df.itertuples():
        contig, start, end = parse_region(row.contig)
        df = load_data(Path(row.path), contig)
        sizes_path = getattr(row, "sizes", None)
        sizes = load_sizes(sizes_path if pd.notna(sizes_path) and sizes_path else None)
        if start is not None:
            df = df[(df.end > start) & (df.start < end)].copy()
            df["start"] = df.start.clip(lower=start)
            df["end"] = df.end.clip(upper=end)
        origin = 0 if start is None else start
        extent = end if end is not None else sizes.get(contig)
        if extent is None and not df.empty:
            extent = int(df.end.max())
        df["start"] -= origin
        df["end"] -= origin
        taxa = choose_taxonomy(df, taxonomy)
        categories.update(taxa.unique())
        prepared.append((row.label, df, taxa, None if extent is None else extent - origin))
    color_map = make_color_map(categories)
    fig, axes = plt.subplots(len(prepared), figsize=(figsize[0], figsize[1] or 2.2 * len(prepared)),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]
    max_extent = 0
    for ax, (label, df, taxa, extent) in zip(axes, prepared):
        bins = bin_intervals_repeat_composition(df, taxa, bin_size, end=extent)
        plot_binned(bins, ax, color_map)
        if df.empty:
            ax.text(0.5, 0.5, "No repeat annotations", transform=ax.transAxes, ha="center")
        ax.set_ylabel(label, rotation=0, ha="right", va="center")
        ax.yaxis.set_label_coords(-0.10, 0.5)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x/1e6:.1f}"))
        ax.xaxis.set_major_locator(mticker.MultipleLocator(5e6))
        ax.xaxis.set_minor_locator(mticker.MultipleLocator(1e6))
        max_extent = max(max_extent, extent or 0)
    if max_extent:
        axes[0].set_xlim(0, max_extent)
    axes[-1].set_xlabel("Position from contig or region start (Mb)")
    add_legend(axes[0], color_map, title=f"Repeat {taxonomy}")
    return fig


def run_from_cli(args):
    fig = plot_multi(load_control_file(Path(args.control)), args.taxonomy, args.bin_size)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
