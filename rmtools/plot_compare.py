"""Compare assemblies using shared taxonomy colors and genomic extent."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from .plot_panel import plot_panel
from .rm_track import load_data, choose_taxonomy, make_color_map
from .depth_track import load_depth, subset_depth
from .agp_track import load_agp, subset_agp
from .universal import parse_region, load_sizes


def _optional_path(row, name):
    value = row.get(name)
    return str(value) if pd.notna(value) and str(value).strip() else None


def run_from_cli(args):
    control = pd.read_csv(args.control, sep="\t", keep_default_na=False)
    if "strain" not in control:
        raise ValueError("Control file missing column: strain")
    if control.empty:
        raise ValueError("Control file has no strains")
    entries = []
    for _, row in control.iterrows():
        paths = {name: _optional_path(row, name) for name in ("rm", "depth", "agp")}
        if not any(paths.values()):
            raise ValueError(f"No tracks supplied for strain {row['strain']}")
        entries.append((row["strain"], paths, load_sizes(_optional_path(row, "sizes"))))

    for region in args.contigs:
        contig, start, end = parse_region(region)
        categories = set()
        extent = end or 0
        for strain, paths, sizes in entries:
            extent = max(extent, sizes.get(contig, 0))
            if paths["rm"]:
                df = load_data(Path(paths["rm"]), contig)
                if not df.empty:
                    extent = max(extent, int(df.end.max()))
                if start is not None:
                    df = df[(df.end > start) & (df.start < end)]
                categories.update(choose_taxonomy(df, args.taxonomy).unique())
            if paths["depth"]:
                df = subset_depth(load_depth(Path(paths["depth"])), contig)
                if not df.empty:
                    extent = max(extent, int(df.pos.max()) + 1)
            if paths["agp"]:
                df = subset_agp(load_agp(Path(paths["agp"])), contig)
                if not df.empty:
                    extent = max(extent, int(df.obj_end.max()))
        color_map = make_color_map(categories)
        fig = plt.figure(figsize=(12, len(entries) * 5))
        fig.suptitle(f"{region} comparison across strains", fontsize=16, fontweight="bold")
        outer = fig.add_gridspec(len(entries), 1, hspace=0.4)
        for i, (strain, paths, sizes) in enumerate(entries):
            heights = [h for name, h in (("rm", 3), ("depth", 2), ("agp", 1)) if paths[name]]
            gs = outer[i].subgridspec(len(heights), 1, height_ratios=heights)
            axes = plot_panel(region, rm_path=paths["rm"], depth_path=paths["depth"],
                              agp_path=paths["agp"], rm_taxonomy=args.taxonomy,
                              rm_bin_size=args.rm_bin, depth_bin_size=args.depth_bin,
                              fig=fig, gs=gs, rm_color_map=color_map,
                              contig_length=sizes.get(contig))
            if start is None and extent:
                axes[0].set_xlim(0, extent)
            axes[0].set_title(strain, loc="left", fontsize=11, fontweight="bold", pad=6)
        out = f"{args.out_prefix}_{region}.pdf"
        print(f"writing {out}")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
