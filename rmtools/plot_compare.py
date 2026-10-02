#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from .plot_panel import plot_panel


def run_from_cli(args):

    control = pd.read_csv(args.control, sep="\t")

    required = ["strain", "rm"]
    for col in required:
        if col not in control.columns:
            raise ValueError(f"Control file missing column: {col}")

    for contig in args.contigs:

        n_strains = len(control)

        fig = plt.figure(figsize=(12, n_strains * 3))

        fig.suptitle(
            f"{contig} comparison across strains",
            fontsize=16,
            fontweight="bold"
        )

        gs_outer = fig.add_gridspec(
            nrows=n_strains,
            ncols=1,
            hspace=0.25
        )

        for i, row in control.iterrows():

            strain = row["strain"]

            rm = row["rm"] if "rm" in row else None
            depth = row["depth"] if "depth" in row else None
            agp = row["agp"] if "agp" in row else None

            tracks = sum([
                rm is not None,
                depth is not None,
                agp is not None
            ])

            panel_spec = gs_outer[i]

            panel_gs = panel_spec.subgridspec(
                tracks,
                1
            )

            axes = plot_panel(
                region=contig,
                rm_path=rm,
                depth_path=depth,
                agp_path=agp,
                rm_taxonomy=args.taxonomy,
                rm_bin_size=args.rm_bin,
                depth_bin_size=args.depth_bin,
                fig=fig,
                gs=panel_gs
            )

            axes[0].set_title(
                strain,
                loc="left",
                fontsize=11,
                fontweight="bold",
                pad=6
            )

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        out = f"{args.out_prefix}_{contig}.pdf"
        print(f"writing {out}")
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()