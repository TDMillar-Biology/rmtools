import argparse
from . import normalize, rm_track, plot_multi, agp_track, depth_track, plot_panel

def main():
    parser = argparse.ArgumentParser(
        prog="rmtools",
        description="RepeatMasker normalization and plotting tools"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    norm = subparsers.add_parser("normalize", help="Normalize RepeatMasker .out")
    norm.add_argument("--rm-out", required=True)
    norm.add_argument("--out", required=True)
    norm.add_argument("--strain", required=True)
    norm.add_argument("--contig", default=None)

    plot = subparsers.add_parser("plot-contig", help="Plot repeats along contig")
    plot.add_argument("--rm", required=True)
    plot.add_argument("--region", required=True)
    plot.add_argument("--taxonomy", choices=["class", "family", "name"], default="class")
    plot.add_argument("--bin-size", type=int, default=None)
    plot.add_argument("--out", required=True)

    parser_multi = subparsers.add_parser("plot-multi")
    parser_multi.add_argument("--control", required=True)
    parser_multi.add_argument("--taxonomy", default="class")
    parser_multi.add_argument("--bin-size", type=int, required=True)
    parser_multi.add_argument("--out", required=True)

    parser_agp = subparsers.add_parser("agp-track")
    parser_agp.add_argument("--agp", required=True)
    parser_agp.add_argument("--region", required=True)
    parser_agp.add_argument("--out", required=True)

    parser_depth = subparsers.add_parser("depth-track")
    parser_depth.add_argument("--depth", required=True)
    parser_depth.add_argument("--region", required=True)
    parser_depth.add_argument("--out", required=True)
    parser_depth.add_argument("--bin-size", type=int, default=10_000)

    parser_panel = subparsers.add_parser("panel")
    parser_panel.add_argument("--rm", required=False)
    parser_panel.add_argument("--depth", required=False)
    parser_panel.add_argument("--agp", required=False)
    parser_panel.add_argument("--region", required=True)
    parser_panel.add_argument("--out", required=True)
    parser_panel.add_argument("--taxonomy", choices=["class", "family", "name"], default="class")
    parser_panel.add_argument("--rm-bin", type=int, default=50_000)
    parser_panel.add_argument("--depth-bin", type=int, default=10_000)
    
    args = parser.parse_args()

    if args.command == "normalize":
        normalize.run_from_cli(args)
    elif args.command == "plot-contig":
        rm_track.run_from_cli(args)
    elif args.command == "plot-multi":
        plot_multi.run_from_cli(args)
    elif args.command == "agp-track":
        agp_track.run_from_cli(args)
    elif args.command == "depth-track":
        depth_track.run_from_cli(args)
    elif args.command == "panel":
        plot_panel.run_from_cli(args)