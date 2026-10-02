import argparse
from .universal import positive_int
from . import normalize, rm_track, plot_multi, agp_track, depth_track, plot_panel, plot_main, size, plot_assembly, plot_compare

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
    plot.add_argument("--bin-size", type=positive_int, default=None)
    plot.add_argument("--out", required=True)

    parser_multi = subparsers.add_parser("plot-multi")
    parser_multi.add_argument("--control", required=True)
    parser_multi.add_argument("--taxonomy", choices=["class", "family", "name"], default="class")
    parser_multi.add_argument("--bin-size", type=positive_int, required=True)
    parser_multi.add_argument("--out", required=True)

    parser_main = subparsers.add_parser("plot-main")
    parser_main.add_argument("--main", nargs='+', required=True)
    parser_main.add_argument("--rm", required=True)
    parser_main.add_argument("--taxonomy", choices=["class", "family", "name"], default="class")
    parser_main.add_argument("--bin-size", default = 50_000, type=positive_int, required=False)
    parser_main.add_argument("--out", required=True)

    parser_agp = subparsers.add_parser("agp-track")
    parser_agp.add_argument("--agp", required=True)
    parser_agp.add_argument("--region", required=True)
    parser_agp.add_argument("--out", required=True)

    parser_depth = subparsers.add_parser("depth-track")
    parser_depth.add_argument("--depth", required=True)
    parser_depth.add_argument("--region", required=True)
    parser_depth.add_argument("--out", required=True)
    parser_depth.add_argument("--bin-size", type=positive_int, default=10_000)

    parser_panel = subparsers.add_parser("panel")
    parser_panel.add_argument("--rm", required=False)
    parser_panel.add_argument("--depth", required=False)
    parser_panel.add_argument("--agp", required=False)
    parser_panel.add_argument("--region", required=True)
    parser_panel.add_argument("--out", required=True)
    parser_panel.add_argument("--taxonomy", choices=["class", "family", "name"], default="class")
    parser_panel.add_argument("--rm-bin", type=positive_int, default=50_000)
    parser_panel.add_argument("--depth-bin", type=positive_int, default=10_000)

    parser_assembly = subparsers.add_parser("plot-assembly")
    parser_assembly.add_argument("--main", nargs="+", required=True)
    parser_assembly.add_argument("--rm")
    parser_assembly.add_argument("--depth")
    parser_assembly.add_argument("--agp")
    parser_assembly.add_argument("--taxonomy", choices=["class", "family", "name"], default="class")
    parser_assembly.add_argument("--rm-bin", type=positive_int, default=50_000)
    parser_assembly.add_argument("--depth-bin", type=positive_int, default=10_000)
    parser_assembly.add_argument("--out", required=True)

    parser_compare = subparsers.add_parser("plot-compare",help="Compare multiple strains across contigs")
    parser_compare.add_argument("--control",required=True,help="TSV listing strain, rm, depth, agp files")
    parser_compare.add_argument("--contigs",nargs="+",required=True,help="Contigs to plot (e.g. chr2L chr2R chr3L)")
    parser_compare.add_argument("--taxonomy",default="class",choices=["class", "family", "name"])
    parser_compare.add_argument("--rm-bin",type=positive_int,default=50_000)
    parser_compare.add_argument("--depth-bin",type=positive_int,default=10_000)
    parser_compare.add_argument("--out-prefix",required=True)

    parser_size = subparsers.add_parser("size")
    parser_size.add_argument("--rm", required=True)
    
    for command_parser in (plot, parser_main, parser_panel, parser_assembly):
        command_parser.add_argument("--sizes", help="Contig lengths TSV or FASTA .fai")

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
    elif args.command == "plot-main":
        plot_main.run_from_cli(args)
    elif args.command == "size":
        size.run_from_cli(args)
    elif args.command == "plot-assembly":
        plot_assembly.run_from_cli(args)
    elif args.command == "plot-compare":
        plot_compare.run_from_cli(args)