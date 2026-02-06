#!/usr/bin/env python3
"""
RepeatMasker normalization module.

Converts RepeatMasker `.out` files into a canonical, query-space,
BED-compatible tabular representation.
"""

from pathlib import Path
import pandas as pd


def split_class_family(class_family: str):
    if "/" in class_family:
        cls, fam = class_family.split("/", 1)
    else:
        cls = class_family
        fam = "NA"
    return cls, fam

def parse_repeatmasker_out(rm_out: Path) -> pd.DataFrame:
    rows = []
    with rm_out.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(("SW", "score", "perc", "----")):
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            rows.append({
                "score": int(parts[0]),
                "perc_div": float(parts[1]),
                "perc_del": float(parts[2]),
                "perc_ins": float(parts[3]),
                "chrom": parts[4],
                "start_1based": int(parts[5]),
                "end_1based": int(parts[6]),
                "strand": parts[8],
                "repeat_name": parts[9],
                "class_family": parts[10],
            })
    return pd.DataFrame(rows)


def normalize_rm(df: pd.DataFrame, strain: str) -> pd.DataFrame:
    ## zero based coordinate system
    df["start"] = df["start_1based"] - 1
    df["end"] = df["end_1based"]

    ## split taxonomy columns
    taxonomy = df["class_family"].apply(
        lambda x: pd.Series(split_class_family(x),
                            index=["repeat_class", "repeat_family"])
    )
    df = pd.concat([df, taxonomy], axis=1)
    df["strain"] = strain

    ## re-encode strandedness from +/C (compliment) to +/- (bed like)
    strand_map = {"+": "+", "C": "-"}
    df['strand'] = df['strand'].apply(lambda x: strand_map[x])

    if df["strand"].isna().any(): # check for failure (indication of malformed input)
        bad = df.loc[df["strand"].isna(), "strand"].unique()
        raise ValueError(f"Unexpected strand values in RM output: {bad}")
    
    norm = df[[
        "chrom",
        "start",
        "end",
        "strand",
        "repeat_name",
        "repeat_class",
        "repeat_family",
        "score",
        "perc_div",
        "perc_del",
        "perc_ins",
        "strain",
    ]].copy()

    return norm


def run_from_cli(args):
    df = parse_repeatmasker_out(Path(args.rm_out))
    if args.contig:
        df = df[df["chrom"] == args.contig]
    norm = normalize_rm(df, args.strain)
    norm.to_csv(args.out, sep="\t", index=False)
