"""Shared coordinate and input helpers."""
from pathlib import Path


def parse_region(region: str):
    """Parse a contig or a zero-based, half-open contig:start-end region."""
    if ":" not in region:
        if not region:
            raise ValueError("Contig must not be empty")
        return region, None, None
    try:
        contig, coords = region.rsplit(":", 1)
        start, end = map(int, coords.split("-"))
        if not contig or not 0 <= start < end:
            raise ValueError
        return contig, start, end
    except ValueError as exc:
        raise ValueError(f"Invalid region '{region}'; expected contig or contig:start-end with 0 <= start < end") from exc


def positive_int(value):
    import argparse
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("Expected a positive integer") from exc
    if result <= 0:
        raise argparse.ArgumentTypeError("Expected a positive integer")
    return result


def load_sizes(path):
    """Read a tab/whitespace-separated chrom,length file or FASTA .fai."""
    if path is None:
        return {}
    sizes = {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, length, *_ = line.split()
            length = int(length)
            if length <= 0 or chrom in sizes:
                raise ValueError(f"Invalid or duplicate contig size: {chrom}")
            sizes[chrom] = length
    return sizes
