#!/usr/bin/env python3
"""
Functions necesary to more than one tool chain in the project
"""

def parse_region(region: str):
    """
    Parse region string:
      - contig
      - contig:start-end

    Returns (contig, start, end) where start/end may be None.
    """
    if ":" not in region:
        return region, None, None

    try:
        contig, coords = region.split(":")
        start, end = coords.split("-")
        return contig, int(start), int(end)
    except Exception:
        raise ValueError(
            f"Invalid region format '{region}'. "
            "Expected contig or contig:start-end"
        )