rmtools

rmtools is a lightweight command-line toolkit for working with RepeatMasker annotations and related genome assembly diagnostics. It provides utilities to normalize RepeatMasker output and generate diagnostic and publication-quality plots that integrate repeat annotations, read depth, and scaffold structure.

The tool is designed primarily for genome assembly curation and validation, particularly in repeat-rich and heterochromatic regions.

---

INSTALLATION

The project is intended to be run from a Python environment with standard scientific packages installed.

Install in editable mode:

    pip install -e .

Dependencies (handled via pyproject.toml) include:
- Python 3.9 or newer
- pandas
- matplotlib

---

COORDINATE CONVENTIONS

All plotting commands use a consistent region syntax:

    CHROM
    CHROM:start-end

Examples:
    X
    X:30000000-35000000

Coordinates are interpreted in assembly coordinate space. Plots are clipped to the requested region, and tracks within a panel share the same x-axis.

---

COMMAND OVERVIEW

General usage:

    rmtools <command> [options]

Available commands:

- normalize
- plot-contig
- plot-multi
- agp-track
- depth-track
- panel

---

NORMALIZE REPEATMASKER OUTPUT

Normalize a RepeatMasker .out file into a tabular TSV format suitable for downstream plotting.

Example:

    rmtools normalize \
      --rm-out input.out \
      --out normalized.tsv \
      --strain BL25211 \
      --contig X

Arguments:
- --rm-out : RepeatMasker .out file
- --out : output TSV
- --strain : strain or assembly label
- --contig : optional contig filter

---

PLOT REPEATS ALONG A CONTIG

Plot RepeatMasker annotations along a contig or genomic region, optionally binned.

Example:

    rmtools plot-contig \
      --rm normalized.tsv \
      --region X:30000000-35000000 \
      --taxonomy class \
      --bin-size 50000 \
      --out repeats.pdf

Arguments:
- --rm : normalized RepeatMasker TSV
- --region : contig or region
- --taxonomy : class, family, or name
- --bin-size : bin size in base pairs (optional)
- --out : output figure

---

PLOT AGP SCAFFOLD STRUCTURE

Visualize scaffold structure from an AGP file.

Example:

    rmtools agp-track \
      --agp ragtag.scaffold.agp \
      --region X_RagTag \
      --out agp.pdf

Each W component is plotted on its own horizontal layer. Gaps are implicit, making scaffold joins and structure easy to inspect.

---

PLOT READ DEPTH

Plot read depth from samtools depth -a output.

Example:

    rmtools depth-track \
      --depth BL25211.depth.tsv \
      --region X:30000000-35000000 \
      --bin-size 10000 \
      --out depth.pdf

This visualization is useful for identifying:
- collapsed repeats (depth spikes)
- over-expanded regions (depth dips)
- scaffold or assembly errors

---

MULTI-TRACK DIAGNOSTIC PANEL (RECOMMENDED)

The panel command combines multiple tracks into a single diagnostic view with a shared x-axis.

Supported tracks:
- RepeatMasker annotations
- Read depth
- AGP scaffold structure

Example:

    rmtools panel \
      --region X:30000000-35000000 \
      --rm BL25211.rm.tsv \
      --depth BL25211.depth.tsv \
      --agp ragtag.scaffold.agp \
      --out panel.pdf

Any subset of tracks may be supplied, but at least one must be provided.

Arguments:
- --region : contig or region
- --rm : RepeatMasker TSV (optional)
- --depth : depth TSV (optional)
- --agp : AGP file (optional)
- --taxonomy : repeat taxonomy level (default: class)
- --rm-bin : repeat bin size (default: 50000 bp)
- --depth-bin : depth bin size (default: 10000 bp)
- --out : output figure

This panel is intended for rapid assembly sanity checking by visualizing sequence composition, read support, and scaffold structure in the same coordinate frame.

---

BATCH PLOTTING

The plot-multi command supports batch plotting using a control file (interface subject to refinement).

Example:

    rmtools plot-multi \
      --control control.tsv \
      --taxonomy class \
      --bin-size 50000 \
      --out batch.pdf

---

TYPICAL WORKFLOW

1. Run RepeatMasker on the assembly
2. Normalize RepeatMasker output using rmtools normalize
3. Map reads to the assembly and compute depth using samtools
4. Generate AGP during scaffolding (e.g. RagTag)
5. Use rmtools panel to inspect regions of interest
6. Identify misassemblies, collapsed repeats, or scaffold artifacts

---

PROJECT STATUS

rmtools is an actively developed research tool. Interfaces and internal APIs may evolve as new diagnostics and plotting modes are added.

The focus is on correctness, interpretability, and reproducibility rather than long-term backward compatibility at this stage.




