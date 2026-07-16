#!/usr/bin/env python3
"""Single-cell periodicity heatmaps from transcriptome BAM.

Replicates the per-cell triplet periodicity heatmaps from Duré et al. 2024
(Figure 1G-I). For each cell, maps P-site positions relative to annotated
start and stop codons across all CDS-containing transcripts, then plots a
heatmap where each row is a single cell and the 3-nt periodicity pattern
is visible as vertical stripes.

Input:
    - Transcriptome dedup BAM with CB tags (from pipeline tx_bam/)
    - riboWaltz annotation CSV (transcript, l_tr, l_utr5, l_cds, l_utr3)
    - P-site offsets CSV (from riboWaltz output)

Output:
    - Heatmap PDFs: periodicity from start and stop codons
    - CSV: per-cell read count summary

Usage:
    python3 single_cell_periodicity.py \\
        --bam results/sample/tx_bam/sample.tx.dedup.bam \\
        --annotation ribowaltz_m25_anno.csv \\
        --psite-offsets results/sample/ribowaltz/psite_offsets.csv \\
        --output-dir results/sample/periodicity
"""

import argparse
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


def load_annotation(path):
    """Load riboWaltz annotation, keep only CDS-containing transcripts.
    Returns dict: transcript_id -> (l_utr5, l_cds, l_utr3, l_tr)
    """
    df = pd.read_csv(path)
    df = df[df["l_cds"] > 0]
    anno = {}
    for _, row in df.iterrows():
        anno[row["transcript"]] = (
            int(row["l_utr5"]),
            int(row["l_cds"]),
            int(row["l_utr3"]),
            int(row["l_tr"]),
        )
    return anno


def load_psite_offsets(path):
    """Load P-site offsets from riboWaltz output.
    Returns dict: read_length -> offset_from_5prime
    """
    df = pd.read_csv(path)
    offsets = {}
    for _, row in df.iterrows():
        length = int(row["length"])
        offset = row["corrected_offset_from_5"]
        if pd.notna(offset):
            offsets[length] = int(offset)
    return offsets


def compute_psite_positions(bam_path, anno, offsets, min_length=26, max_length=34):
    """Extract P-site positions relative to start and stop codons per cell.

    For each read:
    1. Get cell barcode (CB tag) and read length
    2. Look up transcript annotation (UTR5/CDS boundaries)
    3. Compute P-site position = alignment start + offset
    4. Convert to distance from start codon (P-site pos - l_utr5)
       and distance from stop codon (P-site pos - (l_utr5 + l_cds))

    Returns:
        cell_start: dict {cell_barcode: {distance_from_start: count}}
        cell_stop:  dict {cell_barcode: {distance_from_stop: count}}
        cell_total: dict {cell_barcode: total_read_count}
    """
    cell_start = defaultdict(lambda: defaultdict(int))
    cell_stop = defaultdict(lambda: defaultdict(int))
    cell_total = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_reverse:
                continue

            length = read.query_length
            if length is None or length < min_length or length > max_length:
                continue

            # Get cell barcode (CB = corrected, CR = raw; tx BAM may only have CR)
            if read.has_tag("CB"):
                cb = read.get_tag("CB")
            elif read.has_tag("CR"):
                cb = read.get_tag("CR")
            else:
                continue
            if cb == "-":
                continue

            # Get offset for this read length
            if length not in offsets:
                continue

            # Get transcript annotation
            ref_name = read.reference_name
            if ref_name not in anno:
                continue

            l_utr5, l_cds, l_utr3, l_tr = anno[ref_name]

            # P-site position in transcript coordinates
            psite_pos = read.reference_start + offsets[length]

            # Distance from start codon (position 0 = first nt of CDS)
            dist_start = psite_pos - l_utr5

            # Distance from stop codon (position 0 = first nt after CDS)
            dist_stop = psite_pos - (l_utr5 + l_cds)

            cell_start[cb][dist_start] += 1
            cell_stop[cb][dist_stop] += 1
            cell_total[cb] += 1

    return cell_start, cell_stop, cell_total


def build_heatmap_matrix(cell_positions, cell_total, min_pos, max_pos,
                         min_reads=50):
    """Build a matrix for heatmap plotting.

    Rows = cells (sorted by total reads, descending)
    Columns = positions from min_pos to max_pos

    Args:
        cell_positions: dict {cb: {position: count}}
        cell_total: dict {cb: total_reads}
        min_pos: leftmost position to include
        max_pos: rightmost position to include
        min_reads: minimum reads per cell to include

    Returns:
        matrix: numpy array (n_cells x n_positions)
        cell_order: list of cell barcodes in row order
        positions: numpy array of position values
    """
    # Filter cells by minimum reads
    valid_cells = [cb for cb, total in cell_total.items() if total >= min_reads]
    # Sort by total reads (highest on top)
    valid_cells.sort(key=lambda cb: cell_total[cb], reverse=True)

    positions = np.arange(min_pos, max_pos + 1)
    matrix = np.zeros((len(valid_cells), len(positions)))

    for i, cb in enumerate(valid_cells):
        for j, pos in enumerate(positions):
            matrix[i, j] = cell_positions[cb].get(pos, 0)

    return matrix, valid_cells, positions


def plot_heatmap(matrix, positions, output_path, title,
                 landmark_pos=0, landmark_label="Start",
                 log_transform=True, cmap="viridis"):
    """Plot a periodicity heatmap.

    Args:
        matrix: cells x positions array
        positions: position values for x-axis
        output_path: path to save PDF
        title: plot title
        landmark_pos: position of the landmark (start=0 or stop=0)
        landmark_label: label for the landmark
        log_transform: apply log2(x+1) transform
        cmap: colormap
    """
    if log_transform:
        plot_data = np.log2(matrix + 1)
        clabel = "log2(reads + 1)"
    else:
        plot_data = matrix
        clabel = "reads"

    fig, ax = plt.subplots(figsize=(16, max(6, len(matrix) * 0.03)))

    im = ax.imshow(plot_data, aspect="auto", cmap=cmap, interpolation="none",
                   extent=[positions[0] - 0.5, positions[-1] + 0.5,
                           len(matrix) - 0.5, -0.5])

    # Mark the landmark position (start or stop codon)
    landmark_idx = landmark_pos
    ax.axvline(x=landmark_idx, color="red", linewidth=0.8, linestyle="--",
               alpha=0.7)

    # Mark 3-nt periodicity lines in CDS region
    if landmark_label == "Start":
        for p in range(0, positions[-1] + 1, 3):
            ax.axvline(x=p, color="white", linewidth=0.1, alpha=0.3)

    ax.set_xlabel(f"Distance from {landmark_label.lower()} (nt)")
    ax.set_ylabel("Single cells")
    ax.set_title(title)

    # Tick marks every 30 nt
    tick_step = 30
    ticks = np.arange(
        (positions[0] // tick_step) * tick_step,
        positions[-1] + 1,
        tick_step
    )
    ax.set_xticks(ticks)

    cbar = plt.colorbar(im, ax=ax, shrink=0.5, label=clabel)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam", required=True,
                        help="Transcriptome dedup BAM with CB tags")
    parser.add_argument("--annotation", required=True,
                        help="riboWaltz annotation CSV")
    parser.add_argument("--psite-offsets", required=True,
                        help="P-site offsets CSV from riboWaltz")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory")
    parser.add_argument("--min-reads", type=int, default=50,
                        help="Minimum reads per cell to include (default: 50)")
    parser.add_argument("--start-window", type=int, nargs=2, default=[-20, 200],
                        help="Window around start codon (default: -20 200)")
    parser.add_argument("--stop-window", type=int, nargs=2, default=[-200, 20],
                        help="Window around stop codon (default: -200 20)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    print("Loading annotation...")
    anno = load_annotation(args.annotation)
    print(f"  {len(anno)} transcripts with CDS")

    print("Loading P-site offsets...")
    offsets = load_psite_offsets(args.psite_offsets)
    print(f"  Offsets: {offsets}")

    print(f"Processing BAM: {args.bam}")
    cell_start, cell_stop, cell_total = compute_psite_positions(
        args.bam, anno, offsets
    )
    print(f"  {len(cell_total)} cells found")
    print(f"  {sum(cell_total.values()):,} total P-site positions")

    # Filter and report
    n_above = sum(1 for v in cell_total.values() if v >= args.min_reads)
    print(f"  {n_above} cells with >= {args.min_reads} reads")

    if n_above == 0:
        print("ERROR: No cells pass the minimum read filter. Try lowering --min-reads.")
        sys.exit(1)

    # Save per-cell summary
    summary = pd.DataFrame([
        {"cell_barcode": cb, "total_reads": total}
        for cb, total in sorted(cell_total.items(), key=lambda x: -x[1])
    ])
    summary.to_csv(os.path.join(args.output_dir, "cell_read_counts.csv"),
                   index=False)
    print(f"  Median reads/cell: {summary['total_reads'].median():.0f}")

    # Build and plot start codon heatmap
    print("Plotting start codon heatmap...")
    matrix_start, cells_start, pos_start = build_heatmap_matrix(
        cell_start, cell_total,
        min_pos=args.start_window[0], max_pos=args.start_window[1],
        min_reads=args.min_reads
    )
    plot_heatmap(
        matrix_start, pos_start,
        os.path.join(args.output_dir, "periodicity_start.png"),
        f"Distance from start in single cells (26-34 nt footprints)\n"
        f"({len(cells_start)} cells, min {args.min_reads} reads)",
        landmark_pos=0, landmark_label="Start"
    )

    # Build and plot stop codon heatmap
    print("Plotting stop codon heatmap...")
    matrix_stop, cells_stop, pos_stop = build_heatmap_matrix(
        cell_stop, cell_total,
        min_pos=args.stop_window[0], max_pos=args.stop_window[1],
        min_reads=args.min_reads
    )
    plot_heatmap(
        matrix_stop, pos_stop,
        os.path.join(args.output_dir, "periodicity_stop.png"),
        f"Distance from stop in single cells (26-34 nt footprints)\n"
        f"({len(cells_stop)} cells, min {args.min_reads} reads)",
        landmark_pos=0, landmark_label="Stop"
    )

    print(f"Done. Outputs in {args.output_dir}/")


if __name__ == "__main__":
    main()
