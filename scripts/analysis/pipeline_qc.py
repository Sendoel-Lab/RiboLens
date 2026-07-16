#!/usr/bin/env python3
"""Pipeline QC: diagnostic plots from scribo-seq pipeline outputs.

Parses log files and outputs from each pipeline stage to generate a
multi-panel QC report showing read attrition, mapping stats, UMI
distributions, and cell calling.

Usage:
    python3 pipeline_qc.py \\
        --outdir results \\
        --sample young1_S1_L001_R1_001 \\
        --output-dir results/young1_S1_L001_R1_001/qc_report
"""

import argparse
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.io import mmread


# =============================================================================
# Log parsers - extract read counts from each pipeline step
# =============================================================================

def parse_fastp(json_path):
    """Extract read counts before/after trimming from fastp JSON."""
    if not os.path.exists(json_path):
        return None
    with open(json_path) as f:
        data = json.load(f)
    bf = data["summary"]["before_filtering"]
    af = data["summary"]["after_filtering"]
    return {
        "raw_reads": bf["total_reads"],
        "after_trim": af["total_reads"],
    }


def parse_bowtie2(log_path):
    """Extract alignment rate and read counts from bowtie2 log."""
    if not os.path.exists(log_path):
        return None
    with open(log_path) as f:
        text = f.read()
    result = {}
    # Total reads is first line: "12067405 reads; of these:"
    m = re.search(r"(\d+) reads; of these:", text)
    if m:
        result["input_reads"] = int(m.group(1))
    # Unaligned = non-rRNA
    m = re.search(r"(\d+) \([\d.]+%\) aligned 0 times", text)
    if m:
        result["after_rrna"] = int(m.group(1))
    # Alignment rate
    m = re.search(r"([\d.]+)% overall alignment rate", text)
    if m:
        result["rrna_pct"] = float(m.group(1))
    return result


def parse_ribodetector(log_path):
    """Extract non-rRNA read count from ribodetector log."""
    if not os.path.exists(log_path):
        return None
    with open(log_path) as f:
        text = f.read()
    # Strip ANSI color codes (ribodetector uses colored output)
    text = re.sub(r"\x1b\[[0-9;]*m", "", text)
    result = {}
    m = re.search(r"Processed\s+(\d+)\s+sequences in total", text)
    if m:
        result["input_reads"] = int(m.group(1))
    m = re.search(r"Detected\s+(\d+)\s+non-rRNA", text)
    if m:
        result["after_ribodetector"] = int(m.group(1))
    m = re.search(r"Detected.*?(\d+).*?rRNA seq", text)
    if m:
        result["rrna_reads"] = int(m.group(1))
    return result


def parse_star_log(log_path):
    """Extract mapping stats from STAR Log.final.out."""
    if not os.path.exists(log_path):
        return None
    result = {}
    with open(log_path) as f:
        for line in f:
            if "|" in line:
                parts = line.strip().split("|")
                if len(parts) == 2:
                    key = parts[0].strip()
                    val = parts[1].strip()
                    result[key] = val
    return result


def parse_solo_summary(summary_path):
    """Parse STAR Solo Summary.csv."""
    if not os.path.exists(summary_path):
        return None
    result = {}
    with open(summary_path) as f:
        for line in f:
            parts = line.strip().split(",", 1)
            if len(parts) == 2:
                result[parts[0]] = parts[1]
    return result


# =============================================================================
# Plot functions
# =============================================================================

def plot_read_attrition(stages, counts, output_path):
    """Bar chart showing read counts at each pipeline stage.

    Visualises the 'funnel' of reads through the pipeline, showing where
    reads are lost (trimming, rRNA removal, mapping, CB filtering, dedup).
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Bar chart
    colors = plt.cm.Set2(np.linspace(0, 1, len(stages)))
    bars = ax1.barh(range(len(stages) - 1, -1, -1), counts, color=colors)
    ax1.set_yticks(range(len(stages) - 1, -1, -1))
    ax1.set_yticklabels(stages, fontsize=11)
    ax1.set_xlabel("Number of reads", fontsize=12)
    ax1.set_title("Read attrition through pipeline", fontsize=14, fontweight="bold")
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.0f}K"
    ))
    ax1.tick_params(axis="x", labelsize=10)

    # Add count labels on bars
    for bar, count in zip(bars, counts):
        pct = count / counts[0] * 100 if counts[0] > 0 else 0
        ax1.text(bar.get_width() + max(counts) * 0.01,
                 bar.get_y() + bar.get_height() / 2,
                 f"{count:,} ({pct:.1f}%)",
                 va="center", fontsize=10)

    # Step-by-step loss chart
    losses = []
    loss_labels = []
    for i in range(1, len(counts)):
        lost = counts[i - 1] - counts[i]
        losses.append(lost)
        loss_labels.append(f"{stages[i-1]} \u2192 {stages[i]}")

    colors_loss = plt.cm.Reds(np.linspace(0.3, 0.8, len(losses)))
    ax2.barh(range(len(losses) - 1, -1, -1), losses, color=colors_loss)
    ax2.set_yticks(range(len(losses) - 1, -1, -1))
    ax2.set_yticklabels(loss_labels, fontsize=10)
    ax2.set_xlabel("Reads lost", fontsize=12)
    ax2.set_title("Reads lost at each step", fontsize=14, fontweight="bold")
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.0f}K"
    ))
    ax2.tick_params(axis="x", labelsize=10)

    for bar, loss in zip(ax2.patches, losses):
        pct = loss / counts[0] * 100 if counts[0] > 0 else 0
        ax2.text(bar.get_width() + max(losses) * 0.01,
                 bar.get_y() + bar.get_height() / 2,
                 f"{loss:,} ({pct:.1f}%)",
                 va="center", fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_mapping_breakdown(star_log, output_path):
    """Stacked horizontal bar chart of STAR mapping categories."""
    if not star_log:
        return

    total = int(star_log.get("Number of input reads", "0"))
    unique = int(star_log.get("Uniquely mapped reads number", "0"))
    multi = int(star_log.get("Number of reads mapped to multiple loci", "0"))
    too_many = int(star_log.get("Number of reads mapped to too many loci", "0"))
    unmapped_short = int(star_log.get("Number of reads unmapped: too short", "0"))
    unmapped_mm = int(star_log.get("Number of reads unmapped: too many mismatches", "0"))
    unmapped_other = int(star_log.get("Number of reads unmapped: other", "0"))

    sizes = [unique, multi, too_many, unmapped_short, unmapped_mm, unmapped_other]
    labels = ["Unique", "Multi-map", "Too many loci",
              "Unmapped: too short", "Unmapped: mismatch", "Unmapped: other"]
    colors = ["#2ecc71", "#3498db", "#9b59b6", "#e74c3c", "#e67e22", "#95a5a6"]

    # Remove zero-count categories
    nonzero = [(s, l, c) for s, l, c in zip(sizes, labels, colors) if s > 0]
    if not nonzero:
        return
    sizes, labels, colors = zip(*nonzero)

    fig, ax = plt.subplots(figsize=(14, 3))

    # Stacked horizontal bar
    left = 0
    for size, label, color in zip(sizes, labels, colors):
        pct = size / total * 100 if total > 0 else 0
        bar = ax.barh(0, size, left=left, color=color, edgecolor="white",
                       linewidth=0.5, label=f"{label}  {size:,} ({pct:.1f}%)")
        # Add percentage text inside bar if wide enough
        if pct > 5:
            ax.text(left + size / 2, 0, f"{pct:.1f}%",
                    ha="center", va="center", fontsize=11, fontweight="bold",
                    color="white" if pct > 10 else "#333")
        left += size

    ax.set_xlim(0, total)
    ax.set_yticks([])
    ax.set_xlabel("Number of reads", fontsize=12)
    ax.set_title(f"STAR mapping breakdown ({total:,} input reads)",
                 fontsize=14, fontweight="bold")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.0f}K"
    ))
    ax.tick_params(axis="x", labelsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3,
              fontsize=10, frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_umi_knee(solo_base, output_path):
    """Knee plot from STAR Solo raw/filtered matrices."""
    raw_dir = os.path.join(solo_base, "raw")
    filt_dir = os.path.join(solo_base, "filtered")

    if not os.path.isdir(raw_dir) or not os.path.isdir(filt_dir):
        return

    raw_mtx = mmread(os.path.join(raw_dir, "matrix.mtx")).tocsc()
    raw_umis = np.sort(np.array(raw_mtx.sum(axis=0)).flatten())[::-1]

    filt_mtx = mmread(os.path.join(filt_dir, "matrix.mtx")).tocsc()
    filt_umis = np.sort(np.array(filt_mtx.sum(axis=0)).flatten())[::-1]
    filt_genes = np.sort(np.array((filt_mtx > 0).sum(axis=0)).flatten())[::-1]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Knee plot
    ax = axes[0]
    ax.plot(range(1, len(raw_umis) + 1), raw_umis, color="steelblue",
            linewidth=0.5, label="All barcodes")
    ax.plot(range(1, len(filt_umis) + 1), filt_umis, color="red",
            linewidth=2, label=f"Filtered ({len(filt_umis)} cells)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode rank", fontsize=12)
    ax.set_ylabel("UMI count", fontsize=12)
    ax.set_title(f"UMI knee plot\nMedian = {int(np.median(filt_umis)):,} UMI/cell",
                 fontsize=13, fontweight="bold")
    ax.axhline(y=np.median(filt_umis), color="gray", linestyle="--", alpha=0.5)
    ax.legend(fontsize=11)
    ax.tick_params(labelsize=10)

    # UMI per cell histogram
    ax = axes[1]
    ax.hist(filt_umis, bins=30, color="steelblue", edgecolor="white")
    ax.axvline(np.median(filt_umis), color="red", linestyle="--",
               label=f"Median = {int(np.median(filt_umis)):,}")
    ax.set_xlabel("UMI per cell", fontsize=12)
    ax.set_ylabel("Number of cells", fontsize=12)
    ax.set_title("UMI distribution (filtered cells)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11)
    ax.tick_params(labelsize=10)

    # Genes per cell histogram
    ax = axes[2]
    ax.hist(filt_genes, bins=30, color="darkorange", edgecolor="white")
    ax.axvline(np.median(filt_genes), color="red", linestyle="--",
               label=f"Median = {int(np.median(filt_genes)):,}")
    ax.set_xlabel("Genes per cell", fontsize=12)
    ax.set_ylabel("Number of cells", fontsize=12)
    ax.set_title("Gene detection (filtered cells)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11)
    ax.tick_params(labelsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", required=True,
                        help="Pipeline output directory (e.g., results)")
    parser.add_argument("--sample", required=True,
                        help="Sample name")
    parser.add_argument("--output-dir", required=True,
                        help="QC report output directory")
    args = parser.parse_args()

    base = os.path.join(args.outdir, args.sample)
    log_dir = os.path.join(args.outdir, "logs", args.sample)
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Generating QC report for {args.sample}")
    print(f"  Pipeline output: {base}")
    print(f"  Logs: {log_dir}")

    # ---- Parse all logs ----
    fastp = parse_fastp(os.path.join(base, "qc", f"{args.sample}.fastp.json"))
    bowtie2 = parse_bowtie2(os.path.join(log_dir, "bowtie2.log"))
    ribodetector = parse_ribodetector(os.path.join(log_dir, "ribodetector.log"))
    star_log = parse_star_log(os.path.join(base, "star", f"{args.sample}.Log.final.out"))
    solo = parse_solo_summary(
        os.path.join(base, "star", f"{args.sample}.Solo.out", "Gene", "Summary.csv")
    )

    # ---- Build read attrition data ----
    stages = []
    counts = []

    if fastp:
        stages.append("Raw reads")
        counts.append(fastp["raw_reads"])
        stages.append("After fastp\n(trim + 26-34 nt)")
        counts.append(fastp["after_trim"])

    if bowtie2 and "after_rrna" in bowtie2:
        stages.append("After bowtie2\n(rRNA removed)")
        counts.append(bowtie2["after_rrna"])

    if ribodetector and "after_ribodetector" in ribodetector:
        stages.append("After ribodetector\n(2nd rRNA filter)")
        counts.append(ribodetector["after_ribodetector"])

    if star_log:
        unique = int(star_log.get("Uniquely mapped reads number", "0"))
        multi = int(star_log.get("Number of reads mapped to multiple loci", "0"))
        stages.append("STAR mapped\n(unique + multi)")
        counts.append(unique + multi)

    if solo:
        stages.append(f"Cells called\n({solo.get('Estimated Number of Cells', '?')} cells)")
        counts.append(int(solo.get("Unique Reads in Cells Mapped to Gene", "0")))

        stages.append("UMIs in cells\n(deduplicated)")
        counts.append(int(solo.get("UMIs in Cells", "0")))

    if counts:
        print("\n  Read attrition:")
        for stage, count in zip(stages, counts):
            pct = count / counts[0] * 100 if counts[0] > 0 else 0
            print(f"    {stage.replace(chr(10), ' '):40s} {count:>12,} ({pct:.1f}%)")

        plot_read_attrition(
            stages, counts,
            os.path.join(args.output_dir, "read_attrition.png")
        )
        print("  Saved: read_attrition.png")

    # ---- Mapping breakdown ----
    if star_log:
        plot_mapping_breakdown(
            star_log,
            os.path.join(args.output_dir, "mapping_breakdown.png")
        )
        print("  Saved: mapping_breakdown.png")

    # ---- UMI knee plot + cell stats ----
    solo_base = os.path.join(base, "star", f"{args.sample}.Solo.out", "Gene")
    if os.path.isdir(solo_base):
        plot_umi_knee(
            solo_base,
            os.path.join(args.output_dir, "umi_distribution.png")
        )
        print("  Saved: umi_distribution.png")

    print(f"\nQC report complete: {args.output_dir}/")


if __name__ == "__main__":
    main()
