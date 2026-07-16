#!/usr/bin/env python3
"""Generate a combined HTML QC report for multiple scribo-seq samples.

All images are base64-embedded so the report is a single shareable file.
Samples are displayed side-by-side for easy comparison.

Usage:
    python combined_report.py \
        --plates plate1=/mnt/ngs2/finn/results_plate1:414351_1-Plate_1_S1 \
                 plate2=/mnt/ngs2/finn/results_plate2:414351_2-Plate_2_S2 \
        --output /mnt/ngs2/finn/finn_qc_report.html
"""

import argparse
import base64
import json
import os
import re
from datetime import datetime
from pathlib import Path


def img_to_base64(path):
    if not os.path.exists(path):
        return None
    ext = Path(path).suffix.lstrip(".")
    if ext == "jpg":
        ext = "jpeg"
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode()
    return f"data:image/{ext};base64,{data}"


def img_tag(path, width="100%"):
    uri = img_to_base64(path)
    if uri:
        return f'<img src="{uri}" style="max-width:{width}; width:100%;">'
    return '<p class="missing">Not available</p>'


def parse_solo_summary(path):
    if not os.path.exists(path):
        return {}
    stats = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split(",", 1)
            if len(parts) == 2:
                stats[parts[0]] = parts[1]
    return stats


def parse_fastp_json(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        data = json.load(f)
    s = data.get("summary", {})
    bf = s.get("before_filtering", {})
    af = s.get("after_filtering", {})
    return {
        "raw_reads": bf.get("total_reads", 0),
        "after_trim": af.get("total_reads", 0),
        "q30_rate": af.get("q30_rate", 0),
        "gc_content": af.get("gc_content", 0),
    }


def parse_bowtie2_log(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        text = f.read()
    m = re.search(r"([\d.]+)% overall alignment rate", text)
    return m.group(1) + "%" if m else None


def solo_pct(solo, key):
    try:
        return f"{float(solo[key]) * 100:.1f}%"
    except (KeyError, ValueError):
        return "N/A"


def load_plate(outdir, sample):
    """Load all data for a single plate/sample."""
    base = os.path.join(outdir, sample)
    log_dir = os.path.join(outdir, "logs", sample)

    return {
        "sample": sample,
        "outdir": outdir,
        "base": base,
        "solo": parse_solo_summary(
            os.path.join(base, "star", f"{sample}.Solo.out", "Gene", "Summary.csv")
        ),
        "fastp": parse_fastp_json(os.path.join(base, "qc", f"{sample}.fastp.json")),
        "bowtie2_rate": parse_bowtie2_log(os.path.join(log_dir, "bowtie2.log")),
        "qc_report": os.path.join(base, "qc_report"),
        "ribowaltz": os.path.join(base, "ribowaltz"),
        "periodicity": os.path.join(base, "periodicity"),
    }


def build_html(plates, plate_names):
    n = len(plates)
    col_width = f"{100 // n}%"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>scribo-seq Combined QC Report</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
           color: #333; background: #f8f9fa; line-height: 1.6; }}
    .container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
    header {{ background: linear-gradient(135deg, #2c3e50, #3498db); color: white;
             padding: 30px 40px; border-radius: 8px; margin-bottom: 24px; }}
    header h1 {{ font-size: 1.8em; margin-bottom: 4px; }}
    header p {{ opacity: 0.85; font-size: 0.95em; }}
    .section {{ background: white; border-radius: 8px; padding: 24px;
               margin-bottom: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }}
    .section h2 {{ color: #2c3e50; font-size: 1.3em; margin-bottom: 16px;
                  border-bottom: 2px solid #3498db; padding-bottom: 8px; }}
    .section h3 {{ color: #555; font-size: 1.05em; margin: 12px 0 8px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 12px 0; }}
    th {{ background: #f0f4f8; padding: 10px 12px; text-align: left; font-weight: 600;
         border-bottom: 2px solid #ddd; }}
    td {{ padding: 8px 12px; border-bottom: 1px solid #eee; }}
    td:not(:first-child), th:not(:first-child) {{ text-align: right; }}
    .side-by-side {{ display: grid; grid-template-columns: repeat({n}, 1fr); gap: 16px; }}
    .side-by-side img {{ width: 100%; border-radius: 4px; }}
    .plate-label {{ text-align: center; font-weight: 600; color: #2c3e50;
                   margin-bottom: 8px; font-size: 1.05em; }}
    .metric-cards {{ display: grid; grid-template-columns: repeat({n}, 1fr); gap: 16px; margin: 16px 0; }}
    .plate-card {{ background: #f0f4f8; border-radius: 8px; padding: 16px; }}
    .plate-card h4 {{ color: #2c3e50; margin-bottom: 10px; text-align: center; }}
    .plate-card .stats {{ display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }}
    .stat {{ text-align: center; }}
    .stat .value {{ font-size: 1.4em; font-weight: 700; color: #2c3e50; }}
    .stat .label {{ font-size: 0.8em; color: #666; }}
    .missing {{ color: #999; font-style: italic; text-align: center; padding: 40px; }}
    .full-width img {{ width: 100%; border-radius: 4px; }}
    footer {{ text-align: center; color: #999; font-size: 0.85em; padding: 20px; }}
</style>
</head>
<body>
<div class="container">

<header>
    <h1>scribo-seq Combined QC Report</h1>
    <p>Samples: {', '.join(f'<strong>{name}</strong> ({p["sample"]})' for name, p in zip(plate_names, plates))}
    &nbsp;|&nbsp; Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
</header>
"""

    # --- Summary cards ---
    html += '<div class="section"><h2>Summary</h2><div class="metric-cards">'
    for name, p in zip(plate_names, plates):
        solo = p["solo"]
        html += f"""<div class="plate-card"><h4>{name}</h4><div class="stats">
            <div class="stat"><div class="value">{solo.get("Estimated Number of Cells", "N/A")}</div><div class="label">Cells</div></div>
            <div class="stat"><div class="value">{solo.get("Median UMI per Cell", "N/A")}</div><div class="label">Median UMI/cell</div></div>
            <div class="stat"><div class="value">{solo.get("Median Gene per Cell", "N/A")}</div><div class="label">Median genes/cell</div></div>
            <div class="stat"><div class="value">{solo.get("Total Gene Detected", "N/A")}</div><div class="label">Total genes</div></div>
        </div></div>"""
    html += '</div>'

    # --- Comparison table ---
    html += '<table><tr><th>Metric</th>'
    for name in plate_names:
        html += f'<th>{name}</th>'
    html += '</tr>'

    rows = [
        ("Raw reads", lambda p: f"{p['fastp'].get('raw_reads', 0):,}"),
        ("After trimming", lambda p: f"{p['fastp'].get('after_trim', 0):,}"),
        ("rRNA rate (bowtie2)", lambda p: p["bowtie2_rate"] or "N/A"),
        ("Reads to STAR", lambda p: f"{int(p['solo'].get('Number of Reads', 0)):,}"),
        ("Mapped to genome (unique)", lambda p: solo_pct(p["solo"], "Reads Mapped to Genome: Unique")),
        ("Mapped to gene (unique)", lambda p: solo_pct(p["solo"], "Reads Mapped to Gene: Unique Gene")),
        ("Cells detected", lambda p: p["solo"].get("Estimated Number of Cells", "N/A")),
        ("Mean reads/cell", lambda p: p["solo"].get("Mean Reads per Cell", "N/A")),
        ("Median UMI/cell", lambda p: p["solo"].get("Median UMI per Cell", "N/A")),
        ("Median genes/cell", lambda p: p["solo"].get("Median Gene per Cell", "N/A")),
        ("Sequencing saturation", lambda p: solo_pct(p["solo"], "Sequencing Saturation")),
        ("Fraction reads in cells", lambda p: solo_pct(p["solo"], "Fraction of Unique Reads in Cells")),
    ]
    for label, fn in rows:
        html += f'<tr><td>{label}</td>'
        for p in plates:
            html += f'<td>{fn(p)}</td>'
        html += '</tr>'
    html += '</table></div>'

    # --- Side-by-side plot sections ---
    plot_sections = [
        ("Read Attrition", "qc_report", "read_attrition.png",
         "Read count at each pipeline stage showing where reads are lost."),
        ("STAR Mapping Breakdown", "qc_report", "mapping_breakdown.png",
         "Distribution of reads across STAR mapping categories."),
        ("UMI Distribution", "qc_report", "umi_distribution.png",
         "Knee plot, UMI per cell, and genes per cell from STAR Solo filtered matrix."),
    ]

    for title, subdir, filename, desc in plot_sections:
        html += f'<div class="section"><h2>{title}</h2><p>{desc}</p>'
        html += '<div class="side-by-side">'
        for name, p in zip(plate_names, plates):
            path = os.path.join(p[subdir] if subdir in p else os.path.join(p["base"], subdir), filename)
            html += f'<div><div class="plate-label">{name}</div>{img_tag(path)}</div>'
        html += '</div></div>'

    # --- ribowaltz QC plots ---
    rw_qc_plots = [
        ("01_length_dist.png", "Read Length Distribution",
         "Fragment length distribution. Expect a sharp peak at 28-30 nt for ribosome footprints."),
        ("02_length_dist_zoom.png", "Read Length Distribution (20-60 nt)",
         "Zoomed view of the RPF peak region."),
        ("03_psite_region.png", "P-site Region Distribution",
         "Fraction of P-sites in 5'UTR, CDS, and 3'UTR. CDS should dominate."),
        ("04_frame.png", "Reading Frame",
         "P-site counts per reading frame. Frame 0 enrichment confirms true ribosome footprints."),
        ("05_frame_by_length.png", "Reading Frame by Length",
         "Frame analysis per read length. 28-30 nt should show strongest frame 0 bias."),
        # 5' end heatmap replaced by single-cell periodicity below
        ("07_metagene_scaled.png", "Metagene Profile (scaled)",
         "P-site density around start and stop codons, normalised by transcript count."),
        ("08_metagene_raw.png", "Metagene Profile (raw counts)",
         "Raw summed P-site counts around start and stop codons."),
        ("09_metagene_log10.png", "Metagene Profile (log10)",
         "Log10 of raw P-site counts. Makes 3-nt periodicity visible despite start-codon peak."),
    ]

    html += '<div class="section"><h2>riboWaltz QC</h2>'
    html += '<p>CDS-aware ribosome footprint quality metrics from transcriptome-aligned, deduplicated BAMs.</p>'
    for filename, title, desc in rw_qc_plots:
        html += f'<h3>{title}</h3><p style="color:#666;font-size:0.9em;margin-bottom:8px;">{desc}</p>'
        html += '<div class="side-by-side">'
        for name, p in zip(plate_names, plates):
            path = os.path.join(p["ribowaltz"], "qc_plots", filename)
            html += f'<div><div class="plate-label">{name}</div>{img_tag(path)}</div>'
        html += '</div>'
    html += '</div>'

    # --- Single-cell periodicity heatmaps ---
    html += '<div class="section"><h2>Single-Cell Periodicity</h2>'
    html += '<p>Per-cell P-site density relative to start and stop codons (Dur\u00e9 et al. 2024 style). '
    html += 'Each row is one cell, sorted by read count. Vertical stripes at 3-nt intervals confirm '
    html += 'ribosome footprint periodicity at the single-cell level.</p>'

    for filename, title in [
        ("periodicity_start.png", "Distance from start codon"),
        ("periodicity_stop.png", "Distance from stop codon"),
    ]:
        html += f'<h3>{title}</h3><div class="side-by-side">'
        for name, p in zip(plate_names, plates):
            path = os.path.join(p["periodicity"], filename)
            html += f'<div><div class="plate-label">{name}</div>{img_tag(path)}</div>'
        html += '</div>'
    html += '</div>'

    # --- Periodicity plots ---
    html += '<div class="section"><h2>Periodicity (Metagene by Length)</h2>'
    html += '<p>Wide-window metagene profiles. Per-length plots reveal which fragment sizes carry true 3-nt periodicity.</p>'

    per_plots = [
        ("metagene_all_scaled.png", "All lengths (scaled)"),
        ("metagene_all_raw.png", "All lengths (raw)"),
        ("metagene_all_log10.png", "All lengths (log10)"),
    ] + [(f"metagene_{i}nt.png", f"{i} nt") for i in range(20, 35)]

    for filename, title in per_plots:
        html += f'<h3>{title}</h3><div class="side-by-side">'
        for name, p in zip(plate_names, plates):
            path = os.path.join(p["ribowaltz"], "periodicity_plots", filename)
            html += f'<div><div class="plate-label">{name}</div>{img_tag(path)}</div>'
        html += '</div>'
    html += '</div>'

    html += f"""
<footer>
    scribo-seq pipeline combined report &nbsp;|&nbsp; {datetime.now().strftime('%Y-%m-%d')}
</footer>
</div></body></html>"""

    return html


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plates", nargs="+", required=True,
                        help="name=outdir:sample pairs, e.g. Plate1=/path:sample_name")
    parser.add_argument("--output", required=True, help="Output HTML path")
    args = parser.parse_args()

    plates = []
    plate_names = []
    for spec in args.plates:
        name, rest = spec.split("=", 1)
        outdir, sample = rest.split(":", 1)
        plate_names.append(name)
        plates.append(load_plate(outdir, sample))

    html = build_html(plates, plate_names)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as f:
        f.write(html)
    print(f"Report saved: {args.output}")
    print(f"File size: {os.path.getsize(args.output) / 1024:.0f} KB")


if __name__ == "__main__":
    main()
