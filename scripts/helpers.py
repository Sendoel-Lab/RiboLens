"""Pure Python helpers used by the scribo-seq Snakefile.

Kept separate so the Snakefile stays small and the logic is unit-testable
without importing Snakemake.
"""

import glob
import os

_TRUE = {"true", "yes", "on", "1"}
_FALSE = {"false", "no", "off", "0", ""}


def cfg_bool(config, key, default=False):
    """Read a boolean config value, tolerating strings.

    A configfile gives real booleans, but Snakemake's --config passes values
    through as strings, so `--config run_report=false` reads back as the truthy
    string "false" and silently means the opposite of what was typed. Anything
    unrecognised raises rather than being guessed at.
    """
    value = config.get(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        v = value.strip().lower()
        if v in _TRUE:
            return True
        if v in _FALSE:
            return False
        raise ValueError(
            f"config '{key}': cannot read {value!r} as true/false"
        )
    return bool(value)


def parse_sample_sheet(sample_sheet):
    """Parse a bcl2fastq sample sheet CSV and return unique Sample_IDs."""
    samples = set()
    in_data = False
    with open(sample_sheet) as f:
        for line in f:
            line = line.strip()
            if line.startswith("[Data]"):
                in_data = True
                continue
            if (in_data and line
                    and not line.startswith("Lane,")
                    and not line.startswith(",Sample_ID")
                    and not line.startswith("Sample_ID")):
                parts = line.split(",")
                if len(parts) >= 2:
                    sample_id = (parts[1].strip()
                                 if parts[0].strip().isdigit() or parts[0].strip() == ""
                                 else parts[0].strip())
                    if sample_id:
                        samples.add(sample_id)
    return sorted(samples)


def get_fastq_dir(config):
    """Return the directory containing FASTQ files (demux output or user-supplied)."""
    if cfg_bool(config, "run_demux", False):
        return os.path.join(config["output_dir"], "demux")
    return config["fastq_dir"]


def get_samples(config):
    """Return sample list from config, sample sheet, or auto-discover from FASTQ dir."""
    if config.get("samples"):
        return config["samples"]
    if cfg_bool(config, "run_demux", False):
        ss = config.get("sample_sheet", "")
        if ss and os.path.exists(ss):
            return parse_sample_sheet(ss)
        raise ValueError("run_demux is true but sample_sheet not found")
    fq_dir = config["fastq_dir"]
    if not fq_dir:
        raise ValueError("config 'fastq_dir' is required (or provide 'samples' list)")
    fqs = sorted(glob.glob(os.path.join(fq_dir, "*.fastq.gz")))
    if not fqs:
        raise ValueError(f"No *.fastq.gz files found in {fq_dir}")
    return [os.path.basename(f).replace(".fastq.gz", "") for f in fqs]


def get_filtered_fastq(config, outdir, sample):
    """Return path to rRNA-filtered FASTQ (ribodetector output or bowtie2 output)."""
    if cfg_bool(config, "run_ribodetector", False):
        return f"{outdir}/{sample}/ribodetector/{sample}.rd.fastq.gz"
    return f"{outdir}/{sample}/rrna/{sample}.no_rrna.fastq.gz"


def find_sample_fastq(fastq_dir, sample):
    """Find the FASTQ for a sample, handling bcl2fastq naming conventions.

    bcl2fastq produces {Sample_ID}_S{N}_L{Lane}_R1_001.fastq.gz, or
    {Sample_ID}_S{N}_R1_001.fastq.gz with --no-lane-splitting. We also
    support plain {sample}.fastq.gz for pre-demuxed data.
    Returns the resolved path, or a default that Snakemake will flag as missing.
    """
    direct = [
        f"{fastq_dir}/{sample}.fastq.gz",
        f"{fastq_dir}/{sample}_R1_001.fastq.gz",
    ]
    for p in direct:
        if os.path.exists(p):
            return p
    for pattern in [
        f"{fastq_dir}/{sample}_S*_R1_001.fastq.gz",
        f"{fastq_dir}/{sample}_S*_L*_R1_001.fastq.gz",
    ]:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches[0]
    return f"{fastq_dir}/{sample}.fastq.gz"


def get_final_outputs(config, outdir, samples):
    """Build the final target file list from config flags."""
    targets = []
    for s in samples:
        targets.extend([
            f"{outdir}/{s}/genome_bam/{s}.dedup.bam",
            f"{outdir}/{s}/genome_bam/{s}.dedup.bam.bai",
            f"{outdir}/{s}/tx_bam/{s}.tx.dedup.bam",
            f"{outdir}/{s}/tx_bam/{s}.tx.dedup.bam.bai",
            f"{outdir}/{s}/qc/{s}.fastp.json",
        ])
        if cfg_bool(config, "run_bam_split", False):
            targets.append(f"{outdir}/{s}/bams/.split_done")
        if cfg_bool(config, "run_featurecounts", False):
            targets.append(f"{outdir}/{s}/featurecounts/{s}.fc_cds.tsv.gz")
        if cfg_bool(config, "run_salmon", False):
            targets.append(f"{outdir}/{s}/salmon/.salmon_done")
        if cfg_bool(config, "run_ribowaltz", False):
            targets.append(f"{outdir}/{s}/ribowaltz/length_distribution.csv")
        if cfg_bool(config, "run_report", False):
            targets.append(f"{outdir}/{s}/{s}.report.html")
    return targets
