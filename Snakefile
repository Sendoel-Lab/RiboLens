"""
RiboLens - scribo-seq processing pipeline
=========================================
A Snakemake pipeline that takes scribo-seq (single-cell ribosome profiling)
data from BCL or FASTQ through demultiplexing, alignment, deduplication,
optional QC, and a per-sample HTML report.

Usage:
    snakemake -j 16 --configfile my_config.yaml
    snakemake -n --configfile my_config.yaml          # dry run
    snakemake -j 16 --config samples='["sample1"]'    # single sample

    # Run demux only:
    snakemake -j 16 --configfile my_config.yaml -- demux

    # Start from existing FASTQs (set run_demux: false in config):
    snakemake -j 16 --configfile my_config.yaml
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(workflow.snakefile), "scripts"))
from helpers import (
    find_sample_fastq,
    get_fastq_dir,
    get_filtered_fastq as _get_filtered_fastq,
    get_final_outputs,
    get_samples,
)

configfile: "configs/config.yaml"


SAMPLES = get_samples(config)
OUTDIR = config["output_dir"]
FASTQ_DIR = get_fastq_dir(config)


# Snakemake input functions (need wildcards, so wrap the pure helpers).

def get_filtered_fastq(wildcards):
    return _get_filtered_fastq(config, OUTDIR, wildcards.sample)




# =============================================================================
# Validation
# =============================================================================

if config.get("run_featurecounts", False):
    assert config.get("run_bam_split", False), \
        "Config error: run_featurecounts requires run_bam_split to be true"


# =============================================================================
# Wildcard constraints
# =============================================================================

wildcard_constraints:
    sample = "[A-Za-z0-9_.-]+"


# =============================================================================
# Target rule
# =============================================================================

rule all:
    input:
        get_final_outputs(config, OUTDIR, SAMPLES)


# =============================================================================
# Step 0: Demultiplexing (bcl2fastq) -- conditional on run_demux
# =============================================================================

def get_demux_output_fastqs():
    """Return expected FASTQ paths after demux, one per sample."""
    demux_dir = os.path.join(OUTDIR, "demux")
    return [f"{demux_dir}/{s}_R1_001.fastq.gz" for s in SAMPLES]


if config.get("run_demux", False):
    checkpoint demux:
        """Demultiplex BCL files into per-sample FASTQs using bcl2fastq.
        Uses a Snakemake checkpoint so the DAG is re-evaluated after demux
        completes and the actual FASTQ filenames are known.
        """
        output:
            outdir = directory(os.path.join(OUTDIR, "demux")),
        log:
            os.path.join(OUTDIR, "logs", "demux.log")
        benchmark:
            os.path.join(OUTDIR, "benchmarks", "demux.tsv")
        threads: config.get("threads", 16)
        params:
            bcl2fastq    = config.get("bcl2fastq_path", "bcl2fastq"),
            run_dir      = config["bcl_run_dir"],
            sample_sheet = config["sample_sheet"],
            mismatches   = config.get("barcode_mismatches", 1),
            compression  = config.get("fastq_compression", 9),
            extra        = config.get("bcl2fastq_extra_opts", "--ignore-missing-bcls --minimum-trimmed-read-length 5 --mask-short-adapter-reads 5"),
            no_lane_split = "--no-lane-splitting" if config.get("no_lane_splitting", False) else "",
        shell:
            """
            {params.bcl2fastq} -p {threads} \
                -R {params.run_dir} \
                -o {output.outdir} \
                --sample-sheet {params.sample_sheet} \
                --fastq-compression-level {params.compression} \
                --barcode-mismatches {params.mismatches} \
                {params.extra} \
                {params.no_lane_split} \
                2>&1 | tee {log}
            """


def get_sample_fastq(wildcards):
    """Resolve the per-sample FASTQ (waits for demux checkpoint if enabled)."""
    if config.get("run_demux", False):
        fq_dir = checkpoints.demux.get().output.outdir
    else:
        fq_dir = FASTQ_DIR
    return find_sample_fastq(fq_dir, wildcards.sample)


# =============================================================================
# Step 1: Adapter trimming + UMI extraction (fastp)
# =============================================================================

rule fastp:
    input:
        fastq = get_sample_fastq,
    output:
        fastq = f"{OUTDIR}/{{sample}}/trim/{{sample}}.fastp.fastq.gz",
        html  = f"{OUTDIR}/{{sample}}/qc/{{sample}}.fastp.html",
        json  = f"{OUTDIR}/{{sample}}/qc/{{sample}}.fastp.json",
    log:
        f"{OUTDIR}/logs/{{sample}}/fastp.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/fastp.tsv"
    threads: config.get("threads_fastp", 4)
    params:
        adapter = config["adapter_seq"],
        min_len = config.get("min_length", 26),
        max_len = config.get("max_length", 34),
        umi_len = config.get("umi_len", 10),
    shell:
        """
        fastp -Q \
            -i {input.fastq} \
            -o {output.fastq} \
            --umi --umi_loc=read1 --umi_len={params.umi_len} \
            -a {params.adapter} \
            --trim_poly_g --trim_poly_x \
            --dont_eval_duplication \
            -l {params.min_len} --length_limit {params.max_len} \
            -h {output.html} -j {output.json} \
            --thread {threads} \
            2> {log}
        """


# =============================================================================
# Step 2: rRNA removal (bowtie2)
# =============================================================================

rule bowtie2_rrna:
    input:
        fastq = rules.fastp.output.fastq
    output:
        fastq = f"{OUTDIR}/{{sample}}/rrna/{{sample}}.no_rrna.fastq.gz",
    log:
        f"{OUTDIR}/logs/{{sample}}/bowtie2.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/bowtie2.tsv"
    threads: config.get("threads_bowtie2", 16)
    params:
        index = config["rrna_index"],
    shell:
        """
        bowtie2 -p {threads} -D 20 -R 3 -N 1 -L 15 \
            -U {input.fastq} \
            -x {params.index} \
            --un-gz {output.fastq} \
            2> {log} > /dev/null
        """


# =============================================================================
# Step 3: Ribodetector (optional, conditional on run_ribodetector)
# =============================================================================

rule ribodetector:
    input:
        fastq = rules.bowtie2_rrna.output.fastq
    output:
        fastq = f"{OUTDIR}/{{sample}}/ribodetector/{{sample}}.rd.fastq.gz",
    log:
        f"{OUTDIR}/logs/{{sample}}/ribodetector.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/ribodetector.tsv"
    threads: config.get("threads_ribodetector", 12)
    params:
        ribodetector = config.get("ribodetector_path", "ribodetector_cpu"),
        read_len = config.get("ribodetector_len", 29),
        chunk_size = config.get("ribodetector_chunk_size", 6000),
    shell:
        """
        {params.ribodetector} -t {threads} \
            --chunk_size {params.chunk_size} \
            -l {params.read_len} \
            -i {input.fastq} \
            -o {output.fastq} \
            --log {log}
        """


# =============================================================================
# Step 4: Generate synthetic read2 with CB + UMI (cb_umi.py)
# =============================================================================

rule cb_umi:
    input:
        fastq = get_filtered_fastq
    output:
        fastq = f"{OUTDIR}/{{sample}}/cbumi/{{sample}}.cb_umi.fastq.gz",
    log:
        f"{OUTDIR}/logs/{{sample}}/cb_umi.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/cb_umi.tsv"
    params:
        script = os.path.join(config["scripts_dir"], "cb_umi.py"),
    shell:
        """
        ( zcat {input.fastq} | python {params.script} | pigz > {output.fastq} ) 2> {log}
        """


# =============================================================================
# Step 5: STAR alignment with Solo barcode mode
# =============================================================================

rule star_align:
    input:
        cdna  = get_filtered_fastq,
        cbumi = rules.cb_umi.output.fastq,
    output:
        genome_bam = f"{OUTDIR}/{{sample}}/star/{{sample}}.Aligned.sortedByCoord.out.bam",
        tx_bam     = f"{OUTDIR}/{{sample}}/star/{{sample}}.Aligned.toTranscriptome.out.bam",
        log_final  = f"{OUTDIR}/{{sample}}/star/{{sample}}.Log.final.out",
        solo_dir   = directory(f"{OUTDIR}/{{sample}}/star/{{sample}}.Solo.out"),
    log:
        f"{OUTDIR}/logs/{{sample}}/star.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/star.tsv"
    threads: config.get("threads_star", 24)
    params:
        index         = config["star_index"],
        whitelist     = config["whitelist"],
        prefix        = lambda wc: f"{OUTDIR}/{wc.sample}/star/{wc.sample}.",
        multimap_nmax = config.get("star_multimap_nmax", 20),
        mismatch_nmax = config.get("star_mismatch_nmax", 1),
        solo_features = config.get("star_solo_features", "Gene Velocyto"),
        align_ends    = config.get("star_align_ends_type", "EndToEnd"),
        solo_multi    = lambda wc: f"--soloMultiMappers {config['star_solo_multi_mappers']}" if config.get("star_solo_multi_mappers") else "",
        solo_umi_dedup = lambda wc: f"--soloUMIdedup {config['star_solo_umi_dedup']}" if config.get("star_solo_umi_dedup") else "",
        solo_cb_match = lambda wc: f"--soloCBmatchWLtype {config['star_solo_cb_match']}" if config.get("star_solo_cb_match") else "",
        extra         = config.get("star_extra_opts", ""),
    shell:
        """
        rm -rf {params.prefix}_STARtmp

        STAR \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.cdna} {input.cbumi} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 31000000000 \
            --soloCBwhitelist {params.whitelist} \
            --soloType CB_UMI_Simple \
            --soloFeatures {params.solo_features} \
            --soloUMIlen 10 \
            --soloUMIstart 11 \
            --soloCBlen 10 \
            --soloCBstart 1 \
            --outFileNamePrefix {params.prefix} \
            --outSAMattributes UR CR AS NH HI nM MD CB \
            --outFilterMultimapNmax {params.multimap_nmax} \
            --quantMode TranscriptomeSAM \
            --alignEndsType {params.align_ends} \
            --outFilterMismatchNmax {params.mismatch_nmax} \
            {params.solo_multi} \
            {params.solo_umi_dedup} \
            {params.solo_cb_match} \
            {params.extra} \
            2>&1 | tee {log}
        """


# =============================================================================
# Step 6a: Filter reads with missing CB tags (sambamba)
# =============================================================================

rule filter_cb:
    """Remove reads where STAR could not assign a cell barcode (CB='-')."""
    input:
        bam = rules.star_align.output.genome_bam
    output:
        bam = f"{OUTDIR}/{{sample}}/genome_bam/{{sample}}.filtered.bam",
        bai = f"{OUTDIR}/{{sample}}/genome_bam/{{sample}}.filtered.bam.bai",
    log:
        f"{OUTDIR}/logs/{{sample}}/sambamba.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/filter_cb.tsv"
    threads: config.get("threads", 4)
    shell:
        """
        sambamba view -F "[CB] != '-' " -f bam {input.bam} -o {output.bam} 2> {log}
        samtools index {output.bam}
        """


# =============================================================================
# Step 7: UMI deduplication on genome BAM
# =============================================================================

rule dedup_genome:
    input:
        bam = rules.filter_cb.output.bam,
        bai = rules.filter_cb.output.bai,
    output:
        bam = f"{OUTDIR}/{{sample}}/genome_bam/{{sample}}.dedup.bam",
        bai = f"{OUTDIR}/{{sample}}/genome_bam/{{sample}}.dedup.bam.bai",
    log:
        f"{OUTDIR}/logs/{{sample}}/dedup_genome.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/dedup_genome.tsv"
    shell:
        """
        umi_tools dedup \
            -I {input.bam} \
            -S {output.bam} \
            -L {log} \
            --extract-umi-method tag \
            --umi-tag=UR \
            --cell-tag=CB \
            --per-cell
        samtools index {output.bam} 2>> {log}
        """


# =============================================================================
# Step 8a: Sort transcriptome BAM
# =============================================================================

rule sort_tx_bam:
    input:
        bam = rules.star_align.output.tx_bam
    output:
        bam = f"{OUTDIR}/{{sample}}/tx_bam/{{sample}}.tx.sorted.bam",
        bai = f"{OUTDIR}/{{sample}}/tx_bam/{{sample}}.tx.sorted.bam.bai",
    log:
        f"{OUTDIR}/logs/{{sample}}/sort_tx.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/sort_tx.tsv"
    threads: config.get("threads", 4)
    shell:
        """
        samtools sort -@ {threads} {input.bam} -o {output.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


# =============================================================================
# Step 8b: UMI deduplication on transcriptome BAM
# =============================================================================

rule dedup_tx:
    input:
        bam = rules.sort_tx_bam.output.bam,
        bai = rules.sort_tx_bam.output.bai,
        # genome BAM carrying STAR's 1MM-corrected CB (tx BAM only has raw CR)
        cb_source = rules.filter_cb.output.bam,
    output:
        bam = f"{OUTDIR}/{{sample}}/tx_bam/{{sample}}.tx.dedup.bam",
        bai = f"{OUTDIR}/{{sample}}/tx_bam/{{sample}}.tx.dedup.bam.bai",
    log:
        f"{OUTDIR}/logs/{{sample}}/dedup_tx.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/dedup_tx.tsv"
    params:
        cb_bam = lambda wc: f"{OUTDIR}/{wc.sample}/tx_bam/{wc.sample}.tx.cb.bam",
        script = os.path.join(config["scripts_dir"], "transfer_cb_tx.py"),
    shell:
        """
        # Stamp the corrected CB onto the transcriptome reads (by read name) so
        # per-cell dedup/QC group on the whitelist-corrected barcode (capped at
        # the whitelist size) instead of the raw CR (which yields ghost cells).
        conda run -n scribo \
            bash -c '$CONDA_PREFIX/bin/python {params.script} {input.cb_source} {input.bam} {params.cb_bam}' \
            > {log} 2>&1
        samtools index {params.cb_bam} 2>> {log}
        umi_tools dedup \
            -I {params.cb_bam} \
            -S {output.bam} \
            --extract-umi-method tag \
            --umi-tag=UR \
            --cell-tag=CB \
            --per-cell >> {log} 2>&1
        samtools index {output.bam} 2>> {log}
        rm -f {params.cb_bam} {params.cb_bam}.bai
        """


# =============================================================================
# Optional: Split BAM into per-cell BAMs
# =============================================================================

rule split_bam:
    input:
        bam = rules.dedup_genome.output.bam,
        bai = rules.dedup_genome.output.bai,
    output:
        done = touch(f"{OUTDIR}/{{sample}}/bams/.split_done"),
    log:
        f"{OUTDIR}/logs/{{sample}}/split_bam.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/split_bam.tsv"
    params:
        script = os.path.join(config["scripts_dir"], "split_bam_cb.py"),
        prefix = lambda wc: f"{OUTDIR}/{wc.sample}/bams/{wc.sample}",
    shell:
        """
        conda run -n scribo \
            bash -c '$CONDA_PREFIX/bin/python {params.script} {input.bam} CB {params.prefix}' \
            2> {log}
        """


# =============================================================================
# Optional: featureCounts on per-cell BAMs
# =============================================================================

rule featurecounts:
    input:
        done = f"{OUTDIR}/{{sample}}/bams/.split_done",
    output:
        counts = f"{OUTDIR}/{{sample}}/featurecounts/{{sample}}.fc_cds.tsv.gz",
    log:
        f"{OUTDIR}/logs/{{sample}}/featurecounts.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/featurecounts.tsv"
    threads: config.get("threads_featurecounts", 8)
    params:
        gtf          = config["gtf"],
        fc_opts      = config.get("featurecounts_opts", "-M -O --fraction -s 1"),
        feature_type = config.get("featurecounts_feature_type", "CDS"),
        group_by     = config.get("featurecounts_group_by", "gene_name"),
        bam_dir      = lambda wc: f"{OUTDIR}/{wc.sample}/bams",
        tsv          = lambda wc: f"{OUTDIR}/{wc.sample}/featurecounts/{wc.sample}.fc_cds.tsv",
    shell:
        """
        featureCounts {params.fc_opts} \
            -a {params.gtf} \
            -t {params.feature_type} \
            -g {params.group_by} \
            -T {threads} \
            -o {params.tsv} \
            {params.bam_dir}/*.bam \
            2> {log}
        gzip {params.tsv}
        """


# =============================================================================
# Optional: riboWaltz QC (CDS-aware length distribution, P-site, metagene)
# =============================================================================

rule ribowaltz_qc:
    """Run riboWaltz on transcriptome dedup BAM for CDS-level footprint QC.
    Requires the scribo-ribowaltz conda env (see environment_ribowaltz.yaml).
    """
    input:
        bam = rules.dedup_tx.output.bam,
        bai = rules.dedup_tx.output.bai,
    output:
        # qc_plots/ - 9 diagnostic PDFs (Illustrator-compatible vector output)
        qc_length_dist      = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/01_length_dist.pdf",
        qc_length_dist_zoom = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/02_length_dist_zoom.pdf",
        qc_psite_region     = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/03_psite_region.pdf",
        qc_frame            = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/04_frame.pdf",
        qc_frame_by_length  = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/05_frame_by_length.pdf",
        qc_ends_heatmap     = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/06_ends_heatmap.pdf",
        qc_metagene_scaled  = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/07_metagene_scaled.pdf",
        qc_metagene_raw     = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/08_metagene_raw.pdf",
        qc_metagene_log10   = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/09_metagene_log10.pdf",
        # qc_plots/ - matching PNGs (raster copies for HTML report embedding)
        qc_length_dist_png      = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/01_length_dist.png",
        qc_length_dist_zoom_png = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/02_length_dist_zoom.png",
        qc_psite_region_png     = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/03_psite_region.png",
        qc_frame_png            = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/04_frame.png",
        qc_frame_by_length_png  = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/05_frame_by_length.png",
        qc_ends_heatmap_png     = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/06_ends_heatmap.png",
        qc_metagene_scaled_png  = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/07_metagene_scaled.png",
        qc_metagene_raw_png     = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/08_metagene_raw.png",
        qc_metagene_log10_png   = f"{OUTDIR}/{{sample}}/ribowaltz/qc_plots/09_metagene_log10.png",
        # periodicity_plots/ - all-length metagene PDFs
        per_all_scaled      = f"{OUTDIR}/{{sample}}/ribowaltz/periodicity_plots/metagene_all_scaled.pdf",
        per_all_raw         = f"{OUTDIR}/{{sample}}/ribowaltz/periodicity_plots/metagene_all_raw.pdf",
        per_all_log10       = f"{OUTDIR}/{{sample}}/ribowaltz/periodicity_plots/metagene_all_log10.pdf",
        # periodicity_plots/ - all-length metagene PNGs (for HTML report)
        per_all_scaled_png  = f"{OUTDIR}/{{sample}}/ribowaltz/periodicity_plots/metagene_all_scaled.png",
        per_all_raw_png     = f"{OUTDIR}/{{sample}}/ribowaltz/periodicity_plots/metagene_all_raw.png",
        per_all_log10_png   = f"{OUTDIR}/{{sample}}/ribowaltz/periodicity_plots/metagene_all_log10.png",
        # periodicity_plots/ - per-length metagene PDFs (20-34 nt)
        per_by_length = expand(
            f"{OUTDIR}/{{{{sample}}}}/ribowaltz/periodicity_plots/metagene_{{L}}nt.pdf",
            L=range(20, 35),
        ),
        psite_csv  = f"{OUTDIR}/{{sample}}/ribowaltz/psite_offsets.csv",
        length_csv = f"{OUTDIR}/{{sample}}/ribowaltz/length_distribution.csv",
    log:
        f"{OUTDIR}/logs/{{sample}}/ribowaltz.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/ribowaltz.tsv"
    params:
        script     = os.path.join(config["scripts_dir"], "ribowaltz_qc.R"),
        annotation = config["ribowaltz_annotation"],
        bam_dir    = lambda wc: f"{OUTDIR}/{wc.sample}/tx_bam",
        bam_prefix = lambda wc: f"{wc.sample}.tx.dedup",
        output_dir = lambda wc: f"{OUTDIR}/{wc.sample}/ribowaltz",
        cds_only   = "--cds_only" if config.get("ribowaltz_cds_only", True) else "",
        conda_env  = config.get("ribowaltz_conda_env", "scribo-ribowaltz"),
    shell:
        """
        conda run -n {params.conda_env} \
            bash -c '$CONDA_PREFIX/bin/Rscript {params.script} \
            --bam_dir {params.bam_dir} \
            --bam_prefix {params.bam_prefix} \
            --sample {wildcards.sample} \
            --annotation {params.annotation} \
            --output_dir {params.output_dir} \
            {params.cds_only}' \
            2>&1 | tee {log}
        """


# =============================================================================
# Optional: Salmon alevin quantification
# =============================================================================

rule salmon:
    input:
        cdna  = get_filtered_fastq,
        cbumi = rules.cb_umi.output.fastq,
    output:
        done = touch(f"{OUTDIR}/{{sample}}/salmon/.salmon_done"),
    log:
        f"{OUTDIR}/logs/{{sample}}/salmon.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/salmon.tsv"
    threads: config.get("threads_salmon", 8)
    params:
        index     = config.get("salmon_index", ""),
        tgmap     = config.get("salmon_tgmap", ""),
        whitelist = config["whitelist"],
        outdir    = lambda wc: f"{OUTDIR}/{wc.sample}/salmon",
    shell:
        """
        salmon alevin -l ISR -p {threads} \
            -1 {input.cbumi} \
            -2 {input.cdna} \
            -i {params.index} \
            -o {params.outdir} \
            --tgMap {params.tgmap} \
            --bc-geometry 1[1-10] \
            --umi-geometry 1[11-20] \
            --read-geometry 2[1-end] \
            --whitelist {params.whitelist} \
            2> {log}
        """


# =============================================================================
# Visualization / reporting
# =============================================================================
# Three rules that build up a per-sample HTML QC report:
#   pipeline_qc            - read-attrition, mapping breakdown, UMI knee
#   single_cell_periodicity - per-cell P-site heatmaps around start/stop codons
#   report                 - assembles a base64-embedded HTML from the above
# All run in the scribo-analysis conda env.

rule pipeline_qc:
    input:
        fastp_json  = f"{OUTDIR}/{{sample}}/qc/{{sample}}.fastp.json",
        bowtie2_log = f"{OUTDIR}/logs/{{sample}}/bowtie2.log",
        star_log    = f"{OUTDIR}/{{sample}}/star/{{sample}}.Log.final.out",
        solo_dir    = rules.star_align.output.solo_dir,
    output:
        attrition = f"{OUTDIR}/{{sample}}/qc_report/read_attrition.png",
        mapping   = f"{OUTDIR}/{{sample}}/qc_report/mapping_breakdown.png",
        umi       = f"{OUTDIR}/{{sample}}/qc_report/umi_distribution.png",
    log:
        f"{OUTDIR}/logs/{{sample}}/pipeline_qc.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/pipeline_qc.tsv"
    params:
        script    = os.path.join(config["scripts_dir"], "analysis", "pipeline_qc.py"),
        outdir    = OUTDIR,
        report_dir = lambda wc: f"{OUTDIR}/{wc.sample}/qc_report",
        conda_env = config.get("analysis_conda_env", "scribo-analysis"),
    shell:
        """
        conda run -n {params.conda_env} \
            bash -c '$CONDA_PREFIX/bin/python3 {params.script} \
            --outdir {params.outdir} \
            --sample {wildcards.sample} \
            --output-dir {params.report_dir}' \
            2>&1 | tee {log}
        """


rule single_cell_periodicity:
    input:
        bam            = f"{OUTDIR}/{{sample}}/tx_bam/{{sample}}.tx.dedup.bam",
        bai            = f"{OUTDIR}/{{sample}}/tx_bam/{{sample}}.tx.dedup.bam.bai",
        psite_offsets  = f"{OUTDIR}/{{sample}}/ribowaltz/psite_offsets.csv",
    output:
        start = f"{OUTDIR}/{{sample}}/periodicity/periodicity_start.png",
        stop  = f"{OUTDIR}/{{sample}}/periodicity/periodicity_stop.png",
    log:
        f"{OUTDIR}/logs/{{sample}}/single_cell_periodicity.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/single_cell_periodicity.tsv"
    params:
        script     = os.path.join(config["scripts_dir"], "analysis", "single_cell_periodicity.py"),
        annotation = config["ribowaltz_annotation"],
        outdir     = lambda wc: f"{OUTDIR}/{wc.sample}/periodicity",
        min_reads  = config.get("periodicity_min_reads", 50),
        conda_env  = config.get("analysis_conda_env", "scribo-analysis"),
    shell:
        """
        conda run -n {params.conda_env} \
            bash -c '$CONDA_PREFIX/bin/python3 {params.script} \
            --bam {input.bam} \
            --annotation {params.annotation} \
            --psite-offsets {input.psite_offsets} \
            --output-dir {params.outdir} \
            --min-reads {params.min_reads}' \
            2>&1 | tee {log}
        """


rule report:
    input:
        # qc_report
        attrition = rules.pipeline_qc.output.attrition,
        mapping   = rules.pipeline_qc.output.mapping,
        umi       = rules.pipeline_qc.output.umi,
        # periodicity
        per_start = rules.single_cell_periodicity.output.start,
        per_stop  = rules.single_cell_periodicity.output.stop,
        # ribowaltz (used by combined_report for length dist / metagenes)
        ribowaltz_len = f"{OUTDIR}/{{sample}}/ribowaltz/length_distribution.csv",
    output:
        html = f"{OUTDIR}/{{sample}}/{{sample}}.report.html",
    log:
        f"{OUTDIR}/logs/{{sample}}/report.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}/report.tsv"
    params:
        script    = os.path.join(config["scripts_dir"], "analysis", "combined_report.py"),
        outdir    = OUTDIR,
        conda_env = config.get("analysis_conda_env", "scribo-analysis"),
    shell:
        """
        conda run -n {params.conda_env} \
            bash -c '$CONDA_PREFIX/bin/python3 {params.script} \
            --plates {wildcards.sample}={params.outdir}:{wildcards.sample} \
            --output {output.html}' \
            2>&1 | tee {log}
        """
