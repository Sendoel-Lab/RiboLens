#!/usr/bin/env Rscript
# =============================================================================
# riboWaltz QC for scribo-seq
# =============================================================================
# Produces CDS-aware footprint QC from a transcriptome-aligned, deduplicated
# BAM. Replicates the workflow from the original scribo-seq R notebooks
# (kerat/2025/ribowaltz.Rmd).
#
# Outputs (all plots are Illustrator-compatible vector PDFs):
#   qc_plots/                  - 9 diagnostic plots (length dist, P-site region,
#                                 frame, heatmap, metagene) as individual PDFs
#   periodicity_plots/         - wide metagene profiles (all lengths + each
#                                 length 20-34 nt individually) as PDFs
#   psite_offsets.csv           - computed P-site offset per read length
#   length_distribution.csv     - read count per fragment length
#   frame_distribution.csv      - per-length frame-0/1/2 P-site counts + in-frame %
#
# Usage:
#   Rscript ribowaltz_qc.R \
#     --bam_dir    results/young1/tx_bam \
#     --bam_prefix young1.tx.dedup \
#     --sample     young1 \
#     --annotation ribowaltz_m25_anno.csv \
#     --output_dir results/young1/ribowaltz
# =============================================================================

suppressPackageStartupMessages({
  library(riboWaltz)
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(optparse)
})

# =============================================================================
# Command-line arguments
# =============================================================================
option_list <- list(
  make_option("--bam_dir", type = "character",
              help = "Directory containing the transcriptome dedup BAM"),
  make_option("--bam_prefix", type = "character",
              help = "BAM filename prefix (without .bam), e.g. 'young1.tx.dedup'"),
  make_option("--sample", type = "character",
              help = "Sample display name for plots"),
  make_option("--annotation", type = "character",
              help = "Path to riboWaltz annotation CSV (transcript, l_tr, l_utr5, l_cds, l_utr3)"),
  make_option("--output_dir", type = "character",
              help = "Output directory for plots and tables"),
  make_option("--cds_only", type = "logical", default = FALSE,
              action = "store_true",
              help = "Only keep transcripts with CDS (l_cds > 0) [default: use all transcripts]")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("bam_dir", "bam_prefix", "sample", "annotation", "output_dir")) {
  if (is.null(opt[[arg]])) stop(paste0("--", arg, " is required"))
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load transcript annotation
# =============================================================================
# The annotation CSV has one row per transcript with columns:
#   transcript  - ENSEMBL transcript ID (must match BAM reference names)
#   l_tr        - total transcript length
#   l_utr5      - 5' UTR length
#   l_cds       - CDS length
#   l_utr3      - 3' UTR length
#
# By default, all transcripts are kept (including non-coding). This ensures
# the length distribution reflects the full transcriptome-mapped library.
# P-site and metagene analyses internally use only transcripts with CDS.
#
# With --cds_only, only protein-coding transcripts (l_cds > 0) are kept.
# This gives a sharper RPF peak but drops reads on non-coding transcripts.
cat("Loading annotation...\n")
an <- fread(opt$annotation)
cat(sprintf("Annotation: %d transcripts (%d with CDS)\n", nrow(an), sum(an$l_cds > 0)))
if (opt$cds_only) {
  an <- subset(an, l_cds > 0)
  cat(sprintf("CDS-only mode: keeping %d transcripts\n", nrow(an)))
}

# =============================================================================
# 2. Load transcriptome BAM
# =============================================================================
# bamtolist() scans bam_dir for BAM files. The name_samples argument is a
# named vector that acts as both a filter and a rename:
#   names  = BAM filename prefix to match (e.g. "young1.tx.dedup")
#   values = display name to use in all downstream plots (e.g. "young1")
#
# Only reads that map to transcripts present in the annotation are kept.
# Reads on the negative strand are also removed (transcriptome BAMs from
# STAR should have all reads on the + strand; antisense mappings are artefacts).
cat(sprintf("Loading BAM from: %s/%s.bam\n", opt$bam_dir, opt$bam_prefix))
sam <- c(opt$sample)
names(sam) <- opt$bam_prefix
reads_list <- bamtolist(bamfolder = opt$bam_dir, annotation = an, name_samples = sam)

# =============================================================================
# 3. Compute P-site offsets
# =============================================================================
# The P-site (peptidyl-tRNA site) is where the ribosome holds the growing
# peptide chain. Its position within each read depends on the fragment length.
#
# psite() computes the offset empirically for each read length by:
#   1. Collecting all reads whose 5' (or 3') end falls near annotated start
#      codons (within a window of +/- flanking nt around the AUG).
#   2. For each read length, finding the distance from the read end to the
#      start codon that maximises the signal (i.e. the sharpest pile-up).
#   3. That distance is the P-site offset for that read length.
#
# extremity="auto" lets riboWaltz pick whichever end (5' or 3') gives a
# cleaner offset for each length. flanking=6 means it searches +/- 6 nt
# around the start codon.
#
# The output is a table: one row per read length, with the computed offset
# and which extremity was used.
cat("Calculating P-site offsets...\n")
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")
write.csv(psite_offset,
          file.path(opt$output_dir, "psite_offsets.csv"),
          row.names = FALSE)

# psite_info() takes each read and adds a column for its P-site position
# in transcript coordinates, using the per-length offsets computed above.
# It also annotates which transcript region (5'UTR, CDS, 3'UTR) the P-site
# falls in, based on the annotation.
reads_psite_list <- psite_info(reads_list, psite_offset)

# riboWaltz stores plots in named list elements as "plot_{sample_name}"
plot_name <- paste0("plot_", opt$sample)

# =============================================================================
# 3b. In-frame read accounting (numeric table)
# =============================================================================
# The frame plots (04_frame.pdf / 05_frame_by_length.pdf) show the reading-frame
# distribution but riboWaltz does not export the underlying numbers. Compute them
# directly from the P-site table so the in-frame fraction is available as data.
#
# Reading frame is the P-site distance from the CDS start codon modulo 3
# (0 = in-frame / correct codon phase). Only reads whose P-site falls inside the
# CDS are counted, since frame is undefined in the UTRs. The output table has one
# row per read length plus a final length=NA row summing across all lengths:
#   length, n_cds, frame0, frame1, frame2, inframe_pct, sample
cat("Computing frame distribution...\n")
dt   <- reads_psite_list[[opt$sample]]
cdsd <- dt[psite_region == "cds" & !is.na(psite_from_start)]
cdsd[, frame := psite_from_start %% 3]
frame_tab <- cdsd[, .(n_cds  = .N,
                      frame0 = sum(frame == 0),
                      frame1 = sum(frame == 1),
                      frame2 = sum(frame == 2)), by = length][order(length)]
overall <- data.table(length = NA_integer_,
                       n_cds  = sum(frame_tab$n_cds),
                       frame0 = sum(frame_tab$frame0),
                       frame1 = sum(frame_tab$frame1),
                       frame2 = sum(frame_tab$frame2))
frame_out <- rbind(frame_tab, overall)
frame_out[, inframe_pct := round(100 * frame0 / n_cds, 2)]
frame_out[, sample := opt$sample]
setcolorder(frame_out, c("sample", "length", "n_cds",
                         "frame0", "frame1", "frame2", "inframe_pct"))
write.csv(frame_out, file.path(opt$output_dir, "frame_distribution.csv"),
          row.names = FALSE)
cat(sprintf("In-frame (CDS P-sites): %.2f%% of %d reads\n",
            overall[, round(100 * frame0 / n_cds, 2)], overall$n_cds))

# =============================================================================
# 4. QC PDF - 9 diagnostic plots
# =============================================================================
qc_dir <- file.path(opt$output_dir, "qc_plots")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Writing QC plots: %s\n", qc_dir))

# Illustrator-compatible vector PDF: base pdf() device keeps text as editable
# text (not outlined) and useDingbats=FALSE draws points as real vector circles
# rather than glyphs from the Dingbats font (which Illustrator mis-renders).
save_pdf <- function(filename, p, w = 10, h = 6) {
  pdf(file.path(qc_dir, filename), width = w, height = h,
      useDingbats = FALSE, useKerning = FALSE)
  print(p)
  dev.off()
  # also emit a raster PNG (same basename) so the HTML report can embed it;
  # the PDF remains the Illustrator-editable vector version.
  png(file.path(qc_dir, sub("\\.pdf$", ".png", filename)),
      width = w, height = h, units = "in", res = 150)
  print(p)
  dev.off()
}

# ---- Plot 1: Read length distribution (all lengths) ----
# Bar chart of read counts per fragment length. For good ribo-seq data you
# expect a sharp peak at 28-30 nt (the canonical ribosome-protected fragment
# size). A broad or multi-modal distribution suggests contamination with
# non-ribosomal RNA fragments.
length_dist <- rlength_distr(reads_list, sample = opt$sample, colour = "dodgerblue4")
save_pdf("01_length_dist.pdf", length_dist[[plot_name]])

# ---- Plot 2: Read length distribution (20-60 nt zoom) ----
length_dist_filt <- rlength_distr(reads_list, sample = opt$sample,
                                  colour = "dodgerblue4", length_range = 20:60)
save_pdf("02_length_dist_zoom.pdf", length_dist_filt[[plot_name]])

# ---- Plot 3: P-site region distribution (UTR5 / CDS / UTR3) ----
psite_region <- region_psite(reads_psite_list, an, sample = opt$sample,
                             colour = c("darkorange", "darkgreen", "dodgerblue4"))
save_pdf("03_psite_region.pdf", psite_region[["plot"]])

# ---- Plot 4: Reading frame (all lengths combined) ----
fp <- frame_psite(reads_psite_list, an, sample = opt$sample, colour = "dodgerblue4")
save_pdf("04_frame.pdf", fp[[plot_name]])

# ---- Plot 5: Reading frame stratified by read length ----
fl <- frame_psite_length(reads_psite_list, an, sample = opt$sample,
                         cl = 100, colour = "dodgerblue4")
save_pdf("05_frame_by_length.pdf", fl[[plot_name]])

# ---- Plot 6: 5' end heatmap ----
ends_heatmap <- rends_heat(reads_list, an, sample = opt$sample, cl = 100,
                           utr5l = 25, cdsl = 50, utr3l = 25,
                           colour = viridis(20))
save_pdf("06_ends_heatmap.pdf", ends_heatmap[[plot_name]])

# ---- Plot 7: Metagene profile - scaled ----
mp <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                        utr5l = 20, cdsl = 50, utr3l = 20, colour = "dodgerblue4")
save_pdf("07_metagene_scaled.pdf",
         mp[[plot_name]] + geom_line(color = "dodgerblue4", linewidth = 1))

# ---- Plot 8: Metagene profile - unscaled (raw counts) ----
mp_uns <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                            utr5l = 20, cdsl = 50, utr3l = 20,
                            colour = "dodgerblue4", scale_factors = "none")
save_pdf("08_metagene_raw.pdf",
         mp_uns[[plot_name]] + geom_line(color = "dodgerblue4", linewidth = 1))

# ---- Plot 9: Metagene profile - log10 of raw counts ----
mp_log <- mp_uns
mp_log[[plot_name]][["data"]][["mean_scaled_count"]] <-
  log10(mp_log[[plot_name]][["data"]][["mean_scaled_count"]])
mp_log[[plot_name]][["labels"]][["y"]] <- "log10 (# P-sites)"
save_pdf("09_metagene_log10.pdf",
         mp_log[[plot_name]] + geom_line(color = "dodgerblue4", linewidth = 1))

# =============================================================================
# 5. Periodicity PDF - extended metagene profiles
# =============================================================================
# These are wider-window metagene profiles (50 nt UTR, 200 nt CDS) to give a
# broader view of translational dynamics. The per-length plots (26-34 nt) let
# you see which specific fragment lengths carry the clearest 3-nt periodicity,
# helping decide which lengths to include in downstream analyses.
per_dir <- file.path(opt$output_dir, "periodicity_plots")
dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Writing periodicity plots: %s\n", per_dir))

save_per_pdf <- function(filename, p, w = 25, h = 6) {
  pdf(file.path(per_dir, filename), width = w, height = h,
      useDingbats = FALSE, useKerning = FALSE)
  print(p)
  dev.off()
  # also emit a raster PNG (same basename) for HTML report embedding.
  png(file.path(per_dir, sub("\\.pdf$", ".png", filename)),
      width = w, height = h, units = "in", res = 150)
  print(p)
  dev.off()
}

# All lengths - scaled (normalised by transcript count)
mp_all <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                            utr5l = 50, cdsl = 200, utr3l = 50,
                            colour = "dodgerblue4")
save_per_pdf("metagene_all_scaled.pdf",
             mp_all[[plot_name]] + geom_line(color = "dodgerblue4", linewidth = 1))

# All lengths - unscaled (raw summed P-site counts)
mp_all_uns <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                                utr5l = 50, cdsl = 200, utr3l = 50,
                                colour = "dodgerblue4", scale_factors = "none")
save_per_pdf("metagene_all_raw.pdf",
             mp_all_uns[[plot_name]] + geom_line(color = "dodgerblue4", linewidth = 1))

# All lengths - log10 of raw counts
mp_all_log <- mp_all_uns
mp_all_log[[plot_name]][["data"]][["mean_scaled_count"]] <-
  log10(mp_all_log[[plot_name]][["data"]][["mean_scaled_count"]])
mp_all_log[[plot_name]][["labels"]][["y"]] <- "log10 (# P-sites)"
save_per_pdf("metagene_all_log10.pdf",
             mp_all_log[[plot_name]] + geom_line(color = "dodgerblue4", linewidth = 1))

# Per-length metagene profiles (20-34 nt).
# Lengths absent from the trimmed reads (e.g. <26 nt when fastp min_length=26)
# would crash metaprofile_psite, so wrap in tryCatch and emit a placeholder PNG
# to satisfy the declared snakemake outputs.
for (i in 20:34) {
  out_path <- file.path(per_dir, sprintf("metagene_%dnt.pdf", i))
  tryCatch({
    mp_i <- metaprofile_psite(reads_psite_list, an, sample = opt$sample,
                              utr5l = 50, cdsl = 200, utr3l = 50,
                              colour = "dodgerblue4", scale_factors = "none",
                              length_range = i)
    mp_i[[plot_name]][["labels"]][["title"]] <- i
    save_per_pdf(sprintf("metagene_%dnt.pdf", i),
                 mp_i[[plot_name]] + geom_line(color = "dodgerblue4", linewidth = 1))
  }, error = function(e) {
    cat(sprintf("  length %d nt: no plot (%s)\n", i, conditionMessage(e)))
    pdf(out_path, width = 25, height = 6, useDingbats = FALSE, useKerning = FALSE)
    plot.new()
    text(0.5, 0.5, sprintf("No reads at %d nt", i), cex = 3)
    dev.off()
  })
}

# =============================================================================
# 6. Save length distribution table
# =============================================================================
write.csv(length_dist[["count_dt"]],
          file.path(opt$output_dir, "length_distribution.csv"),
          row.names = FALSE)

cat("riboWaltz QC complete.\n")
