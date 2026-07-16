#!/usr/bin/env python3
"""Transfer STAR's whitelist-corrected CB tag onto the transcriptome BAM.

STAR writes the corrected cell barcode (CB) only to the genome alignments, not to
the --quantMode TranscriptomeSAM output, which carries only the raw CR. Per-cell
analyses run on the transcriptome BAM therefore group on the uncorrected CR and
report more than `whitelist size` cells, because 1-2 bp barcode sequencing errors
appear as extra "ghost cells".

This script builds a read-name -> CB map from the genome BAM (where CB is the
1MM-corrected whitelist barcode) and stamps that CB onto every read of the
transcriptome BAM. Transcriptome reads whose name has no valid CB in the genome
BAM (CB == '-' or absent) are dropped, so downstream tx-level dedup/QC can group
on the corrected CB and is capped at the whitelist size.

Usage:
    transfer_cb_tx.py <genome_bam_with_CB> <tx_bam_in> <tx_bam_out>
"""
import argparse
import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("genome_bam", help="genome BAM carrying corrected CB tags")
    ap.add_argument("tx_bam", help="input transcriptome BAM (raw CR only)")
    ap.add_argument("out_bam", help="output transcriptome BAM with CB tags")
    args = ap.parse_args()

    # read-name -> corrected CB (multimapper lines share name+CB, so last wins == same)
    cb_by_read = {}
    with pysam.AlignmentFile(args.genome_bam, "rb") as gb:
        for r in gb:
            if not r.has_tag("CB"):
                continue
            cb = r.get_tag("CB")
            if cb and cb != "-":
                cb_by_read[r.query_name] = cb

    kept = dropped = 0
    with pysam.AlignmentFile(args.tx_bam, "rb") as tin, \
            pysam.AlignmentFile(args.out_bam, "wb", template=tin) as tout:
        for r in tin:
            cb = cb_by_read.get(r.query_name)
            if cb is None:
                dropped += 1
                continue
            r.set_tag("CB", cb, value_type="Z")
            tout.write(r)
            kept += 1

    print(f"transfer_cb_tx: kept={kept} dropped_no_CB={dropped} "
          f"distinct_cells={len(set(cb_by_read.values()))}")


if __name__ == "__main__":
    main()
