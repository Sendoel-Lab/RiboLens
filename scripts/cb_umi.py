#!/usr/bin/env python3
"""Build the synthetic CB+UMI read that STARsolo consumes as read 2.

Reads a FASTQ on stdin and writes, for every record, a 20nt read made of the
10nt cell barcode followed by the 10nt UMI, matching the STAR geometry in the
Snakefile (--soloCBstart 1 --soloCBlen 10 --soloUMIstart 11 --soloUMIlen 10):

    @LH00289:304:233HLJLT3:1:1101:8757:1112:CGNCTGGCCG 1:N:0:ATCACG+TAAGACAGCA
                                            ^^^^^^^^^^        ^^^^^^ ^^^^^^^^^^
                                            UMI (fastp)       i7     i5 = CB
    ->
    @LH00289:304:233HLJLT3:1:1101:8757:1112:CGNCTGGCCG 1:N:0:ATCACG+TAAGACAGCA
    TAAGACAGCACGNCTGGCCG
    +
    EEEEEEEEEEEEEEEEEEEE

Both parts come out of the header: fastp (--umi --umi_loc=read1) moves the UMI
from the 5' end of read 1 into the read name, and the cell barcode is the i5
index. Quality is synthetic; STAR's --soloCBmatchWLtype 1MM ignores it.

The header is parsed by field rather than by offset from the end of the line, so
an i7 of any length -- or no i7 at all -- is handled. A header carrying no usable
barcode is an error rather than a guess: that is archive-stripped input (see the
GSE280255 appendix in the README), and taking the last 10 characters regardless
would turn `length=151` into a "barcode" and carry on silently.

Usage:
    zcat input.fastq.gz | python cb_umi.py | pigz > output.cb_umi.fastq.gz
"""

import argparse
import sys

_STRIPPED_HINT = (
    "If this FASTQ came from SRA/GEO, the archive dropped the index field that "
    "carries the cell barcode. Re-dump it with a defline that restores the index:\n"
    "  fastq-dump -Z --defline-seq '@$sn 1:N:0:$sg' --defline-qual '+' <SRR>\n"
    "See the GSE280255 appendix in the README."
)


def parse_header(header, cb_len, umi_len):
    """Return (cb, umi) from an Illumina header carrying a fastp UMI suffix."""
    parts = header.split(None, 1)
    name = parts[0]
    if len(parts) < 2:
        raise ValueError(
            f"no index field (nothing follows the read name): {header!r}\n{_STRIPPED_HINT}"
        )
    index = parts[1].rsplit(":", 1)[-1]
    cb = index.rsplit("+", 1)[-1]
    if len(cb) != cb_len:
        raise ValueError(
            f"expected a {cb_len}nt cell barcode in the index field, got {cb!r} "
            f"({len(cb)}nt): {header!r}\n{_STRIPPED_HINT}"
        )
    umi = name[-umi_len:]
    if len(umi) != umi_len:
        raise ValueError(
            f"expected a {umi_len}nt UMI at the end of the read name, got {umi!r}: "
            f"{header!r}\nIs this fastp output (fastp --umi --umi_loc=read1)?"
        )
    return cb, umi


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--cb-len", type=int, default=10, help="cell barcode length (default: 10)")
    ap.add_argument("--umi-len", type=int, default=10, help="UMI length (default: 10)")
    args = ap.parse_args()

    qual = "E" * (args.cb_len + args.umi_len)
    out = sys.stdout
    n = 0
    for i, line in enumerate(sys.stdin):
        if i % 4 != 0:
            continue
        header = line.rstrip("\n")
        cb, umi = parse_header(header, args.cb_len, args.umi_len)
        out.write(f"{header}\n{cb}{umi}\n+\n{qual}\n")
        n += 1
    print(f"cb_umi: wrote {n} CB+UMI reads", file=sys.stderr)


if __name__ == "__main__":
    main()
