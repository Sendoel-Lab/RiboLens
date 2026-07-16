# Extracts cell barcode (CB, 10nt) and UMI (10nt) from FASTQ reads via stdin/stdout pipe.
# Usage: zcat input.fastq.gz | python cb_umi.py | pigz > output.cb_umi.fastq.gz
# CB is extracted from positions -35:-25, UMI from -11:-1 of the sequence line.
# Output is a synthetic FASTQ with 20nt (CB+UMI) as the read, quality all 'E'.

from sys import stdin

n=0

for line in stdin:
    n+=1
    if n%4==1:
        # 10nt cb + 10nt umi
        print(line + line[-11:-1]+ line[-35:-25] +'\n+\nEEEEEEEEEEEEEEEEEEEE')
