from Bio import SeqIO
import math
import os
from pathlib import Path
import math


def determine_number_reads(input_fastq):
    # determine the number of reads in fastq file without reading into memory
    print('coutning lines')
    with open(input_fastq) as f:
        num_lines = sum(1 for line in f)
    return num_lines / 4

input_fastq = snakemake.input[0]

num_reads = determine_number_reads(input_fastq) - len(snakemake.output[0])

reads = SeqIO.parse(input_fastq, 'fastq')
buffer = []

# determine number of reads per file
reads_per_file = math.floor(num_reads / len(snakemake.output))

if reads_per_file < 1:
    # more files than reads
    reads_per_file = 1

current_file = 0

# write 1 record to each file to ensure all files are make to make snakemake
# happy
for each_file in snakemake.output:
    with open(each_file, 'w') as handle:
        SeqIO.write([next(reads)], handle, 'fastq')


while current_file < len(snakemake.output):
    try:
        buffer.append(next(reads))
        if len(buffer) >= reads_per_file:
            output_path = snakemake.output[current_file]
            SeqIO.write(buffer, output_path, 'fastq')
            buffer = []
            current_file += 1
            print(current_file)
    except StopIteration:
        break
    
# try to keep the iterator going so we don't miss any records
try:
    with open(snakemake.output[-1], 'a+') as handle:
        for r in reads:
            SeqIO.write([r], handle, 'fastq')
except StopIteration:
    print('Done')
    