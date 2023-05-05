# splits one fastq with many reads into many fastq with
# few reads
from Bio import SeqIO
from pathlib import Path
import os


def determine_number_reads(input_fastq):
    # determine the number of reads in fastq file without reading into memory
    reads = 0
    records = SeqIO.parse(input_fastq, 'fastq')
    for _ in records:
        reads += 1
    return reads


input_fastq = snakemake.input[0]

num_reads = determine_number_reads(input_fastq)
num_outputs = snakemake.params['num_outputs']

reads = SeqIO.parse(input_fastq, 'fastq')
buffer = []

# determine number of reads per file
reads_per_file = num_reads / num_outputs

if reads_per_file < 1:
    # more files than reads
    reads_per_file = 1

current_file = 0


while current_file < num_outputs:
    try:
        buffer.append(next(reads))
        if len(buffer) >= reads_per_file:
            output_path = snakemake.input[current_file]
            SeqIO.write(buffer, output_path, 'fastq')
            buffer = []
            current_file += 1
    except StopIteration:
        break