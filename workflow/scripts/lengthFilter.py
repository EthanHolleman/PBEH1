# Read fastq file, write reads that are in specific length range to the
# specified output file
from Bio import SeqIO

input_fastq = snakemake.input[0]
output_fastq = snakemake.output[0]

lower_bound = snakemake.params['lower']
upper_bound = snakemake.params['upper']

with open(output_fastq, 'a+') as output_handle:
    records = SeqIO.parse(input_fastq, 'fastq')
    for each_record in records:
        seq_len = len(each_record.seq)
        if seq_len > lower_bound and seq_len < upper_bound:
            SeqIO.write([each_record], output_handle, 'fastq')


