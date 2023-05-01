# counts the number of unique reads in fastq file
# does not count fwd / rev reads as seperate reads
# write file name and read count to not header tsv
# file for easy concat


from Bio import SeqIO
import pandas as pd


fastq_path = snakemake.input[0]
records = SeqIO.parse(fastq_path, 'fastq')
flow_cell = snakemake.params['flow_cell']



def parse_read_name(read_name):
    return read_name.split('/')[1].strip()


def write_read_lengths(unique_read_lengths, output_path):
    with open(output_path, 'w') as handle:
        for name, read_len in unique_read_lengths:
            handle.write(
                f'{name}\t{read_len}\n'
            )


read_ids = set([])
read_count_unique = 0
total_reads = 0
unique_read_lengths = []

for r in records:
    r_id = parse_read_name(r.description)
    if r_id not in read_ids:
        read_ids.add(r_id)
        read_count_unique += 1
        unique_read_lengths.append(
            (r.description, len(r.seq))
        )
    total_reads += 1


count_df = pd.DataFrame([
    {
        'filepath': snakemake.input[0],
        'uniqueReadCount': read_count_unique,
        'totalReads': total_reads,
        'flow_cell': flow_cell
    }
])

count_df.to_csv(snakemake.output['read_counts'], sep='\t', header=None, index=None)
write_read_lengths(unique_read_lengths, snakemake.output['read_lengths'])

