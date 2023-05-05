from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

fasta_dir = Path(snakemake.params['genome_dir'])
Path(snakemake.output['shifted_dir']).mkdir(exist_ok=True, parents=True)


def move_start_position_downstream(sequence):
    # remove any whitespace from the sequence
    sequence = sequence.replace('\n', '')

    # calculate the length of the sequence
    seq_length = len(sequence)

    # calculate the new start position
    new_start = (seq_length - 1000) % seq_length

    # construct the modified sequence
    new_sequence = sequence[new_start:] + sequence[:new_start]

    return new_sequence


def generate_bed6_entry(record):

    return [
        record.id,
        0,
        len(record.seq),
        record.id,
        0,
        '+'
    ]


def get_shifted_output_path(output_dir, ref_path):
    return Path(output_dir).joinpath(f'{ref_path.stem}.shift.fasta')

shifted_records = []
bed_entries = []

for each_ref in fasta_dir.iterdir():
    if each_ref.suffix == '.fasta':
        record = SeqIO.read(each_ref, 'fasta')
        shifted_seq = move_start_position_downstream(record.seq)
        record.seq = shifted_seq
        shifted_records.append(record)
        bed_entries.append(generate_bed6_entry(record))
        
        shifted_record_path = get_shifted_output_path(
            snakemake.output['shifted_dir'], each_ref
        )
        SeqIO.write([record], shifted_record_path, 'fasta')


SeqIO.write(shifted_records, snakemake.output['genome'], 'fasta')
pd.DataFrame(bed_entries).to_csv(
    snakemake.output['bed6'], index=False, header=None, sep='\t'
    )


        
        