from Bio import SeqIO
from pathlib import Path
import pandas as pd

def genbank_to_fasta(genbank_path):
    record = SeqIO.read(str(genbank_path), 'genbank')
    
def generate_bed6_entry(record):

    return [
        record.id,
        0,
        len(record.seq),
        record.id,
        0,
        '+'
    ]

def main():


    ref_dir = snakemake.params['genome_dir']
    records, bed6 = [], []
    for each_path in Path(ref_dir).iterdir():
        if each_path.suffix == '.gb':
            records.append(
                SeqIO.read(str(each_path), 'genbank')
            )
            records[-1].id = each_path.stem
            records[-1].description = ''
            bed6.append(generate_bed6_entry(records[-1]))
            
    SeqIO.write(records, snakemake.output['genome'], 'fasta')
    pd.DataFrame(bed6).to_csv(
        snakemake.output['bed6'], index=False, header=None, sep='\t'
        )


if __name__ == '__main__':
    main()
