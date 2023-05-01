from Bio import SeqIO
import pandas as pd

id_plasmid = pd.read_csv(
    snakemake.input['plasmid_id'], sep='\t'
)

records = SeqIO.parse(snakemake.input['reads'], 'fastq')

def get_read_ids(records):
    read_ids = set([])
    for each_read in records:
        read_ids.add(each_read.description)
    return read_ids


def get_un_id_records(id_df, all_record_ids):
    id_set = set(id_df['read_name'])
    all_set = set(all_record_ids)
    
    return all_set - id_set



def write_un_id_reads(un_id_reads, records_path, un_id_path):

    with open(un_id_path, 'a+') as handle:
        for each_read in SeqIO.parse(records_path, 'fastq'):
            if each_read.description in un_id_reads:
                SeqIO.write([each_read], handle, 'fastq')


read_ids = get_read_ids(records)
un_id_reads = get_un_id_records(id_plasmid, read_ids)

# write reads with no assinged plasmid to fastq file 
write_un_id_reads(
    un_id_reads, snakemake.input['reads'], 
    snakemake.output[0]
    )


