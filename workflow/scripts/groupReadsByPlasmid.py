import pandas as pd
from Bio import SeqIO
from pathlib import Path



def assign_groups(plasmid_id_df):
    group_dict = {}
    for i, row, in plasmid_id_df.iterrows():
        group_dict[row.read_name] = row.ref_name
    return group_dict


def generate_output_path(output_dir, group):
    num = Path(output_dir).parts[-1]
    if isinstance(group, str):
        group = group.replace('-', '_')
    print(group)
    return Path(output_dir).joinpath(f'{group}.group.fastq')


def write_to_fastq(record, group, output_dir):
    fastq_name = generate_output_path(output_dir, group)
    with open(fastq_name, 'a+') as handle:
        SeqIO.write(record, handle, 'fastq')
    return fastq_name


def initialize_output(output_dir, plasmid_id_df):
    for each_plasmid in list(plasmid_id_df['ref_name']):
        fastq_name = generate_output_path(output_dir, each_plasmid)
        with open(fastq_name, 'w') as handle:
            pass
    


def write_reads_to_groups(group_dict, reads_path, output_dir, plasmid_id_df):

    #initialize_output(output_dir, plasmid_id_df)
    locations = []
    records = SeqIO.parse(reads_path, 'fastq')
    for each_record in records:
        try:
            group = group_dict[each_record.description]
        except KeyError:
            group = 'None'
        loc = write_to_fastq(each_record, group, output_dir)
        locations.append(loc)
    return locations


def write_location_sheet(locations, output_path):
    df = pd.DataFrame(locations)
    df.columns = ['paths']
    df.to_csv(output_path, sep='\t', index=None)



def main():

    output_dir = snakemake.output[0]
    table_path = Path(snakemake.output[0]).joinpath('groups.tsv')
    if not Path(output_dir).is_dir():
        Path(output_dir).mkdir()

    plasmid_id_path = snakemake.input['plasmid_id_df']
    reads = snakemake.input['reads']

    plasmid_df = pd.read_csv(plasmid_id_path, sep='\t')
    groups = assign_groups(plasmid_df)
    locations = write_reads_to_groups(
        groups, reads, output_dir, plasmid_df)
    write_location_sheet(locations, table_path)


if __name__ == '__main__':
    main()

    