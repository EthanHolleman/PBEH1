import pandas as pd
from pathlib import Path


def read_blast_output(filepath):
    blast_headers = [
        'Query ID',
        'Subject ID',
        'Percent Identity',
        'Alignment Length',
        'Mismatches',
        'Gap Openings',
        'Query Start',
        'Query End',
        'Subject Start',
        'Subject End',
        'E-value',
        'Bit Score'
    ]
    df = pd.read_csv(filepath, header=None, sep='\t')
    df.columns = blast_headers
    return df


def get_alignments(blast_df):
    subject_ids = set(blast_df['Subject ID'])

    align_dict = {}
    for each_id in subject_ids:
        blast_df.loc[blast_df['Subject ID'] == each_id]
        align = blast_df.loc[
            (blast_df['Subject ID'] ==each_id) & (blast_df['Percent Identity'] >= 90)
            & (blast_df['Alignment Length'] >= 50)
            ]
        if len(align):
            align_dict[each_id] = list(align['Query ID'])
    
    return align_dict


def identify_plasmids(id_df, align_dict):
    matches = []
    for each_read in align_dict:
        mask = (id_df[align_dict[each_read]] == 1)  # get rows where same columns specified == 1
        matching_rows = id_df[mask.all(axis=1)]  # use this to get those rows
        row_sum_match = matching_rows.loc[:, matching_rows.columns != 'ref_name']
        row_sum = row_sum_match.sum(axis=1)  # calculate row sum
        sum_match = matching_rows[row_sum == len(align_dict[each_read])]
        # only return rows where the row sum = the number of columns that way plasmids that require
        # 2 alignments do not show up when a read from a plasmid that requires a single alignment
        # appears
        try:
            reference = list(sum_match['ref_name'])[0]
        except IndexError:
            reference = None
        matches.append({
            'ref_name': reference,
            'read_name': each_read
        })

    return pd.DataFrame(matches)



def add_ref_rel_path(ref_dir, plasmid_ids):
    plasmid_ids['ref_path'] = plasmid_ids.apply(
        lambda row: str(Path(ref_dir).joinpath(str(row['ref_path']))), axis=1
    )
    return plasmid_ids



def main():

    blast_tsv = snakemake.input['blast_tsv']
    ref_dir = snakemake.params['ref_dir']
    id_df_path = snakemake.params['id_df']

    id_df = pd.read_csv(id_df_path, sep='\t')
    print(id_df.columns)
    print('==='*20)
    blast_df = read_blast_output(blast_tsv)
    alignments = get_alignments(blast_df)
    plasmid_ids = identify_plasmids(id_df, alignments)
    print('='*20)
    print(plasmid_ids)
    #plasmid_ids_rel = add_ref_rel_path(ref_dir, plasmid_ids)

    plasmid_ids.to_csv(
        snakemake.output[0], sep='\t', index=None
    )


if __name__ == '__main__':
    main()

        
