import pandas as pd


def read_bowtie_tsv(filepath):
    df = pd.read_csv(filepath, sep='\t', names=[
        'query', 'subject', 'loc', 'seq'
    ])
    return df


def get_alignments(bowtie_tsb):
    subject_ids = set(bowtie_tsb['subject'])

    align_dict = {}
    for each_id in subject_ids:
        align = bowtie_tsb.loc[bowtie_tsb['subject'] == each_id]

        if len(align):
            align_dict[each_id] = list(align['query'])
    
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



def main():

    bowtie_tsv = snakemake.input['bowtie_tsv']
    ref_dir = snakemake.params['ref_dir']
    id_df_path = snakemake.params['id_df']

    id_df = pd.read_csv(id_df_path, sep='\t')

    bowtie_df = read_bowtie_tsv(bowtie_tsv)
    alignments = get_alignments(bowtie_df)
    plasmid_ids = identify_plasmids(id_df, alignments)
    
    #plasmid_ids_rel = add_ref_rel_path(ref_dir, plasmid_ids)

    plasmid_ids.to_csv(
        snakemake.output[0], sep='\t', index=None
    )


if __name__ == '__main__':
    main()