import pandas as pd


plasmid_df = pd.read_csv(snakemake.input['plasmid_id'], sep='\t')

length_df = pd.read_csv(snakemake.input['lengths'], sep='\t')
length_df.columns = ['read_name', 'length']

merge = plasmid_df.merge(length_df, on='read_name')
merge.to_csv(snakemake.output[0], sep='\t', index=None, header=None)