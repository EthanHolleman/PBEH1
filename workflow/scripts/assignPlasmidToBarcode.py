# This script merges a read's barcode and plasmid assignment so they
# are easily available in one file. This makes plotting R-loop frequency
# per plasmid and sample easier because read counts per sample per plasmid
# can be easily calculated from this merged tsv file


import pandas as pd


def read_plasmid_tsv(filepath):
    df = pd.read_csv(filepath, sep='\t')
    df.columns = ['plasmid', 'read', 'length']
    return df


flow_cell_barcodes = pd.read_csv(snakemake.input['flow_cell_barcodes'])
flow_cell_plasmids =read_plasmid_tsv(snakemake.input['flow_cell_plasmids'])

merge_df = flow_cell_plasmids.merge(flow_cell_barcodes, on='read')

pd.write_csv(merge_df, snakemake.output[0])

