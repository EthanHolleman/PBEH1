# read primers from sample assignment and convert to a fasta file that
# can be used for blast search against reads

import pandas as pd

samples_df = pd.read_csv(
    snakemake.input['sample_df'], sep='\t'
)

output_fasta = snakemake.output[0]

# Get rows relevant to primer info

df = samples_df[[
    'FWDPrimer', 'REVPrimer', 'FwdPrimerSeq', 'RevPrimerSeq'
]]

records = []
names = set([])

for i, each_row in df.iterrows():
    fwd_name = each_row['FWDPrimer'].strip()
    fwd_seq = each_row['FwdPrimerSeq'].strip()
    rev_name = each_row['REVPrimer'].strip()
    rev_seq = each_row['RevPrimerSeq'].strip()

    if fwd_name not in names:
        records.append((fwd_name, fwd_seq))
        names.add(fwd_name)
    
    if rev_name not in names:
        records.append((rev_name, rev_seq))
        names.add(rev_name)

with open(snakemake.output[0], 'w') as handle:
    for name, seq in records:
        handle.write(f'>{name}\n')
        handle.write(f'{seq}\n')






