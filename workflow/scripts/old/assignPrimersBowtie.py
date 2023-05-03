import pandas as pd


def read_blast_tsv(filepath):
    df = pd.read_csv(filepath, sep='\t', header=None)
    df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    return df


def find_paired_barcodes(blast_df, sample_table, fwd_primer_names, rev_primer_names, min_length=46):
    # Identify which barcodes were assigned to which reads and produce a table
    # with this information. If a read was not assinged a fwd and or rev barcode
    # the value will be "none". Min length gives the min alignment length
    # to count a primer alignment towards a read

    primer_table = []

    read_names = set(blast_df['sseqid'])
    for each_read in read_names:
        read_alignments = blast_df.loc[(blast_df['sseqid'] == each_read) & (blast_df['length'] >= min_length)]
        aligned_primers = list(read_alignments['qseqid'])
        fwd_primer, rev_primer = 'none', 'none'

        # determine which is forward and which is reverse primer
        for each_primer in aligned_primers:
            if each_primer in fwd_primer_names:
                fwd_primer = each_primer
            if each_primer in rev_primer_names:
                rev_primer = each_primer
        

        # pull out some info about the location of primer alignments if they do align
        if fwd_primer != 'none':
            fwd_alignment = read_alignments.loc[read_alignments['qseqid'] == fwd_primer]

            fwd_alignment_dict = {
                'fwd_alignment_start': list(fwd_alignment['sstart'])[0],
                'fwd_alignment_end': list(fwd_alignment['send'])[0],
                'fwd_alignment_ident': list(fwd_alignment['pident'])[0],
                'fwd_alignment_len': list(fwd_alignment['length'])[0],
                'fwd_alignment_mismatch': list(fwd_alignment['mismatch'])[0]
            }
        else:
            fwd_alignment_dict = {
                'fwd_alignment_start': None,
                'fwd_alignment_end': None,
                'fwd_alignment_ident': None,
                'fwd_alignment_len': None,
                'fwd_alignment_mismatch': None
            }
        
        if rev_primer != 'none':
            rev_alignment = read_alignments.loc[read_alignments['qseqid'] == rev_primer]

            rev_alignment_dict = {
            'rev_alignment_start': list(rev_alignment['sstart'])[0],
            'rev_alignment_end': list(rev_alignment['send'])[0],
            'rev_alignment_ident': list(rev_alignment['pident'])[0],
            'rev_alignment_len': list(rev_alignment['length'])[0],
            'rev_alignment_mismatch': list(rev_alignment['mismatch'])[0]
        }
        else:
            rev_alignment_dict = {
                'rev_alignment_start': None,
                'rev_alignment_end': None,
                'rev_alignment_ident': None,
                'rev_alignment_len': None,
                'rev_alignment_mismatch': None
        }

        # try and find a sample that matches this primer assignment (forward and reverse primer match sample)
        sample = sample_table.loc[
            (sample_table['FWDPrimer'] == fwd_primer) 
            & (sample_table['REVPrimer'] == rev_primer)
            ]
        
        if len(sample) > 0:
            sample = list(sample['Sample-ID'])
            if len(sample) == 1:
                sample = sample[0]
        else:
            sample = -1
        read_entry = {
            'read': each_read,
            'fwd_primer': fwd_primer,
            'rev_primer': rev_primer,
            'sample_id': sample

        }
        read_entry.update(fwd_alignment_dict)
        read_entry.update(rev_alignment_dict)

        primer_table.append(read_entry)

    return pd.DataFrame(primer_table)


blast_tsv = snakemake.input['blast']
sample_table = snakemake.input['sample_table']

blast_df = read_blast_tsv(blast_tsv)
sample_df = pd.read_csv(sample_table, sep='\t')

fwd_primers = set(sample_df['FWDPrimer'])
rev_primers = set(sample_df['REVPrimer'])


# do the thing

primer_table = find_paired_barcodes(
    blast_df, sample_df, fwd_primers, rev_primers
)

primer_table.to_csv(
    snakemake.output[0], sep='\t', index=None
)