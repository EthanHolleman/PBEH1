from Bio import SeqIO


import pandas as pd



def hamming_distance(string1, string2, max_distance):
    if len(string1) != len(string2):
        return -1  # error, strings are of different length
    distance = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            distance += 1
            if distance > max_distance:
                return distance
    return distance


def windowed_hamming(big_string, small_string, max_dist):
    for i in range(len(big_string)-len(small_string)):
        substring = big_string[i:i+len(small_string)]
        dist = hamming_distance(substring, small_string, max_dist)
        if dist < max_dist:
            return i, dist
    
    return -1, dist


def read_blast_tsv(filename):
    column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sstrand']
    blast_df = pd.read_csv(filename, sep='\t', header=None, names=column_names)
    
    # add extra column that tells you what primer type aligned
    blast_df['primer_dir'] = blast_df.apply(
        lambda row: row['qseqid'].split('-')[0],
        axis=1
    )
    
    return blast_df


def get_fwd_primer_alignment(read_name, blast_df):

    read_df = blast_df.loc[(blast_df['sseqid'] == read_name) & 
                           (blast_df['primer_dir'] == 'FWD')]
    if len(read_df) > 0:
        return read_df.iloc[0]
    else:
        return None
    


def get_best_reverse_primer(read_name, blast_df):
    '''Select the "best" reverse primer to use for a given read based on BLAST
    alignment. This needs to be done because while there is only one annealing
    sequence for forward primers there are multible for reverse primers. In
    many cases plasmids have binding sites for all reverse primers which
    can result in a read that was amplified with one reverse primer containing
    the sequence of another because the primer used for amplification was 
    downstream of the other primer.
    
    To deal with this, this function just selects the reverse primer that
    is closest to the end of the read and assumes that was the one that
    was actually used for amplification. 

    Args:
        read_name (str): Read description (fasta)
        blast_df (DataFrame): Pandas dataframe of blast results
    '''
    read_df = blast_df.loc[(blast_df['sseqid'] == read_name) & 
                           (blast_df['primer_dir'] == 'REV')]
    # get only reverse primer alignments
    if len(read_df) > 0:
        
        # alignments should all have same orrientation
        orr = set(read_df['sstrand'])
        assert len(orr) == 1
        orr = orr.pop()
        if orr == 'FWD':
            read_df = df.sort_values(by=['sstart'], ascending=True)
        elif orr == 'REV':
            read_df = df.sort_values(by=['sstart'], ascending=False)
        else:
            raise Exception(f'Orrientation of the read was wack:{orr}')
        
        return read_df.iloc[0]

    else:
        return None


def search_for_barcode(primer_df, read, fwd_barcodes, rev_barcodes):
    # if the strand is minus we want to look for the barcode at the position
    # with the higher value (end of alignment)
    # however if it is a forward orrientation alignment then we want to look
    # at the end with a lower value
    barcode_hits = []
    
    if primer_df != None:
    
        for i, each_primer in primer_df.iterrows():
            if each_primer.primer_dir == 'FWD':
                
                barcodes = fwd_barcodes
                
                if each_primer.sstrand == 'plus':
                    search_start, search_end = 0, each_primer.sstart
                else:
                    search_start, search_end = each_primer.sstart, len(read.seq)
            else:
                
                barcodes = rev_barcodes
                
                if each_primer.sstrand == 'plus':
                    search_start, search_end = each_primer.send, len(read.seq)
                else:
                    search_start, search_end = each_primer.sstart, len(read.seq)
            
            search_seq = read.seq[search_start, search_end]
            print(search_seq)
            
            for each_barcode in barcodes:
                loc, dist = windowed_hamming(search_seq, each_barcode, 1)
                if loc != -1:  # only not negative one when a match occurs
                    barcode_hits.append(
                        {
                            'barcode': each_barcode,
                            'primer': each_primer.primer_dir,
                            'read_name': read.description,
                            'barcode_start': loc,
                            'dist': dist
                        }
                    )
    
    return pd.DataFrame(barcode_hits)


def filter_for_best_barcode(barcode_hits):
    if len(barcode_hits) > 0:
        # determine the min score (best barcode match) value in the df
        min_dist = min(barcode_hits['dist'])
        best_hits = barcode_hits.loc[barcode_hits['dist'] == min_dist]
        if len(best_hits) == 1:
            return best_hits.iloc[0]
    return None
        
                
def concat_series_to_df(series_list):
    # Create empty list to store valid Series
    valid_series = []
    
    # Iterate over series in series_list and append to valid_series if type is Series
    for s in series_list:
        if isinstance(s, pd.Series):
            valid_series.append(s)
    
    # Concatenate valid_series into a DataFrame
    if valid_series:
        return pd.concat(valid_series, axis=1)
    else:
        return None


def get_read_barcodes(read, blast_df, samples_df):
    
    read_name = read.description
    
    fwd_primer = get_fwd_primer_alignment(read_name, blast_df)
    rev_primer = get_best_reverse_primer(read_name, blast_df)
    
    #print(fwd_primer)
    #print(rev_primer)

        
    fwd_barcodes = set(samples_df['fwdBarcode'])
    rev_barcodes = set(samples_df['revBarcode'])

    best_fwd, best_rev = None, None
    if fwd_primer != None:
        fwd_search = search_for_barcode(fwd_primer, read, fwd_barcodes, rev_barcodes)
        best_fwd = filter_for_best_barcode(fwd_search)
    if rev_primer != None:
        rev_search = search_for_barcode(rev_barcodes, read, fwd_barcodes, rev_barcodes)
        best_rev = filter_for_best_barcode(rev_search)
    
    return concat_series_to_df([best_fwd, best_rev])


def assign_barcodes_all_reads(reads, blast_df, samples_df):
    all_codes = []
    for each_read in reads:
        barcodes = get_read_barcodes(each_read, blast_df, samples_df)
        all_codes.append(barcodes)
    
    return pd.concat(all_codes)

        
        
        
    



# reads_path = snakemake.input['reads']
# blast_tsv = snakemake.input['blast']
# samples = snakemake.input['samples']

# samples_df = pd.read_csv(samples, sep='\t')
# blast_df = read_blast_tsv(blast_tsv)
# reads = SeqIO.parse(reads_path, 'fastq')

# barcoded_reads = assign_barcodes_all_reads(reads, blast_df, samples_df)
# barcoded_reads.to_csv(
#     snakemake.output[0], sep='\t', index=None
# )


reads_path = '/home/ethollem/workflows/PBEH1/workflow/output/reads/seperatedTrim/PBEH1-2/13/sep.fastq'
blast_tsv = '/home/ethollem/workflows/PBEH1/workflow/output/barcode/primerBLAST/PBEH1-2/13.blastn.search.tsv'
samples = '/home/ethollem/workflows/PBEH1/resources/PBEH-1-samples.tsv'

samples_df = pd.read_csv(samples, sep='\t')
blast_df = read_blast_tsv(blast_tsv)
reads = SeqIO.parse(reads_path, 'fastq')

print(blast_df.head())

barcoded_reads = assign_barcodes_all_reads(reads, blast_df, samples_df)
barcoded_reads.to_csv(
    snakemake.output[0], sep='\t', index=None
)
    
    
    

    

# next step is to just find the barcode that is next to the aligned primer

    
    
    



