import re
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

# Parse alignments of primer sequences to reads produced by bowtie program
# Since primers are very similar (many only differ in their barcodes)
# and we run bowtie with -a param each read will align many of the same primer
# type (fwd primer 1, 2, 3 etc) and we need to find which primer (if any)
# both aligns (bowtie figures this out) and has a matching barcode
# To determine if there is a matching barcode this program will read an alignment
# from a sam file, lookup the sequnce from the read at the aligned position
# and then search for a match to the barcode that should be contained in the 
# aligned primer with hamming distance based metric.
# 
# The barcode is considered a match if it meets two criteria. The hamming
# distance must be below a user defined value (I set as 2) and there cannot
# be multible barcodes that result in a hamming distance that meets this
# threshold. For example if barcode 1 and barcode 2 both have a hamming
# distance of 2 for a given read then the read is not assigned a barcode
# and considered ambigious

def read_sam_file_to_dataframe(sam_file_path):
    col_names = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL',
                 'TAG1', 'VALUE1', 'TAG2', 'VALUE2', 'TAG3', 'VALUE3', 'TAG4', 'VALUE4', 'TAG5']
    sam_dataframe = pd.read_csv(sam_file_path, sep='\t', 
                                comment='@', names=col_names, 
                                header=None, low_memory=False)
    return sam_dataframe

def read_barcodes(samples_df):
    return set(samples_df['fwdBarcodes']), set(samples_df['revBarcodes'])


def hamming(search, target):
    assert len(target) == len(search)
    dist = 0
    for i in range(len(target)):
        if search[i] != target[i]:
            dist += 1
    
    return dist

def find_barcode(barcode, seq, hams_allowed=0):
    barcode_rc = str(Seq(barcode).reverse_complement())
    best_pos, best_score = -1, len(barcode)
    for i in range(len(seq)-len(barcode)):
        search_seq = seq[i:i+len(barcode)]
        dist = hamming(search_seq, barcode)
        if dist == 0:  # cant get better
            return i, 0
        if dist < best_score:
            best_score = dist
            best_pos = i
        
    return best_pos, best_score


def get_alignment_end(pos, cigar):
    end_pos = pos - 1
    num_bases = int(''.join(filter(str.isdigit, cigar)))
    for op in re.findall('[A-Z]', cigar):
        if op == 'M' or op == 'D' or op == 'N':
            end_pos += num_bases
        if op == 'M' or op == 'I' or op == 'S':
            num_bases = 0
    return end_pos


def score_alignment(read_seq, barcode, sam_row, hams_allowed=2):
    align_end = get_alignment_end(int(sam_row.POS), sam_row.CIGAR)
    search_seq = read_seq[int(sam_row.POS)-1: align_end-1]
    return find_barcode(barcode, search_seq, hams_allowed)
    
    
    
def process_sam_entry(sam_row, read_seq, sample_df):
    # process a primer alignment entry. Each entry represents a valid
    # bowtie alignment of a barcoded primer to a read but we need to
    # verify that the actual barcode portion of the read matches the aligned
    # primer
    primer_name = sam_row.QNAME
    primer_type = None
    
    if primer_name in set(sample_df['FWDPrimer']):
        barcode = sample_df.loc[sample_df['FWDPrimer'] == primer_name].iloc[0]['fwdBarcode']
        primer_type = 'fwd'
    else:
        barcode = sample_df.loc[sample_df['REVPrimer'] == primer_name].iloc[0]['revBarcode']
        primer_type = 'rev'
    
    barcode_rc = str(Seq(barcode).reverse_complement())
    score = score_alignment(read_seq, barcode, sam_row)
    score_rc = score_alignment(read_seq, barcode_rc, sam_row)
    if score[1] > score_rc[1]:  # best forward is worse than best rc
        score = score_rc

    return {
        'read': sam_row.RNAME,
        'primer_name': primer_name,
        'pos': score[0],
        'dist': score[1],
        'cigar': sam_row.CIGAR,
        'primerType': primer_type
    }
    

def row_is_valid(sam_row):
    if sam_row.RNAME != '*':
        return True
    else:
        return False


def process_sam_file(sam_df, sample_df, reads_dict):
    barcode_dicts = []
    for i, each_row in sam_df.iterrows():
        if row_is_valid(each_row):
            read = reads_dict[each_row.RNAME]
            barcode_dicts.append(process_sam_entry(each_row, read, sample_df))
    return pd.DataFrame(barcode_dicts)



def assign_primers_to_read(ops_df, read_name, max_ops=2):
    all_aligns = ops_df.loc[ops_df['read'] == read_name]
    
    def get_best_primer(primer_type):
        aligns = all_aligns.loc[all_aligns['primerType'] == primer_type]
        if len(aligns) == 0:
            return -1   # no aligned reads :(
        best_ops_val = min(aligns['dist'])
        if best_ops_val > max_ops:
            return -1
        else:
            # get all rows that have this min ops
            min_ops_rows = aligns.loc[aligns['dist'] == best_ops_val]
            if len(min_ops_rows) > 1:
                # more than one alignment has meet the min ops threshold
                # so for now this read is ambigious. Might be able to
                # make a better determination by seeing if the other primer
                # has a good barcode because this will reduce the number of
                # barcodes that have to be considered for the other primer
                return -1
            else:
                return min_ops_rows.iloc[0]
    
    fwd_primer = get_best_primer('fwd')
    rev_primer = get_best_primer('rev')
    
    primer_dict = {
        'fwd_primer': 'none',
        'rev_primer': 'none',
        'read': read_name,
        'fwd_dist': -1,
        'rev_dist': -1
    }

    if isinstance(fwd_primer, pd.Series):
        primer_dict['fwd_primer'] = fwd_primer['primer_name']
        primer_dict['fwd_dist'] = fwd_primer['dist']
    
    if isinstance(rev_primer, pd.Series):
        primer_dict['rev_primer'] = rev_primer['primer_name']
        primer_dict['rev_dist'] = rev_primer['dist']
    
    return primer_dict


def assign_primers(ops_df, max_ops=3):
    reads = set(ops_df['read'])
    assignments = []
    for each_read in reads:
        f = assign_primers_to_read(ops_df, each_read, max_ops)
        assignments.append(f)
    return pd.DataFrame(assignments)



def get_sample_id(primers_df, samples_df):
    
    def assign_sample(row):
        sample = samples_df.loc[(samples_df['FWDPrimer'] == row.fwd_primer) 
                       & (samples_df['REVPrimer'] == row.rev_primer)]
        if len(sample) > 0:
            return list(sample['Sample-ID'])[0]
        else:
            return -1
    
    primer_df['sample_id'] = primer_df.apply(lambda row: assign_sample(row), axis=1)
    
    return primer_df
        
    
    
        
    

reads = snakemake.input['reads']
sam =  snakemake.input['sam']
samples =  snakemake.input['samples']

sam_df = read_sam_file_to_dataframe(sam)
samples = pd.read_csv(samples, sep='\t')
reads_dict = {r.description: str(r.seq) for r in SeqIO.parse(reads, 'fastq')}


ops_df = process_sam_file(
    sam_df, samples, reads_dict
)

primer_df = assign_primers(ops_df)
primer_df = get_sample_id(primer_df, samples)

primer_df.to_csv(
    snakemake.output[0], sep='\t', index=None
    )
    
        
    
    