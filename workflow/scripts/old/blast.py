from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import os

num_cores = os.cpu_count()

def hamming_distance(string1, string2, max_distance):
    if len(string1) != len(string2):
        return -1  # error, strings are of different length
    distance = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            distance += 1
            if distance > max_distance:
                return False
    return True


def windowed_hamming(big_string, small_string, max_dist):
    for i in range(len(big_string)-len(small_string)):
        substring = big_string[i:i+len(small_string)]
        dist = hamming_distance(substring, small_string, max_dist)
        if dist:
            return i
    
    return -1


def get_primer_barcodes(sample_df):
    
    fwd_primers = list(sample_df['FwdPrimerSeq'])
    fwd_codes = list(sample_df['fwdBarcode'])
    
    rev_primers = list(sample_df['RevPrimerSeq'])
    rev_codes = list(sample_df['revBarcode'])
    
    return dict(zip(fwd_primers, fwd_codes)), dict(zip(rev_primers, rev_codes))



def primer_search(read, read_name, primer_barcode_dict):
    
    # first identify the primer species (forward or reverse) by hamming distance
    # search, barcodes allow up to 10% mismatches
    primer_hits = {}
    for each_primer in primer_barcode_dict:
        search = windowed_hamming(read, each_primer, 1-int(0.9*len(each_primer)))
        if search != -1:
            primer_hits[each_primer] = search
            print('hit foward')
    
    if not primer_hits:  # no hits? Try the reverse complement
        print('RC')
        read = str(Seq(read).reverse_complement())
        for each_primer in primer_barcode_dict:
            search = windowed_hamming(read, each_primer, 1-int(0.9*len(each_primer)))
            if search != -1:
                primer_hits[each_primer] = search
        
    
    # found some number of primers that have meet mismatch threshold
    barcode_hits = []
    if primer_hits:
        for each_hit in primer_hits:
            
            hit_loc = primer_hits[each_hit]
            read_hit_substring = read[hit_loc:hit_loc+len(each_hit)]
            
            barcode = primer_barcode_dict[each_hit]
            
            barcode_search = windowed_hamming(
                read_hit_substring, barcode, 1)
            
            if barcode_search != -1:

                # found a barcode that meets the thresehold so can
                # identify this barcoded primer
                barcode_hits.append(
                    {
                        'primer_seq': each_hit,
                        'read_name': read_name,
                        'barcode': barcode,
                        'actual_seq': read_hit_substring[barcode_search:barcode_search+len(barcode)]
                    }
                )
    return pd.DataFrame(barcode_hits)



def identify_read_primers(record, sample_df):
    fwd_primers, rev_primers = get_primer_barcodes(sample_df)
    fwd_search = primer_search(str(record.seq), record.description, fwd_primers)
    rev_search = primer_search(str(record.seq), record.description, rev_primers)
    
    id = 'ambig'
    
    if len(fwd_search) > 0 and len(rev_search) > 0:
        
        print('primers found')
        
        fwd_seqs, rev_seqs = set(fwd_search['primer_seq']), set(rev_search['primer_seq'])
        sample = sample_df.loc[(sample_df['FwdPrimerSeq'].isin(fwd_seqs)) & (sample_df['FwdPrimerSeq'].isin(rev_seqs))]
    
        if len(sample) == 1:
            id =  list(sample['Sample-ID'])[0]
        elif len(sample) > 1:
            id =  'ambig'
        else:
            id = 'none'
        
    return {
        'read_name': record.description,
        'sample': id
    }
    


reads = snakemake.input['reads']
samples = snakemake.input['samples']

samples_df = pd.read_csv(
    samples, sep='\t'
)

read_ids = []
reads = SeqIO.parse(reads, 'fastq')

ids = []

with ThreadPoolExecutor(max_workers=num_cores) as executor:
    # Map the function over the argument tuples using the thread pool     
    args_list = [(r, samples_df) for r in reads]
    results = executor.map(lambda args: identify_read_primers(*args), args_list)
    for future in concurrent.futures.as_completed(results):
        result = future.result()
        ids.append(result)


id_df = pd.DataFrame(read_ids)
id_df.to_csv(
    snakemake.output[0],
    sep='\t', index=None, header=None
)
