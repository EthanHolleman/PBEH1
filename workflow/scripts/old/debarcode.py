from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
import gzip


def limited_hamming(search, target, hams_allowed):
    assert len(target) == len(search)
    dist = 0
    for i in range(len(target)):
        if search[i] != target[i]:
            dist += 1
            if dist > hams_allowed:
                return False
    
    return True

def find_barcode(barcode, seq, hams_allowed=0):
    barcode_rc = str(Seq(barcode).reverse_complement())
    for i in range(len(seq)-len(barcode)):
        search = seq[i:i+len(barcode)]
        if limited_hamming(search, barcode, hams_allowed):
            return i, False
        if limited_hamming(search, barcode_rc, hams_allowed):
            return i, True
    return None, -1


def search(read, barcode_list, search_depth=200, left=True, hams_allowed=0):
    if len(read) < search_depth:
        search_string = str(read)
    else:
        if left:
            # search the left side of the read (fwd part of read)
            search_string = read[:search_depth]
        else:
            # otherwise search the end of the read
            search_string = read[-200:]
    
    for each_barcode in barcode_list:
        # look for each barcode in the search string
        result, rc = find_barcode(each_barcode, search_string, hams_allowed)
        if result:
            # if found the barcode return the barcode and where it was found
            # if this is from the end of the read need to adjust position
            if not left:
                return each_barcode, len(read) - result, rc
            else:
                return each_barcode, result, rc

    # went through all barcodes but did not find any return none
    return None, -1, -1


def front_back_barcode_search(read, barcodes, *args):
    front = search(read, barcodes, *args)
    back = search(read, barcodes, left=False, *args)
    return front, back


def extract_barcodes_from_df(sample_df):
    codes = set(sample_df['fwdBarcode']).union(set(sample_df['revBarcode']))
    barcodes = []
    for c in codes:
        if isinstance(c, str):
            barcodes.append(c)
    return barcodes


def read_fastq(filepath):
    
    return SeqIO.parse(filepath, 'fastq')



def main():

    # read file
    reads =read_fastq(snakemake.input['fastq'])
    sample_df = snakemake.params['sample_df']
    barcodes = extract_barcodes_from_df(sample_df)
    
    barcode_list = []
    first_write = False

    for each_read in reads:
        codes = front_back_barcode_search(each_read, barcodes)
        code_dict = {
            'fwdBarcode': codes[0][0],
            'fwdBarcodeLoc': codes[0][1],
            'fwdRC': codes[0][2],
            'revBarcode': codes[1][0],
            'revBarcodeLoc': codes[1][1],
            'revRC': codes[1][2],
            'read': each_read.description,
            'readLength': len(each_read)
        }
        barcode_list.append(code_dict)
        if len(barcode_list) % 250 == 0:
            # write barcodes every 250 entries to avoid storing
            # everything in memory
            df = pd.DataFrame(barcode_list)
            if first_write:
                df.to_csv(
                    snakemake.output[0],
                    index=False, header=True, sep='\t'
                )
                first_write = False
            else:
                df.to_csv(
                    snakemake.output[0],
                    index=False, header=False, mode='a', sep='\t'
                )
            barcode_list = []



        
    
if __name__ == '__main__':
    main()

