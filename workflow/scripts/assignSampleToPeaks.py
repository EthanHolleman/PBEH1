from pathlib import Path
import pandas as pd



def read_peak_bed(bed_path):
    # read peak bed file and add extra column that is just
    # the read name
    df = pd.read_csv(bed_path, sep='\t', header=None)
    print(df.head())
    df.columns = ['peakID', 'start', 'end']
    
    reads = [row['peakID'].split('.')[-1] for _, row in df.iterrows()]
    df['read'] = reads
    return df


def read_sample_assignments(tsv_path):
    return pd.read_csv(tsv_path, sep='\t')



def get_peak_files(mapping_dir):
    peak_dir = Path(mapping_dir).joinpath('PEAKS_LOCAL')
    if not peak_dir.is_dir():
        return []
    return list(peak_dir.iterdir())


def assign_strand_to_peak(peak_file):
    peak_file_stem = Path(peak_file).stem
    chunks = peak_file_stem.split('_')
    for each_chunk in chunks:
        if each_chunk == 'Neg':
            return '-'
        if each_chunk == 'Pos':
            return '+'


def determine_call_type(peak_file):
    peak_file_stem = Path(peak_file).stem
    call_type = peak_file_stem.split('.PEAK')[0].split('_')[-1]
    return call_type
    

def merge_peaks_and_sample_assignment(peak_file, sample_df, plasmid_name):
    peak_df = read_peak_bed(peak_file)
    merge_df = peak_df.merge(sample_df, on='read')
    strand = assign_strand_to_peak(str(peak_file))
    call_type = determine_call_type(peak_file)
    merge_df['strand'] = strand
    merge_df['call_type'] = call_type
    merge_df['plasmid_name'] = plasmid_name
    return merge_df


peak_dir = snakemake.input['peaks']
output_path = snakemake.output[0]

sample_df = read_sample_assignments(snakemake.input['samples'])

peak_files = get_peak_files(peak_dir)
merged_peaks = []

plasmid_name = snakemake.params['plasmid_name']


# Assign all peaks to a sample and annotate based on the file names
# which include info on call type (not exactly sure what that is)
# and the strand the call was made on

if peak_files:  # peaks were called 

    for each_file in peak_files:
        merged_peaks.append(
            merge_peaks_and_sample_assignment(
                each_file, sample_df, plasmid_name
                )
        )
    all_merged = pd.concat(merged_peaks)
    all_merged.to_csv(output_path, sep='\t', index=None)
    print(f'{len(all_merged)}', 'peaks assigned to samples')

else:  # no peaks called for this input write empty file
    print('NO PEAKS CALLED')
    with open(output_path, 'w') as handle:
        pass
    
    
    





