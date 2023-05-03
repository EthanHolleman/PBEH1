
from pathlib import Path
import os
import subprocess
import pandas as pd



# read in input from snakemake rule

grouped_reads_dir = snakemake.input['reads']
output_dir = snakemake.output['output_dir']
valid_refs = snakemake.params['plasmid_id']
index = snakemake.input['index']
label = snakemake.params['label']
genome_dir = snakemake.params['genome_dir']


refs = list(pd.read_csv(valid_refs, sep='\t')['ref_name'])


def get_reads(grouped_reads_dir, refs):
    reads_paths = []
    for each_path in Path(grouped_reads_dir).iterdir():
        if each_path.stem.split('.group')[0] in refs:
            reads_paths.append(each_path)
    return reads_paths


def run_footloop(reads, output_dir, label, genome, index):
    
    target = Path(output_dir).joinpath(reads.stem.split('.group')[0])
    target.mkdir(parents=True)
    
    cmd = f'submodules/footLoop/footLoop.pl -r {reads} \
            -n {target} -l {label} -g {genome} -i {index}'
    
    
    if os.stat(reads).st_size > 0:
        #output = subprocess.run(cmd.split(' '))
        os.system(cmd)


def determine_genome(validated_read_path, genome_dir):
    genome_name = Path(validated_read_path).name.split('.')[0]
    return Path(genome_dir).joinpath(genome_name).with_suffix('.fasta')


# these are files grouped by plasmid that we want to map
read_files = get_reads(grouped_reads_dir, refs)

for each_read_file in read_files:
    
    # determine what file we should use as the reference
    genome = determine_genome(each_read_file, genome_dir)
    # check to make sure these are not read files without any plasmids 
    # mapped to them
    if 'nan' not in Path(each_read_file).name and 'None' not in Path(each_read_file).name:
        # run footloop program
        run_footloop(each_read_file, output_dir, label, genome, index)
        break









 