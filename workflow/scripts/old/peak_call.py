
from pathlib import Path
import os
import subprocess
import pandas as pd



# read in input from snakemake rule
mapping_dir = snakemake.input[0]
output_dir = snakemake.output[0]



def get_mapping_dirs(mapping_dir):
    # collect all the plasmid dirs that are stored in the mapping directory
    return [p for p in Path(mapping_dir).iterdir() if p.is_dir()]


def make_output_plasmid_dir(plasmid_dir, output_dir):
    # name and create a path to store footpeak results for a specific plasmid
    plasmid_name = Path(plasmid_dir).name
    full_path = Path(output_dir).joinpath(plasmid_name)
    full_path.mkdir(exist_ok=True, parents=True)
    return full_path


def run_footpeak(plasmid_dir, output_path):
    target_dir = make_output_plasmid_dir(plasmid_dir, output_path)
    cmd = f'submodules/footLoop/footPeak.pl -n {str(plasmid_dir)} -o {str(target_dir)}'
    os.system(cmd)



# Get all directories within the mapping dir that we need to call peaks on
footpeak_dirs = get_mapping_dirs(mapping_dir)

# run footpeak on each of these directories
for each_plasmid_dir in footpeak_dirs:
    run_footpeak(each_plasmid_dir, output_dir)









 