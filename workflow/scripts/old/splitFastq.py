# splits one fastq with many reads into many fastq with
# few reads
from Bio import SeqIO
from pathlib import Path
import gzip
import shutil
import os


def determine_number_reads(input_fastq):
    # determine the number of reads in fastq file without reading into memory
    reads = 0
    records = SeqIO.parse(input_fastq, 'fastq')
    for _ in records:
        reads += 1
    return reads


def gzip_file(file_path):
    """
    Compresses a file using gzip.

    Args:
        file_path (str): The path of the file to be compressed.

    Returns:
        str: The path of the gzipped file.
    """
    # Create the gzipped file path
    gzipped_file_path = str(file_path) + '.gz'

    # Open the original file in binary mode
    with open(file_path, 'rb') as file_in:
        # Open the gzipped file in binary write mode
        with gzip.open(gzipped_file_path, 'wb') as file_out:
            # Copy the contents of the original file to the gzipped file
            shutil.copyfileobj(file_in, file_out)
    os.remove(file_path)

    return gzipped_file_path



input_fastq = snakemake.input[0]

num_reads = determine_number_reads(input_fastq)
num_outputs = snakemake.params['num_outputs']

reads = SeqIO.parse(input_fastq, 'fastq')
buffer = []

# determine number of reads per file
reads_per_file = num_reads / num_outputs

if reads_per_file < 1:
    # more files than reads
    reads_per_file = 1

current_file = 0

if not Path(snakemake.output['output_dir']).is_dir():
    Path(snakemake.output['output_dir']).mkdir()


while current_file < num_outputs:
    try:
        buffer.append(next(reads))
        if len(buffer) >= reads_per_file:
            output_path = Path(
                snakemake.output['output_dir']).joinpath(
                f'{snakemake.params["flow_cell"]}_{current_file}.fastq'
            )
            SeqIO.write(buffer, output_path, 'fastq')
            gzip_file(output_path)
            buffer = []
            current_file += 1
    except StopIteration:
        break


