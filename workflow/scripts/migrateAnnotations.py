# Move annotations from a genbank file to a fasta file
# that has had its start position shifted by some amount
# Produces a new genbank file based on the sequence
# of the fasta file but with as many annotations from the original
# genbank file transfered as possible. 

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import ExactPosition
from Bio import SeqIO
from pathlib import Path


def transfer_feature(feature, genbank, record):
    feature_seq = str(feature.extract(genbank).seq)
    start = str(record.seq).find(feature_seq)
    if start != -1:
        end = start + len(feature_seq)
        new_loc = FeatureLocation(start, end)
        feature.location = new_loc
        record.features.append(feature)
    
    return record


def transfer_all_features(genbank, record):
    record.name = genbank.name
    record.id = genbank.id
    record.description = genbank.description
    record.annotations["molecule_type"] = "DNA"
    for each_feature in genbank.features:
        record = transfer_feature(each_feature, genbank, record)
    return record


def get_output_path(record_path, output_path):
    return Path(output_path).joinpath(
        f'{Path(record_path).stem}.shift.gb'
    )



genbank_files = [f for f in Path(snakemake.params['genbank_files']).iterdir() if f.suffix == '.gb']
fasta_files = [f for f in Path(snakemake.input).iterdir() if f.suffix == '.fasta']

fasta_lookup = {f.stem: f for f in fasta_files}

# match genbank file to their shifted fasta version
file_pairs = []
for each_genbank in genbank_files:
    try:
        file_pairs.append(
            (each_genbank, fasta_lookup[each_genbank.stem])
        )
    except KeyError:
        continue

# transfer the annotations from the orignal genbank using the
# shifted fasta to create an annotated version of the shifted
# genbank file
shifted_gb = []

for gb, fa in file_pairs:
    gb_record = SeqIO.read(gb, 'genbank')
    fa_record = SeqIO.read(fa, 'fasta')
    shifted = SeqRecord(fa_record.seq)
    shifted_gb.append(
        transfer_all_features(gb_record, shifted)
    )
    output_path = get_output_path(
        gb_record, snakemake.output[0]
    )
    SeqIO.write(shifted_gb, output_path, 'genbank')
    


    



fa = SeqIO.read(shifted_fa, 'fasta')
record = SeqRecord(fa.seq)
record = transfer_all_features(gb, record)

SeqIO.write(record, snakemake.output[0], 'genbank')
