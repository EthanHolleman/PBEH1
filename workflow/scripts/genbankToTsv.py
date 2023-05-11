# Takes a list of genbank files and produces a tsv file that contains
# the features (minus and plus ones I want to add or remove for plotting
# purposes). Having features in tsv format makes things a lot easier for
# ggplot functions. It should be notes that since the start positions of
# genbank files are shifted before alignment that the features need to
# be extract from a genbank file that compensates for this otherwise things
# will not line up correctly. See migrateAnnotations.py script. 


from Bio import SeqIO
import pandas as pd
from pathlib import Path


def get_name(feature):
    """
    Returns the name of a feature from a Biopython GenBank record's feature object.
    """
    if "label" in feature.qualifiers:
        return feature.qualifiers["label"][0]
    elif "product" in feature.qualifiers:
        return feature.qualifiers["product"][0]
    else:
        return feature.type


def genbank_to_tsv(genbank_path):
    feats = []
    record = SeqIO.read(genbank_path, 'genbank')

    for each_feat in record.features:
        
        name = get_name(each_feat)
        if 'prime' not in name and 'Anchor' not in name:
            # remove some features that I don't really add value to plotting
            feats.append({
                'start': int(each_feat.location.start),
                'end':   int(each_feat.location.end),
                'plasmid': Path(genbank_path).stem,
                'strand': each_feat.location.strand,
                'name': get_name(each_feat)
            }
            )
            if 'Variable region' in name and 'T7' in Path(genbank_path).stem:
                # if feature and genbank meet these conditions this should
                # be a T7 Init VR plasmid. For plotting we also want to annotate
                # the SNRPN region which extends 531 bp downstream of the end of
                # the variable region
                feats.append(
                    {
                    'start': int(each_feat.location.end),
                    'end':   int(each_feat.location.end+531),
                    'plasmid': Path(genbank_path).stem,
                    'strand': each_feat.location.strand,
                    'name': 'SNRPN'
                    }
                )
                
    feats_df = pd.DataFrame(feats)
    
    return feats_df



genbank_paths = snakemake.input
features_dfs = [genbank_to_tsv(each_path) for each_path in genbank_paths]
all_plasmids_df = pd.concat(features_dfs)
all_plasmids_df.to_csv(tsv_path, sep='\t', index=None)