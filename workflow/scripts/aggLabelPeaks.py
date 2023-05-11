# After calling all peaks for a given flowcell, merge all of these files
# back into one large file so can access all peaks at once for plotting
# I'm using pandas here instead of cat as a way to avoid having
# header rows in the middle of files; sue me. 

import pandas as pd

peak_paths = snakemake.input
peak_dfs = [pd.read_csv(path, sep='\t') for path in peak_paths]
big_df = pd.concat(peak_dfs)
big_df.to_csv(snakemake.output[0], sep='\t', index=None)
