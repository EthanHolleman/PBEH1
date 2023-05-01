from pathlib import Path
import sys
import shutil
import pandas as pd

target_dir = sys.argv[1]
filename = sys.argv[2]

files = Path(target_dir).iterdir()

records = []

for f in files:
    dir_path = f.parent.parent.joinpath(f.stem)
    dir_path.mkdir(parents=True, exist_ok=True)
    filepath = dir_path.joinpath(filename)
    shutil.copyfile(f, filepath)
    

# df = pd.DataFrame(records)
# df.columns = 'split_files'
# df.to_csv(
#     snakemake.output[0], index=False
# )


