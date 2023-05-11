#!/usr/bin/python

from pathlib import Path
import os

err = sorted(
    [f for f in Path('.').iterdir() if f.suffix=='.err'],
    key=lambda f: int(f.stem.split('j')[-1])
)

os.system(f'tail -n 50 {err[-1]}')