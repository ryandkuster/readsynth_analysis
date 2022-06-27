import os
import pandas as pd
import sys


"""
sys.argv[1] is the samtools depth file with header (-H) argument
sys.argv[2] is the filename to print to stdout
"""

df = pd.read_csv(sys.argv[1], sep="\t")
colname = df.columns[2]
df = df[df[colname] > 0]
print(f'{sys.argv[2]}\t{df[colname].sum() / len(df)}\t{df[colname].median()}')
