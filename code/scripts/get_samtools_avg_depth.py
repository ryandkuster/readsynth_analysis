import os
import pandas as pd
import sys


"""
sys.argv[1] is the samtools depth file with header (-H) argument
"""

df = pd.read_csv(sys.argv[1], sep="\t")
colname = df.columns[2]
print(df[colname].sum() / len(df))
