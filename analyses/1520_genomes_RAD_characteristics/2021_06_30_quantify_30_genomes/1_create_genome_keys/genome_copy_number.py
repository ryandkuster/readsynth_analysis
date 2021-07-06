import pandas as pd
import random


df = pd.read_csv('sampled_files.txt', header=None)

copy_ls = []
for i in range(30):
    copy_ls.append(random.randint(1,100))

df['copies'] = copy_ls
df.to_csv('sampled_files_key.txt', header=None, sep='\t', index=None)
