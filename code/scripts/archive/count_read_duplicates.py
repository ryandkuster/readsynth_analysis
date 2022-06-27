import numpy as np
import pandas as pd
import sys


'''
sys.argv[1] is the sam file to open
sys.argv[2] is the number of lines to skip in sam file
sys.argv[3] is the minimum number of times a pos must appear to be kept
'''

output_file = 'duplicate_counts.csv'
with open(sys.argv[1]) as f, open(output_file, 'w') as o:
    for i, line in enumerate(f):
        if i > int(sys.argv[2]):
            o.write('\t'.join(line.split('\t')[:9]) + '\n')


#df = pd.read_csv(output_file, header=None, sep ='\t')
df = pd.read_csv(output_file, header=None, delimiter=r"\s+")

forward_pos_df = df[3].value_counts().to_frame()
reverse_pos_df = df[7].value_counts().to_frame()

for_dt = forward_pos_df.to_dict()[3]
rev_dt = reverse_pos_df.to_dict()[7]

df['forward_count'] = df[3].map(for_dt)
df['reverse_count'] = df[7].map(rev_dt)

df = df[(df['forward_count'] >= int(sys.argv[3])) & (df['reverse_count'] >= int(sys.argv[3]))]

df.to_csv(output_file)

#TODO
'''
create a way to count the number of each read, while also tracking the read length
remove the reverse strand (it will create double entries)
'''
