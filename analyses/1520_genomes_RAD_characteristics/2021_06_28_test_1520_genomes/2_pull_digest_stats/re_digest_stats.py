import os
import pandas as pd
import sys 


"""
sys.argv[1] is the genome_stats.csv file
sys.argv[2] is the directory with raw_digest...csv files
sys.argv[3] is the directory with .fna genomes
"""

df = pd.read_csv(sys.argv[1])

frags_total = []
frags_1_to_200 = []
frags_201_to_400 = []
frags_401_to_600 = []
total_gc_ls = []
total_n_ls = []
total_all_ls = []
seq = ''


def count_gc(total_gc, seq):
    for base in seq:
        if base in ['G', 'C']:
            total_gc += 1

    return total_gc


for idx, row in df.iterrows():
    print(f'processing genome {idx}')
    genome = row['accession_id']
    total_gc = 0 
    total_n = 0 
    total_all = 0 
    seq = ''
    for filename in os.listdir(sys.argv[2]):
        if filename.startswith('raw_digest_' + genome):
            count_df = pd.read_csv(os.path.join(sys.argv[2], filename))
            template_df = count_df.loc[count_df['strand'] == '+']
            frags_total.append(template_df.shape[0])
            frags_1_to_200.append(template_df['length'].between(1, 200, inclusive=True).sum())
            frags_201_to_400.append(template_df['length'].between(201, 400, inclusive=True).sum())
            frags_401_to_600.append(template_df['length'].between(401, 600, inclusive=True).sum())

    for filename in os.listdir(sys.argv[3]):
        if filename.startswith(genome):
            with open(os.path.join(sys.argv[3], filename)) as f:
                for line in f:
                    if line.startswith('>') and seq:
                        total_gc = count_gc(total_gc, seq)
                        total_n += seq.count('N')
                        total_all += len(seq)
                        seq = ''
                    elif line.startswith('>'):
                        continue
                    else:
                        seq += line.rstrip().upper()
            
                total_gc = count_gc(total_gc, seq)
                total_n += seq.count('N')
                total_all += len(seq)

    total_gc_ls.append(total_gc)
    total_n_ls.append(total_n)
    total_all_ls.append(total_all)

df['frags_total'] = frags_total
df['frags_1_to_200'] = frags_1_to_200
df['frags_201_to_400'] = frags_201_to_400
df['frags_401_to_600'] = frags_401_to_600
df['total_gc_ls'] = total_gc_ls
df['total_n_ls'] = total_n_ls
df['total_all_ls'] = total_all_ls

df.to_csv('digested_' + sys.argv[1], index=None)
