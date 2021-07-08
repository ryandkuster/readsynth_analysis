import pandas as pd
import sys


"""
sys.argv[1] is the header-only file from the mock fastq community
sys.argv[2] is the sampled_files_key.txt file from the first step
sys.argv[3] is the digested_genome_stats.csv file from the 1520 EDA
"""

# open header file and count the occurances of each genome name
# this represents the raw fastq reads resulting from the simulation
# and sampling

count_dt = {}
with open(sys.argv[1]) as f:
    for line in f:
        GCA_id = line.split(':')[6].split()[0].split('_ASM')[0]
        if GCA_id in count_dt:
            count_dt[GCA_id] += 1
        else:
            count_dt[GCA_id] = 1
read_total = sum(count_dt.values())

# opent the digested_genome_stats file to assist in normalizing the abundance
# profiles for fair comparison

df_stats = pd.read_csv(sys.argv[3])
df = pd.read_csv(sys.argv[2], sep='\t', header=None, names=['file', 'copies'])
copy_total = df['copies'].sum()

copy_ratio_ls = []
read_ls = []
read_ratio_ls = []

final_df = pd.DataFrame()
for idx, row in df.iterrows():
    copy_ratio_ls.append(row['copies']/copy_total)
    for k, v in count_dt.items():
        if row['file'].startswith(k):
            read_ls.append(v)
            read_ratio_ls.append(v/read_total)
            tmp_df = df_stats.loc[df_stats['accession_id'] == k]
            final_df = pd.concat([final_df, tmp_df])

final_df['copies'] = df['copies']
final_df['copy_ratio'] = copy_ratio_ls
final_df['reads'] = read_ls
final_df['read_ratio'] = read_ratio_ls

final_df.to_csv('sampled_genome_stats.csv', index=None)
