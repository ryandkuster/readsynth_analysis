import pandas as pd
import sys


"""
sys.argv[1] is the bracken output (e.g. '..bracken_genuses.kreport')
sys.argv[2] is the desired level (e.g. 'G' for genus)
sys.argv[3] is the sampled_genome_stats.csv file
"""

# pull only the relevant levels and sort descending by abundance

df = pd.DataFrame()
percent_ls = []
read_ls = []
id_ls = []
name_ls = []

with open(sys.argv[1]) as f:
    for line in f:
        line = line.rstrip()
        if line.split()[3] == sys.argv[2]:
            percent_ls.append(float(line.split()[0]))
            read_ls.append(line.split()[1])
            id_ls.append(line.split()[4])
            name_ls.append(line.split()[5])

df['percent'] = percent_ls
df['reads'] = read_ls
df['taxid'] = id_ls
df['name'] = name_ls
df = df.sort_values('percent', ascending=False)
            
#df.to_csv('bracken_' + sys.argv[2] + '_counts.csv', index=None)            

# open the sampled_genome_stats file and search for taxids that were
# either false positives or false negatives

true_pos_reads = []
true_pos_percent = []
true_pos_names = []
genus_names = []

df_stats = pd.read_csv(sys.argv[3])

for idx1, row1 in df_stats.iterrows():
    for idx2, row2 in df.iterrows():
        if row2['name'] in row1['taxonomy']:
            genus_names.append(row2['name'])
            true_pos_names.append(row1['taxonomy'])
            true_pos_reads.append(row2['reads'])
            true_pos_percent.append(row2['percent'])
    if row1['taxonomy'] not in true_pos_names:
        genus_names.append(' ')
        true_pos_reads.append(0)
        true_pos_percent.append(0)

df_stats['bracken_G'] = genus_names
df_stats['bracken_G_reads'] = true_pos_reads
df_stats['bracken_G_percent'] = true_pos_percent

df_stats = df_stats.sort_values('bracken_G')
df_stats.to_csv('profiled_genome_stats.csv', index=None)
