import itertools
import numpy as np
import os
import pandas as pd
import sys


'''
sys.argv[1] is directory with readsynth output
sys.argv[2] is kmer size (k)
sys.argv[3] is a prefix for the output name
'''


def kmer_combos(k):
    kmers = []
    k_list = ['ACGTN' for i in range(k)]
    kmer_ls = list(itertools.product(*k_list))
    for i in kmer_ls:
        kmers.append(''.join(i))
    return kmers        


def count_kmers(seq, k, kmers):
    kmers_dt = {kmer:0 for kmer in kmers}
    for i in range(len(seq)-k+1):
        kmers_dt[seq[i:i+k]] += 1

    return kmers_dt


k = int(sys.argv[2])
kmers = kmer_combos(k)
final_df = pd.DataFrame(columns = kmers)


for i in os.listdir(sys.argv[1]):
    if i.startswith('counts'):
        genome = i[len('counts_'):-len('.gz.csv')]
        print(genome)
        df = pd.read_csv(i)
        df = df[df['internal'] == 0]
        total_kmers_dt = {j:0 for j in kmers}
        for seq in df['seq']:
            kmers_dt = count_kmers(seq, k, kmers)
            total_kmers_dt = {kmer: total_kmers_dt[kmer] + kmers_dt[kmer] for kmer in kmers}
        final_df.loc[genome] = total_kmers_dt.values()


kmers_final = [i for i in kmers if 'N' not in i]
gen_ls = final_df.index.values.tolist()
std_ls = [final_df[i].std() for i in kmers_final]
kmers_final = [i for i, j in zip(kmers_final, std_ls) if j > 0]
std_ls = [i for i in std_ls if i > 0]

final_df_e = final_df[kmers_final].copy()
final_df_m = final_df[kmers_final].copy()
final_df_e[gen_ls] = 0
final_df_m[gen_ls] = 0


for idx, gen1 in enumerate(gen_ls):
    for gen2 in gen_ls:
        e_ls = final_df.loc[gen1][kmers_final] - final_df.loc[gen2][kmers_final]
        e_ls = e_ls ** 2
        euc_sum = sum(e_ls.to_list())
        final_df_e.at[gen1, gen2] = euc_sum
        final_df_e.at[gen2, gen1] = euc_sum

        m_ls = (final_df.loc[gen1][kmers_final] / std_ls) - (final_df.loc[gen2][kmers_final] /  std_ls)
        m_ls = m_ls ** 2
        mah_sum = sum(m_ls.to_list())
        final_df_m.at[gen1, gen2] = mah_sum
        final_df_m.at[gen2, gen1] = mah_sum
    print(idx)
final_df_e.drop(kmers_final, axis=1, inplace=True)
final_df_e.to_csv(f'{sys.argv[3]}_{str(k)}_euc_distance.tsv', sep='\t')

final_df_m.drop(kmers_final, axis=1, inplace=True)
final_df_m.to_csv(f'{sys.argv[3]}_{str(k)}_mah_distance.tsv', sep='\t')



