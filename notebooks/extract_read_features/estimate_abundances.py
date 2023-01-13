import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys

from scipy.stats.stats import pearsonr


# open and preprocess recreated fragments file
current_file = sys.argv[1]
df = pd.read_csv(current_file, index_col=0)
kmers_ls = df.columns.to_list()[15:-1]
kmers_ls = [i for i in kmers_ls if 'N' not in i]
print(kmers_ls)
df = df[df['internal']==0]
p_tot_frag = int(df['observed'].sum())
print(f'{p_tot_frag} total fragments observed')
df.sort_values('rel_abund', inplace=True)
p_uniq_frag = df.shape[0]
print(f'{p_uniq_frag} unique fragments')
p_taxa_no = len(df['genome'].unique())
print(f'{p_taxa_no} taxa observed')
gen_ls = list(df['genome'].unique())

# create the ground truth list in order of genomes list
e_ls = []
for gen in gen_ls:
    e = df[df['genome']==gen]['rel_abund'].unique()[0]
    e_ls.append(float(e))
e_ls_sum = sum(e_ls)
e_ls = [i/e_ls_sum for i in e_ls]
print(len(e_ls))
print(sum(e_ls))


def process_ratios(tmp_df, gen_ls):
    '''
    for a given length of fragment, within each genome x, get the average
    observed count as avg, then divide each observed count for every
    genome by avg this average and save as a column named after x
    '''
    for gen in gen_ls:
        avg = tmp_df[tmp_df['genome']==gen]['observed'].mean()
        tmp_df[gen] = tmp_df['observed'] / avg
    return tmp_df


def scale_ratios(np_arr):
    '''
    this approach assumes the first column, first row is a reliable
    representation of the real count data...
    '''
    rel_base = np_arr[0,0]
    for idx, i in enumerate(np_arr[0,:]):
        col_scale = rel_base/i
        np_arr[:,idx] = np_arr[:,idx]*col_scale
    return np_arr


def get_max_idx(ratios_df):
    max_ratios = ratios_df.count().max()
    for idx, i in enumerate(ratios_df.columns.to_list()):
        if ratios_df[i].count() == max_ratios:
            max_genome = i
            max_idx = idx
    return max_idx


def scale_ratios_to_max(np_ratios, max_idx):
    np_arr = np.copy(np_ratios)
    manp = np_arr[:,max_idx]

    for idx, i in enumerate(np_arr.T):
        col_scale = manp/i
        np_arr[:,idx] = np.nanmean(col_scale)*i
    return np_arr


def average_over_columns(np_arr):
    avg_ls = np.nanmean(np_arr, axis=1).tolist()
    return avg_ls


def return_rel_abund(o_ls):
    rel_ls = []
    for i in o_ls:
        rel_ls.append(i/sum(o_ls))
    return rel_ls


test_df = df.copy()
try:
    test_df.drop(gen_ls, inplace=True, axis=1)
except KeyError:
    pass
test_df = test_df.reindex(columns = test_df.columns.tolist() + gen_ls)

final_df = pd.DataFrame()

for j in range(0, test_df['length'].max()+1):
    if j % 100 == 0:
        print(f'processing fragments of {j}bp')
    tmp_df = test_df[test_df['length']==j].copy()
    if tmp_df.shape[0] > 0:
        tmp_df = process_ratios(tmp_df, gen_ls)
        final_df = pd.concat([final_df, tmp_df])

ratios_df = pd.DataFrame(0, index=gen_ls, columns=gen_ls)

for gena in gen_ls:
    for genb in gen_ls:
        ratios_df.loc[gena, genb] = final_df[final_df['genome']==gena][genb].mean()

np_ratios = np.array(ratios_df)
max_idx = get_max_idx(ratios_df)
scaled_arr = scale_ratios_to_max(np_ratios, max_idx)
o_ls = average_over_columns(scaled_arr)
o_ls = return_rel_abund(o_ls)

ratios_df.to_csv(os.path.join('./ratios',f'{os.path.basename(current_file)[:-4]}_ratios.csv'))

def count_to_rel(count_ls):
    rel_ls = []
    for i in count_ls:
        rel_ls.append(i/sum(count_ls))
    return rel_ls

mean_ls = []
median_ls = []

for gen in gen_ls:
    mean_ls.append(df[df['genome'] == gen]['observed'].mean())
    median_ls.append(df[df['genome'] == gen]['observed'].median())

mean_ls = count_to_rel(mean_ls)
median_ls = count_to_rel(median_ls)

abundances_df = pd.DataFrame()
abundances_df['genome'] = gen_ls
abundances_df['expected'] = e_ls
abundances_df['ratio'] = o_ls
abundances_df['mean'] = mean_ls
abundances_df['median'] = median_ls
abundances_df.to_csv(os.path.join('./ratios',f'{os.path.basename(current_file)[:-4]}_abundances.csv'))

ticks = [i for i in range(len(o_ls))]
labels = [i for i in gen_ls]

a = "{:.3f}".format(pearsonr(e_ls,o_ls)[0])
b = "{:.3f}".format(pearsonr(e_ls,mean_ls)[0])
c = "{:.3f}".format(pearsonr(e_ls,median_ls)[0])

plt.figure(figsize=(20,10))
plt.scatter(ticks, e_ls, c='black', marker='_', s=40, alpha=0.45,label='ground truth               pearson r')
plt.scatter(ticks, o_ls, c='green', s=20, alpha=0.25,     label=f'ratio estimation          {a}')
plt.scatter(ticks, mean_ls, c='orange', s=20, alpha=0.25, label=f'mean depth                {b}')
plt.scatter(ticks, median_ls, c='red', s=20, alpha=0.25,  label=f'median depth             {c}')
plt.legend()
plt.ylabel('relative abundance')
plt.xlabel('taxa')
plt.savefig(os.path.join('./figures', f'{os.path.basename(current_file)[:-4]}_abundance.tif'), dpi=350)
plt.savefig(os.path.join('./figures', f'{os.path.basename(current_file)[:-4]}_abundance.png'), dpi=350)
plt.close()


ticks = [i for i in range(len(o_ls))]
labels = [i for i in gen_ls]

a = "{:.3f}".format(pearsonr(e_ls,o_ls)[0])
b = "{:.3f}".format(pearsonr(e_ls,mean_ls)[0])
c = "{:.3f}".format(pearsonr(e_ls,median_ls)[0])

plt.figure(figsize=(20,10))
plt.scatter(ticks, e_ls, c='black', marker='_', s=40, alpha=0.45,label='ground truth               pearson r')
plt.scatter(ticks, o_ls, c='green', s=20, alpha=0.25,     label=f'ratio estimation          {a}')
plt.scatter(ticks, mean_ls, c='orange', s=20, alpha=0.25, label=f'mean depth                {b}')
plt.scatter(ticks, median_ls, c='red', s=20, alpha=0.25,  label=f'median depth             {c}')
plt.legend()
plt.yscale('log')
plt.ylabel('relative abundance (log10)')
plt.xlabel('taxa')
plt.savefig(os.path.join('./figures', f'log_{os.path.basename(current_file)[:-4]}_abundance.tif'), dpi=350)
plt.savefig(os.path.join('./figures', f'log_{os.path.basename(current_file)[:-4]}_abundance.png'), dpi=350)
plt.close()

