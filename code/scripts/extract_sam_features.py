import argparse
import itertools
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statistics
import sys


'''
NOTE: this is designed specifically for ECORI and MSEI libraries

-i is a abundance table (key) similar to the readsynth input, except
that instead of the path to the reference genomes, it will give the path to each

    mapped sam file for that genome
        column1 = path to sam file
        column2 = abundance per each genome
        column3 = path to raw_digest file from simulation

-o is the name of the output csv to write

for each sam file, load all non-header rows as pandas df
use a readsynth raw_digest file (column 3 in key) to pull out the rows that
match the expected start and stop sites and save the number of read pairs
that do match (basically add count data to the raw digest dataframe)

key column 2 is a label to apply to the whole dataframe (e.g. read depth or
abundance that is expected)
'''


def main():
    args = parse_user_input()

    final_df = pd.DataFrame()
    key = pd.read_csv(args.i, header=None)
    key_sum = sum(key[1].to_list())
    for sam, label, cut_sites in zip(key[0].to_list(), key[1].to_list(), key[2].to_list()):
        genome_name = cut_sites.split('/')[-1][11:-8]
        print(genome_name)
        df = process_sam(sam, args)
        cut_sites = pd.read_csv(cut_sites)
        df, cut_sites = remove_unmatched(df, cut_sites)
        cut_sites = count_fragments(df, cut_sites)
        cut_sites['copy_number'] = label
        cut_sites['rel_abund'] = label/key_sum
        cut_sites['genome'] = genome_name
        final_df = pd.concat([final_df, cut_sites], axis=0)

    if args.k:
        kmers = kmer_combos(args.k)
        final_df[kmers] = 0
        final_df[kmers] = final_df['seq'].apply(lambda x: pd.Series(count_kmers(x, args.k, kmers)))

    final_df.reset_index(inplace=True)
    final_df.to_csv(f'{args.o}.csv',index=None)
    final_dt = final_df.groupby('length')['observed'].sum().to_dict()
    write_to_json(f'{args.o}.json', final_dt)

    save_histogram(final_df, f'{args.o}.pdf')

    if args.p:
        final_df = remove_highly_similar(final_df, args)
        final_df.reset_index(inplace=True,drop=True)
        final_df.to_csv(f'{args.o}_no_dupes.csv',index=None)
        save_histogram(final_df, f'{args.o}_no_dupes.pdf')
        final_dt = final_df.groupby('length')['observed'].sum().to_dict()
        write_to_json(f'{args.o}_no_dupes.json', final_dt)

    calculate_cut_efficiency(final_df, args)

    print('\nwriting sam files with unique hits that match motifs')

    for sam, label, cut_sites in zip(key[0].to_list(), key[1].to_list(), key[2].to_list()):
        genome_name = cut_sites.split('/')[-1][11:-8]
        print(genome_name)
        tmp_df = final_df[final_df['genome'] == genome_name].copy()
        q_ls = [j for i in tmp_df['sam_queries'].to_list() for j in i]
        print(len(q_ls))
        write_modified_sam(sam, q_ls)


def parse_user_input():
    parser = argparse.ArgumentParser(description='simulate RAD libary')

    parser.add_argument('-i', type=str, required=True,
                        help='path to key csv (sam,abundance,raw_digest')

    parser.add_argument('-o', type=str, required=True,
                        help='prefix name of files to store output')

    parser.add_argument('-p', type=float, required=False,
                        help='percent of read that must have unique hamming distance')

    parser.add_argument('-k', type=int, required=False,
                        help='kmer size, if kmer summary is desired')

    parser.add_argument('-s', type=int, required=True,
                        help='when adjusting for imperfect adapter removal, soft clips to allow')

    args = parser.parse_args()

    return args


def process_sam(sam, args):
    with open(sam) as f:
        for idx, line in enumerate(f):
            if line.startswith('@'):
                pass
            else:
                sam_header = idx
                break

    #skip the header rows we just learned from the loop above and open the file
    df = pd.read_csv(sam, skiprows=sam_header, usecols=range(11), sep='\s+', header=None)
    cols = 'query flag ref start qual cigar pair pair_pos length seq score'.split()
    df.columns = cols

    # get the position where the mapped read ends
    df['end'] = df['start'] + df['seq'].str.len() - 1

    # apply boolean based on bitwise flag to reflect if first/second in pair
    df['forward'] = df['flag'].apply(lambda x: check_first_in_pair(x))

    # drop the unneeded sam columns
    df = df.drop(['flag', 'qual', 'pair', 'pair_pos', 'score'], axis=1)
    df = adjust_clipping(df, args)

    # create two separate dataframes for the first/second (R1/R2) reads
    df_for = df[df['forward'] == True].copy()
    df_rev = df[df['forward'] == False].copy()

    # the following two lines should fix the issue with imperfect adapter removal affecting motif start/ends
    df_for.loc[:, 'end'] = np.where(df_for['length']<0, df_for['end']-df_for['l_clip'], df_for['end'])
    df_rev.loc[:, 'end'] = np.where(df_rev['length']<0, df_rev['end']-df_rev['l_clip'], df_rev['end'])

    # merge the first/second read dataframes together so each row is a pair
    df = df_for.merge(df_rev, on="query")

    # calculate where the ecori or msei cut site should be located, adjusted both
    # for overhang and directionality
    # np.where is using a true/false condition for the last two elements
    df['m1_pos'] = np.where(df['length_x']>0, df['start_x']-2, df['end_x']-5)
    df['m2_pos'] = np.where(df['length_x']<0, df['start_y']-2, df['end_y']-3)
    return df


def adjust_clipping(df, args):
    l_clip = []
    cigars_ls = df['cigar'].to_list()
    for i in cigars_ls:
        if i.count('S') == 1:
            if i.endswith('S'):
                try:
                    l_clip.append(0)
                except ValueError:
                    l_clip.append(0)
            else:
                try:
                    soft_clips = int(i.split('S')[0])
                    if soft_clips <= args.s:
                        l_clip.append(soft_clips)
                    else:
                        l_clip.append(0)
                except ValueError:
                    l_clip.append(0)
        else:
            l_clip.append(0)
    df['l_clip'] = l_clip
    return df


def remove_unmatched(df, cut_sites):
    '''
    increase processing speed by removing non-relevant sites from sam reads and sim fragments
    '''
    sam_sites = df['m1_pos'].unique().tolist() + df['m2_pos'].unique().tolist()
    sim_sites = cut_sites['start'].unique().tolist() + cut_sites['end'].unique().tolist()
    sim_sites = list(set(sim_sites))

    df = df[df['m1_pos'].isin(sim_sites)]
    df = df[df['m2_pos'].isin(sim_sites)]
    cut_sites = cut_sites[cut_sites['start'].isin(sam_sites)]
    cut_sites = cut_sites[cut_sites['end'].isin(sam_sites)]

    return df, cut_sites


def check_first_in_pair(flag):
    '''
    check sam flag to see if read is first in pair (returns true)
    '''
    return bin(flag)[-7] == '1'


def count_fragments(df, cut_sites):
    '''
    for each simulated fragment, find the number of occurances from the sam file
    '''
    cut_sites.reset_index(inplace=True, drop=True)
    count_ls = []
    queries_ls = []
    np_cut_sites = np.array(cut_sites)
    for i in np_cut_sites:
        chrom = i[0]
        c_start = i[2]
        c_end = i[3]
        c_m1 = i[4]
        c_m2 = i[5]
        if c_m1 == 'GAATTC':
            tmp_df_top = df[(df['m1_pos']==c_start) & (df['m2_pos']==c_end) & (df['ref_x']==chrom)]
            count_ls.append(tmp_df_top.shape[0])
            if tmp_df_top.shape[0] == 0:
                queries_ls.append([])
            else:
                queries_ls.append(tmp_df_top['query'].to_list())
        elif c_m1 == 'TTAA':
            tmp_df_bot = df[(df['m1_pos']==c_end) & (df['m2_pos']==c_start) & (df['ref_x']==chrom)]
            count_ls.append(tmp_df_bot.shape[0])
            if tmp_df_bot.shape[0] == 0:
                queries_ls.append([])
            else:
                queries_ls.append(tmp_df_bot['query'].to_list())

    cut_sites['observed'] = count_ls
    cut_sites['sam_queries'] = 0
    cut_sites['sam_queries'].astype('object')
    cut_sites['sam_queries'] = queries_ls
    cut_sites = cut_sites[cut_sites['observed']>0]

    return cut_sites


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

    return list(kmers_dt.values())


def write_to_json(filename, obj):
    with open(filename, "w") as outfile:
        json.dump(obj, outfile)


def save_histogram(final_df, title):
    try:
        ax = sns.histplot(data=final_df, x='length', hue='genome',
                          weights=final_df['observed'], multiple="stack",
                          binwidth=8, element="step")
    except IndexError:
        print('singular read lengths, cannot produce histogram')
        return

    old_legend = ax.legend_
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left',
              borderaxespad=0)
    plt.savefig(title, bbox_inches='tight')
    plt.close()


def remove_highly_similar(df, args):
    df_dupes = pd.DataFrame()

    for i in range(1, df['length'].max()+1):
        tmp_df = df[df['length']==i].copy()
        seq_ls = tmp_df['seq'].to_list()
        gen_ls = tmp_df['genome'].to_list()

        hamm_ls = []
        for idxa, (seqa, gena) in enumerate(zip(seq_ls, gen_ls)):
            close_ones = 0
            for idxb, (seqb, genb) in enumerate(zip(seq_ls, gen_ls)):
                if idxa == idxb:
                    pass
                else:
                    hamm = get_hamming_distance(seqa, seqb)
                    if hamm <= args.p:
                        close_ones += 1
            hamm_ls.append(close_ones)
        tmp_df.loc[:,'neighbors'] = hamm_ls
        df_dupes = pd.concat([df_dupes, tmp_df])

    final_df = df_dupes[df_dupes['neighbors'] == 0].copy()

    return final_df


def get_hamming_distance(a, b):
    hamm = 0
    for i, j in zip(a, b):
        if i != j:
            hamm += 1
    return hamm/len(a)


def calculate_cut_efficiency(final_df, args):
    cut_prob_dt = {}
    tmp_df = final_df[final_df['internal']==1].copy()
    for idx, i in tmp_df.iterrows():
        s = i['start']
        e = i['end']
        c = i['chrom']
        g = i['genome']
        l_o = i['length']
        outer = i['observed']
        inner = final_df[(final_df['internal'] == 0) &
                         (final_df['start'] >= s) &
                         (final_df['end'] <= e) &
                         (final_df['chrom'] == c) &
                         (final_df['genome'] == g)].copy()
        if inner.shape[0] == 1:
            l_i = inner['length'].mean()
            inner = inner['observed'].mean()
            ratio = inner/outer
            cut_prob = (ratio - 1)/(ratio)
            if (l_o, l_i) in cut_prob_dt:
                cut_prob_dt[(l_o, l_i)].append(cut_prob)
            else:
                cut_prob_dt[(l_o, l_i)] = [cut_prob]

    cut_prob_ls = []
    for k, v in cut_prob_dt.items():
        if k[0] < 450 and k[1] > 100:
            cut_prob_ls += v

    with open(f'{args.o}_cut_prob.txt', 'w') as o:
        for k, v in cut_prob_dt.items():
            o.write(f'{k} : {v}\n')
        o.write('\n')
        o.write(f'median cut probability: {statistics.median(cut_prob_ls)}')


def write_modified_sam(sam, q_ls):
    '''
    using the limited reads in df, write a sam formatted file to capture only those
    reads that match this simulated sites
    '''
    with open(sam) as f, open(os.path.join(os.path.dirname(sam), f'motif_only_{os.path.basename(sam)}'), 'w') as o:
        for line in f:
            if line.startswith('@'):
                o.write(line)
            elif line.split()[0] in q_ls:
                o.write(line)
            else:
                pass


if __name__ == '__main__':
    main()
