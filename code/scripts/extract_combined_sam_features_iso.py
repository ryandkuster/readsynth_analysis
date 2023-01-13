import argparse
import itertools
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import seaborn as sns
import statistics
import sys


'''
this script is meant to find the chromosomes simulated for iso-length enzymes only
fragment recreation is not performed (as fragments should be identical in length)
-o is the name of the output csv to write
'''


def main():
    args = parse_user_input()
    df = process_sam(args)
    final_df = pd.DataFrame()
    key = pd.read_csv(args.i, header=None)
    key_sum = sum(key[0].to_list())

    for label, cut_sites in zip(key[0].to_list(), key[1].to_list()):
        if cut_sites.endswith('.fna.csv'):
            genome_name = cut_sites.split('/')[-1][11:-8]
        elif cut_sites.endswith('.fna.gz.csv'):
            genome_name = cut_sites.split('/')[-1][11:-11]
        else:
            sys.exit(f'genome name suffix not known')
        print(genome_name)

        try:
            cut_sites = pd.read_csv(cut_sites)
        except FileNotFoundError:
            print(f'{cut_sites} not found, continuing')
            continue

        unique_seq = list(cut_sites['chrom'].unique())
        if df is False:
            continue
        subset_df = df[df['ref_x'].isin(unique_seq)]
        subset_df = pd.DataFrame({'observed' : subset_df.groupby(["seq_x"]).size()}).reset_index()
        subset_df['copy_number'] = label
        subset_df['rel_abund'] = label/key_sum
        subset_df['genome'] = genome_name
        final_df = pd.concat([final_df, subset_df], axis=0)

    print('sam files processed')

    final_df.reset_index(inplace=True)
    final_df.to_csv(f'{args.o}.csv',index=None)


def parse_user_input():
    parser = argparse.ArgumentParser(description='simulate RAD libary')

    parser.add_argument('-i', type=str, required=True,
                        help='path to key csv (sam,abundance,raw_digest')

    parser.add_argument('-o', type=str, required=True,
                        help='prefix name of files to store output')

    parser.add_argument('-sam', type=str, required=True,
                        help='combined sam file')

    args = parser.parse_args()

    return args


def process_sam(args):
    print('processing combined sam file')
    #skip the header rows we just learned from the loop above and open the file
    cols = 'query flag ref start qual cigar pair pair_pos length seq score d12 d13 d14 d15 d16 d17 d18 d19 d20'.split()
    df = pd.read_csv(args.sam, names=cols, sep='\s+', header=None)
    cols = 'query flag ref start qual cigar pair pair_pos length seq score'.split()
    df = df[cols]

    #df = pd.read_csv(args.sam, usecols=range(11), sep='\s+', header=None)
    #cols = 'query flag ref start qual cigar pair pair_pos length seq score'.split()
    #df.columns = cols

    # drop reads with 0 scores (e.g. multiple-mappings)
    df = df[df['qual'] > 0]

    # apply boolean based on bitwise flag to reflect if first/second in pair
    df['forward'] = df['flag'].apply(lambda x: check_first_in_pair(x))

    # drop the unneeded sam columns
    df = df.drop(['flag', 'qual', 'pair', 'pair_pos', 'score'], axis=1)

    # create two separate dataframes for the first/second (R1/R2) reads
    df_for = df[df['forward'] == True].copy()
    df_rev = df[df['forward'] == False].copy()

    # merge the first/second read dataframes together so each row is a pair
    df = df_for.merge(df_rev, on="query")

    return df


def check_first_in_pair(flag):
    '''
    check sam flag to see if read is first in pair (returns true)
    '''
    return bin(flag)[-7] == '1'


def count_matches(cigar):
    '''
    check cigar for numbers preceding M, sum them
    '''
    count_m = 0
    m_split = [i for i in cigar.split('M') if i != '']

    for i in m_split:
        if i[-1].isnumeric():
            count_m += int(re.findall(r"[^\W\d_]+|\d+", i)[-1])
    return count_m


def bad_cigar(cigar):
    '''
    check cigar for characters other than M or S
    '''
    status = 0
    #fail_ls = ['I', 'D', 'N', 'H', 'X']
    fail_ls = ['I', 'D']
    for i in cigar:
        if i in fail_ls:
            status = 1
    return status


def count_fragments_nosam(df, cut_sites, args):
    '''
    for each simulated fragment, find the number of occurances from the sam file
    '''
    cut_sites.reset_index(inplace=True, drop=True)
    count_ls = []
    np_cut_sites = np.array(cut_sites)
    for i in np_cut_sites:
        chrom = i[0]
        c_start = i[2]
        c_end = i[3]
        c_m1 = i[4]
        c_m2 = i[5]
        if c_m1 == args.m1:
            tmp_df_top = df[(df['m1_pos']==c_start) & (df['m2_pos']==c_end) & (df['ref_x']==chrom)]
            count_ls.append(tmp_df_top.shape[0])
        elif c_m1 == args.m2:
            tmp_df_bot = df[(df['m1_pos']==c_end) & (df['m2_pos']==c_start) & (df['ref_x']==chrom)]
            count_ls.append(tmp_df_bot.shape[0])

    cut_sites['observed'] = count_ls
    cut_sites = cut_sites[cut_sites['observed']>0]

    return cut_sites


if __name__ == '__main__':
    main()
