import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys

'''
sys.argv[1] is sam alignment summary
sys.argv[2] is cut sites summary
'''

def main():
    # open sam summary and enzyme cut sites summary files
    df = pd.read_csv(sys.argv[1], sep='\t')
    cut_sites = pd.read_csv(sys.argv[2])

    # Use sam flags to identify forward vs. reverse reads.
    df['forward'] = df['flag'].apply(lambda x: check_first_in_pair(x))


    '''
    EcoRI

       c r1
       . .
    NNNNGAATTCNNNN
    NNNNCTTAAGNNNN
       .    .
       c    r2
    '''

    # Keep only the forward reads (starging with EcoRI cut site).
    df_for = df[df['forward'] == True]

    # Create 'motif_pos' column to capture where EcoRI motif actually begins before cut is made.
    # (note: the position in the cut_sites dataframe are based on 0 indexing)

    # Use the 'start' position for these:
    top_strand = df_for[(df_for['length'] > 0)].copy()
    top_strand['motif_pos'] = top_strand['start'] - 2

    # Use the 'end' position for these:
    bot_strand = df_for[(df_for['length'] < 0)].copy()
    bot_strand['motif_pos'] = bot_strand['end'] - 5

    df_for = pd.concat([top_strand, bot_strand])
    df_for['motif_aligned'] = check_in_cut_sites(cut_sites[cut_sites['enzyme']==1], df_for)

    '''
    MseI

       c r1
       . .
    NNNNTTAANNNN
    NNNNAATTNNNN
       .  .
       c  r2
    '''

    # Keep only the reverse reads  (starging with MseI cut site).
    df_rev = df[df['forward'] == False]

    # Create 'motif_pos' column to capture where EcoRI motif actually begins before cut is made.  
    # (note: the position in the cut_sites dataframe are based on 0 indexing)

    # Use the 'start' position for these:
    top_strand = df_rev[(df_rev['length'] > 0)].copy()
    top_strand['motif_pos'] = top_strand['start'] - 2

    # Use the 'end' position for these:
    bot_strand = df_rev[(df_rev['length'] < 0)].copy()
    bot_strand['motif_pos'] = bot_strand['end'] - 3

    df_rev = pd.concat([top_strand, bot_strand])
    df_rev['motif_aligned'] = check_in_cut_sites(cut_sites[cut_sites['enzyme']==2], df_rev)

    df = pd.concat([df_for, df_rev])
    df = df[df['motif_aligned'] == True]
    df.reset_index(inplace=True, drop=True)
    df['abs_length'] = df['length'].abs()
    df = df[df['abs_length'] < 10000]
    df['counts'] = df.groupby(['query'])['motif_aligned'].transform('count')
    df = df[df['counts'] == 2]



    df.sort_values(by=['motif_pos'], ascending = (True), inplace=True)
    df_for = df[df['forward'] == True]
    df_rev = df[df['forward'] == False]
    df = df_for.merge(df_rev, on='query', how='left')
    df.reset_index(inplace=True, drop=True)

    final_dt = df['abs_length_x'].value_counts().to_dict()
    write_to_json('fragments_' + os.path.basename(sys.argv[1]) + '.json', final_dt)

    plt.figure(figsize=(10,8))
    plt.hist(final_dt.keys(), weights=final_dt.values(), bins=100)
    plt.savefig('hist_' + os.path.basename(sys.argv[1]) + '.png')
    plt.close()

    df = df[['query', 'abs_length_x', 'motif_pos_x', 'motif_pos_y', 'length_x']]
    np_array = np.array(df)
    cut_loci = {}

    for i in range(np_array.shape[0]):
        if np_array[i,:][4] > 0:
            tmp_tup = (np_array[i,:][2], np_array[i,:][3])
        else:
            tmp_tup = (np_array[i,:][3], np_array[i,:][2])
        if tmp_tup in cut_loci:
            cut_loci[tmp_tup] += 1
        else:
            cut_loci[tmp_tup] = 1

    all_cut_sites = []

    for i in cut_loci.keys():
        if i[0] not in all_cut_sites:
            all_cut_sites.append(i[0])
        if i[1] not in all_cut_sites:
            all_cut_sites.append(i[1])

    all_cut_sites.sort()

    fragments = []

    for i, count in cut_loci.items():
        try:
            start = all_cut_sites.index(i[0])
            end = all_cut_sites.index(i[1])
            fragments.append([i[1]-i[0], end - start - 1, count])
        except ValueError:
            pass

    fragments = np.array(fragments)
    max_internals = np.max(fragments[:,1])
    last_target, last_internals, last_count = (0, 0), 0, 0
    ratios = []

    for idx, i in enumerate(all_cut_sites):
        for j in range(max_internals+1):
            try:
                target = (all_cut_sites[idx], all_cut_sites[idx+j+1])
            except IndexError:
                break
            if target in cut_loci:
                if target[0] == last_target[0]:
                    if j - last_internals == 1:
                        ratios.append([last_count, cut_loci[target]])
                last_internals = j
                last_target = target
                last_count = cut_loci[target]

    ratios = pd.DataFrame(np.array(ratios), columns=['inner', 'outer'])
    ratios['ratio'] = ratios['inner']/ratios['outer']

    prob = get_prob(ratios['ratio'].mean())
    write_to_json('prob_' + os.path.basename(sys.argv[1]) + '.json', prob)

def check_first_in_pair(flag):
    return bin(flag)[-7] == '1'


def check_in_cut_sites(cut_sites, df_for):
    site_dt = cut_sites.groupby('chr')['pos'].apply(list).to_dict()
    results = []
    np_for = np.array(df_for)

    for i in np_for:
        chr = i[2]
        try:
            results.append(i[9] in site_dt[chr])
        except KeyError:
            results.append(False)

    return results


def get_prob(ratio):
    return (1-ratio) / -ratio


def write_to_json(filename, obj):
    with open(filename, "w") as outfile:
        json.dump(obj, outfile)


if __name__ == '__main__':
    main()
