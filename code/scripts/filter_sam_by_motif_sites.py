import numpy as np
import pandas as pd
import sys


def main():
    final_file = sys.argv[1][:-4] + '_motif_matching.sam'
    with open(sys.argv[1]) as f, open(final_file, 'w') as o:
        for idx, line in enumerate(f):
            if line.startswith('@'):
                o.write(line)
            else:
                sam_header = idx
                break

    df = pd.read_csv(sys.argv[1], skiprows=sam_header, usecols=range(11), sep='\s+', header=None)
    cut_sites = pd.read_csv(sys.argv[2])

    cols = 'query flag ref start qual cigar pair pair_pos length seq score'.split()
    df.columns = cols
    df['end'] = df['start'] + df['seq'].str.len() - 1
    df['forward'] = df['flag'].apply(lambda x: check_first_in_pair(x))

    # Keep only the forward reads (starging with EcoRI cut site).
    df_for = df[df['forward'] == True]
    # Create 'motif_pos' column to capture where EcoRI motif actually begins before cut is made.

    '''
    EcoRI

       c r1
       . .
    NNNNGAATTCNNNN
    NNNNCTTAAGNNNN
       .    .
       c    r2
    '''

   # (note: the position in the cut_sites dataframe are based on 0 indexing)

    # Use the 'start' position for these:
    top_strand = df_for[(df_for['length'] > 0)].copy()
    top_strand['motif_pos'] = top_strand['start'] - 2 

    # Use the 'end' position for these:
    bot_strand = df_for[(df_for['length'] < 0)].copy()
    bot_strand['motif_pos'] = bot_strand['end'] - 5 

    df_for = pd.concat([top_strand, bot_strand])
    site_dt = cut_sites[cut_sites['enzyme']==1].groupby('chr')['pos'].apply(list).to_dict()
    df_for['motif_aligned'] = check_in_cut_sites(site_dt, df_for)

    # Keep only the reverse reads  (starging with MseI cut site).
    df_rev = df[df['forward'] == False]

    # Create 'motif_pos' column to capture where EcoRI motif actually begins before cut is made.  

    '''
    MseI

       c r1
       . .
    NNNNTTAANNNN
    NNNNAATTNNNN
       .  .
       c  r2
    '''

   # (note: the position in the cut_sites dataframe are based on 0 indexing)

    # Use the 'start' position for these:
    top_strand = df_rev[(df_rev['length'] > 0)].copy()
    top_strand['motif_pos'] = top_strand['start'] - 2

    # Use the 'end' position for these:
    bot_strand = df_rev[(df_rev['length'] < 0)].copy()
    bot_strand['motif_pos'] = bot_strand['end'] - 3

    df_rev = pd.concat([top_strand, bot_strand])
    site_dt = cut_sites[cut_sites['enzyme']==2].groupby('chr')['pos'].apply(list).to_dict()
    df_rev['motif_aligned'] = check_in_cut_sites(site_dt, df_rev)

    df = pd.concat([df_for, df_rev])
    df = df[df['motif_aligned'] == True]
    df.reset_index(inplace=True, drop=True)
    df['abs_length'] = df['length'].abs()
    df = df[df['abs_length'] < 10000]
    df['counts'] = df.groupby(['query'])['motif_aligned'].transform('count')
    df = df[df['counts'] == 2]
    df_file = sys.argv[1][:-4] + '_motif_matching.tsv'
    df.to_csv(df_file, columns = cols, sep='\t', index=None, header=None)

    with open(df_file) as f, open(final_file, 'a') as o:
        for line in f:
            o.write(line)


def check_first_in_pair(flag):
    return bin(flag)[-7] == '1'


def check_in_cut_sites(site_dt, df_for):
    results = []
    np_for = np.array(df_for)

    for i in np_for:
        chr = i[2]
        try:
            results.append(i[13] in site_dt[chr])
        except KeyError:
            results.append(False)

    return results


if __name__ == '__main__':
    main()
