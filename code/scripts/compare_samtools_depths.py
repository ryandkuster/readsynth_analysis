import argparse
import pandas as pd
import sys


def parse_user_input():
    parser = argparse.ArgumentParser(description='compare samtools depths')

    parser.add_argument('-d', type=str, required=True, nargs='+',
                        help='list of depth files to combine')

    parser.add_argument('-n', type=str, required=True, nargs='+',
                        help='for each file, give a unique name for column header')

    args = parser.parse_args()

    return args


def get_first_df(args):
    df = pd.read_csv(args.d[0], sep='\s+')
    df.rename(columns={df.columns[2] : args.n[0]}, inplace=True)
    return df


def combine_depths(args):
    for depth, name in zip(args.d[1:], args.n[1:]):
        print(f'adding {name}')
        tmp_df = pd.read_csv(depth, sep='\s+')
        tmp_df.rename(columns={tmp_df.columns[2] : name}, inplace=True)
        args.df = args.df.merge(tmp_df, how='outer', on=['#CHROM', 'POS'])
    return args.df


def bc_distance(args):
    args.df['abs_diff'] = abs(args.df[args.n[0]] - args.df[args.n[1]])
    args.df['sum'] = args.df[args.n[0]] + args.df[args.n[1]]
    print(args.df['abs_diff'].sum()/args.df['sum'].sum())


if __name__ == '__main__':
    args = parse_user_input()
    if len(args.d) != len(args.n):
        sys.exit('inputs -d and -n differ in length')
    args.df = get_first_df(args)
    args.df = combine_depths(args)
    args.df.drop('#CHROM', axis=1, inplace=True)
    args.df.drop('POS', axis=1, inplace=True)
    args.df.fillna(0, inplace=True)

    print('calculating pearson')
    pearson_df = args.df.corr(method='pearson')
    pearson_df.to_csv('pearson_corr.csv')

    print('calculating spearman')
    spearman_df = args.df.corr(method='spearman')
    spearman_df.to_csv('spearman_corr.csv')

    if len(args.d) == 2:
        bc_distance(args)
    


