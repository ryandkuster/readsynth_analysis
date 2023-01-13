import argparse
import Levenshtein
import pandas as pd


def main():
    args = parse_user_input()
    df = pd.read_csv(args.i)
    print(df)


def parse_user_input():
    parser = argparse.ArgumentParser(description='simulate RAD libary')

    parser.add_argument('-i', type=str, required=True,
                        help='path to key csv (sam,abundance,raw_digest')

    parser.add_argument('-p', type=float, required=False,
                        help='percent of read that must have unique hamming distance')

    args = parser.parse_args()
 
    return args


def remove_highly_similar(df, args):
    gen_ls = list(df['genome'].unique())
    remove_seqs = []

    for idx, gen in enumerate(gen_ls[:-1]):
        print(f'processing {gen}')
        query_ls = df[df['genome']==gen]['seq'].to_list()
        query_ls = list(set(query_ls))
        ref_ls = df[(df['genome'].isin(gen_ls[idx+1:]))]['seq'].to_list()
        ref_ls = list(set(ref_ls))
        seq_dt = make_seq_dict(ref_ls)
        for seq in query_ls:
            remove_seqs = get_lev_ratio_dict(seq, seq_dt, remove_seqs, args)

    remove_seqs = list(set(remove_seqs))
    df = df[~df['seq'].isin(remove_seqs)]
    return df


def make_seq_dict(ref_ls):
    seq_dt = {}
    for seq in ref_ls:
        if len(seq) in seq_dt:
            seq_dt[len(seq)].append(seq)
        else:
            seq_dt[len(seq)] = [seq]
    return seq_dt


def get_lev_ratio_dict(seq, seq_dt, remove_seqs, args):
    lower, upper = get_lev_ratio_limits(seq, args)
    for seq_len in range(int(round(lower)), int(round(upper)+1)):
        try:
            for ref_seq in seq_dt[seq_len]:
                if Levenshtein.ratio(seq, ref_seq) > args.p:
                    remove_seqs += [seq, ref_seq]
        except KeyError:
            pass
    return remove_seqs


def get_lev_ratio_limits(seq, args):
    lower = (args.p*len(seq))/((1-args.p)+1)
    upper = (len(seq)+((1-args.p)*len(seq)))/args.p
    return lower, upper


def get_levenshtein_ratio(a, b):
    return Levenshtein.ratio(a, b)


if __name__ == '__main__':
    main()
