import itertools
import pandas as pd
import sys


'''
sys.argv[1] is the file containing sequences to open
sys.argv[2] is the kmer size
sys.argv[3] is the name of the output file
'''


def kmer_combos(k):
    kmers = []
    #k_list = ['ACGT' for i in range(k)]
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

print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])

#df = pd.read_csv(sys.argv[1], index_col=0)
df = pd.read_csv(sys.argv[1])
k = int(sys.argv[2])
kmers = kmer_combos(k)
df[kmers] = 0
df[kmers] = df['seq'].apply(lambda x: pd.Series(count_kmers(x, k, kmers)))
df.to_csv(sys.argv[3])

