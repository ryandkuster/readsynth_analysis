import os
import sys


"""

sys.argv[1] is a newline-separated file containing only the
forward reads from a sorted sam file.

sys.argv[2] is a newline-separated file containing only the
forward reads from a sorted sam file.

sys.argv[3] is the k-mer size

reads in the sam file should already conform to the '+' template
strand

analyses and python code taken from:
https://sourmash.readthedocs.io/en/latest/kmers-and-minhash.html

"""

def jaccard_similarity(a, b):
    a = set(a)
    b = set(b)

    intersection = len(a.intersection(b))
    union = len(a.union(b))

    return intersection / union


def jaccard_similarity_duplicates(a, b):
    intersection = 0
    for k, v in a.items():
        if k in b:
            intersection += min(a[k], b[k])
    union = sum(a.values()) + sum(b.values())

    return intersection / union


def jaccard_containment(a, b):
    a = set(a)
    b = set(b)

    intersection = len(a.intersection(b))

    return intersection / len(a)


def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers


def kmer_dictionary(all_kmers, kmers):
    for i in kmers:
        if i in all_kmers:
            all_kmers[i] += 1
        else:
            all_kmers[i] = 1

    return all_kmers


def read_kmers_from_file(filename, ksize):
    all_kmers = {}
    with open(filename) as f:
        for sequence in f:

            kmers = build_kmers(sequence.rstrip(), ksize)
            all_kmers = kmer_dictionary(all_kmers, kmers)

    return all_kmers


#def read_kmers_from_file(filename, ksize):
#    all_kmers = []
#    with open(filename) as f:
#        for sequence in f:
#
#            kmers = build_kmers(sequence.rstrip(), ksize)
#            all_kmers += kmers
#
#    return all_kmers


if __name__ == '__main__':
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    ksize = int(sys.argv[3])
    print('getting kmers for first file')
    all_kmers_1 = read_kmers_from_file(file1, ksize)
    print('getting kmers for second file')
    all_kmers_2 = read_kmers_from_file(file2, ksize)
    #print(os.path.basename(file1) + ' vs ' + os.path.basename(file2), jaccard_similarity(all_kmers_1, all_kmers_2))
    print(os.path.basename(file1) + ' vs ' + os.path.basename(file2), jaccard_similarity_duplicates(all_kmers_1, all_kmers_2))

