#!/usr/bin/env python

import os
import subprocess
import sys

accessions = sys.argv[1]


def get_ftp_link(id):
    '''
    this is totally hacky, but NCBI links aren't named consistently
    '''
    genome = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/' + \
    id[4:7] + '/' + id[7:10] + '/' + id[10:13] + '/' + \
    'GCA_' + id[4:] + '_ASM' + id[6:12] + \
    'v1/GCA_' + id[4:] + '_ASM' + id[6:12] + 'v1_genomic.fna.gz'
    genome_file = genome[-42:]
    md5 = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/' + \
    id[4:7] + '/' + id[7:10] + '/' + id[10:13] + '/' + \
    'GCA_' + id[4:] + '_ASM' + id[6:12] + \
    'v1/md5checksums.txt'

    if id == 'GCA_003439895.1':
        genome = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/439/895/GCA_003439895.1_OM02-37/GCA_003439895.1_OM02-37_genomic.fna.gz'
        genome_file = 'GCA_003439895.1_OM02-37_genomic.fna.gz'
        md5 = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/439/895/GCA_003439895.1_OM02-37/md5checksums.txt'

    if not os.path.isfile(os.path.join('./genomes', genome_file)):
        subprocess.call(['wget', genome, '-P', './genomes'])
        subprocess.call(['wget', md5, '-P', './genomes'])

    return genome_file


def master_md5_file(genome_file, o):
    with open(os.path.join('./genomes', 'md5checksums.txt')) as f2:
        for line in f2:
            if genome_file in line:
                o.write(line)
    subprocess.call(['rm', os.path.join('./genomes', 'md5checksums.txt')])


with open(accessions) as f, \
     open(os.path.join('./genomes', 'all_md5.txt'), 'a') as o:
    for i, line in enumerate(f):
        if i > 1:
            fields = line.rstrip().split()
            id = fields[0]
            genome_file = get_ftp_link(id)
            master_md5_file(genome_file, o)

