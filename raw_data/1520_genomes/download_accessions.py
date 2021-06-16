#!/usr/bin/env python

import os
import pandas as pd
import subprocess
import sys

accesions = sys.argv[1]

column_names = ['accession_id',
                'assembly',
                'level',
                'biosample',
                'strain',
                'taxonomy',
                'taxid',
                'coverage',
                'platform',
                'total-length',
                'spanned-gaps',
                'unspanned-gaps',
                'region-count',
                'scaffold-count',
                'scaffold-N50',
                'scaffold-L50',
                'scaffold-N75',
                'scaffold-N90',
                'contig-count',
                'contig-N50',
                'contig-L50',
                'total-gap-length',
                'molecule-count',
                'top-level-count',
                'component-count']

df = pd.DataFrame([], columns=column_names)


def get_stats_file(id):
    '''
    this is totally hacky, but NCBI links aren't named consistently
    '''

    stats = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/' + \
    id[4:7] + '/' + id[7:10] + '/' + id[10:13] + '/' + \
    'GCA_' + id[4:] + '_ASM' + id[6:12] + \
    'v1/GCA_' + id[4:] + '_ASM' + id[6:12] + 'v1_assembly_stats.txt'

    stats_file = stats[-46:]

    if id == 'GCA_003439895.1':
        stats = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/439/895/GCA_003439895.1_OM02-37/GCA_003439895.1_OM02-37_assembly_stats.txt'
        stats_file = 'GCA_003439895.1_OM02-37_assembly_stats.txt'

    if not os.path.isfile(os.path.join('./assembly_stats', stats_file)):
        subprocess.call(['wget', stats, '-P', './assembly_stats'])

    return stats_file


def get_assembly_stats(stats_file, stats_dt):
    with open(os.path.join('./assembly_stats', stats_file)) as f2:
        for line in f2:
            if line.startswith('# Taxid:'):
                stats_dt['taxid'] = (line.rstrip().split()[2])
            if line.startswith('# Genome coverage:'):
                stats_dt['coverage'] = (line.rstrip().split()[3])
            if line.startswith('# Sequencing technology:'):
                stats_dt['platform'] = (' '.join(line.rstrip().split()[3:]))
            if line.startswith('all'):
                line = line.rstrip().split()
                stats_dt[line[4]] = line[5]

    return stats_dt


with open(accesions) as f:
    for i, line in enumerate(f):
        stats_dt = {}
        if i > 1:
            fields = line.rstrip().split()
            stats_dt['accession_id'] = fields[0]
            stats_dt['assembly'] = fields[1]
            stats_dt['level'] = fields[2]
            stats_dt['biosample'] = fields[3]
            stats_dt['strain'] = fields[4]
            stats_dt['taxonomy'] = ' '.join(fields[5:])
            id = fields[0]
            stats_file = get_stats_file(id)
            stats_dt = get_assembly_stats(stats_file, stats_dt)
            df = df.append(stats_dt, ignore_index = True)

    print(df)
    df.to_csv('genome_stats.csv', index=False)
