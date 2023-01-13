import os
import pandas as pd
import subprocess
import sys


'''
this script requires an up-to-date 'assembly_summary_genbank.txt' file as sys.argv[1]

wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

sys.argv[2] is a tab_separated file from bracken output
sys.argv[3] is the location to store the downloaded genome
'''


def main():
    outfile = open('summary.txt', 'w')
    failfile = open('fail_summary.txt', 'w')
    outfile.write(f'ncbi_name\tname\tftp\tmean\tcat\tlev\trelease\trep_gen\tdate\tbest_assembly\tunclassified\n')
    if not os.path.exists(sys.argv[3]):
        sys.exit('outdir does not exist')
    tax_df = open_taxids()
    df = pd.read_csv(sys.argv[1], skiprows=1, sep='\t')
    df['seq_rel_date']= pd.to_datetime(df['seq_rel_date'])
    df['organism_name'] = df['organism_name'].str.lower()
    accession_ls = []
    u_df = pd.DataFrame({'genus': [], 'species': [], 'mean': []})

    for name, binom in tax_df.groupby(['genus', 'species']):
        tmp_tax_df = tax_df[(tax_df['genus']==name[0]) & (tax_df['species']==name[1])]
        targets = tmp_tax_df.shape[0] 
        g, s, unclassified = check_unclassified(name[0], name[1])
        tmp_tax_df.loc[:,'genus'] = g
        tmp_tax_df.loc[:,'species'] = s

        if unclassified:
            u_df = pd.concat([u_df, tmp_tax_df])
            continue
        else:
            tmp_df = check_refseq_category_s(df, g, s, accession_ls, targets)
            if tmp_df.shape[0] >= targets:
                tmp_df = choose_best_assemblies(tmp_df, targets)
                accession_ls += tmp_df['# assembly_accession'].to_list()
                for idx, m in enumerate(tmp_tax_df['mean'].to_list()):
                    best_assembly = tmp_df['ftp_path'].to_list()[idx]
                    best_assembly = best_assembly + '/' + best_assembly.split('/')[-1] + '_genomic.fna.gz'
                    org = tmp_df['organism_name'].to_list()[idx].replace(' ', '_')
                    cat = tmp_df['refseq_category'].to_list()[idx]
                    lev = tmp_df['assembly_level'].to_list()[idx]
                    rel = tmp_df['release_type'].to_list()[idx]
                    rep = tmp_df['genome_rep'].to_list()[idx]
                    date = tmp_df['seq_rel_date'].to_list()[idx]
                    outfile.write(f'{org}\t{g}_{s}\t{best_assembly}\t{m}\t{cat}\t{lev}\t{rel}\t{rep}\t{date}\t{best_assembly}\t{unclassified}\n')
                    subprocess.call(['wget', best_assembly, '-P', './genomes'])
            else:
                u_df = pd.concat([u_df, tmp_tax_df])
                continue
    
    print('\nprocessing unclassified...\n')
    for g in u_df['genus'].unique():
        targets = u_df[u_df['genus']==g].shape[0]
        tmp_df = check_refseq_category_g(df, g, accession_ls, targets)
        if tmp_df.shape[0] > 0:
            tmp_df = choose_best_assemblies(tmp_df, targets)
            accession_ls += tmp_df['# assembly_accession'].to_list()
            tmp_u_df = u_df[u_df['genus']==g]
            for idx, (s, m) in enumerate(zip(tmp_u_df['species'].to_list(), tmp_u_df['mean'].to_list())):
                try:
                    best_assembly = tmp_df['ftp_path'].to_list()[idx]
                    best_assembly = best_assembly + '/' + best_assembly.split('/')[-1] + '_genomic.fna.gz'
                    org = tmp_df['organism_name'].to_list()[idx].replace(' ', '_')
                    cat = tmp_df['refseq_category'].to_list()[idx]
                    lev = tmp_df['assembly_level'].to_list()[idx]
                    rel = tmp_df['release_type'].to_list()[idx]
                    rep = tmp_df['genome_rep'].to_list()[idx]
                    date = tmp_df['seq_rel_date'].to_list()[idx]
                    outfile.write(f'{org}\t{g}_{s}\t{best_assembly}\t{m}\t{cat}\t{lev}\t{rel}\t{rep}\t{date}\t{best_assembly}\t{unclassified}\n')
                    subprocess.call(['wget', best_assembly, '-P', './genomes'])
                except IndexError:
                    failfile.write(f'na\t{s}_{g}\tna\t{m}\tna\tna\tna\tna\tna\tna\tna\n')
        else:
            for m in u_df[u_df['genus']==g]['mean'].to_list():
                failfile.write(f'na\t{s}_{g}\tna\t{m}\tna\tna\tna\tna\tna\tna\tna\n')
    outfile.close()


def open_taxids():
    tax_df = pd.read_csv(sys.argv[2])
    tax_df = tax_df[['genus', 'species', 'mean']]
    tax_df['genus'] = tax_df['genus'].str.lower()
    tax_df['species'] = tax_df['species'].str.lower()
    return tax_df


def check_unclassified(g, s):
    unclassified = False
    if 'unclassified' in g:
        g = g.split('_')[0]
        unclassified = True
    if 'unclassified' in s:
        s = s.split('_')[0]
        unclassified = True
    return g, s, unclassified


def check_refseq_category_s(df, g, s, accession_ls, targets):
    tmp_df = pd.DataFrame()
    refseq_ls = ['representative genome', 'reference genome']
    assembly_ls = ['Complete Genome', 'Chromosome', 'Scaffold']

    for i in refseq_ls:
        tmp_tmp_df = df.loc[(df['organism_name'].str.contains(s)) & (df['organism_name'].str.contains(g)) & (df['refseq_category'] == i)]
        tmp_df = pd.concat([tmp_df, tmp_tmp_df])
        tmp_df = tmp_df[~tmp_df['# assembly_accession'].isin(accession_ls)]
        if tmp_df.shape[0] >= targets:
            return tmp_df

    for i in assembly_ls:
        tmp_tmp_df = df.loc[(df['organism_name'].str.contains(s)) & (df['organism_name'].str.contains(g)) & (df['assembly_level'] == i)]
        tmp_df = pd.concat([tmp_df, tmp_tmp_df])
        tmp_df.drop_duplicates(inplace=True)
        tmp_df = tmp_df[~tmp_df['# assembly_accession'].isin(accession_ls)]
        if tmp_df.shape[0] >= targets:
            return tmp_df
    return tmp_df 


def check_refseq_category_g(df, g, accession_ls, targets):
    tmp_df = pd.DataFrame()
    refseq_ls = ['representative genome', 'reference genome']
    assembly_ls = ['Complete Genome', 'Chromosome', 'Scaffold']

    for i in refseq_ls:
        tmp_tmp_df = df.loc[(df['organism_name'].str.startswith(g)) & (df['refseq_category'] == i)]
        tmp_df = pd.concat([tmp_df, tmp_tmp_df])
        tmp_df = tmp_df[~tmp_df['# assembly_accession'].isin(accession_ls)]
        if tmp_df.shape[0] >= targets:
            return tmp_df

    for i in assembly_ls:
        tmp_tmp_df = df.loc[(df['organism_name'].str.startswith(g)) & (df['assembly_level'] == i)]
        tmp_df = pd.concat([tmp_df, tmp_tmp_df])
        tmp_df.drop_duplicates(inplace=True)
        tmp_df = tmp_df[~tmp_df['# assembly_accession'].isin(accession_ls)]
        if tmp_df.shape[0] >= targets:
            return tmp_df
    return tmp_df 


def choose_best_assembly(tmp_df):
    '''
    choose assemblies with:
        major release type
        full genome rep
        newest release date
    '''
    tmp_df = tmp_df.loc[(tmp_df['release_type'] == 'Major') & \
                               (tmp_df['genome_rep'] == 'Full')].copy()
    tmp_df.sort_values(by='seq_rel_date', inplace=True)
    tmp_df = tmp_df[-1:]
    return tmp_df


def choose_best_assemblies(tmp_df, targets):
    '''
    choose assemblies with:
        major release type
        full genome rep
        newest release date
    '''
    tmp_df = tmp_df.loc[(tmp_df['release_type'] == 'Major') & \
                               (tmp_df['genome_rep'] == 'Full')].copy()
    tmp_df.sort_values(by='seq_rel_date', inplace=True)
    tmp_df = tmp_df[-targets:]
    return tmp_df


if __name__ == '__main__':
    main()
