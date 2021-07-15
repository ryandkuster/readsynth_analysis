import sys

"""

sys.argv[1] is the R2 file with naming discrepencies
due to sequential numbering in original SRA file that could not be
split using the fastq-dump command

e.g.:
@SRR5298272.2 length=75
CCTTCTTCTCTCCTTTCCCCCATTTCTCCTCCTCCTTCCCCCTCCTTCTCTCCTCCACCCTTCCCCTTTTCTTTC
+SRR5298272.2 length=75
8,8-8;;,;;,,,,,66;,,++,;<,;,,;;;;,,,;,,,;77,86;;,,,,,;,,,;686,;6,,6,;6;,,,,

header discrepency blocks bwa mem... so here this script will replace
each sequential number with -1

"""

with open(sys.argv[1]) as f, open('modified_' + sys.argv[1], 'w') as o:
    for idx, line in enumerate(f):
        if idx % 2 == 0:
            line_ls = line.split(' ')
            line_ls_0 = line_ls[0].split('.')
            line_ls[0] = line_ls_0[0] + '.' + str(int(line_ls_0[1])-1)
            line = ' '.join(line_ls)
        o.write(line)
