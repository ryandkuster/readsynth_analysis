import os
import sys

"""

produces a fastq file with subset of the kraken2 output containing the
target taxids only

sys.argv[1] is the resulting kraken2 reads file (e.g. cseqs_1.fq)
sys.argv[2] is the desired taxid to subset

traverse fastq output from Kraken2 with header formatted as follows:
@SRR5298272.7 length=76 kraken:taxid|657322

"""

krak_f = sys.argv[1]
target = sys.argv[2]
krak_o = 'taxid_' + target + '_' + os.path.basename(krak_f)

with open(krak_f) as f, open(krak_o, 'w') as o:
    for idx, line in enumerate(f):
        if idx % 4 == 0:
            new_read = ''
            taxid = line.rstrip().split('kraken:taxid|')[1]
        new_read += line
        if taxid == target and (idx+1) % 4 == 0:
            o.write(new_read)

