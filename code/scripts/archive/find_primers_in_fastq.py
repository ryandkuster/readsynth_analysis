import sys


'''
sys.argv[1] is the fastq file to assess
sys.argv[2] is the kraken2 taxid expected to end the header line
sys.argv[3] is the length of the sequence to assess for primer content
'''

query = sys.argv[2]
primer_len = int(sys.argv[3])
primer_dt = {}


def update_dt(primer_dt, seq):
    if seq in primer_dt:
        primer_dt[seq] = primer_dt[seq] + 1
    else:
        primer_dt[seq] = 1

    return primer_dt


with open(sys.argv[1]) as f:
    line_no = 0
    for line in f:
        line_no += 1
        if line_no == 1:
            header = line.rstrip()
        if header.endswith(query) and line_no == 2:
            seq = line[:primer_len+1]
            primer_dt = update_dt(primer_dt, seq)
        if line_no == 4:
            line_no = 0


for k, v in primer_dt.items():
    if v > 10:
        print(k)
