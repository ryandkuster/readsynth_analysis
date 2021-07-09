import pandas as pd
import sys


"""
sys.argv[1] is the kracken output (e.g. 'cseqs.fq')
"""

df = pd.DataFrame()

kraken_dt = {}
with open(sys.argv[1]) as f:
    line_no = 0
    for line in f:
        line_no += 1
        if line_no == 1:
            taxid = line.rstrip().split('kraken:taxid|')[1]
            if taxid in kraken_dt:
                kraken_dt[taxid] += 1
            else: 
                kraken_dt[taxid] = 1
        if line_no == 4:
            line_no = 0

id_ls = []
count_ls = []

for k, v in kraken_dt.items():
    id_ls.append(k)
    count_ls.append(v)

df['kraken_taxid'] = id_ls
df['count'] = count_ls

df = df.sort_values('count', ascending=False)
df.to_csv('kraken_' + sys.argv[1] + '_counts.csv', index=None)            
