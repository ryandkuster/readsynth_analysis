import pandas as pd
import sys

'''
1 is mean
2 is sd
3 is counts
'''
mean_ls = []
sd_ls = []
count_ls = []

with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2, open(sys.argv[3]) as f3:
    for i, j, k in zip(f1, f2, f3):
        mean_ls.append(float(i.split()[4]))
        sd_ls.append(float(j.split()[5]))
        count_ls.append(float(k.split()[3]))

df = pd.DataFrame({'mean': mean_ls, 'sd': sd_ls, 'count': count_ls})
print(df)

df['mean_w'] = df['mean'] * df['count']
df['sd_w'] = df['sd'] * df['count']
print(df['mean_w'].sum() / df['count'].sum())
print(df['sd_w'].sum() / df['count'].sum())
