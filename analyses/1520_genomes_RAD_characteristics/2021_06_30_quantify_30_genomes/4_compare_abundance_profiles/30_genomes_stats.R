library(corrplot)
library(ggplot2)

setwd("~/github/read_simulation/analyses/1520_genomes_RAD_characteristics/2021_06_30_quantify_30_genomes/4_compare_abundance_profiles")

df = read.csv(file='sampled_genome_stats.csv')

df$gc_ratio = df$total_gc / (df$total_all_ls - df$total_n_ls)
df$normalized_reads = df$reads / df$total.length
df$normalized_reads_ratio = df$normalized_reads / sum(df$normalized_reads)

mydata = df[c("copy_ratio", "read_ratio", "normalized_reads_ratio", "gc_ratio")]
mydata.cor = cor(mydata, use = "pairwise.complete.obs")
corrplot(diag = F, method = "pie", mydata.cor)

