library(corrplot)
library(ggplot2)
library(plotly)

setwd("~/github/read_simulation/analyses/1520_genomes_RAD_characteristics/2021_06_30_quantify_30_genomes/4_compare_abundance_profiles")

df = read.csv(file='sampled_genome_stats.csv')

df$gc_ratio = df$total_gc / (df$total_all_ls - df$total_n_ls)
df$normalized_reads = df$reads / df$total.length
df$normalized_reads_ratio = df$normalized_reads / sum(df$normalized_reads)

ggplot(df, aes(x=copies, y=gc_ratio)) +
  geom_point(size=2, shape=19, alpha = .3, color="black") +
  theme_classic()

mydata = df[c("copy_ratio", "read_ratio", "normalized_reads_ratio", "gc_ratio")]
mydata.cor = cor(mydata, use = "pairwise.complete.obs")
corrplot(diag = F, method = "pie", mydata.cor)


# after pulling bracken abundance profiles, compress into genus levels
df = read.csv(file='profiled_genome_stats.csv')
df = df[c(26:39)]
df = aggregate(.~bracken_G+bracken_G_percent+bracken_G_reads,data=df,FUN=sum)
df$copy_ratio = df$copies / sum(df$copies)
df$read_ratio = df$reads / sum(df$reads)
df$gc_ratio = df$total_gc / (df$total_all_ls - df$total_n_ls)
df$normalized_reads = df$reads / df$total_all_ls
df$normalized_reads_ratio = df$normalized_reads / sum(df$normalized_reads)


mydata = df[c(2:17)]
mydata.cor = cor(mydata, use = "pairwise.complete.obs")
corrplot(diag = F, method = "color", mydata.cor)

ggplot(df, aes(x=reads, y=copy_ratio)) +
  geom_point(size=3, shape=19, color="black") +
  theme_classic()

ggplot(df, aes(x=bracken_G_reads, y=reads)) +
  geom_point(size=2, shape=19, alpha = .3, color="black") +
  theme_classic() +
  geom_text(aes(label=bracken_G), nudge_x = 12000, nudge_y = 7000, check_overlap = T)

ggplot(df, aes(x=bracken_G_reads, y=gc_ratio)) +
  geom_point(size=2, shape=19, alpha = .3, color="black") +
  theme_classic() +
  geom_text(aes(label=bracken_G))

ggplot(df, aes(x=copies, y=gc_ratio)) +
  geom_point(size=2, shape=19, alpha = .3, color="black") +
  theme_classic() +
  geom_text(aes(label=bracken_G))

fig = plot_ly(x=df$reads, y=df$bracken_G_reads, z=df$gc_ratio, type="scatter3d", mode="markers", color=df$gc_ratio)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'expected reads'),
                                   yaxis = list(title = 'bracken reads'),
                                   zaxis = list(title = 'gc ratio')))
fig




