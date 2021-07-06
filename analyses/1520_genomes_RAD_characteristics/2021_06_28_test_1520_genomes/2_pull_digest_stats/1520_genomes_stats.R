library(corrplot)
library(ggplot2)

df = read.csv(file='digested_genome_stats.csv')
df$gc_ratio = df$total_gc / (df$total_all_ls - df$total_n_ls)
df$normalized_frags_total = df$frags_total / df$total.length
df$normalized_frags_1_to_200 = df$frags_1_to_200 / df$total.length
df$normalized_frags_201_to_400 = df$frags_201_to_400 / df$total.length
df$normalized_frags_401_to_600 = df$frags_401_to_600 / df$total.length

# produce correlation object to visualize all relations in dataframe
mydata = df[c(10:37)]

# remove columns with missing data
mydata = mydata[ , -which(names(mydata) %in% c("unspanned.gaps","region.count","molecule.count"))]
mydata.cor = cor(mydata, use = "pairwise.complete.obs")
corrplot(mydata.cor)

ggplot(mydata, aes(x=gc_ratio, y=normalized_frags_1_to_200)) +
  geom_point(size=3, shape=15, alpha=0.15, color="black") +
  geom_smooth(method=lm, color="seagreen3") +
  theme_classic()

                                                                                        