library(tidyverse)

setwd("~/github/read_simulation/analyses/snipen_2021_RMS/2021_07_26_test_mock_data/7_compare_bracken_profiles/")

df_50 = read.csv(file='50_percent_bracken.txt', header=F, sep="\t")
df_50$group = rep("50_read_abundance", length(df_50$V1))

df_75 = read.csv(file='75_percent_bracken.txt', header=F, sep="\t")
df_75$group = rep("75_read_abundance", length(df_75$V1))

df_95 = read.csv(file='95_percent_bracken.txt', header=F, sep="\t")
df_95$group = rep("95_read_abundance", length(df_95$V1))

df_s1 = read.csv(file='SRR10199716_bracken.txt', header=F, sep="\t")
df_s1$group = rep("s1_read_abundance", length(df_s1$V1))

df_s2 = read.csv(file='SRR10199724_bracken.txt', header=F, sep="\t")
df_s2$group = rep("s2_read_abundance", length(df_s2$V1))

df_s3 = read.csv(file='SRR10199725_bracken.txt', header=F, sep="\t")
df_s3$group = rep("s3_read_abundance", length(df_s3$V1))

df = rbind(df_50, df_75, df_95, df_s1, df_s2, df_s3)

# pull species levels only (all 20 here are species) and remove spaces in names
df = df[ which(df$V4=="S"), ]

df = as.data.frame(apply(df,2,function(x)gsub('\\s+', '',x)))
percent_df = df[c(1,6,7)]

# get the true abundance profiles
df_true = read.csv(file='named_sampled_files_key.txt', header=F, sep="")
names = rep("true", length(df_true$V1))
df_true$group = names
df_true$percents = (df_true$V3/sum(df_true$V3))*100
name_key = df_true[c(1,4)]

df_true = df_true[c(6,4,5)]
df_true = rename(df_true, V1 = percents, V6 = V4)
percent_df = rbind(percent_df, df_true)

# get the abundances as depths
df_50_depth = read.csv(file='50_percent_avg_depths.txt', header=F, sep="")
df_50_depth$group = rep("50_depth", length(df_50_depth$V1))
df_50_depth$percents = (df_50_depth$V2/sum(df_50_depth$V2))*100
df_50_depth = merge(df_50_depth[, c(1, 3, 4)], name_key[, c(1, 2)])

df_75_depth = read.csv(file='75_percent_avg_depths.txt', header=F, sep="")
df_75_depth$group = rep("75_depth", length(df_75_depth$V1))
df_75_depth$percents = (df_75_depth$V2/sum(df_75_depth$V2))*100
df_75_depth = merge(df_75_depth[, c(1, 3, 4)], name_key[, c(1, 2)])

df_95_depth = read.csv(file='95_percent_avg_depths.txt', header=F, sep="")
df_95_depth$group = rep("95_depth", length(df_95_depth$V1))
df_95_depth$percents = (df_95_depth$V2/sum(df_95_depth$V2))*100
df_95_depth = merge(df_95_depth[, c(1, 3, 4)], name_key[, c(1, 2)])

df_s1_depth = read.csv(file='SRR10199716_avg_depths.txt', header=F, sep="")
df_s1_depth$group = rep("s1_depth", length(df_s1_depth$V1))
df_s1_depth$percents = (df_s1_depth$V2/sum(df_s1_depth$V2))*100
df_s1_depth = merge(df_s1_depth[, c(1, 3, 4)], name_key[, c(1, 2)])

df_s2_depth = read.csv(file='SRR10199724_avg_depths.txt', header=F, sep="")
df_s2_depth$group = rep("s2_depth", length(df_s2_depth$V1))
df_s2_depth$percents = (df_s2_depth$V2/sum(df_s2_depth$V2))*100
df_s2_depth = merge(df_s2_depth[, c(1, 3, 4)], name_key[, c(1, 2)])

df_s3_depth = read.csv(file='SRR10199725_avg_depths.txt', header=F, sep="")
df_s3_depth$group = rep("s3_depth", length(df_s3_depth$V1))
df_s3_depth$percents = (df_s3_depth$V2/sum(df_s3_depth$V2))*100
df_s3_depth = merge(df_s3_depth[, c(1, 3, 4)], name_key[, c(1, 2)])

depth_df = rbind(df_50_depth, df_75_depth, df_95_depth, df_s1_depth, df_s2_depth, df_s3_depth)

depth_df = depth_df[c(3,4,2)]
depth_df = rename(depth_df, V1 = percents, V6 = V4, )
percent_df = rbind(percent_df, depth_df)
percent_df = percent_df[order(percent_df$V6),]

rm(list=setdiff(ls(), c("depth_df", "percent_df")))

# "unmelt" the data so it can be used in a correlation
df = pivot_wider(percent_df, names_from = group, values_from = V1)
#write.csv(df, "all_data.csv")

library(corrplot)
df = df[c(2:14)]
df <- data.frame(sapply(df, function(x) as.numeric(as.character(x))))
df.cor = cor(df, use = "pairwise.complete.obs")
corrplot(diag = F, method = "color", df.cor)

####################################################
####################################################
####################################################

# the following reproduces all the above depth samples without calculating percentages
rm(list=ls())

# get the true abundance profiles
df_true = read.csv(file='named_sampled_files_key.txt', header=F, sep="")
names = rep("true", length(df_true$V1))
df_true$group = names
df_true$percents = (df_true$V3/sum(df_true$V3))*100
name_key = df_true[c(1,4)]

df_true = df_true[c(3,4,5)]
df_true = rename(df_true, V1 = V3, V6 = V4)

# get the abundances as depths
df_50_depth = read.csv(file='50_percent_avg_depths.txt', header=F, sep="")
df_50_depth$group = rep("50_depth", length(df_50_depth$V1))
df_50_depth = merge(df_50_depth[, c(1, 2, 3)], name_key[, c(1, 2)])

df_75_depth = read.csv(file='75_percent_avg_depths.txt', header=F, sep="")
df_75_depth$group = rep("75_depth", length(df_75_depth$V1))
df_75_depth = merge(df_75_depth[, c(1, 2, 3)], name_key[, c(1, 2)])

df_95_depth = read.csv(file='95_percent_avg_depths.txt', header=F, sep="")
df_95_depth$group = rep("95_depth", length(df_95_depth$V1))
df_95_depth = merge(df_95_depth[, c(1, 2, 3)], name_key[, c(1, 2)])

df_s1_depth = read.csv(file='SRR10199716_avg_depths.txt', header=F, sep="")
df_s1_depth$group = rep("s1_depth", length(df_s1_depth$V1))
df_s1_depth = merge(df_s1_depth[, c(1, 2, 3)], name_key[, c(1, 2)])

df_s2_depth = read.csv(file='SRR10199724_avg_depths.txt', header=F, sep="")
df_s2_depth$group = rep("s2_depth", length(df_s2_depth$V1))
df_s2_depth = merge(df_s2_depth[, c(1, 2, 3)], name_key[, c(1, 2)])

df_s3_depth = read.csv(file='SRR10199725_avg_depths.txt', header=F, sep="")
df_s3_depth$group = rep("s3_depth", length(df_s3_depth$V1))
df_s3_depth = merge(df_s3_depth[, c(1, 2, 3)], name_key[, c(1, 2)])

depth_df = rbind(df_50_depth, df_75_depth, df_95_depth, df_s1_depth, df_s2_depth, df_s3_depth)

depth_df = depth_df[c(2,4,3)]
depth_df = rename(depth_df, V1 = V2, V6 = V4, )
raw_depth_df = rbind(df_true, depth_df)
raw_depth_df = raw_depth_df[order(raw_depth_df$V6),]

rm(list=setdiff(ls(), "raw_depth_df"))

# "unmelt" the data so it can be used in a correlation
df = pivot_wider(raw_depth_df, names_from = group, values_from = V1)
#write.csv(df, "raw_depth_data.csv")


library(corrplot)
df = df[c(2:8)]
df <- data.frame(sapply(df, function(x) as.numeric(as.character(x))))
df.cor = cor(df, use = "pairwise.complete.obs")
corrplot(diag = F, method = "color", df.cor)
