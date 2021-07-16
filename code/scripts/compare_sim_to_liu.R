library(corrplot)

setwd("~/github/read_simulation/analyses/liu_2017_RMS/3_compare_rms_to_sim_mapping/")

df = read.csv(file='compare_all.depth.txt', sep=)


mydata = df[c(3:7)]
mydata.cor = cor(mydata, use = "pairwise.complete.obs")
corrplot(diag = F, method = "color", mydata.cor)

