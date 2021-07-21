library(corrplot)

setwd("~/github/read_simulation/code/scripts")

df = read.csv(file='../../analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/compare_all.depth.txt', sep="\t")


mydata = df[c(3:7)]
mydata.cor = cor(mydata, use = "pairwise.complete.obs")
corrplot(diag = F, method = "color", mydata.cor)

