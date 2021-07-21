library(corrplot)

setwd("~/github/read_simulation/code/scripts")

df = read.csv(file='../../analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/compare_all.depth.txt', sep="\t")


mydata = df[c(3:7)]
mydata.cor = cor(mydata, use = "pairwise.complete.obs")
corrplot(diag = F, method = "color", mydata.cor)

df = read.table(file='../../analyses/liu_2017_RMS/2021_07_20_updated_copy_number/2_map_to_ref/compare_depths.txt', header=T)
mydata = df[c(3,7,8,4,5,6)]

mydata.cor = cor(mydata, use = "pairwise.complete.obs")
corrplot(diag = F, method = "color", mydata.cor)
