library(ConsensusClusterPlusChris)
library(pryr) #memory monitoring

sink(file = "/cbcb/project-scratch/ZCL/consensus_plus/output.txt")

print("libraries loaded")

filename <- "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/cpm_hd_lcm_zh.csv"

d <- read.csv(filename, header = TRUE)

rownames(d) <- d$X

d$X <- NULL
d$X.1 <- NULL

dt <- t(d)

dt <- dt[, apply(dt, 2, var, na.rm=TRUE) > 0.05]

variance <- apply(dt, 2, var, na.rm=TRUE)

topvar <- names(rev(sort(variance))[1:100])

sub_dt <- dt[, topvar]

print("making distnace object")
dist_dt <- as.dist(1 - cor(dt, method = "pearson"))

rm(d)
#rm(dt)

setwd("/cbcb/project-scratch/ZCL/consensus_plus/")

print("begin clustering")
print(mem_used())
res <- ConsensusClusterPlus(dt, maxK = 15, minNumK = 17, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "hc", title = "cluster", plot = "pngBMP", writeTable = TRUE, verbose = TRUE)

exit(1)

print(mem_used())
print("finish clustering")

calcICL(res,title="calc", maxK = 15, minNumK = 17,plot="pngBMP",writeTable=TRUE)

sink()


