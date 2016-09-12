library(ConsensusClusterPlus)
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

topvar <- names(rev(sort(variance))[1:10000])

sub_dt <- dt[, topvar]

print("making distnace object")
dist_dt <- as.dist(1 - cor(sub_dt, method = "pearson"))

rm(d)
rm(dt)

print("begin clustering")
print(mem_used())
res <- ConsensusClusterPlus(dist_dt, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "hc", title = "/cbcb/project-scratch/ZCL/consensus_plus", plot = "png", writeTable = TRUE, verbose = TRUE)

print(mem_used())
print("finish clustering")

calcICL(res,title="/cbcb/project-scratch/ZCL/consensus_plus",plot="png",writeTable=TRUE)

sink()


