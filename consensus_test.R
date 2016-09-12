library(ConsensusClusteringPlus)
library(pryr) #memory monitoring

sink(file = "/cbcb/project-scratch/ZCL/consensus_plus/output.txt")

filename <- "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/cpm_hd_lcm_zh.csv"

d <- read.csv(filename, header = TRUE)

rownames(d) <- d$X

d$X <- NULL
d$X.1 <- NULL

dt <- t(d)

dt <- dt[, apply(dt, 2, var, na.rm=TRUE) > 0.05]

dist_dt <- as.dist(1 - cor(dt, method = "pearson"))

rm(d)

print(mem_used())
res <- ConsensusClusteringPlus(dist_dt, maxK = 3, pItem = 0.8, pfeature = 1, clusterAlg = "hc", title = "/cbcb/project-scratch/ZCL/consensus_plus", plot = "png", verbose = TRUE)

print(mem_used())
