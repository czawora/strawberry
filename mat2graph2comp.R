
args <- commandArgs(TRUE)

library(data.table)
library(igraph)


###some inputs including adjmat file

#adjmat
exp_name <- as.character(args[1])
adjmat_name <- as.character(args[2])
cutoff <- as.numeric(args[3])


sink(file = paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name, "/cutoffs/", cutoff, ".mgclog", sep = ""))

print(args)

adjmat <- fread(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name,"/", adjmat_name, sep = ""), header= TRUE, sep = ",")

m_adjmat <- as.matrix(adjmat)

print("as.matrix")

rm(adjmat)

m_adjmat[m_adjmat >= cutoff] <- 1
m_adjmat[m_adjmat < cutoff] <- 0

print("1,0 matrix")

g <- graph_from_adjacency_matrix(m_adjmat)

write_graph(g,file = paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name, "/cutoffs/", cutoff, ".graphml", sep = ""), format = "graphml")

