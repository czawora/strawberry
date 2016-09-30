
args <- commandArgs(TRUE)

norm <- as.character(args[1])

cluster_dir <- "/cbcb/project-scratch/ZCL/wgcna/"

all_nets <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/config_clusters.csv", header = TRUE)

norm_sub <- all_nets[which(all_nets$norm == norm), c("config", "merge")]

print(norm_sub)

print(paste(nrow(norm_sub), " networks in the norm subset", sep = ""))

write(as.character(norm_sub$config), file = paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/", norm, "_networks.txt", sep =""))

