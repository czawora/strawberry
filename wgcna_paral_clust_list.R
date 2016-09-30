##this script reads a networks cluster info and outputs an
##adjacency matrix of clusters

library(pryr)
require(utils)

args <- commandArgs(TRUE)

cluster_dir <- "/cbcb/project-scratch/ZCL/wgcna/"

name <- as.character(args[1])
merge <- as.numeric(args[2])

file.create(paste("/cbcb/project-scratch/ZCL/wgcna_consensus/complete/", name, sep =""))

print("read cluster")
if(merge == 1){
	clust <- read.csv(paste(cluster_dir, name, "/merged_clusters.csv", sep = "" ), header = TRUE)
} else {
	clust <- read.csv(paste(cluster_dir, name, "/original_clusters.csv", sep =""), header = TRUE)
}

clust_count <- 1

all_pairs <- data.frame(pairs1 = character(), pairs2 = character())

print("cluster read in")
print(mem_used())


mat <- matrix(0, nrow = length(unique(clust$genes)), ncol = length(unique(clust$genes)))

rownames(mat) <- as.character(unique(clust$genes))
colnames(mat) <- as.character(unique(clust$genes))

print("matrix made")


for(m in unique(clust$group)){

	print(paste(clust_count, " / ", length(unique(clust$group)), sep = ""))

	if(m == "grey"){
		next
	}

	mod_genes <- as.character(clust[which(clust$group == m), "genes"])

	print(length(mod_genes))

	combos <- expand.grid(mod1 = mod_genes, mod2 = mod_genes)

	# combos$both <- paste(combos$mod1,"AAAAA", combos$mod2, sep = "")

	all_pairs <- rbind(all_pairs, data.frame(pairs1 = combos$mod1, pairs2 = combos$mod2))

	# print(length(combos$both))

	# write(combos$both, file = paste("/cbcb/project-scratch/ZCL/wgcna_consensus/", name, "_hd_clust_list.txt", sep = ""), append = TRUE)

	print(mem_used())
	clust_count <- clust_count + 1
}

mat[ all_pairs$pairs1, all_pairs$pairs2 ] = 1

write(mat, file = paste("/cbcb/project-scratch/ZCL/wgcna_consensus/lcpm_mats/", name, "_hd_clust_mat.txt", sep = ""))

file.remove(paste("/cbcb/project-scratch/ZCL/wgcna_consensus/complete/", name, sep =""))
