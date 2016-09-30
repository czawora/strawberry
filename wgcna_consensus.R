suppressMessages(library("WGCNA"))
library(pryr)

args <- commandArgs(TRUE)

norm <- as.character(args[1])

cluster_dir <- "/cbcb/project-scratch/ZCL/wgcna/"

all_nets <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/config_clusters.csv", header = TRUE)

#outf <- paste("/cbcb/project-scratch/ZCL/wgcna_consensus/", norm, "_hd_clust_list.csv", sep = "")
#write.csv("network,one,two", file = outf, append = TRUE)

norm_sub <- all_nets[which(all_nets$norm == norm), c("config", "merge")]

print(paste(nrow(norm_sub), " networks in the norm subset", sep = ""))

consensus_mat <- "x"
#clust_df <- data.frame(network = character(), one = character, two = character())

for(n in nrow(norm_sub)){

	name <- as.character(norm_sub[n,"config"])
	merge <- as.numeric(norm_sub[n, "merge"])
	print(name)
	print(mem_used())

	print("read cluster")
	if(merge == 1){
		clust <- read.csv(paste(cluster_dir, name, "/merged_clusters.csv", sep = "" ), header = TRUE)
	} else {
		clust <- read.csv(paste(cluster_dir, name, "/original_clusters.csv", sep =""), header = TRUE)
	}

	if(consensus_mat == "x"){

		print("make matrix")

		consensus_mat <- matrix(0 , nrow = length(unique(as.character(clust$genes))), ncol = length(unique(as.character(clust$genes))))

		rownames(consensus_mat) <- unique(as.character(clust$genes))
		colnames(consensus_mat) <- rownames(consensus_mat)
	}


	print(head(consensus_mat, 1))

	print(paste("running through ", length(unique(clust$group)), " clusters", sep = ""))

	count <- 1
	for(m in unique(clust$group)){

		if(m == "grey"){
			next
		}

		mod_genes <- as.character(clust[which(clust$group == m), "genes"])

		print(paste(count, " / ",length(unique(clust$group)) , "   ", m ," genes: ", length(mod_genes), sep = "" ))

		consensus_mat[mod_genes, mod_genes] = consensus_mat[mod_genes, mod_genes] + 1

		print(mem_used())
	
		count <- count + 1	
	}

}

consensus_mat[rownames(consensus_mat), colnames(consensus_mat)] <- 1

consensus_mat <- consensus_mat / nrow(norm_sub)

write.csv(consensus_mat, file = paste("/cbcb/project-scratch/ZCL/wgcna_consensus/", norm, "_hd_clust_list.csv", sep = ""))
