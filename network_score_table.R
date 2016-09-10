
setwd("/cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/scores/")

print("-----------read net names")
net_names <- read.delim("/cbcb/project-scratch/ZCL/go_enrich_results/exon_hd.txt", header = FALSE)
net_names <- as.character(net_names$V1)

### get the list of networks to create the table for

df <- data.frame(network = character(), norm = character(), power = double(), minModuleSize = double(), merge = double(), num_clusters = double(), uniqueGO = double(), pval_score = double())

for(n in net_names){

	print(paste("-----------", n, sep = ""))
	###are there header?

	net <- read.delim(paste(n, "_uniqueGO_sumP.txt", sep =""), header = TRUE, sep = " ")
	u <- as.numeric(net$uniqueGO_score)
	p <- as.numeric(net$sum_neg_log_p)

	##parse the network name for the info

	print("--------splitting")
	sub_split <- strsplit(n, "sub_")[[1]][2]
	x_split <- strsplit(sub_split, "_")[[1]]

	n_norm <- x_split[1]
	n_mms <- as.numeric(x_split[5])
	n_merge <- as.numeric(x_split[9])
	n_power <- as.numeric(x_split[3])

	#get the number of clusters for that network
	
	if(n_merge == 1){
		clust <- read.delim(paste("/cbcb/project-scratch/ZCL/wgcna/", n , "/merged_clusters.csv",sep = ""), sep= ",", header = TRUE)
	} else {
		clust <- read.delim(paste("/cbcb/project-scratch/ZCL/wgcna/", n , "/original_clusters.csv",sep = ""), sep= ",", header = TRUE)
	}

	moduleColors <- clust$group
	moduleColors <- moduleColors[moduleColors != "grey"]
	n_clust_num <- length(unique(moduleColors))

	df <- rbind(df, data.frame(network = n, norm = n_norm, power = n_power, minModuleSize = n_mms, merge = n_merge, num_clusters = n_clust_num, uniqueGO = u, pval_score = p))

}


write.csv(df, file = "hd_summary_scores.csv")

print("-------------complete")