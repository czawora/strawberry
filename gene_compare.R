
compare_clust <- function(c1, c2){

	print("in compare_clust")

	c1 <- c1[which(c1$group != "grey"), ]
	c2 <- c2[which(c2$group != "grey"), ]

	rownames(c1) <- c1$genes
	rownames(c2) <- c2$genes

	c1_genes <- as.character(c1$genes)
	c2_genes <- as.character(c2$genes)

	gene_union <- intersect(c1_genes, c2_genes)
	print(gene_union)
	l_gene_union <- length(gene_union)
	print(l_gene_union)

	sim_scores <- c()

	for(g in gene_union){

		print("current gene")
		print(g)

		c1_group_name <- c1[g, "group"]
		c2_group_name <- c2[g, "group"]

		if(is.null(c1_group_name) == TRUE){
			score_res = 0
		} else if(is.null(c2_group_name) == TRUE){
			score_res = 0
		} else{

			#two cluster sets
			genes_in_group1 <- as.character(c1[which(c1$group == c1_group_name), "genes"])
			genes_in_group2 <- as.character(c2[which(c2$group == c2_group_name), "genes"])

			print("group 1 genes")
			print(genes_in_group1)
			print("group 2 genes")
			print(genes_in_group2)

			l2 <- as.double(length(genes_in_group2))
			print(paste("l2_genes ", l2, sep = ""))
			l1 <- as.double(length(genes_in_group1))
			print(paste("l1_genes ", l1, sep = ""))

			intersect_buddies <- intersect(genes_in_group1, genes_in_group2)
			print(paste("intersect_buddies: ", intersect_buddies, sep = ""))

			intersect_l <- as.double(length(intersect_buddies))

			if(l2 < l1){

				score_res <- intersect_l/l2
			} else {
				score_res <- intersect_l/l1
			}

			if(intersect_l == 1){

				score_res = 0
			}

			if(is.na(score_res)==TRUE){
				score_res = 0
			}

		}
		print(g)
		print(score_res)
		sim_scores <- c(sim_scores, score_res)

		print(paste(length(sim_scores), "/", length(gene_union), sep =""))


	}

	print("c1 mag")
	print(length(c1_genes))
	print("c2 mag")
	print(length(c2_genes))
	print(sum(sim_scores)/l_gene_union)
	print(sim_scores)
	print(sum(sim_scores))
	return(sum(sim_scores)/l_gene_union)
}

# already_compared <- function(n1, n2, check_df){

# 	return_val <- -1

# 	prev_comp1 <- check_df[which(check_df$x == n1 & check_df$y == n2), ]
# 	prev_comp2 <- check_df[which(check_df$x == n2 & check_df$y == n1), ]

# 	if(nrow(prev_comp1) > 0){

# 		return_val <- prev_comp1$score

# 	} else if(nrow(prev_comp2) > 0) {

# 		return_val <- prev_comp2$score
# 	}

# 	return(return_val)
# }

args<-commandArgs(TRUE)

n1 <- as.character(args[1])
n2 <- as.character(args[2])

print(n1)
print(n2)

n1 <- gsub("\n", "", n1)
n2 <- gsub("\n", "", n2)

print(n1)
print(n2)


comparisons <- data.frame(x = character(), y = character(), score = double())

# clust1 <- read.delim(paste("/Users/Chris/Documents/strawberry/ZCL/", n1, sep = ""), sep = ",", header = TRUE)
# print(clust1)
# clust2 <- read.delim(paste("/Users/Chris/Documents/strawberry/ZCL/", n2, sep = ""), sep = ",", header = TRUE)
# print(clust2)

network_dir1 <- paste("/cbcb/project-scratch/ZCL/wgcna/", n1, sep ="")

#get network info
sub_split1 <- strsplit(n1, "sub_")[[1]][2]
set_split1 <- strsplit(n1, "sub_")[[1]][1]
x_split1 <- strsplit(sub_split1, "_")[[1]]

n_set1 <- strsplit(set_split1, "_")[[1]][1]
n_norm1 <- x_split1[1]
n_mms1 <- as.numeric(x_split1[5])
n_merge1 <- as.numeric(x_split1[9])
n_power1 <- as.numeric(x_split1[3])


print(paste(network_dir1, "/merged_clusters.csv", sep =""))
print(paste(network_dir1, "/original_clusters.csv", sep =""))
#read network clusters
print("reading clusters")
if(n_merge1 == 1){
	clust1 <- read.delim(paste(network_dir1, "/merged_clusters.csv", sep =""), sep =",", header = TRUE)
} else {
	clust1 <- read.delim(paste(network_dir1, "/original_clusters.csv", sep =""), sep =",", header = TRUE)
}

network_dir2 <- paste("/cbcb/project-scratch/ZCL/wgcna/", n2, sep ="")

#get network info
sub_split2 <- strsplit(n2, "sub_")[[1]][2]
set_split2 <- strsplit(n2, "sub_")[[1]][1]
x_split2 <- strsplit(sub_split2, "_")[[1]]

n_set2 <- strsplit(set_split2, "_")[[1]][1]
n_norm2 <- x_split2[1]
n_mms2 <- as.numeric(x_split2[5])
n_merge2 <- as.numeric(x_split2[9])
n_power2 <- as.numeric(x_split2[3])

#read network clusters
print("reading clusters")
if(n_merge2 == 1){
	clust2 <- read.delim(paste(network_dir2, "/merged_clusters.csv", sep =""), sep =",", header = TRUE)
} else {
	clust2 <- read.delim(paste(network_dir2, "/original_clusters.csv", sep =""), sep =",", header = TRUE)
}

result <- compare_clust(clust1, clust2)



comparisons <- rbind(comparisons, data.frame(x = n1, y = n2, score = result))
print(comparisons)
write.csv(comparisons, file = paste("/cbcb/project-scratch/ZCL/comparison_mats/gene_comp/results/", n1, "____", n2, ".csv", sep =""))


