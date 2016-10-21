
#####
#####This script is for creating 1,0 adjacency matrices for the clusters produced by subsamp_wgcna.R
####


library("data.table")
library(pryr)

args <- commandArgs(TRUE)

exp_name <- as.character(args[1])
runNum <- as.character(args[2])

##create matrix from gene list
gene_f <- read.delim("/cbcb/lab/smount/ZCL/gene_list.txt", sep = "\t", header = FALSE)

genes <- as.character(gene_f$V1)

#############################################################
###construct matrix as data.table
#############################################################

#make matrix
print("make matrix")
#con_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
con_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
print(dim(con_mat))


#convert to data frame
print("make data frame")
con_df <- as.data.frame(con_mat)
rm(con_mat)

#set colnames
colnames(con_df) <- genes
con_df$g <- genes

print("make data table")
#convert to data table
con_dt <- data.table(con_df)
rm(con_df)

setkey(con_dt, g)

#############################################################
#############################################################
#############################################################


exp_dir <- "/cbcb/project-scratch/ZCL/wgcna/consensus/"


#open cluster list
clust_name <- paste(exp_dir, exp_name, "/clusters/",  runNum, ".csv", sep = "")

clust <- read.csv(clust_name)

for(m in unique(clust$group)){

	print(mem_used())
	if(m == "grey"){
		next
	}

	mod_genes <- as.character(clust[which(clust$group == m), "genes"])

	print(length(mod_genes))

	con_dt[mod_genes[1:length(mod_genes)], mod_genes[1:length(mod_genes)] := 1]

	# combos <- expand.grid(mod1 = mod_genes, mod2 = mod_genes)

}

# print(sum(con_dt[, lapply(.SD, sum, na.rm=TRUE), .SDcols=genes]))

fwrite(con_dt, file = paste(exp_dir, exp_name, "/adjmat/", runNum, ".csv" ,sep = ""))



