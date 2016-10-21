##probably need 120GB run


library("data.table")

args <- commandArgs(TRUE)

exp_name <- as.character(args[1])
lst_name <- as.character(args[2])


lst <- as.character(read.table(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name, "/", lst_name, sep = ""), sep = "\n", header = FALSE)$V1)


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

###############################################
###############################################
###############################################

for(i in 1:length(lst)){

	sublist_name <- lst[i]

	print(i)

	sublist <- as.character(read.delim(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name, "/sub_genes/", sublist_name, sep = ""), sep = "\n", header = FALSE)$V1)


	count <- 0
	for(s in sublist){
		
		print(count)
		b <- con_dt[sublist, s, with = FALSE]
		con_dt[sublist, s := b + 1, with = FALSE]

		count <- count + 1
	}



}

fwrite(con_dt, paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name, "/indicator_mat.csv", sep = ""))









