

args <- commandArgs(TRUE)

library(data.table)

###some inputs including adjmat file

#adjmat
exp_name <- as.character(args[1])
adjmat_name <- as.character(args[2])
cutoff <- as.numeric(args[3])

print(args)

adjmat <- fread(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name,"/", adjmat_name, sep = ""), header= TRUE, sep = ",")

##can use arbitrary cutoff
#cutoff
col <- colors()

lidx <- 0
clist <- list()

adjmat <- as.data.frame(adjmat)
##index by rownames
rownames(adjmat) <- colnames(adjmat)
cnames <- colnames(adjmat)

gene_convered <- 0

for(r in rownames(adjmat)){

	if(gene_convered < nrow(adjmat)){

		print(r)
		print(gene_convered)
		
		##if rownames in existing list

		success <- FALSE
		if(lidx > 0){
			
			for(i in 1:lidx){

				print(i)
				
				s <- sum(r %in% clist[[i]])

				if(s > 0){

					gene_convered <- gene_convered - length(clist[[i]])
					row <- adjmat[r,]
					clist[[i]] <- c(clist[[i]], cnames[row >= cutoff])
					clist[[i]] <- unique(clist[[i]])
					gene_convered <- gene_convered + length(clist[[lidx]])
					success <- TRUE

				} 
			}

			if(success == FALSE){
				
				lidx <- lidx + 1

				row <- adjmat[r,]
				clist[[lidx]] <- cnames[row >= cutoff]
				gene_convered <- gene_convered + length(clist[[lidx]])

			}

		} else {
			
			lidx <- lidx + 1

			row <- adjmat[r,]
			clist[[lidx]] <- cnames[row >= cutoff]
			gene_convered <- gene_convered + length(clist[[lidx]])


		}
	}
}


##write out statistics

total_vec <- c()

for(l in 1:lidx){

	total_vec <- c(total_vec, clist[[l]])
}

write(paste(exp_name, cutoff, "\n", sep = "\t"), file = paste("/cbcb/project-scratch/wgcna/consensus/", exp_name, "/cutoffs/", cutoff, ".info", sep = ""), append = TRUE)

write(paste("cluster number", lidx, "\n", sep = "\t"), file = paste("/cbcb/project-scratch/wgcna/consensus/", exp_name, "/cutoffs/", cutoff, ".info", sep = ""), append = TRUE)

write(paste("genes above cutoff", length(total_vec), "\n", sep = "\t"), file = paste("/cbcb/project-scratch/wgcna/consensus/", exp_name, "/cutoffs/", cutoff, ".info", sep = ""), append = TRUE)


#######################

##write out clusters in 2 column format

write(paste("gene", "cluster\n", sep = ","),file = paste("/cbcb/project-scratch/wgcna/consensus/", exp_name, "/cutoffs/", cutoff, ".csv", sep = "") ,append = TRUE)

for(l in 1:lidx){


	clust_color <- col[l]

	out_df <- data.frame(gene = clist[[l]], color = clust_color)

	write.csv(out_df, file = paste("/cbcb/project-scratch/wgcna/consensus/", exp_name, "/cutoffs/", cutoff, ".csv", sep = ""), row.names = FALSE, col.names = FALSE, append = TRUE)
}



