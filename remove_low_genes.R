selection <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/hd_lcm_zh.csv")
# fruit <- read.delim("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/fruit_samples.txt", header = FALSE, sep="\n")
# vegetative <- read.delim("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/vegetative_samples.txt", header = FALSE, sep="\n")

#get subset of columns
# columns <- c("X", as.character(fruit$V1), as.character(vegetative$V1))
# selection <- all_counts[, columns]

#reassign rownames
rownames(selection) <- selection[, 1]
selection[, 1] <- NULL

selection <- as.matrix(selection)


low_genes <- c()
low_gene_names <-c()
# remove low expressed genes (if expression is not above 5 95% of the time and there is never expression above 20)

colNum <- ncol(selection)
for(r in 1:nrow(selection)){

	low_count <- 0
	above20 <- FALSE
	for(c in colnames(selection)){

		if(selection[r,c] <= 10){
			low_count <- low_count + 1
		}
		if(selection[r,c] > 20){
			above20<-TRUE
			break
		}

	}

	if(above20 == FALSE && (low_count/colNum >= .95)){
		
		low_genes <- c(low_genes, r)
		low_gene_names <- c(low_gene_names, rownames(selection)[r])
	}
}

print(length(low_genes))
# write.csv(low_genes, file = "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/low_genes.csv")
selection <- selection[-(low_genes), ]

write.csv(selection, file = "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/high_genes.csv")
# write.csv(low_gene_names, file = "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/low_gene_names.csv")


#read in annotation
annotation <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/HD_annotation.csv")

keep_rows <- c()
gene_names <- rownames(selection)

for(i in 1:nrow(annotation)){

	if(annotation[i, "GeneID"] %in% gene_names){
		keep_rows <- c(keep_rows, i)
	}

}

annotation <- annotation[keep_rows,]
write.csv(annotation, file = "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/high_annotation.csv")
