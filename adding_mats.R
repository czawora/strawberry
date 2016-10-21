

library(data.table)
library(pryr)

args <- commandArgs(TRUE)

exp_name <- as.character(args[1])

mat1 <- as.character(args[2])
mat2 <- as.character(args[3])

out_name <- as.character(args[4])

sort <- as.character(args[5])

out_file <- paste("/cbcb/project-scratch/ZCL/wgcna/consensus/",exp_name, "/adjmat/", out_name, ".csv" , sep = "")


print(args)
print(out_file)

x <- fread(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name, "/adjmat/", mat1, sep = ""), sep = ",", header = TRUE)

if("g" %in% colnames(x) == TRUE){
	
	setkey(x, g)

	setkey(x, NULL)

	x[, g := NULL]

}

mx <- as.matrix(x)
## do columns names stick around after matrix transform
rm(x)


y <- fread(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name, "/adjmat/", mat2, sep = ""), sep = ",", header = TRUE)

if("g" %in% colnames(y) == TRUE){

	setkey(y, g)

	setkey(y, NULL)

	y[, g := NULL]

}

my <- as.matrix(y)
## do columns names stick around after matrix transform
rm(y)


mz <- mx + my


fwrite(as.data.frame(mz), out_file)






