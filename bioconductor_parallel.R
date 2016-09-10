library("Rsamtools")
library("BiocParallel")
library("BatchJobs")

reg <- makeRegistry(id="bioconductor_p")
print(reg)

library("WGCNA")
library("methods")
library("edgeR")
library(gplots)
library(RColorBrewer)

powers <- c(1,2,4,8,12)
module_size <- c(30,45,60,75)
TOM <- c(TRUE, FALSE)

param <- as.data.frame(c())

for(p in powers){

  for(m in module_size){

    for(t in TOM){

      param <- rbind(param, c(p,m,t))
    }
  }
}

colnames(param) <- c("power", "module", "tom")

 
#making outputs directories for different paramters
#dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

print("importing data")
selection <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/rpkm_hd_lcm_zh.csv", header=TRUE)
# fruit <- read.delim("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/fruit_samples.txt", header = FALSE, sep="\n")
# vegetative <- read.delim("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/vegetative_samples.txt", header = FALSE, sep="\n")

#get subset of columns
# columns <- c("X", as.character(fruit$V1), as.character(vegetative$V1))
# selection <- all_counts[, colnames(all_counts)]

#transpose with only number to keep numeric
transpose <- t(selection[,2:ncol(selection)])
#add colnames
colnames(transpose) <- selection[,"X"]

# #reassign rownames
# rownames(selection) <- selection[, "X"]
# selection[, "X"] <- NULL

# #read in annotation
# annotation <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/HD_annotation.csv")

# #normalize with rpkm no log
# cpm_vals <- rpkm(selection, annotation[, "Length"])

# dim(cpm_vals)

#remove non variance columsn
print("removing columns with no variance")
transpose <- transpose[,apply(transpose, 2, var, na.rm=TRUE) != 0]

run_tests <- function(param_row, data_df){

	library("WGCNA")
	library("methods")
	library("edgeR")
	library(gplots)
	library(RColorBrewer)


  print(param_row)

  class(param_row)
  softPower <- param_row[1]
  minModuleSize <- param_row[2]
  TOM_on <- param_row[3]

  print("making adjacency matrix")
  adjacency_res <- adjacency(data_df, power = softPower)

  dim(adjacency_res)

  if (TOM_on == 1){
    # Turn adjacency into topological overlap
    print("making TOMsimilarity")
    TOM <- TOMsimilarity(adjacency_res)
    distanceMat <- 1-TOM
  
  } else {

    distanceMat <- 1 - adjacency_res
  }

  print("building gene Tree")
  # Call the hierarchical clustering function
  geneTree <- hclust(as.dist(distanceMat), method = "average")

  # Module identification using dynamic tree cut:
  print("cutting modules")
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = distanceMat, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, verbose= 4)

  table(dynamicMods)

  # Convert numeric lables into colors
  dynamicColors <- labels2colors(dynamicMods)
  table(dynamicColors)

  print("plotting dendrogram")
  png(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/graphs/", softPower, "_", minModuleSize, "_", TOM_on, "_", "rpkm", "_dendro.png", sep = ""),width=1500,height=1500)
  
  # Plot the dendrogram and colors underneath
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,addGuide = TRUE, hang=0.3, guideHang = 0.05, main = "Gene dendrogram and module colors")

  dev.off()

}

#apply the run tests function to every row of param, also passing in the data for analysis
# apply(param, 1, run_tests, transpose)

# png(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/graphs/", "rpkm_veg_fruit", "_cluster.png", sep = ""),width=500,height=500)
# # Plot the resulting clustering tree (dendrogram)
# plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
# labels = FALSE, hang = 0.04);

# dev.off()

# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=30, by=2))
# # Call the network topology analysis function
# sft = pickSoftThreshold(rpkm_vals, powerVector = powers, verbose = 5)


# # Plot the results:
# cex1 = 0.9
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2))

# png(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/graphs/", "cpm_veg_fruit", "_modelfit.png", sep = ""),width=500,height=500)

# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
# xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
# main = paste("Scale independence"))
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
# labels=powers,cex=cex1,col="red")
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")

# dev.off()


res <- batchMap(reg, run_tests, param, transpose)

print(reg)

done <- submitJobs(reg)






