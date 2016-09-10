
#setup
args<-commandArgs(TRUE)

#gather command line arguments
if (length(args) < 2){
  # softPower <- 1
  # minModuleSize <- 50
  # TOM_on <- 0
  stop()
} else {
 softPower <- as.numeric(args[1])
 minModuleSize <- as.numeric(args[2])
 TOM_on <- as.numeric(args[3])
 filename <- as.character(args[5])
 merge_eigengenes <- as.numeric(args[4])
}

#create directory for output if doesnt exist
dir_name <- paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna/", strsplit(args[5], "\\.")[1], "_power_", args[1] , "_modulesize_", args[2] , "_TOM_" , args[3] , "_merge_", args[4], "/", sep = "")

dir.create(dir_name)

setwd(dir_name)

write.table(args, file = "config.txt")

#redirect output
sink(file = "output.txt")


suppressMessages(library("WGCNA"))
suppressMessages(library("methods"))
suppressMessages(library("edgeR"))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))


# stop()

print("importing data")
selection <- read.csv(filename, header=TRUE)

#subset into 


#transpose with only number to keep numeric
transpose <- t(selection[,2:ncol(selection)])
#add colnames
colnames(transpose) <- selection[,"X"]


#remove non variance columsn
print("removing columns with no variance")
transpose <- transpose[,apply(transpose, 2, var, na.rm=TRUE) != 0]


print("making adjacency matrix")
adjacency_res <- adjacency(transpose, power = softPower)

dim(adjacency_res)

if (TOM_on == 1){
  # Turn adjacency into topological overlap
  print("making TOMsimilarity")
  TOM <- TOMsimilarity(adjacency_res, verbose =2)
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

#stop output redirection
sink()

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

print("plotting dendrogram")
# png(paste(strsplit(args[5], "\\.")[1], "_power_", args[1] , "_modulesize_", args[2] , "_TOM_" , args[3] , "_merge_", args[4], "_dendro.png", sep = ""),width=1500,height=1500)

png("x.png", ,width=1500,height=1500)
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,addGuide = TRUE, hang=0.3, guideHang = 0.05, main = "Gene dendrogram and module colors")

dev.off()





