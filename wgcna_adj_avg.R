
suppressMessages(library("WGCNA"))
suppressMessages(library("methods"))
suppressMessages(library("edgeR"))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
nthr <- enableWGCNAThreads()

args<-commandArgs(TRUE)

f <- as.character(args[1])
traits <- as.character(args[2])

powers <- c(1,2,4,8,12,16)
minModuleSize <- c(40, 60, 90, 120, 150, 180, 210)

selection <- read.csv(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/", f, sep = ""), header = TRUE)

#set rownames
rownames(selection) <- selection[,"X"]
#remove useless column
selection[,"X"] <- NULL
selection[,"X.1"] <- NULL


print(dim(selection))
#transpose with only number to keep numeric
transpose <- t(selection)
#selection <- NULL

#"gene35204-v1.0-hybrid"

#remove non variance columsn
print("removing columns with no variance")
transpose <- transpose[,apply(transpose, 2, var, na.rm=TRUE) > 0.05]
print(dim(transpose))


master <- matrix(data = 0, nrow = ncol(transpose), ncol = ncol(transpose))


for(p in powers){

	print(paste("power: ", p, sep = ""))
	master = master + adjacency(transpose, power = p, type = "signed")
}

master <- master/6

write.csv(master, file = "/cbcb/project-scratch/ZCL/wgcna_consensus/cpm_avg.csv")

distanceMat <- 1 - master

geneTree <- hclust(as.dist(distanceMat), method = "average")


print("cutting modules")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = distanceMat, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 100, verbose= 4)

table(dynamicMods)


dynamicColors <- labels2colors(dynamicMods)
l <- length(unique(dynamicColors))
l
table(dynamicColors)

if(l > 200){
	print("too many modules")
	quit(save = "no")
}

#write out clusters
gene_names <- rownames(adjacency_res)
clusters <- as.data.frame(cbind(genes= gene_names, group = dynamicColors))
clusters <- clusters[order(clusters$group), ]

write.csv(clusters, file = "/cbcb/project-scratch/ZCL/wgcna_consensus/original_clusters.csv")

print("plotting dendrogram")
png("/cbcb/project-scratch/ZCL/wgcna_consensus/orignal_dendro.png",width=1500,height=1500)

# png("x.png", ,width=1500,height=1500)
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,addGuide = TRUE, hang=0.3, guideHang = 0.05, main = "Gene dendrogram and module colors")

dev.off()

moduleColors <- dynamicColors

MEList <- moduleEigengenes(transpose, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result

png("/cbcb/project-scratch/ZCL/wgcna_consensus/clustered_me.png", width =1500, height=1500)
plot(METree, main = "Clustering of module eigengenes",
   xlab = "", sub = "")

MEDissThres <- 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

dev.off()


merge <- mergeCloseModules(transpose, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
moduleColors <- mergedColors
print(table(moduleColors))
print(length(unique(moduleColors)))
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

png("/cbcb/project-scratch/ZCL/wgcna_consensus/merged_dendro.png",width=1500,height=1500)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                  c("Dynamic Tree Cut", "Merged dynamic"),
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05)
dev.off()


#write out merged clusters
gene_names <- rownames(adjacency_res)
mclusters <- as.data.frame(cbind(genes= gene_names, group = mergedColors))
mclusters <- mclusters[order(mclusters$group), ]
write.csv(mclusters , file = "merged_clusters.csv")



if(traits == "lcm"){

  trait_file <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/lcm_traits.csv", header = TRUE)

} else {

  trait_file <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/hd_zh_traits.csv", header = TRUE)
}

trait_file[,"X"] <- NULL
trait_file[, "x"] <- NULL

nGenes <- ncol(transpose)
nSamples <- nrow(transpose)
#tissue names
tissues <- rownames(transpose)

#using previously calculated eigengenes

MEs = moduleEigengenes(transpose, moduleColors)$eigengenes
orderedMEs <- orderMEs(MEs)

moduleTraitCor <- cor(orderedMEs, trait_file, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

png("/cbcb/project-scratch/ZCL/wgcna_consensus/tissue_correlation.png", height = 2000, width = 3000)

labeledHeatmap(Matrix = moduleTraitCor, colorLabels = FALSE, xLabels = names(trait_file), yLabels = paste(names(orderedMEs), ".", sep =""), colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1,1), main = paste("Module-trait relationships"))

dev.off()

print("complete")

