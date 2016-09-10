
#setup
args<-commandArgs(TRUE)

#gather command line arguments
if (length(args) < 4){
   softPower <- 1
   minModuleSize <- 50
   TOM_on <- 0
   filename <- "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/hd_zh_sub_lrpkm.csv"
   traits <- "hd_zh"
   merge_eigengenes <- 1
  #stop()
} else {
  softPower <- as.numeric(args[1])
  minModuleSize <- as.numeric(args[2])
  TOM_on <- as.numeric(args[3])
  traits <- as.character(args[6])
  filename <- as.character(args[5])
  merge_eigengenes <- as.numeric(args[4])
}

#create directory for output if doesnt exist
new_folder <- strsplit(strsplit(filename, "\\/")[[1]][9], "\\.")[[1]][1]
dir_name <- paste("/cbcb/project-scratch/ZCL/wgcna/", new_folder, "_power_", softPower , "_modulesize_", minModuleSize , "_TOM_" , TOM_on , "_merge_", merge_eigengenes, "/", sep = "")

dir.create(dir_name)

#fail marker
file.create(paste("/cbcb/project-scratch/ZCL/wgcna/failed/", new_folder, "_power_", softPower , "_modulesize_", minModuleSize , "_TOM_" , TOM_on , "_merge_", merge_eigengenes, sep =""))

setwd(dir_name)

write.table(args, file = "config.txt")

#redirect output
sink(file = "output.txt")

suppressMessages(library("WGCNA"))
suppressMessages(library("methods"))
suppressMessages(library("edgeR"))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
nthr <- enableWGCNAThreads()
print(nthr)

print("importing data")
selection <- read.csv(filename, header=TRUE)

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

print("making adjacency matrix")
adjacency_res <- adjacency(transpose, power = softPower, type="signed")

dim(adjacency_res)

if (TOM_on == 1){
  # Turn adjacency into topological overlap
  print("making TOMsimilarity")
  TOM <- TOMsimilarity(adjacency_res, TOMType ="signed", verbose =2)
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
l <- length(unique(dynamicColors))
l
table(dynamicColors)

if(l > 200){

        file.remove(paste("/cbcb/project-scratch/ZCL/wgcna/failed/", new_folder, "_power_", softPower , "_modulesize_", minModuleSize , "_TOM_" , TOM_on , "_merge_", merge_eigengenes, sep =""))
        print("too many clusters")
        file.create(paste("/cbcb/project-scratch/ZCL/wgcna/2many/",  new_folder, "_power_", softPower , "_modulesize_", minModuleSize , "_TOM_" , TOM_on , "_merge_", merge_eigengenes, sep =""))

        quit(save = "no")
}

#write out clusters
gene_names <- rownames(adjacency_res)
clusters <- as.data.frame(cbind(genes= gene_names, group = dynamicColors))
clusters <- clusters[order(clusters$group), ]

write.csv(clusters , file = "original_clusters.csv")

#plot dendogram
print("plotting dendrogram")
png(paste(new_folder, "_power_", softPower , "_modulesize_", minModuleSize , "_TOM_" , TOM_on , "_merge_", merge_eigengenes, "_dendro.png", sep = ""),width=1500,height=1500)

# png("x.png", ,width=1500,height=1500)
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,addGuide = TRUE, hang=0.3, guideHang = 0.05, main = "Gene dendrogram and module colors")

dev.off()

moduleColors <- dynamicColors
#if we want to merge cluster based on eigengenes
if (merge_eigengenes == 1){

  # Calculate eigengenes
  MEList <- moduleEigengenes(transpose, colors = dynamicColors)
  MEs <- MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss <- 1-cor(MEs)
  # Cluster module eigengenes
  METree <- hclust(as.dist(MEDiss), method = "average")
  # Plot the result

  png("clustered_me.png", width =1500, height=1500)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")

  MEDissThres <- 0.25
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")

  dev.off()

  # Call an automatic merging function
  merge <- mergeCloseModules(transpose, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  moduleColors <- mergedColors
  print(table(moduleColors))
  print(length(unique(moduleColors)))
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;

  png(paste(new_folder, "_power_", softPower , "_modulesize_", minModuleSize , "_TOM_" , TOM_on , "_merge_", merge_eigengenes, "_merged_dendro.png", sep = ""),width=1500,height=1500)
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

}

##identifying eigengene correlation with tissue

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

png("tissue_correlation.png", height = 2000, width = 3000)

labeledHeatmap(Matrix = moduleTraitCor, colorLabels = FALSE, xLabels = names(trait_file), yLabels = paste(names(orderedMEs), ".", sep =""), colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1,1), main = paste("Module-trait relationships"))

dev.off()

# geneModuleMembership <- as.data.frame(cor(transpose, orderedMEs, user = "p"))
# MMPvalue <- as.data.frame(corPvalueStudents(as.matrix(geneModuleMembership), nSamples))

print("complete")

#stop output redirection
sink()

file.remove(paste("/cbcb/project-scratch/ZCL/wgcna/failed/", new_folder, "_power_", softPower , "_modulesize_", minModuleSize , "_TOM_" , TOM_on , "_merge_", merge_eigengenes, sep =""))
file.create(paste("/cbcb/project-scratch/ZCL/wgcna/success/", new_folder, "_power_", softPower , "_modulesize_", minModuleSize , "_TOM_" , TOM_on , "_merge_", merge_eigengenes, sep =""))

