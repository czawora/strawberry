suppressMessages(library(WGCNA))
nthr <- enableWGCNAThreads()

args<-commandArgs(TRUE)
n1_name <- as.character(args[1])
n1_merge <- as.character(args[2])
n2_name <- as.character(args[3])
n2_merge <- as.character(args[4])
out_dir <- as.character(args[5])

file.create(paste(out_dir, "/failed/", n1_name, "$$$$$", n2_name, sep=""))

if(n1_merge == "1"){
  n1 <- read.delim(paste("/cbcb/project-scratch/ZCL/wgcna/", n1_name, "/merged_clusters.csv", sep= ""), sep =",", header = TRUE)
}else{
  n1 <- read.delim(paste("/cbcb/project-scratch/ZCL/wgcna/", n1_name, "/original_clusters.csv", sep= ""), sep =",", header = TRUE)
}

if(n2_merge == "1"){
  n2 <- read.delim(paste("/cbcb/project-scratch/ZCL/wgcna/", n2_name, "/merged_clusters.csv", sep= ""), sep =",", header = TRUE)
}else{
  n2 <- read.delim(paste("/cbcb/project-scratch/ZCL/wgcna/", n2_name, "/original_clusters.csv", sep= ""), sep =",", header = TRUE)
}


# n1 <- read.delim("/cbcb/project-scratch/ZCL/wgcna/hd_zh_sub_lcpm_power_12_modulesize_120_TOM_0_merge_0/original_clusters.csv", sep = ",", header=TRUE)
# n2 <- read.delim("/cbcb/project-scratch/ZCL/wgcna/hd_zh_sub_lcpm_power_12_modulesize_120_TOM_0_merge_0/original_clusters.csv", sep = ",", header=TRUE)
#
# out_dir <- "/cbcb/project-scratch/ZCL/comparison_mats/2016_08_08__03_24_26/results"

n1 <- n1[order(n1$genes), ]
n2 <- n2[order(n2$genes), ]

n1 <- n1[which(n1$group != "grey"), ]
n2 <- n2[which(n2$group != "grey"), ]

rownames(n1) <- n1$genes
rownames(n2) <- n2$genes

common_genes <- intersect(rownames(n1), rownames( n2))
common_n1 <- n1[common_genes, ]
common_n2 <- n2[common_genes, ]

n1_clusters <- unique(common_n1$group)
print(length(n1_clusters))

n2_clusters <- unique(common_n2$group)
print(length(n2_clusters))

print("making overlapTable")
o <- overlapTable(common_n1$group, common_n2$group)$countTable

rowMaxes <- apply(o, 1, max)

#get cluster conserved percentages
for(clust in n1_clusters){
  l <- nrow(common_n1[which(common_n1$group == clust), ])
  rowMaxes[clust] <- rowMaxes[clust]/as.double(l)
}

print("print rowMaxes")
print(rowMaxes)

network1_conservation <- sum(rowMaxes)/length(rowMaxes)
print("net1 conservation score")
print(network1_conservation)



rowMaxes <- NULL
#flip the countTable and calculate conservation for other network
o <- t(o)

rowMaxes <- apply(o, 1, max)

#get cluster conserved percentages
for(clust in n2_clusters){
  l <- nrow(common_n2[which(common_n2$group == clust), ])
  rowMaxes[clust] <- rowMaxes[clust]/as.double(l)
}

print("print rowMaxes")
print(rowMaxes)

network2_conservation <- sum(rowMaxes)/length(rowMaxes)

if(is.na(network2_conservation) == TRUE){
  if(is.na(network1_conservation) == TRUE){
  quit(save = "no")
  }
}
print("net2 conservation score")
print(network2_conservation)

df1 <- data.frame(checking_network = n1_name, conservation_in = n2_name, score = network1_conservation)
df2 <- data.frame(checking_network = n2_name, conservation_in = n1_name, score = network2_conservation)

out_df <- rbind(df1, df2)

print(out_dir)
#print(paste(out_dir,"/", n1, "_", n2, ".csv", sep=""))
write.csv(out_df, file = paste(out_dir,"/", n1_name, "_", n2_name, ".csv", sep=""))


#this will procuce the result of 1st argument labels in rows and the intersection counts against the 2nd argument labels in the columns. The key here is to measure how much the sub-clusters of the larger cluster network are convserved in the more numerous cluster network. This will effectively give the similarity of the two networks. So the first argument should be the network with more total networks
#It would not make sense to see if the more fractured clusters became subclusters of the larger clusters. This is to be expected in some combination by inducing larger clusters. What is more informative between two networks is, as a large cluster network is "split" into greater, smaller clusters, do the sub-networks of the large cluster stay together?
file.remove(paste(out_dir, "/failed/", n1_name, "$$$$$", n2_name, sep=""))


file.create(paste(out_dir, "/success/", n1_name, "$$$$$", n2_name, sep=""))





