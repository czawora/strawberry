## automated go enrichment tests over each cluster of every output specified

suppressMessages(library(topGO))
suppressMessages(library(WGCNA))

 allowWGCNAThreads()

 #setup
args<-commandArgs(TRUE)

run <- as.character(args[2])
mod <- as.character(args[1])
set <- as.character(args[3])
norm <- as.character(args[4])
merged <- as.numeric(args[5])

# run <- "hd_zh_sub_lrpkm_power_8_modulesize_90_TOM_0_merge_0"
# mod <- "blue"
# set <- "hd"
# norm <- "lrpkm"

file.create(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/go_enrich/failed/",  mod, "_", run  , sep =""))

############################
#read in count data to reconstruct transpose again

if(set == "hd"){
	hd_data <- read.csv(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/hd_zh_sub_", norm, ".csv", sep=""), header=TRUE)
} else {
	hd_data <- read.csv(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/lcm_sub_", norm, ".csv", sep=""), header=TRUE)
}
#set rownames
rownames(hd_data) <- hd_data[,"X"]
#remove useless column
hd_data[,"X"] <- NULL
hd_data[,"X.1"] <- NULL

#transpose with only number to keep numeric
transpose <- t(hd_data)

#remove non variance columsn
transpose <- transpose[,apply(transpose, 2, var, na.rm=TRUE) > 0.05]
############################


###read in trait file
if(set == "hd"){
	trait_file <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/hd_zh_traits.csv", header = TRUE)
} else {
	trait_file <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/lcm_traits.csv", header = TRUE)
}
trait_file[,"X"] <- NULL
trait_file[, "x"] <- NULL

setwd("/cbcb/project-scratch/ZCL/wgcna/")

##check if merged modules or not
# merged <- strsplit(run, "_")[[1]][length(strsplit(run, "_")[[1]])]

print("merge determination")
print(merged)

#read in cluster assignemnts
if(merged == 1){

	clust <- read.delim(paste(run, "/merged_clusters.csv", sep=""), sep= ",", header = TRUE)

} else {

	clust <- read.delim(paste(run, "/original_clusters.csv", sep=""), sep= ",", header = TRUE)

}

##set to output directory
setwd("/cbcb/project-scratch/ZCL/go_enrich_results/")

dir.create(run)

setwd(paste("/cbcb/project-scratch/ZCL/go_enrich_results/", run, sep=""))

sink(file = paste(mod, "_output.txt", sep =""))

print("clusters loaded")

################################
##calculate module eigengene correlations with traits
moduleColors <- clust$group
MEs = moduleEigengenes(transpose, moduleColors)$eigengenes
orderedMEs <- orderMEs(MEs)
moduleTraitCor <- cor(orderedMEs, trait_file, use = "p")
###############################

if(!file.exists("cluster_counts.txt")){
	write.table(rev(sort(table(clust$group))), file = "cluster_counts.txt")
}

clust$X <- NULL

##read in gene names
##this is crucial, the gene universe cannot be obtained from the GO annotations since 
##only half the genes have annotations, also it cannot be obtained from the cluster files
##because each run of the network pruned out difference sets of genes based on variance
##this list is the definitive list of all gene names in Fragaria Vesca

print("getting gene universe")
gene_n <- read.delim("/cbcb/lab/smount/ZCL/gene_names.txt", header = FALSE)
gene_u <- as.character(gene_n$V1)

print("reading GO annots")
#read in go mappings
gene2GO <- readMappings(file = "/cbcb/lab/smount/ZCL/go_annots/gene2GO.txt")

#read in gene names
gene_name_conver <- read.delim("/cbcb/lab/smount/ZCL/id_conversion.fve.csv", sep=";", header = TRUE)
rownames(gene_name_conver) <- gene_name_conver$id

print(mod)

##write out traits highly correlated with this module
trait_cor_mod_ME <- moduleTraitCor[paste("ME", mod ,sep=""), ]
print("calculated trait correlations")
trait_cor_mod_ME <- trait_cor_mod_ME[order(trait_cor_mod_ME)]
l <- length(trait_cor_mod_ME)
one_minus_length <- l - 1
most_neg <- trait_cor_mod_ME[1:2]
most_pos <- trait_cor_mod_ME[one_minus_length:l]

print("write correlations")

write.table(c(most_neg, most_pos), file = paste(mod, "_cor_traits.txt", sep=""))

#select the genes for the chosen cluster
print("select interesting genes")
interestGenes <- clust[which(clust$group ==  mod), "genes"]

renamed <- c()
#rename genes in the interest set for output use only
for(g in as.character(interestGenes)){

	spl <- strsplit(g, "-")[[1]][1]
	renamed <- c(renamed,spl)
}

print("write out converted names")
conver_names <- as.character(gene_name_conver[renamed, "gene_id"])
write.table(conver_names, file = paste(mod, "_renamed_genes.txt", sep =""), col.names = FALSE, row.names=FALSE)

print("create gList")
#create a vector with 1 for "in cluster" and 0 for "not in cluster"
gList <- factor(as.integer(gene_u %in% interestGenes))

names(gList) <- gene_u

print("create stat test")
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

##try running for each ontology type
for(onto in c("MF", "BP", "CC")){

	print(onto)
	Godata <- new("topGOdata", ontology = onto, allGenes = gList, annot = annFUN.gene2GO, gene2GO = gene2GO)
	print(Godata)

	resultFisher <- getSigGroups(Godata, test.stat)

	print(resultFisher)

	print("allRes")
	allRes <- GenTable(Godata, classic= resultFisher, orderBy= "classic", ranksOf = "classic")

	print("adjust")
	allRes$adj_classic <- p.adjust(allRes$classic, method = "fdr")

	print("subset")
	stat_sig <- allRes[which(allRes$adj_classic <= 0.05), ]

	print("write")
	write.table(stat_sig[order(stat_sig$adj_classic), ], file = paste(mod, "_", onto, "_fisher.txt", sep = ""), sep="\t", row.names=FALSE)

}

print("complete")
#stop output redirection
sink()

file.remove(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/go_enrich/failed/", mod, "_", run  , sep =""))
file.create(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/go_enrich/success/", mod, "_", run  , sep =""))





