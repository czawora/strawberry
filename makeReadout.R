args<-commandArgs(TRUE)
network_name <- as.character(args[1])

##creating browsable csv file for network
suppressMessages(library(topGO))
suppressMessages(library(WGCNA))

#network_name <- "hd_zh_sub_log_power_8_modulesize_40_TOM_0_merge_0"
network_dir <- paste("/cbcb/project-scratch/ZCL/wgcna/", network_name, sep ="")
out_name <- paste("OUT_", network_name, ".csv", sep = "")

#get network info
sub_split <- strsplit(network_name, "sub_")[[1]][2]
set_split <- strsplit(network_name, "sub_")[[1]][1]
x_split <- strsplit(sub_split, "_")[[1]]

n_set <- strsplit(set_split, "_")[[1]][1]
n_norm <- x_split[1]
n_mms <- as.numeric(x_split[5])
n_merge <- as.numeric(x_split[9])
n_power <- as.numeric(x_split[3])


print("reading transpose")
if(n_set == "hd"){
	hd_data <- read.csv(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/hd_zh_sub_", n_norm, ".csv", sep=""), header=TRUE)
} else {
	hd_data <- read.csv(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/lcm_sub_", n_norm, ".csv", sep=""), header=TRUE)
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

#read network clusters
print("reading clusters")
if(n_merge == 1){
	clust <- read.delim(paste(network_dir, "/merged_clusters.csv", sep =""), sep =",", header = TRUE)
} else {
	clust <- read.delim(paste(network_dir, "/original_clusters.csv", sep =""), sep =",", header = TRUE)
}

clust$X <- NULL
clust$renamed_genes <- "x"
clust$renamed_genes <- as.character(clust$renamed_genes)
clust$genes <- as.character(clust$genes)


for(i in seq(1,nrow(clust))){
	newname <- strsplit(clust[i, "genes"], "-v1.0")[[1]][1]
	clust[i, "renamed_genes"] <- newname 
}

rownames(clust) <- clust$genes

#split gene names on this ->>>>  -v1.0

#read in traits
print("reading traits")
if(n_set == "hd"){
	trait_file <- read.delim("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/hd_zh_traits.csv",sep=",", header = TRUE)
} else {
	trait_file <- read.delim("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/lcm_traits.csv", sep=",", header = TRUE)
}
trait_file[,"X"] <- NULL
trait_file[, "x"] <- NULL

convert <- read.delim("/cbcb/lab/smount/ZCL/FVE_ATH.csv", header = TRUE, sep = ",")
rownames(convert) <- convert$gene_id


clust$FVE_gene <- "x"
clust$FVE_gene <- convert[clust$renamed_genes, "FV_gene"]

clust$orthos <- "x"
clust$orthos <- convert[clust$renamed_genes, "orthos"]

clust$orthos_names <- "x"
clust$orthos_names <- convert[clust$renamed_genes, "orthos_name"]


###now get to work
moduleColors <- clust$group
MEs = moduleEigengenes(transpose, moduleColors)$eigengenes
orderedMEs <- orderMEs(MEs)
moduleTraitCor <- cor(orderedMEs, trait_file, use = "p")

print("reading GO annots")
#read in go mappings
gene2GO <- readMappings(file = "/cbcb/lab/smount/ZCL/go_annots/gene2GO.txt")

GO2geneID <- inverseList(gene2GO)


#read in FVE gene descriptions
FVE_desc <- read.delim("/cbcb/lab/smount/ZCL/gene_products.csv", header = TRUE, sep = ",")
rownames(FVE_desc) <- FVE_desc$fve_id


out_df <- data.frame(cluster = character(), size = double(), traitCor = character(), topGO = character(), topGO_description = character(), topGOgenes = character(), topGOgenes_product = character(), topGOgenes_orthologue= character(), hubGenes = character(), hubGenes_products = character(), hubGenes_orthologues = character(), clusterGenes = character(), clusterProducts = character(), arabidopsis = character())

#run through each module
for(mod in unique(clust$group)){

	if(mod != "grey"){

		GOterms <- data.frame(GO.ID = character(), Term = character(), pval= double())

		print(mod)

		clust_size <- nrow(clust[which(clust$group == mod), ])

		print(paste("cluster size = ", clust_size,sep =""))

		this_mod <- clust[which(clust$group == mod), ]

		##get the 2 highest and 2 lower correlated traits

		traitsCors <- moduleTraitCor[paste("ME", mod, sep=""), ]
		traitsCors <- traitsCors[order(traitsCors)]

		bottom2 <- traitsCors[1:2]
		bottom2_names <- names(bottom2)

		l_traitsCors <- length(traitsCors)
		one_minus_length <- l_traitsCors - 1
		top2 <- traitsCors[one_minus_length:l_traitsCors]
		top2_names <- names(top2)

		trait_str <- paste(top2_names[2], ":", top2[2], "  ", top2_names[1], ":", top2[1], "  ", bottom2_names[2], ":", bottom2[2], "  ", bottom2_names[1], ":", bottom2[1], sep = "")

		#get gene products for this cluster
		this_mod$product <- FVE_desc[this_mod$FVE_gene, "desc"]
		this_mod$name_product <- paste(this_mod$renamed_genes, this_mod$product, sep = ":")
		print(colnames(this_mod))

		##get the GO terms for the module
		BP <-  paste("/cbcb/project-scratch/ZCL/go_enrich_results/", network_name, "/",mod, "_BP_fisher.txt", sep = "")
		MF <-  paste("/cbcb/project-scratch/ZCL/go_enrich_results/", network_name, "/", mod, "_MF_fisher.txt", sep = "")
		CC <-  paste("/cbcb/project-scratch/ZCL/go_enrich_results/", network_name, "/", mod, "_CC_fisher.txt", sep = "")

		res.BP <- read.delim(BP, sep ="\t", header=TRUE)
		res.MF <- read.delim(MF, sep ="\t", header=TRUE)
		res.CC <- read.delim(CC, sep ="\t", header=TRUE)

		GOterms <- rbind(GOterms, data.frame(GO.ID = res.BP$GO.ID, Term=res.BP$Term, pval = res.BP$adj_classic))
		GOterms <- rbind(GOterms, data.frame(GO.ID = res.MF$GO.ID, Term=res.MF$Term, pval = res.MF$adj_classic))
		GOterms <- rbind(GOterms, data.frame(GO.ID = res.CC$GO.ID, Term=res.CC$Term, pval = res.CC$adj_classic))

		GOterms <- GOterms[order(GOterms$pval), ]

		topGO_term <- as.character(GOterms[1,"GO.ID"])
		topGO_desc <- as.character(GOterms[1,"Term"])

		#print(topGO_term)
		topGenes <- NULL
		topGenes <- GO2geneID[[topGO_term]]

		topGenes_r <- NULL
		topGenes_r <- this_mod[topGenes, c("genes", "renamed_genes", "name_product")]
		topGenes_r <- topGenes_r[which(is.na(topGenes_r$genes) != TRUE), "renamed_genes"]


		GOcount = 2

		while(is.null(topGenes_r )== TRUE || length(topGenes_r) == 0){
			
			topGO_term <- as.character(GOterms[GOcount,"GO.ID"])
			topGO_desc <- as.character(GOterms[GOcount,"Term"])

			topGenes <- GO2geneID[[topGO_term]]

			GOcount = GOcount + 1

			print("topGO term genes")
			#get genes with top GO term
			topGenes_r <- this_mod[topGenes, c("genes", "renamed_genes", "name_product")]
			topGenesFVE <- topGenes_r$name_product
			topGenes_r <- topGenes_r[which(is.na(topGenes_r$genes) != TRUE), "renamed_genes"]

		}

		if(is.null(topGenes_r )== FALSE){
			topGO_arab <- this_mod[topGenes_r, "orthos_names"]

			if(length(topGO_arab) != 0){
				topGO_arab <- topGO_arab[which(topGO_arab != "")]
				topGO_arab_str <- paste(topGO_arab, collapse=" ")
			#however not all those genes are in this cluster
			#topGenes <- topGenes[!is.na(topGenes)]
			} else{
				topGO_arab_str <- ""
			}

			topGenes_str <- paste(topGenes_r,  collapse=" ")

		} else {
			topGenes_str <- ""
			topGO_arab_str <- ""

		}

		#topGenesFVE string
		topGenesFVE_string <- paste(as.character(topGenesFVE[!is.na(topGenesFVE)]), collapse = " ")
		
		##got the top genes and GO term for this module

		##now need the hub genes

		#slice transpose

		clust_genes <- clust[which(clust$group == mod), "genes"]
		sub_transpose <- transpose[, clust_genes]

		print("here")
		connect <- softConnectivity(sub_transpose, type="signed", verbose=1)
		print("here2")
		names(connect) <- colnames(sub_transpose)

		rev_connect <- rev(connect[order(connect)])

		if(length(rev_connect) > 20){
			rev_connect <- rev_connect[1:20]
		}

		#get top hub gene orthologues
		print("here4")
		connect_orthos <- clust[names(rev_connect), "orthos_names"]

		connect_orthos_str <- paste(connect_orthos, collapse = " ")

		hubs_str <- paste(clust[names(rev_connect[1]), "renamed_genes"], ":", rev_connect[1], sep ="")
		count = 2
		for(thing in clust[names(rev_connect[2:length(rev_connect)]), "renamed_genes"]){

			hubs_str <- paste(hubs_str, "  ", thing, ":", rev_connect[count], sep ="")
			count = count + 1
		}

		connect_products <- this_mod[names(rev_connect), "name_product"]
		connect_products_str <- paste(as.character(connect_products), collapse = " ")

		#got hub gene ranking

		#now stringify the arabidopsis orthologues

		clust_orthos <- clust[which(clust$group == mod),  "orthos_names"]
		#clust_orthos <- clust_orthos[which(clust_orthos != "")]

		orthos_str <- paste(clust_orthos, collapse=";")

		print("here5")
		#all genes in cluster with their product 
		print(this_mod$FVE_gene)
		
		prods <- paste(as.character(this_mod$name_product) , collapse = " ")
		print(prods)

		#all genes in cluster
		cG <- paste(as.character(this_mod[, "renamed_genes"]), collapse = " ")

		#out_df <- data.frame(cluster = character(), traitCor = character(), topGO = character(), topGOgenes = character(), hubGenes = character(), arabidopsis = character())

		print("here6")
		out_df <- rbind(out_df, data.frame(cluster = mod,size = clust_size, traitCor = trait_str, topGO = topGO_term, topGO_description = topGO_desc, topGOgenes = topGenes_str, topGOgenes_product = topGenesFVE_string,  topGOgenes_orthologue = topGO_arab_str, hubGenes = hubs_str, hubGenes_products = connect_products_str, hubGenes_orthologues = connect_orthos_str, clusterGenes = cG, clusterProducts = prods, arabidopsis = orthos_str))
		print("here7")
	}
}

write.csv(out_df, file = paste("/cbcb/project-scratch/ZCL/csv_out/", network_name, ".csv", sep = ""), row.names = FALSE)


