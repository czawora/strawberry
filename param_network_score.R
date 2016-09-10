
args<-commandArgs(TRUE)

#read network to score
network <- as.character(args[1])
#network <- "hd_zh_sub_rld_power_8_modulesize_40_TOM_0_merge_1"

#create fail marker
file.create(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/failed/",  network , sep =""))

setwd(paste("/cbcb/project-scratch/ZCL/go_enrich_results/", network, sep = ""))

sink(file = paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/scores/", network, "_output.txt", sep =""))
print("start")


print("-----------read clusters")
clusters <- read.delim("cluster_counts.txt", header = TRUE, sep =" ")

cluster_names <- rownames(clusters)

GOterms <- data.frame(GO.ID = character(), Term = character())
pvals <- c()

print("-----------loop through clusters")
for(cn in cluster_names){

	print(cn)
	if(cn != "grey"){

		BP <-  paste(cn, "_BP_fisher.txt", sep = "")
		MF <-  paste(cn, "_MF_fisher.txt", sep = "")
		CC <-  paste(cn, "_CC_fisher.txt", sep = "")

		res.BP <- read.delim(BP, sep ="\t", header=TRUE)
		res.MF <- read.delim(MF, sep ="\t", header=TRUE)
		res.CC <- read.delim(CC, sep ="\t", header=TRUE)

		GOterms <- rbind(GOterms, data.frame(GO.ID = res.BP$GO.ID, Term=res.BP$Term))
		GOterms <- rbind(GOterms, data.frame(GO.ID = res.MF$GO.ID, Term=res.MF$Term))
		GOterms <- rbind(GOterms, data.frame(GO.ID = res.CC$GO.ID, Term=res.CC$Term))

		pvals <- c(pvals, res.BP$adj_classic,res.MF$adj_classic,res.CC$adj_classic )


	}

}

print("------------create scores")
uniqueGO <- length(as.character(unique(GOterms$GO.ID)))
sumP <- sum(-1 * log10(pvals))

print("-------------write out")
out_df <- data.frame(net = network, uniqueGO_score = uniqueGO, sum_neg_log_p = sumP)

setwd("/cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/scores")

write.table(out_df, file = paste(network, "_uniqueGO_sumP.txt", sep =""))

print("-------------complete")

sink()
#remove fail marker
file.remove(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/failed/",  network , sep =""))
#create success marker
file.create(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/net_scores/success/",  network , sep =""))
