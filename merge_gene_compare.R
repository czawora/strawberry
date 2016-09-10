
library("pheatmap")
library("RColorBrewer")


networks <- as.character(read.delim("/cbcb/project-scratch/ZCL/comparison_mats/top50.txt", 'r', header=FALSE)$V1)

network_num <- length(networks)

comparison_table <- as.data.frame(matrix(nrow = network_num, ncol = network_num))

rownames(comparison_table) <- networks
colnames(comparison_table) <- networks

out_file <- as.character(read.delim("/cbcb/project-scratch/ZCL/comparison_mats/gene_comp/results/out_files.txt", header = FALSE, sep="\n")$V1)


for(o1 in out_file){

	df <- read.csv(paste("/cbcb/project-scratch/ZCL/comparison_mats/gene_comp/results/", o1, sep = ""), header = TRUE)

	x <- as.character(df$x)
	y <- as.character(df$y)
	score <- as.double(df$score)

	print("-------")
	print(x)
	print(y)
	print(score)

	comparison_table[x,y] <- score
	comparison_table[y,x] <- score
}

write.csv(comparison_table, file = "/cbcb/project-scratch/ZCL/comparison_mats/gene_comp/results/comp_table.csv")

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png(file = "/cbcb/project-scratch/ZCL/comparison_mats/gene_comp/results/comp_heat.png", height= 1000, width=1000)
pheatmap(as.matrix(comparison_table), clustering_distance_rows="correlation", clustering_distance_cols="correlation",col=rev(colors))
dev.off()