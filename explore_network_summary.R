
library('ggplot2')

net <- read.delim("/Users/Chris/Documents/strawberry/ZCL/summary_scores.csv", sep = ",", header = TRUE)
net$X <- NULL

colnames(net)

ggplot() + aes(net$uniqueGO ) + geom_histogram(binwidth=5*max(net$uniqueGO)/100, colour="black", fill="blue")  + geom_density() + labs(x = "# of unique GO terms", y = "# of networks with unique GO term score" )
5*max(net$uniqueGO)/100
ggplot() + aes(net$pval_score)+ geom_histogram(binwidth=5*max(net$pval_score)/100, colour="black", fill="orange")  + geom_density() + labs(x = "p-value score", y = "# of networks with p-value score" )
5*max(net$pval_score)/100

unique_order <- net[rev(order(net$uniqueGO)), ]
pval_order <- net[rev(order(net$pval_score)), ]

m <- cbind(net$uniqueGO, net$pval_score)
colnames(m) <- c("uniqueGO_rank", "pval_rank")
cor(m[, "uniqueGO_rank"], m[, "pval_rank"], method="kendall", use="pairwise")

#top uGO networks
unique_order[which(unique_order$merge == 1),]
#top pvalnetworks

x <- unique_order$norm
x
x[x != "lcpm"]

cor(net_hd$merge, net_hd$uniqueGO)
cor(net_hd$merge, net_hd$pval_score)


####--------------------------------------------------------
####--------------------------------------------------------
####--------------------------------------------------------
####--------------------------------------------------------

rownames(net) <- net$network




