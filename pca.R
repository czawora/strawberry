
library(ggbiplot)
library(stats)

input <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/lrpkm_hd_lcm_zh.csv")
groups <- c("hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd","hd","hd","hd","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","zh","zh","zh","zh","zh","zh","zh","zh")

dates <- read.delim(file = "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/sample_dates.txt", sep="\t")
dates$combo <- dates$date#paste(dates$date, "_", dates$lane, sep = "")

d<-as.data.frame(t(input[,2:ncol(input)]))

#reassign rownames

colnames(d) <- input[,"X"]

d$groups <- groups

#mtch batch info with row names
matches <- match(rownames(d), dates$sampleID, nomatch=nrow(dates))
#add batch info to data
d$batch <- dates[matches,"combo"]
#rename
#coerce to character
d$batch <- as.character(d$batch)

FUN <- function(row_val, desired_label){

		if(row_val == desired_label){

			return(row_val)
		}
		else{
			return("other")
		}
}

d <- d[ which(d$groups=="hd" | d$groups=="zh"),  ]
print(d$batch)

lookups <- c("July_7_2011", "Sept_30_2011", "Nov_2_2011", "Nov_30_2011", "Jan_4_2012")

# for (lookup in lookups){
	b <- d
# (final_groups <- as.character(sapply(b$batch, FUN, lookup)))
# final_groups <- paste(d$groups, "_", d$batch, sep = "")
final_groups <- b$groups
b$groups <- NULL
b$batch <- NULL

dim(b)

#remove no variance columsn
b <- b[,apply(b, 2, var, na.rm=TRUE) != 0]

pca <- prcomp(b,
                 center = TRUE,
                 scale. = TRUE) 

# print(pca)

# png("/cbcb/lab/smount/ZCL/bioconductor_scripts/pca/components.png")


# plot(pca, type = "l")

# dev.off()

# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

for(i in c(1)){
	for(j in 1:10){

		if (i != j){

			png(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/pca/", lookup,"_pca",i, "_",j, "_labels.png", sep=""), width = 1500, height = 1500)

 			g <- ggbiplot(pca, choices = c(i,j), obs.scale = 1, var.scale = 1, 
              groups = final_groups, ellipse = TRUE, labels = rownames(b),
              circle = TRUE, varname.size=1, var.axes=FALSE)
			g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
			plot(g)

			dev.off()


		}
	}
}
#}
