
library(stats)
suppressMessages(library(ggbiplot))
suppressMessages(library(cbcbSEQ))
# suppressMessages(library("WGCNA"))
# suppressMessages(library("methods"))
# suppressMessages(library("edgeR"))
# suppressMessages(library(gplots))
# suppressMessages(library(RColorBrewer))

setwd("/Users/Chris/Documents/strawberry/ZCL/counts/")

# count_files <- c("rpkm_hd_lcm_zh.csv", "rld_hd_lcm_zh.csv", "lrpkm_hd_lcm_zh.csv",
#                  "log_raw_hd_lcm_zh.csv", "lcpm_hd_lcm_zh.csv", "cpm_hd_lcm_zh.csv",
#                  "hd_lcm_zh.csv")

#count_files <- c( "rld_hd_lcm_zh.csv", "lrpkm_hd_lcm_zh.csv",
                 #"log_raw_hd_lcm_zh.csv", "lcpm_hd_lcm_zh.csv")

count_files <- c("rpkm_hd_lcm_zh.csv")
groups <- c("hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd","hd","hd","hd","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","zh","zh","zh","zh","zh","zh","zh","zh")

dates <- read.delim(file = "/Users/Chris/Documents/strawberry/ZCL/counts/sample_dates.txt", sep="\t")

for(read_file in count_files){

  print(read_file)

 # read_file <- paste("intron_", read_file, sep ="")
  input <- read.csv(read_file)

  if (read_file == "hd_lcm_zh.csv"){
    input$X.1 <- NULL
  }

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

   #d <- d[ which(d$groups=="zh" | d$groups == "hd"),  ]
  print(d$batch)

#   lookups <- c("July_7_2011", "Sept_30_2011", "Nov_2_2011", "Nov_30_2011", "Jan_4_2012")

#   for (lookup in lookups){
#     print(lookup)
  b <- d
#   (final_groups <- as.character(sapply(b$batch, FUN, lookup)))
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

  print("plotting")
  for(i in c(1)){
    for(j in c(2)){#1:10){

      if (i != j){

        #png(paste("/Users/Chris/Documents/strawberry/ZCL/graphs/pca/", read_file,"/","pca",i, "_",j, "_labels.png", sep=""), width = 1500, height = 1500)

        g <- ggbiplot(pca, choices = c(i,j), obs.scale = 1, var.scale = 1,
                      groups = final_groups, ellipse = TRUE, labels = rownames(b),
                      circle = TRUE, varname.size=1, labels.size = 5, var.axes=FALSE)
        g <- g + theme(legend.direction = 'horizontal',
                       legend.position = 'top')
        plot(g)

        #dev.off()


      }
    }
  }
  #}
}

