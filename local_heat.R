library("pheatmap")
library("RColorBrewer")


setwd("/Users/Chris/Documents/strawberry/ZCL/counts/")

count_files <- c("rpkm_hd_lcm_zh.csv", "rld_hd_lcm_zh.csv", "lrpkm_hd_lcm_zh.csv",
                 "log_raw_hd_lcm_zh.csv", "lcpm_hd_lcm_zh.csv", "cpm_hd_lcm_zh.csv",
                 "hd_lcm_zh.csv")

count_files <- c("rld_hd_lcm_zh.csv",  "lrpkm_hd_lcm_zh.csv", "log_raw_hd_lcm_zh.csv", "lcpm_hd_lcm_zh.csv")

for(read_file in count_files){

  read_file <- "log_raw_hd_lcm_zh.csv"
  print(read_file)

  intron_read_file <- paste("intron_", read_file, sep ="")
  introns <- read.csv(intron_read_file)
  input <- read.csv(read_file)

  if (read_file == "hd_lcm_zh.csv"){
    input$X.1 <- NULL
  }

  #change rownames to first colomn
  rownames(input) <- input$X
  input$X <- NULL

  rownames(introns) <- introns$X
  introns$X <- NULL
  rownames(introns) <- paste("I_", rownames(introns), sep="")

  #now transpose briefly
  t_input <- as.matrix(t(input))
  t_introns <- as.matrix(t(introns))

  t_input <- merge(t_input, t_introns, by.x="row.names", by.y="row.names", sort=FALSE)

  rownames(t_input) <- t_input$Row.names
  t_input$Row.names <- NULL

#  rownames(t_input) <- NULL

  print(t_input[1:5,1:5])

  dim(t_input)

  #calculate correlation
   c <- cor(t_input)


#   #calcualte distance
  d <- dist(t_input)
  tree <- hclust(d, method = "average")
  png(paste(read_file, "_sample_hclust.png", sep = ""), width =1500, height=1500)

  plot(tree)

  dev.off()
  dMat <- as.matrix(d)

  #get colors
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

  png(paste("/Users/Chris/Documents/strawberry/ZCL/graphs/heat/",read_file,"_cor.png", sep=""), width=2000, height = 2000)

  pheatmap(c,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         col=rev(colors))

  dev.off()

  png(paste("/Users/Chris/Documents/strawberry/ZCL/graphs/heat/",read_file,"_dist.png", sep=""), width=2000, height = 2000)

  pheatmap(dMat,
         clustering_distance_rows=d,
         clustering_distance_cols=d,
         col=colors)

  dev.off()

}
