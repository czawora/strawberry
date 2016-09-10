setwd("/Users/Chris/Documents/strawberry/ZCL/counts")

files <- c( "rld_hd_lcm_zh.csv", "lrpkm_hd_lcm_zh.csv", "log_raw_hd_lcm_zh.csv", "lcpm_hd_lcm_zh.csv")

for(name in files){

  print(name)
  exons <- read.delim(name, sep = ",", header = TRUE)
  introns <- read.delim(paste("intron_", name, sep=""), sep = ",", header =TRUE)

  rownames(exons) <- exons$X
  rownames(introns) <- introns$X

  exons$X <- NULL
  introns$X <- NULL

  sub_exons <- exons[rownames(introns), ]
  ratio <- sub_exons/introns


  exon_names <- rownames(exons)
  exon_names <- setdiff(exon_names,rownames(introns))

  diff_exons <- as.data.frame(matrix(1000000, nrow=length(exon_names), ncol = 92))
  colnames(diff_exons) <- colnames(exons)
  rownames(diff_exons) <- exon_names

  ratio <- rbind(ratio, diff_exons)
  rownames(ratio) <- paste("r_", rownames(ratio), sep ="")
  rownames(exons) <- paste("e_", rownames(exons), sep ="")

  combined <- rbind(exons, ratio)

  write.csv(combined, file = paste("exon_ei_ratio_", name, sep=""))
}