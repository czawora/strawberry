hd_names <- as.character(read.csv("/Users/Chris/Documents/strawberry/ZCL/hd_zh_names.txt", header=FALSE)$V1)
lcm_names <- as.character(read.csv("/Users/Chris/Documents/strawberry/ZCL/lcm_names.txt", header= FALSE)$V1)


count_files <- c("exon_intron_lcpm_hd_lcm_zh.csv", "exon_intron_log_raw_hd_lcm_zh.csv", "exon_intron_lrpkm_hd_lcm_zh.csv", "exon_intron_rld_hd_lcm_zh.csv")

for(file in count_files){

  print(file)
  setwd("/Users/Chris/Documents/strawberry/ZCL/counts/")

  d <- read.csv(file)

  #create new filename
  splits <- strsplit(file, "_")
  norm_prefix <- splits[[1]][3]
  hd_filename = paste("ei_hd_zh_sub_", norm_prefix, ".csv", sep = "")
  lcm_filename =  paste("ei_lcm_sub_", norm_prefix, ".csv", sep = "")

#   if (file == "hd_lcm_zh.csv"){
# #     d$X.1 <- NULL
#     hd_filename = paste("intron_hd_zh_sub", "_raw", ".csv", sep = "")
#     lcm_filename = paste("intron_lcm_sub", "_raw", ".csv", sep = "")
#   }

  #subset
  hd_zh_sub <- d[,c("X",hd_names)]
  lcm_sub <- d[,c("X",lcm_names)]

  setwd("/Users/Chris/Documents/strawberry/ZCL/counts/subsets/")

  write.csv(hd_zh_sub , file = hd_filename)
  write.csv(lcm_sub, file = lcm_filename)

}