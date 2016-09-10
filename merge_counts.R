


hd <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/HD_counts.csv")
lcm <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/LCM_counts.csv")
zh <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/Zhongchi3383_counts.csv")

tmp <- merge(hd, lcm, by="X")
final <- merge(tmp, zh, by="X")

write.csv(final, file = "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/hd_lcm_zh.csv")