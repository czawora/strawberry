###computing topGO functional entihcment analysis for some hd exon only result sets

setwd("/cbcb/lab/smount/ZCL/bioconductor_scripts/wgcna/")
selec <- read.delim("hd_zh_sub_lrpkm_power_1_modulesize_90_TOM_0_merge_1/merged_clusters.csv", sep=",", header =TRUE)

selec$X <- NULL
