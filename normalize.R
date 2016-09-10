
library("methods")
library("edgeR")
library("DESeq2")

selection <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/hd_lcm_zh.csv")

#reassign rownames
rownames(selection) <- selection[, "X"]
selection[, "X"] <- NULL
selection[,"X.1"] <- NULL

head(selection)


#read in annotation
annotation <- read.csv("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/ordered_annots.csv")

# print("here1")
rld <- varianceStabilizingTransformation(as.matrix(selection), blind=FALSE)

l <- log2(selection + 1)
print("here2")
cpm_vals <- cpm(selection, annotation[, "Length"])
print("here3")
cpm_log_vals <- cpm(selection, annotation[, "Length"], log=TRUE)
print("here4")
rpkm_vals <- rpkm(selection, annotation[, "Length"])
print("here5")
rpkm_log_vals <- rpkm(selection, annotation[, "Length"], log=TRUE)

write.csv(l, file= "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/log_raw_hd_lcm_zh.csv")
write.csv(rld, file= "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/rld_hd_lcm_zh.csv")
write.csv(cpm_vals, file= "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/cpm_hd_lcm_zh.csv")
write.csv(cpm_log_vals, file= "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/lcpm_hd_lcm_zh.csv")
write.csv(rpkm_vals, file= "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/rpkm_hd_lcm_zh.csv")
write.csv(rpkm_log_vals, file= "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/lrpkm_hd_lcm_zh.csv")


