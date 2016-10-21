# library("Rsamtools")
# library("GenomicFeatures")
# library("GenomicAlignments")
# library("BiocParallel")

# .libPaths("/cbcb/lab/smount/programs/R-3.1.2/library")
library("Rsubread")
library("limma")
library("edgeR")
library("methods")
library("WGCNA")

args <- commandArgs(TRUE)
allowWGCNAThreads()

bams <- read.delim(paste("/cbcb/project-scratch/ZCL/", args[1], "_mapped_reads_tophat/bam_paths.txt", sep = "") ,  sep = "\n", header = FALSE)
filenames <- as.character(bams$V1)

fc <- featureCounts(files = filenames, annot.ext="/cbcb/lab/smount/ZCL/Fragaria_vesca_v2.0.a1.transcripts.gtf", nthreads=4, isGTFAnnotationFile=TRUE, useMetaFeatures=TRUE, GTF.attrType="gene_id", GTF.featureType="CDS")

#(fc$annotation)

# write.csv(fc$counts, file = paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/", args[1], "_counts.csv", sep = ""))
# write.csv(fc$counts_junction, file = paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/", args[1], "_counts_junction.csv", sep = ""))
# write.csv(fc$annotation, file = paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/", args[1], "_annotation.csv", sep = ""))
# write.csv(fc$targets, file = paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/", args[1], "_targets.csv", sep = ""))
# write.csv(fc$stat, file = paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/", args[1], "_stat.csv", sep = ""))


# x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

# x_rpkm <- rpkm(x, x$genes$Length)

# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=30, by=2))
# # Call the network topology analysis function
# sft = pickSoftThreshold(x_rpkm, powerVector = powers, verbose = 5)


# png(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/graphs/", args[1], "_modelfit.png", sep = ""),width=500,height=500)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2))
# cex1 = 0.9


# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
# xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
# main = paste("Scale independence"))
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
# labels=powers,cex=cex1,col="red")
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")

# dev.off()

# png(paste("/cbcb/lab/smount/ZCL/bioconductor_scripts/graphs/", args[1], "_connectivity.png", sep = ""),width=500,height=500)

# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
# xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
# main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# dev.off()

# write.csv(fc$counts, file = "/cbcb/lab/smount/ZCL/bioconductor_scripts/out.csv")
# write.csv(x_rpkm, file = "/cbcb/lab/smount/ZCL/bioconductor_scripts/out_rpkm.csv")


# bamfiles <- BamFileList("/cbcb/project-scratch/ZCL/HD_mapped_reads_tophat/anther10-1_5bp_trim/anther10-1_5bp_trim.bam" , index = character(), yieldSize = 650000000)
# seqinfo(bamfiles[1])

# gfffile <- "/cbcb/lab/smount/ZCL/Fragaria_vesca_v2.0.a1.transcripts.gff3"

# (txdb <- makeTxDbFromGFF(gfffile, format="gff3", circ_seqs=character()))
# (ebg <- exonsBy(txdb, by="gene"))

# register(SerialParam())

# se <- summarizeOverlaps(features=ebg, reads=bamfiles,
#                         mode="Union",
#                         singleEnd=TRUE,
#                         ignore.strand=TRUE)

# countdata <- assay(se)

# head(countdata, 20)

# colSums(countdata)


#write.csv(assay(se), file = "/cbcb/project-scratch/ZCL/HD_mapped_reads_tophat/count_matrix.csv")



