# library("Rsamtools")
# library("GenomicFeatures")
# library("GenomicAlignments")
# library("BiocParallel")
library("Rsubread")
library("limma")
library("edgeR")
library("methods")
library("WGCNA")

bams <- read.delim("/cbcb/project-scratch/ZCL/Zhongchi3383_mapped_reads_tophat/bam_paths.txt" ,  sep = "\n", header = FALSE)
filenames <- as.character(bams$V1)

fc <- featureCounts(files = filenames, annot.ext="/cbcb/lab/smount/ZCL/Fragaria_vesca_v2.0.a1.transcripts.gtf", nthreads=4, isGTFAnnotationFile=TRUE, useMetaFeatures=TRUE, GTF.attrType="gene_id", GTF.featureType="CDS")

x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

x_rpkm <- rpkm(x, x$genes$Length)

pdf("/cbcb/lab/smount/ZCL/bioconductor_scripts/graphs/rpkm.pdf",width=6,height=4,paper='special')
dev.off()

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



