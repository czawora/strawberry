
library(GenomicRanges)
library(rtracklayer)
library(Rsubread)
library(methods)
library(RIPSeeker)


gffRangedData<-import.gff("/cbcb/lab/smount/ZCL/Fragaria_vesca_v2.0.a1.transcripts_cds_as_introns.gff3")
# print(colnames(myGranges))
myGranges<-as(gffRangedData, "GRanges")

# print(colnames(myGranges))
# class()
# colnames(myGranges) <- c("seqnames", "start", "end" ,"width" ,"strand", "source" ,  "type", "score", "phase", "id")

# exportGRanges(myGranges, "/cbcb/lab/smount/ZCL/bioconductor_scripts/granges.txt", "txt")

# txDB <- makeTxDbFromGRanges(myGranges)
annotation <- createAnnotationFile(myGranges)

print(head(annotation, 10))

