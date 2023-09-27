#####################################################################
#                             Fig6K                                 #
#####################################################################
#                             Fig7K                                 #
#####################################################################


library(Rsamtools)
library(GenomicRanges)

read.BAM <- function(fn, ext=0) {
  what=c("rname", "strand", "pos", "qwidth")
  TSS.counts=NULL
  param=ScanBamParam(what = what)
  bam=scanBam(fn, param=param)[[1]]
  ix= !is.na(bam$rname) & !is.na(bam$pos)
  if(ext==0) {
    qwidth=bam$qwidth[ix]
    IRange.reads=GRanges(seqnames=Rle(bam$rname[ix]),
                         ranges=IRanges(bam$pos[ix], width=bam$qwidth[ix]))
  } else {
    idx.minus=bam$strand=="-"
    idx.minus=idx.minus[ix]
    ss=bam$pos[ix]; ee=bam$pos[ix]+bam$qwidth[ix]
    ss2=ss
    ss2[idx.minus]=ee[idx.minus]-ext
    IRange.reads=GRanges(seqnames=Rle(bam$rname[ix]),
                         ranges=IRanges(ss2, width=ext))
  }
  IRange.reads
}

tmp = c("Pu.1_cWT_combined.bam","Pu.1_cKO_combined.bam",
        "H3K4me3_cWT_combined.bam","H3K4me3_cKO_combined.bam",
        "H3K27ac_cWT_combined.bam","H3K27ac_cKO_combined.bam",)
up_DGEs <- read.table(file = "../data/up_tss_2000.bed",sep="\t")
down_DEGs <- read.table(file = "../data/down_tss_2000.bed",sep="\t")

totalcounts = c(1:length(tmp))
for(i in 1:length(tmp)) {
  cat("totalcounts_", i, ",")
  bamfile = file.path(tmp[i])
  reads.ranges=read.BAM(bamfile)
  totalcounts[i]=length(reads.ranges)  
}

peak <- up_DGEs
#peak <- down_DEGs
peak$V2 <- as.numeric(peak$V2)
peak$V3 <- as.numeric(peak$V3)
ran=IRanges(start=peak$V2, end=peak$V3)
peak.range=GRanges(seqnames=Rle(peak$V1), ranges=ran)

result=matrix(0, nrow=length(peak.range), ncol=length(tmp)+1)
result[,length(tmp)+1] <- peak$V4

for(j in 1:length(tmp)) {
  cat("count", j, ",")
  bamfile = file.path(tmp[j])
  reads.ranges=read.BAM(bamfile)
  result[,j] = countOverlaps(peak.range, reads.ranges)/(totalcounts[j]/10000000)
}

colnames(result) <- c('Pu.1_cWT','Pu.1_cKO','H3K4me3_cWT','H3K4me3_cKO','H3K27ac_cWT','H3K27ac_cKO','gene_id')
result <- as.data.frame(result)
write.csv(result, file="up_tss_2000_count.csv", quote=F,row.names=F)
#write.csv(result, file="down_tss_2000_count.csv", quote=F,row.names=F)
