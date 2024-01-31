#####################################################################
#                             Fig7A                                 #
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

tmp = c("cWT-1.bam","cWT-2.bam","cWT-3.bam",
        "cKO-1.bam","cKO-2.bam","cKO-3.bam")
# Jdp_site <- c("chr12",85597916,85600916,"Jdp2")

totalcounts = c(1:length(tmp))
for(i in 1:length(tmp)) {
  cat("totalcounts_", i, ",")
  bamfile = file.path(tmp[i])
  reads.ranges=read.BAM(bamfile)
  totalcounts[i]=length(reads.ranges)  
}

ran=IRanges(start=as.numeric(85597916), end=as.numeric(85600916))
peak.range=GRanges(seqnames=Rle('chr12'), ranges=ran)
result=matrix(0, nrow=length(peak.range), ncol=length(tmp)+1)
result[,length(tmp)+1] <- 'Jdp2'

for(j in 1:length(tmp)) {
  cat("count", j, ",")
  bamfile = file.path(tmp[j])
  reads.ranges=read.BAM(bamfile)
  result[,j] = countOverlaps(peak.range, reads.ranges)/(totalcounts[j]/10000000)
}

result <- as.data.frame(result)
colnames(result) <- c("cWT-1","cWT-2","cWT-3","cKO-1","cKO-2","cKO-3")
# plot by GraphPad Prism

#####################################################################
#                             Fig7I                                 #
#####################################################################
tmp = c("Pu.1_cWT_combined.bam","Pu.1_cKO_combined.bam",
        "H3K4me3_cWT_combined.bam","H3K4me3_cKO_combined.bam",
        "H3K27ac_cWT_combined.bam","H3K27ac_cKO_combined.bam")
down_DGEs <- read.table(file = "../data/down_tss_2000.bed",sep="\t")

totalcounts = c(1:length(tmp))
for(i in 1:length(tmp)) {
  cat("totalcounts_", i, ",")
  bamfile = file.path(tmp[i])
  reads.ranges=read.BAM(bamfile)
  totalcounts[i]=length(reads.ranges)  
}

peak <- down_DGEs
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

result <- as.data.frame(result)
colnames(result) <- c('Pu.1_cKO','Pu.1_cWT','H3K27ac_cKO','H3K27ac_cWT','H3K4me3_cKO','H3K4me3_cWT','gene_id')
write.csv(result, file="down_gene_tss_count.csv", quote=F,row.names=F)

down_result = melt(result,measure.vars = c('Pu.1_cKO','Pu.1_cWT','H3K27ac_cKO','H3K27ac_cWT','H3K4me3_cKO','H3K4me3_cWT'))
down_result$grodown = unlist(lapply(strsplit(as.character(down_result$variable),'_'),'[',1))
down_result$genotype = unlist(lapply(strsplit(as.character(down_result$variable),'_'),'[',2))
down_result$genotype = factor(down_result$genotype,levels = c('cWT','cKO'))
down_result$value = as.numeric(down_result$value)

pu1_down_plot = ggplot(filter(down_result,grodown == 'Pu.1'),aes(x=genotype,y=value,color=genotype))+
  geom_jitter(size=1,width =0)+
  geom_boxplot(fill="#FFFFFF")+
  scale_color_manual(values = c(cWT ='black',cKO='red'))+
  stat_compare_means(comparisons = list(c('cWT','cKO')),method = "t.test",
                     label = 'p.signif')+
  scale_y_continuous(expand = c(0,0),limits = c(0,120),breaks = c(0,30,60,90,120))+
  theme_classic()+
  labs(title = 'Pu.1',y="Normalized counts",x="")+
  theme(axis.text.y = element_text(color='black',size=8),
        axis.text.x = element_text(color='black',size=12),
        axis.title.y = element_text(color='black',size=16),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')
ggsave('Figure7I-1.pdf',pu1_down_plot,width = 3,height = 3.8)

H3K4me3_down_plot = ggplot(filter(down_result,grodown == 'H3K4me3'),aes(x=genotype,y=value,color=genotype))+
  geom_jitter(size=1,width =0)+
  geom_boxplot(fill="#FFFFFF")+
  scale_color_manual(values = c(cWT ='black',cKO='red'))+
  stat_compare_means(comparisons = list(c('cWT','cKO')),method = "t.test",
                     label = 'p.signif')+
  scale_y_continuous(expand = c(0,0),limits = c(0,300),breaks = c(0,100,200,300))+
  theme_classic()+
  labs(title = 'H3K4me3',y="",x="")+
  theme(axis.text.y = element_text(color='black',size=8),
        axis.text.x = element_text(color='black',size=12),
        axis.title.y = element_text(color='black',size=16),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')
ggsave('Figure7I-2.pdf',H3K4me3_down_plot,width = 2.8,height = 3.8)


H3K27ac_down_plot = ggplot(filter(down_result,grodown == 'H3K27ac'),aes(x=genotype,y=value,color=genotype))+
  geom_jitter(size=1,width =0)+
  geom_boxplot(fill="#FFFFFF")+
  scale_color_manual(values = c(cWT ='black',cKO='red'))+
  stat_compare_means(comparisons = list(c('cWT','cKO')),method = "t.test",
                     label = 'p.signif')+
  scale_y_continuous(expand = c(0,0),limits = c(0,500),breaks = c(0,100,200,300,400,500))+
  theme_classic()+
  labs(title = 'H3K27ac',y="",x="")+
  theme(axis.text.y = element_text(color='black',size=8),
        axis.text.x = element_text(color='black',size=12),
        axis.title.y = element_text(color='black',size=16),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')
ggsave('Figure7I-3.pdf',H3K27ac_down_plot,width = 2.8,height = 3.8)
