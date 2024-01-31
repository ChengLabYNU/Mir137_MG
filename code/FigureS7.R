#####################################################################
#                             FigS7A                                #
#####################################################################
library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)

gene = read.csv('../../H3K27ac_cWT_specific_peak_and_down_gene.txt',header = F)$V1

# GO enrich
go_enrich = enrichGO(gene,OrgDb=org.Mm.eg.db,keyType ='SYMBOL',
                     ont = 'ALL', pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",qvalueCutoff = 0.2)

go_data = head(arrange(go_enrich@result,pvalue),30)
go_data$log10pvalue = -log10(go_data$pvalue)
go_data$log10pvalue = factor(go_data$log10pvalue,levels = go_data$log10pvalue)

p_go = ggplot(data = go_data,
              aes(x= log10pvalue,y=reorder(Description,log10pvalue),fill=ONTOLOGY))+
  geom_col()+
  geom_text(aes(x=-log10(pvalue)+max(-log10(pvalue))*0.02,
                y=Description,
                label=Count),size=2.8)+
  theme_test()+
  facet_grid(ONTOLOGY~.,scales = "free", space = "free")+
  theme(axis.title = element_text(color='black',size=10),
        axis.text = element_text(color = 'black',size=10),
        strip.background = element_rect(fill='#ffffff'),
        strip.text = element_text(size=10),
        legend.position='none')+
  labs(x=expression(-log[10](pvalue)),y='Description')+
  scale_fill_manual(values = c(BP='#bc3c29',
                               CC='#0072b5',
                               MF='#e18727'))+
  scale_x_continuous(expand = c(0,0),limits = c(0, 12))
ggsave('FigureS7A-1.pdf',
       p_go,height = 6,width = 7)

en_id = bitr(gene,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db)$ENTREZID
kegg_enrich = enrichKEGG(en_id,organism ='mmu',
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.2)

p_kegg = ggplot(data = head(arrange(kegg_enrich@result,pvalue),19),
                aes(x=-log10(pvalue),y=reorder(sub(" - Mus musculus \\(house mouse\\)",'',Description),pvalue)))+
  geom_col(fill='#20854e')+
  geom_text(aes(x=-log10(pvalue)+max(-log10(pvalue))*0.02,
                label=Count),size=2.8)+
  theme_test()+
  labs(x=expression(-log[10](pvalue)),y='Description')+
  theme(axis.title = element_text(color='black',size=10),
        axis.text = element_text(color = 'black',size=10))+
  scale_x_continuous(expand = c(0,0),limits = c(0, 3.5))
ggsave('FigureS7A-2.pdf',
       p_kegg,height = 4,width = 7)


#####################################################################
#                             FigS7B                                #
#####################################################################
gene = c('Tgfb1','Asap1','Rock2','Eng')
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

#获取totalcounts用于后续标准化
totalcounts = c(1:length(tmp))
for(i in 1:length(tmp)) {
  cat("totalcounts_", i, ",")
  bamfile = file.path(tmp[i])
  reads.ranges=read.BAM(bamfile)
  totalcounts[i]=length(reads.ranges)  
}

#读取bam计算count并绘图

peak$V2 <- as.numeric(peak$V2)
peak$V3 <- as.numeric(peak$V3)
ran=IRanges(start=peak$V2, end=peak$V3)
peak.range=GRanges(seqnames=Rle(peak$V1), ranges=ran)

#设置result矩阵
result=matrix(0, nrow=length(peak.range), ncol=length(tmp)+1)
result[,length(tmp)+1] <- peak$V4

#读取bam文件并计数
for(j in 1:length(tmp)) {
  cat("count", j, ",")
  bamfile = file.path(tmp[j])
  reads.ranges=read.BAM(bamfile)
  result[,j] = countOverlaps(peak.range, reads.ranges)/(totalcounts[j]/10000000)
}

#设置行名，转换data.frame，并存储为csv文件
result <- as.data.frame(result)
colnames(result) <- c('cWT-1','cWT-2','cWT-3','cKO-1','cKO-2','cKO-4','gene_id')

for (i in gene){
  tmp_data = filter(gene,gene_id == i)
  tmp_data = melt(tmp_data,id.vars = 'gene_id')
  tmp_data$genotype = unlist(lapply(strsplit(as.character(tmp_data$variable),'-'),'[',1))
  top_limit = max(tmp_data$value)
  
  assign(paste0('p_',i),ggplot(tmp_data,aes(x=genotype,y=value,color=genotype))+
           geom_bar(stat="summary",fun=mean,width = 0.5,fill='#FFFFFF')+
           stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",
                        width=0.2)+
           geom_jitter(size=2, width=0)+
           scale_color_manual(values = c(cWT ='black',cKO='red'))+
           stat_compare_means(comparisons = list(c('cWT','cKO')),method = "t.test",
                              label = 'p.signif', lable.y = top_limit*1.1)+
           scale_y_continuous(expand = c(0,0),limits = c(0,top_limit*1.2))+
           theme_classic()+
           ggtitle(paste0(i,' TSS 1k H3K27ac count'))+
           theme(axis.text = element_text(color='black',size=8),
                 plot.title = element_text(hjust = 0.5)))
}

ggsave('FigureS7B-1',p_Tgfb1, width = 3,height = 2.7)
ggsave('FigureS7B-2',p_Asap1, width = 3,height = 2.7)
ggsave('FigureS7B-3',p_Rock2, width = 3,height = 2.7)
ggsave('FigureS7B-4',p_Eng, width = 3,height = 2.7)