#####################################################################
#                            FigS6A                                 #
#####################################################################

library(tidyverse)
library(eulerr)

peak_annotation = read.csv('ATAC-seq_peak_annotation.txt',header = T,sep = '\t')

WT_peak = filter(peak_annotation,group %in% 'WT Specific')
KO_peak = filter(peak_annotation,group %in% 'KO Specific')
com_peak = filter(peak_annotation,group %in% 'Common')

WT_peak_site = paste0(WT_peak[,3],WT_peak[,4],WT_peak[,5])
KO_peak_site = paste0(KO_peak[,3],KO_peak[,4],KO_peak[,5])
com_peak_site = paste0(com_peak[,3],com_peak[,4],com_peak[,5])


venn_list=c("WT peaks" = length(WT_peak_site),
            "KO peaks" = length(KO_peak_site),
            "WT peaks&KO peaks" = length(com_peak_site))
p_venn = plot(euler(venn_list),fills = list(fill=c("black","red"),alpha=0.5),
              quantities = c(length(WT_peak_site),length(KO_peak_site),length(com_peak_site)))
ggsave('Figure6A.pdf',p_venn,width = 4,height = 3)

#####################################################################
#                            FigS6B                                 #
#####################################################################
library(ggrepel)
library(ggsci)

peak_annotation$anno = unlist(lapply(strsplit(peak_annotation$Annotation,' \\('),'[',1))
peak_annotation$anno = factor(peak_annotation$anno,
                              levels = c("3' UTR","5' UTR","exon","Intergenic","intron","non-coding","promoter-TSS","TTS"))

p_bar = ggplot(peak_annotation,aes(x=group,fill=anno))+
  geom_bar(position = "fill") +
  geom_text_repel(aes(label = ..count..), stat = "count", 
                  position = 'fill',size=2.8)+
  labs(title = 'ATAC-seq peaks distribution',y='Percent')+
  theme_classic()+
  scale_fill_nejm()+
  scale_y_continuous(expand = c(0,0),limits = c(0,1),labels = scales::percent)+
  theme(plot.title = element_text( hjust = 0.5),
        axis.title = element_text(color = 'black',size=12),
        axis.text = element_text(color = 'black',size=8))

ggsave('FigureS6B.pdf',p_bar,width = 4.2,height = 3)