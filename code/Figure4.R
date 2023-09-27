library(Seurat)
library(ggplot2)
library(MySeuratWrappers)
library(dplyr)
library(pheatmap)

load('data/cKW_0.6_final.rda')
load(file="data/all.RData")
#####################################################################
#                             Fig4A                                 #
#####################################################################


deg.data <- read.csv(file="data/deseq2.csv", header= T, sep = ",")
deg.data$logP <- -log10(deg.data$padj)

deg.data$Group = "Not-significant"
deg.data$Group[which((deg.data$padj < 0.05) & (deg.data$log2FoldChange > 1))] = "Up-regulated"
deg.data$Group[which((deg.data$padj < 0.05) & (deg.data$log2FoldChange < -1))] = "Down-regulated"
deg.data$Label = ""


deg.data <- deg.data[order(deg.data$r), ]
up.genes <- head(deg.data$X[which(deg.data$Group == "Up-regulated")],9)
down.genes <- head(deg.data$X[which(deg.data$Group == "Down-regulated")],5)
deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
deg.data$Label[match(deg.top10.genes, deg.data$X)] <- deg.top10.genes

set.seed(1234)
p1 = ggscatter(deg.data, x = "log2FoldChange", 
               y = "logP", color = "Group", palette = c("#2f5688", "#BBBBBB","#000000","#CC0000"), 
               size = 1.5,
               xlab = "log2FoldChange",
               ylab = "-log10(padj)", ylim=c(0,350),
               xlim=c(-10,12),) + theme_base() + 	 
  geom_text_repel(aes(label=deg.data$Label), size=4, #family="Arial", fontface="bold", 
                  box.padding=unit(1, "lines"), point.padding=unit(0.5, "lines"), segment.color = "Black", 
                  segment.size = 0.5, force = 1, max.overlaps = Inf,
                  nudge_x = c(-0.5, 0.5),nudge_y = c(-0.5, 0.5),) + 
  theme(text = element_text(size = 20), axis.text = element_text(size = 18))
p1
ggsave(filename = "Fig4A.pdf",plot = p1,width =12,height =10,dpi = 800)


#####################################################################
#                             Fig4B                                 #
#####################################################################

dat = apply(test, 2, function(x) { x[which(is.na(x))] = 0; x } )
dat <- as.matrix(dat)
dat = t(scale(t(dat)))
dat = dat[match(tt$geneid, rownames(dat)),]
dat = dat[apply(dat, 1, function(x) { !all(is.na(x))}),]
colAnnot = data.frame("Group"=Group$Group)
row.names(colAnnot) = colnames(dat)
annotation_row = data.frame( Geneset = tt$geneSet)
row.names(annotation_row) = rownames(dat)
p= pheatmap(dat, annotation_col = colAnnot, annotation_row = annotation_row, cluster_cols=FALSE, border_color = NA, 
            cluster_rows = F, show_colnames = FALSE, fontsize_row = 6,treeheight_row = 0, show_rownames = T)
ggsave(filename = "Fig4B.pdf",plot = p,width = 10,height =12,dpi = 800)

#####################################################################
#                             Fig4C                                 #
#####################################################################

pv <- VlnPlot(scRNA1, features = c('Cx3cr1','Hexb','Tmem119','Olfml3','Fcrls','Nav2','Nav3','Pik3ip1','Ccl3','Ccl4','Ifit2','Ifit3','Pf4','Cxcr2','Plac8','H2-Eb1','Nkg7'),
              cols=c('#339966','#CC9900','#993300','#6633CC','#CC6699','#CC00CC','#CC0033','#006666','#FF9999'),
              group.by= "group", 
              stacked=T,
              pt.size=0)+
  theme(axis.ticks.y= element_blank())+
  theme(axis.text.y= element_blank())
ggsave(plot=pv,file='Fig4C.pdf',width=12,height=11)

detach ("package:MySeuratWrappers")

#####################################################################
#                             Fig4D                                 #
#####################################################################

pu <- DimPlot(scRNA1,pt.size=0.8,group.by='group',cols=c('#339966','#CC9900','#993300','#6633CC','#CC6699','#CC00CC','#CC0033','#006666','#FF9999'))
ggsave(plot=pu,file='Fig4D.pdf',width=12,height=11)

#####################################################################
#                             Fig4E                                 #
#####################################################################

put <- DimPlot(scRNA1,pt.size=0.8,group.by='TYPE',cols=c('black','red'),reduction='umap')
ggsave(plot=put,file='Fig4E.pdf',width=12,height=11)


#####################################################################
#                             Fig4F                                 #
#####################################################################

pm1 <- FeaturePlot(scRNA1,reduction='umap',feature='Tmem119')
pm2 <- FeaturePlot(scRNA1,reduction='umap',feature='Cx3cr1')
pm3 <- FeaturePlot(scRNA1,reduction='umap',feature='Hexb')
pm <- pm1/pm2/pm3

pW1 <- FeaturePlot(scRNA1,reduction='umap',feature='Olfml3',cols=c('lightgrey','black'),max.cutoff=3)
pW2 <- FeaturePlot(scRNA1,reduction='umap',feature='Fcrls',cols=c('lightgrey','black'),max.cutoff=3)
pW3 <- FeaturePlot(scRNA1,reduction='umap',feature='Ccr5',cols=c('lightgrey','black'),max.cutoff=3)
pW <- pW1/pW2/pW3

pK1 <- FeaturePlot(scRNA1,reduction='umap',feature='Pik3ip1',cols=c('lightgrey','red'),max.cutoff=3)
pK2 <- FeaturePlot(scRNA1,reduction='umap',feature='Ccl3',cols=c('lightgrey','red'),max.cutoff=3)
pK3 <- FeaturePlot(scRNA1,reduction='umap',feature='Ccl4',cols=c('lightgrey','red'),max.cutoff=3)
pK <- pK1/pK2/pK3

ggsave(plot=pm,file='Fig4F_MG.pdf',width=3.5,height=9)
ggsave(plot=pW,file='Fig4F_cWT.pdf',width=3.5,height=9)
ggsave(plot=pK,file='Fig4F_cKO.pdf',width=3.5,height=9)


#####################################################################
#                             Fig4G                                 #
#####################################################################

bardata <- read.csv('data/Forbarplot.csv')
bardata$group <- factor(bardata$group,levels = c('T/NK','DC','MNC','NEUT','PVM','MG-4','MG-3','MG-2','MG-1'))
bardata$TYPE <- factor(bardata$TYPE,levels = c('cWT','cKO'))
pb <- ggplot(bardata,aes(x=TYPE,weight=Proportion,fill=group))+
  geom_bar(position = 'stack')+scale_fill_manual(values = c('#FF9999','#006666','#CC0033','#CC00CC','#CC6699','#6633CC','#993300','#CC9900','#339966'))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave(plot=pb,file='Fig4G.pdf',width=6,height=8)