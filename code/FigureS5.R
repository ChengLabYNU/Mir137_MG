# 单细胞MG的EISA数据处理
library(eisaR)
library(GenomicFeatures)
library(reshape2)
library(tidyverse)

Rex <- read.csv('../data/sc_MG_un_exon_sum.csv',row.names=1,check.names=F)
Rin <- read.csv('../data/sc_MG_un_intron_sum.csv',row.names=1,check.names=F)
Rall <- Rex + Rin
fracIn <- colSums(Rin)/colSums(Rall)
summary(fracIn)

Nex <- t(t(Rex) / colSums(Rex) * mean(colSums(Rex)))
Nin <- t(t(Rin) / colSums(Rin) * mean(colSums(Rin)))

NLex <- log2(Nex + 8)
NLin <- log2(Nin + 8)

Dex <- NLex[,c("cKO-1","cKO-2")] - NLex[,c("cWT-1","cWT-2")]
Din <- NLin[,c("cKO-1","cKO-2")] - NLin[,c("cWT-1","cWT-2")]
Dex.Din <- Dex - Din

quantGenes <- rownames(Rex)[ rowMeans(NLex) > 5.0 & rowMeans(NLin) > 5.0 ]
length(quantGenes)

cor(Dex[quantGenes,1], Dex[quantGenes,2])
cor(Din[quantGenes,1], Din[quantGenes,2])
cor(Dex.Din[quantGenes,1], Dex.Din[quantGenes,2])

####Statistical analysis####
library(edgeR)
cnt <- data.frame(Ex = Rex, In = Rin)
y <- DGEList(counts = cnt, genes = data.frame(ENTREZID = rownames(cnt)))
y <- y[quantGenes, ]
y <- calcNormFactors(y)
region <- factor(c("ex","ex","ex","ex","in","in","in","in"), levels = c("in", "ex"))
cond <- rep(factor(c("cKO","cKO","cWT","cWT")), 2)
design <- model.matrix(~ region * cond)
rownames(design) <- colnames(cnt)
design
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
tt <- topTags(lrt, n = nrow(y), sort.by = "none")
head(tt$table[order(tt$table$FDR, decreasing = FALSE), ])
sig <- tt$table$FDR < 0.05
sum(sig)

sig.dir <- sign(tt$table$logFC[sig])
cols <- ifelse(sig, "#E41A1C", "#22222244")

tmp <- as.data.frame(cbind(rowMeans(Dex),rowMeans(Din),cols))
colnames(tmp) <- c('Dexon','Dintron','cols')
plot_data <- tmp[quantGenes,]
text_tmp <- paste0('R=',round(cor.test(plot_data$Dexon,plot_data$Dintron)$estimate, 4))
cor.test(plot_data$Dexon,plot_data$Dintron, method="pearson")

p1 <- ggplot(plot_data,aes(x=Dexon,y=Dintron))+
  geom_point(color = cols)+
  xlab(expression(paste(Delta,"exon"))) +
  ylab(expression(paste(Delta,"intron")))+
  ylim(floor(min(plot_data$Dintron)),ceiling(max(plot_data$Dintron)))+
  xlim(floor(min(plot_data$Dexon)),ceiling(max(plot_data$Dexon)))+
  annotate("text", label = text_tmp, x = 0.9*min(plot_data$Dexon), y = 0.9*max(plot_data$Dintron), size = 3)+  theme_bw()
ggsave(plot=p1,file='FigS5B.pdf',width=8,height=8)

