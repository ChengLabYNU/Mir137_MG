library(Seurat)
library(ggplot2)
library(ggcor)
load('data/cKW_0.6_final.rda')



norm.data <- as.data.frame(scRNA1@assays$RNA@data)
barcode <- colnames(norm.data)
WT1.sc.data <- norm.data[,colnames(norm.data) %in% grep("^cWT_1",barcode,value=T)]
KO1.sc.data <- norm.data[,colnames(norm.data) %in% grep("^cKO_1",barcode,value=T)]
WT1.norm.sc.data <- data.frame(submitted_id=rownames(WT1.sc.data),W1norm_exp=rowSums(WT1.sc.data))
KO1.norm.sc.data <- data.frame(submitted_id=rownames(KO1.sc.data),K1norm_exp=rowSums(KO1.sc.data))

WT2.sc.data <- norm.data[,colnames(norm.data) %in% grep("^cWT_2",barcode,value=T)]
KO2.sc.data <- norm.data[,colnames(norm.data) %in% grep("^cKO_2",barcode,value=T)]
WT2.norm.sc.data <- data.frame(submitted_id=rownames(WT2.sc.data),W2norm_exp=rowSums(WT2.sc.data))
KO2.norm.sc.data <- data.frame(submitted_id=rownames(KO2.sc.data),K2norm_exp=rowSums(KO2.sc.data))


table(scRNA1$orig.ident)

#cKO_1 cKO_2 cWT_1 cWT_2 
# 4197  3034  3101  2304

WT1.norm.sc.data$W1mean_exp <- WT1.norm.sc.data$W1norm_exp/3101
WT2.norm.sc.data$W2mean_exp <- WT2.norm.sc.data$W2norm_exp/2304
KO1.norm.sc.data$K1mean_exp <- KO1.norm.sc.data$K1norm_exp/4197
KO2.norm.sc.data$K2mean_exp <- KO2.norm.sc.data$K2norm_exp/3034

B2.norm.data <- full_join(WT2.norm.sc.data,KO2.norm.sc.data)
B1.norm.data <- full_join(WT1.norm.sc.data,KO1.norm.sc.data)
B2B1.data <- full_join(B2.norm.data,B1.norm.data,"submitted_id")
B2B1.data$scKO_mean <- (B2B1.data$K1norm_exp + B2B1.data$K2norm_exp)/2
B2B1.data$scWT_mean <- (B2B1.data$W1norm_exp + B2B1.data$W2norm_exp)/2
B2B1.data [is.na (B2B1.data)] <- 0


summary(B2B1.data$W1mean_exp)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00051 0.09464 0.06560 5.62308

summary(B2B1.data$W2mean_exp)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00066 0.09551 0.06813 5.97092

summary(B2B1.data$K1mean_exp)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.000000 0.000531 0.092705 0.063712 5.485085

summary(B2B1.data$K2mean_exp)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.000000 0.000662 0.090839 0.059481 5.891204

WT.B2B1 <- B2B1.data %>%  filter(W2mean_exp>0& W2mean_exp<0.07)%>%  filter(W1mean_exp>0& W1mean_exp<0.07)
KO.B2B1 <- B2B1.data %>%  filter(K2mean_exp>0& K2mean_exp<0.07)%>%  filter(K1mean_exp>0& K1mean_exp<0.07)


WT <- select(WT.B2B1,submitted_id,W2mean_exp,W1mean_exp)
KO <- select(KO.B2B1,submitted_id,K2mean_exp,K1mean_exp)
alldata1 <- full_join(WT,KO,"submitted_id") 
alldata1 [is.na (alldata1)] <- 0
row.names(alldata1) <- alldata1$submitted_id
alldata1 <- subset(alldata1,select = -c(submitted_id))

#####################################################################
#                             FigS4B                                #
#####################################################################

pdf("FigS4B.pdf")
quickcor(alldata1, cluster = TRUE,cor.test = TRUE) +
geom_colour() +
geom_mark(size=5,color="black",fontface=5)+
scale_fill_gradientn(colours = c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"))+ geom_panel_grid(colour = "white",size = 1)
dev.off()

#####################################################################
#                             FigS4C                                #
#####################################################################

sc_data <- select(B2B1.data,submitted_id,scKO_mean,scWT_mean)
colnames(sc_data) <- c('geneid','scKO_mean','scWT_mean')
Inde_Bulk <- read.csv('data/bulk_cKW.csv')
Inde_Bulk <- select(Inde_Bulk,geneid,bKO_mean,bWT_mean)
Isb_data <- full_join(Inde_Bulk,sc_data,'geneid')
Isb_data [is.na (Isb_data)] <- 0

library(ggpubr)
psbk <-ggplot(Isb_data, aes(x = log2(scKO_mean),y=log2(bKO_mean))) + 
  geom_point(size=.3,color="#F8766D",alpha=0.2) +
  geom_smooth(method=lm ,size=1.5,color = "black") +
  stat_cor(label.y = 15, size =5) +
  stat_regline_equation(label.y = 14, size =5) +
  coord_fixed() +
  theme_gray() +
  theme_classic()+
  ylim(0,16)+
  theme(text = element_text(size=15))

psbw <-ggplot(Isb_data, aes(x = log2(scWT_mean),y=log2(bWT_mean))) + 
  geom_point(size=.3,color="#1f78b4",alpha=0.2) +
  geom_smooth(method=lm ,size=1.5,color = "black") +
  stat_cor(label.y = 15, size =5) +
  stat_regline_equation(label.y = 14, size =5) +
  coord_fixed() +
  theme_gray() +
  theme_classic()+
  ylim(0,16)+
  theme(text = element_text(size=15))

ggsave(filename = "FigS4C_KO.pdf",plot = psbk,width = 6,height = 4,dpi = 600)
ggsave(filename = "FigS4C_WT.pdf",plot = psbw,width = 6,height = 4,dpi = 600)



#####################################################################
#                             FigS4E                                #
#####################################################################

genelist <- read.csv('data/markers_list2.csv')
pdp <- DotPlot(scRNA1,features = genelist$gene,group.by = 'group')+ theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
ggsave (plot=pdp,file='FigS4E.pdf',width=12,height=6)

#####################################################################
#                             FigS4F                                #
#####################################################################

vc <- VlnPlot(scRNA1,group.by='group',split.by='TYPE',cols=c('black','red'),feature='Cx3cr1',pt.size=0)+theme(axis.title.x = element_blank())
vh <- VlnPlot(scRNA1,group.by='group',split.by='TYPE',cols=c('black','red'),feature='Hexb',pt.size=0)+theme(axis.title.x = element_blank())
vt <- VlnPlot(scRNA1,group.by='group',split.by='TYPE',cols=c('black','red'),feature='Tmem119',pt.size=0)+theme(axis.title.x = element_blank())
va <- vc/vh/vt

ggsave(plot=va,file='FigS4F.pdf')

#####################################################################
#                             FigS4G                                 #
#####################################################################

bardata <- read.csv('data/Forbarplot.csv')
bardata$group <- factor(bardata$group,levels = c('T/NK','DC','MNC','NEUT','PVM','MG-4','MG-3','MG-2','MG-1'))
bardata$TYPE <- factor(bardata$TYPE,levels = c('cWT','cKO'))
pb <- ggplot(bardata,aes(x=TYPE,weight=Proportion,fill=group))+
  geom_bar(position = 'stack')+scale_fill_manual(values = c('#FF9999','#006666','#CC0033','#CC00CC','#CC6699','#6633CC','#993300','#CC9900','#339966'))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave(plot=pb,file='Fig4G.pdf',width=6,height=8)
