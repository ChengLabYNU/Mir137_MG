library(Seurat)
library(ggplot2)
load('data/cKW_0.6_final.rda')


#####################################################################
#                             Fig5B                                #
#####################################################################

pfj <- FeaturePlot(scRNA1,reduction='umap',cols=c('lightgrey','darkred'),max.cutoff=2,pt.size=0.8,features = 'Jdp2')
pfm <- FeaturePlot(scRNA1,reduction='umap',cols=c('lightgrey','darkred'),max.cutoff=2,pt.size=0.8,features = 'Mtus2')

ggsave(plot=pfj,file='Fig5B_Jdp2.pdf',width=12,height=11,dpi=540)
ggsave(plot=pfm,file='Fig5B_Mtus2.pdf',width=12,height=11,dpi=540)

#####################################################################
#                             Fig5G                                #
#####################################################################

library(ggsci)
library(scales)
colors <- pal_futurama("planetexpress", alpha = 0.7)(5)

sylamer_plot <- function(data){
  # Any particular words that you want highlighted
  chosenOligos <- c("gcaata","caataa","agcaat","gcaataa","agcaata","agcaataa","aagcaata")
  chosenOligos <- toupper(chosenOligos)
  chosenOligos <- gsub("U","T",chosenOligos)
  chosenOligos <- gsub("[^AGCT]+","",chosenOligos)
  chosenOligos <- unique(chosenOligos)
  sylOutput <- data
  if (!file.exists(sylOutput)) {
    stop(paste("Can't find the expected output file: ",sylOutput,sep=""), call.=FALSE)
  }
  
  sylTable <- read.table(sylOutput, sep="\t", row.names=1, header=T, check.names=F)
  sylTable <- cbind("0"=0,sylTable)   # To add an initial column of 0s
  xVals <- as.numeric(colnames(sylTable))
  
  yMin   <- min(sylTable)
  yMax   <- max(sylTable)
  yRange <- yMax - yMin
  minYlim <- round(yMin-yRange/10)
  maxYlim <- round(yMax+yRange/10)
  if (minYlim > -8) {minYlim = -8}
  if (maxYlim < 8) {maxYlim = 8}
  
  plot(NULL, xlab="Sorted sequences", ylab="log10(enrichment P-value)", axes=T, 
       main=paste("Sylamer landscape using words of length: ",nchar(row.names(sylTable)[1]),sep=""),
       ylim=c(minYlim,maxYlim), xlim=range(xVals))
  
  maxPlot <- 1000
  
  if (maxPlot > 0) {
    sylTable <- sylTable[order(apply(sylTable,1,function(x) max(abs(x))),decreasing=TRUE),]
    if (maxPlot > nrow(sylTable)) {
      maxPlot <- nrow(sylTable)
    }
  } else {
    maxPlot <- nrow(sylTable)
  }
  
  for (i in 1:maxPlot) {
    lines(xVals,sylTable[i,], col='grey')
  }
  
  abline(h=0)
  abline(h=log10(nrow(sylTable))+2,lty=2)
  abline(h=-log10(nrow(sylTable))-2,lty=2)
  abline(v=0,lty=2)
  abline(v=10000,lty=2)
  
  oligosChosen <- c()
  oligosChosen <- chosenOligos[chosenOligos %in% rownames(sylTable)]
  
  names(colors) <- oligosChosen
  for (i in rev(seq_along(oligosChosen))) {
    lines(xVals, sylTable[oligosChosen[i],], col=colors[oligosChosen[i]], lwd=2)
    index <- which.max(abs(sylTable[oligosChosen[i],]))
    text(names(sylTable[oligosChosen[i],][index]),
         sylTable[oligosChosen[i],][index],
         labels = rownames(sylTable[oligosChosen[i],]),
         col = colors[oligosChosen[i]])
  }
}

imageOut  <- "Fig5G1.pdf"
pdf(file=imageOut, width=10, height=11)
par(mfcol=c(3,3))
for (i in c('6','7','8')){
  plot_data <- paste0("../data/bulk_top1w_",i,"_utr3.txt")
  print(plot_data)
  sylamer_plot(plot_data)
}
dev.off()

imageOut  <- "Fig5G2.pdf"
pdf(file=imageOut, width=10, height=11)
par(mfcol=c(3,3))
for (i in c('6','7','8')){
  plot_data <- paste0("../data/sc_top1w_",i,"_utr3.txt")
  print(plot_data)
  sylamer_plot(plot_data)
}
dev.off()