library('grid')
library('VennDiagram')
args <- commandArgs(TRUE)
numA<- as.numeric(args[1])
numB<- as.numeric(args[2])
numAB<- as.numeric(args[3])
sampleA <- args[4]
sampleB <- args[5]
#outputPngFile<-args[6]
#outputPngFile<-args[6]
#outputPdfFile<-args[7]
outputFilePrefix <- args[6]
outputPdfFile <- args[7]
cols<-c('#5da5da','#decf3f') 
#draw.pairwise.venn(area1=numA,area2=numB,cross.area=numAB,category=c(sampleA,sampleB),lwd=rep(1,1),lty=rep(2,2),col=cols,fill=cols,cat.col=cols)
#png(outputPngFile)
#dev.off()
#png(outputPngFile)
#pdf(outputPdfFile)
#draw.pairwise.venn(area1=numA,area2=numB,cross.area=numAB,category=c(sampleA,sampleB),lwd=rep(1,1),lty=rep(2,2),col=cols,fill=cols,cat.col=cols )
#dev.off()
#dev.off()
draw<-function() {draw.pairwise.venn(area1=numA,area2=numB,cross.area=numAB,
                   category=c(sampleA,sampleB),lwd=rep(1,1),lty=rep(2,2),
                   col=cols,fill=cols,cat.col=cols )
}
png(outputFilePrefix)
draw()
dev.off()

pdf(outputPdfFile)
draw()
dev.off()
