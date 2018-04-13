library("gplots") 

args <- commandArgs(TRUE)
dataFile<-args[1]
outputCorrelationFile<-args[2]
outputPngFile<-args[3]
outputPdfFile<-args[4]

#pngFile<-"d:/r/heatmaps/test.png"
#dataFile<-"d:/r/heatmaps/chr_normalized2.txt"

data<-read.table(dataFile,sep="\t",header=TRUE,row.names=1,comment.char = "!")
#data <- gsub("[#]", "", data)
head (data)
ncol(data)
data = data[rowSums(data)>=0.5*ncol(data),]

cs=colSums(data)
data = t( t(data)/cs) * 1e6;
data = log2(data+1)

#rownames(data)
correlationData<-cor(data,method='pearson',use='pairwise.complete.obs')
#colnames(correlationData) = c("",col(data)
write.table(correlationData,outputCorrelationFile, sep="\t", quote=FALSE,col.names=NA)
#colnames(correlationData) = cname
correlationData
draw<-function() {
	heatmap.2(correlationData,col=colorRampPalette(c("#0072bc", "#ffffff","#cc2424")),
		       margins=c(15,15),
		       srtCol=NULL,trace="none",density.info=c("none"),key=T,keysize=1.2)
}

png(outputPngFile,width=400,height=400,units="mm",res=300)
draw()
dev.off()

pdf(outputPdfFile)
draw()
dev.off()

