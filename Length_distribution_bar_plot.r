args <- commandArgs(TRUE)
png(args[2])
test<-read.table(args[1],sep="\t",header=TRUE)
barplot(test[,2],xlab=colnames(test)[1],ylab=colnames(test)[2],col="#0072bc",names.arg=test[,1])
dev.off()

