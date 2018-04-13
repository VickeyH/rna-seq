library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
args <- commandArgs(TRUE)

infile <- args[1]
outPng <- args[2]
outPdf <- args[3]
data <- read.table(infile,sep="\t",header=FALSE)
type=data$V2
type=factor(type, levels=c('UP','DOWN'))
pic <- ggplot(data,aes(x=data$V1,y=data$V3,fill=type)) + geom_bar(stat="identity", position=position_dodge(width = 0.9)) + scale_fill_manual(values=c("#cc2424", "#00a651")) + xlab("VS") + ylab("Numbers") + ggtitle("DEG stat") + geom_text(aes(label = data$V3, vjust = 0.5, hjust = 0, angle=90), show.legend = FALSE , position=position_dodge(width = 0.9)) + ylim(min(data$V3, 0)*1.1, max(data$V3)*1.1) + theme_bw()+theme(panel.background=element_rect(fill='transparent', color='black'),panel.border=element_rect(fill='transparent', color='transparent'),panel.grid=element_blank(),axis.text.x=element_text(angle=30,vjust=0.5))

ggsave(pic,file=outPng)
ggsave(pic,file=outPdf)
