#!/usr/bin/env Rscript

# Rscript run_DEG.R -e exprSet.txt -i group.txt -c 'groupA-groupB' -s counts -m DESeq2

####################### group.file ######
# group   samples
# group1  sample1
# group1  sample2
# group2  sample3
# group2  sample4
######################
## please make sure that the header must be include and the first column is the group information 
## Also please make sure that the order of exprSet and samples are identical.
##############################################################################


library("optparse")
 
option_list = list(
	make_option(c("-e", "--expression"), type="character", default=NULL, 
              help="expression matrix file name", metavar="character"),
	make_option(c("-i", "--group"), type="character", default=NULL, 
              help="group information file name.", metavar="character"),
	make_option(c("-c", "--contrast"), type="character", default=NULL, 
              help="How to compare them,such as case-control", metavar="character"),
	make_option(c("-s", "--style"), type="character", default='signals', 
              help="What kind of data,read counts or signals?[default = signals]", metavar="character"),
	make_option(c("-m", "--method"), type="character", default="limma", 
              help="which method to use for the DEG[default = limma].For signals, you can choose limma or t.test.For read counts , you can choose DESeq2 or edgeR", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$expression) | is.null(opt$group) ){
  print_help(opt_parser)
  stop("expression matrix file and group information file must be supplied", call.=FALSE)
}
if (is.null(opt$contrast)  ){
  print_help(opt_parser)
  stop("The comparison information must be supplied", call.=FALSE)
}

# program...

exprSet=read.table(opt$expression,stringsAsFactors = F,header = T,row.names=1)
group_info=read.table(opt$group,stringsAsFactors = F,header = T,sep="\t")

rownames(group_info) <- group_info[,2]
group_list=group_info[,1]

run_limma <- function(exprSet,group_list,contrast){
	suppressPackageStartupMessages(library(limma))
	design=model.matrix(~0+factor(group_list))
	colnames(design)=levels(factor(group_list))
	fit=lmFit(exprSet,design)
	cont.matrix=makeContrasts(contrasts=contrast ,levels = design)
	fit2=contrasts.fit(fit,cont.matrix)
	fit2=eBayes(fit2)
	results=topTable(fit2,adjust='BH',n=Inf)
	write.table(results,paste0("limma_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)

}

run_t.test <- function(exprSet,group_list,contrast){

	keep = group_list   %in%  strsplit(contrast,'-')[[1]]
	dat=exprSet[,keep]
	group_list=group_list[keep]

	group1 = which(group_list == strsplit(contrast,'-')[[1]][1] )
	group2 = which(group_list == strsplit(contrast,'-')[[1]][2] )

	dat1 = dat[, group1]
	dat2 = dat[, group2]
	dat = cbind(dat1, dat2)   
	library(pi0)
	pvals = matrix.t.test(dat, 1, 
		   length(group1), length(group2))
	p.adj = p.adjust(pvals, method = "BH")
	avg_1 = rowMeans(dat1)
	avg_2 = rowMeans(dat2)
	FC = avg_2/avg_1
	results = cbind(avg_1, avg_2, FC, pvals, p.adj)
	colnames(results) = c("avg_1", "avg_2", "logFC", "P.Value", "adj.P.Val") 
	DEG_t_test=results[order(results[,4]),]
	write.table(DEG_t_test,paste0("t_test_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)

}

run_DESeq2 <- function(exprSet,group_list,contrast){
	suppressPackageStartupMessages(library(DESeq2))
    
    group1 = which(group_list == strsplit(contrast,'-')[[1]][1] )
    group2 = which(group_list == strsplit(contrast,'-')[[1]][2] )

	geneLists=row.names(exprSet)
	keepGene=rowSums(edgeR::cpm(exprSet)>0) > 0
	table(keepGene);dim(exprSet)
	dim(exprSet[keepGene,])
	exprSet=exprSet[keepGene,]
	rownames(exprSet)=geneLists[keepGene]

	(colData <- data.frame(row.names=colnames(exprSet), group_list=group_list) )
	dds <- DESeqDataSetFromMatrix(countData = exprSet,
								  colData = colData,
								  design = ~ group_list)

    dds2 <- DESeq(dds)
    dds <- estimateSizeFactors(dds)
    norcounts <- counts(dds, normalized=T)
    
	rld <- rlogTransformation(dds)
	exprMatrix_rlog=assay(rld) 

	exprMatrix_rpm=as.data.frame(norcounts) 
	head(exprMatrix_rpm) 
    
    diff_gene_deseq2 <- results(dds2,contrast = c("group_list",strsplit(contrast,'-')[[1]] ))
    diff_gene_deseq2 <- diff_gene_deseq2[order(diff_gene_deseq2$padj),]
    
    up_down <- c()
    for (n in seq(1,nrow(diff_gene_deseq2))){
        if (is.na(diff_gene_deseq2$padj[n])){
            up_down <- c(up_down,"no")
            next
        }
        if (diff_gene_deseq2$padj[n] <= 0.05 & diff_gene_deseq2$log2FoldChange[n] >= log(2.000000,2))
            up_down <- c(up_down,"up")
        else if (diff_gene_deseq2$padj[n] <= 0.05 & diff_gene_deseq2$log2FoldChange[n] <= -log(2.000000,2))
            up_down <- c(up_down,"down")
        else
            up_down <- c(up_down,"no")
    }

    DESeq2_DEG <- data.frame(gene_id=rownames(diff_gene_deseq2),baseMean=diff_gene_deseq2$baseMean,group1=rowMeans(norcounts[rownames(diff_gene_deseq2),group1]),group2=rowMeans(norcounts[rownames(diff_gene_deseq2),group2]),Log2FC=diff_gene_deseq2$log2FoldChange,P_value=diff_gene_deseq2$pvalue,Q_value=diff_gene_deseq2$padj,up_down=up_down,row.names=rownames(diff_gene_deseq2))
    colnames(DESeq2_DEG) <- c("gene_id","baseMean",strsplit(contrast,'-')[[1]],"Log2FC","P_value","Q_value","up_down")
    
	write.table(cbind(data.frame(gene_id=rownames(norcounts)),norcounts),file="DESeq2.exprMatrix_rpm.txt",quote=F,sep = "\t",row.names=F)
    write.table(cbind(data.frame(gene_id=rownames(exprMatrix_rlog)),exprMatrix_rlog),file='DESeq2.exprMatrix_rlog2.txt',quote=F,sep = "\t",row.names=F)
    write.table(DESeq2_DEG,file=paste0("DESeq2_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = F)
    
}


run_edgeR <- function(exprSet,group_list,contrast){
		suppressPackageStartupMessages(library(edgeR)) 
        
        group1 = which(group_list == strsplit(contrast,'-')[[1]][1] )
        group2 = which(group_list == strsplit(contrast,'-')[[1]][2] )

		keep = group_list   %in%  strsplit(contrast,'-')[[1]]
		exprSet=exprSet[,keep]
		group_list=group_list[keep]

		d <- DGEList(counts=exprSet,group=factor(group_list))
        k_feature <- rowSums(cpm(d)>0) > 0
        d <- d[k_feature,]
		#d$samples$lib.size <- colSums(d$counts)
		d <- calcNormFactors(d)
		d$samples 
		dge=d

		design <- model.matrix(~factor(group_list))
		rownames(design)<-colnames(dge)
		colnames(design)<-levels(factor(group_list))
		dge <- estimateDisp(dge,design)

		# #To perform quasi-likelihood F-tests:
		# fit <- glmQLFit(dge,design)
		# qlf <- glmQLFTest(fit,coef=2)
        # lrt$FDR <- p.adjust(qlf$table$PValue, method='BH')
        
		# To perform likelihood ratio tests:
		fit <- glmFit(dge,design)
		lrt <- glmLRT(fit,coef=2)
        lrt$updown <- decideTestsDGE(lrt)
        lrt$table$updown[lrt$updown==1] <- "down"
        lrt$table$updown[lrt$updown==-1] <- "up"
        lrt$table$updown[lrt$updown==0] <- "no"
        lrt$table$group1mean <- rowMeans(cpm(d)[,group1])
        lrt$table$group2mean <- rowMeans(cpm(d)[,group2])
        lrt$table$baseMean <- rowMeans(cpm(d))

		nrDEG=topTags(lrt, n=nrow(d))
		nrDEG=as.data.frame(nrDEG)
        nrDEG=nrDEG[c("baseMean","group1mean","group2mean","logFC","logCPM","LR","PValue","FDR","updown")]
        colnames(nrDEG) <- c("baseMean",strsplit(contrast,'-')[[1]],"logFC","logCPM","LR","PValue","FDR","updown")
		edgeR_lrt_DEG=nrDEG

		# nrDEG=topTags(qlf, n=nrow(exprSet))
		# nrDEG=as.data.frame(nrDEG)
		# head(nrDEG)
		# edgeR_qlf_DEG=nrDEG

		write.table(edgeR_lrt_DEG,file=paste0("edgeR_lrt_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)
		#write.table(edgeR_qlf_DEG,file=paste0("edgeR_qlf_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)

		
}




if (opt$style == "counts") {

	if (opt$method == "DESeq2") {
	    
			suppressPackageStartupMessages(library(DESeq2))
			run_DESeq2(exprSet,group_list,opt$contrast)
	    }

	if (opt$method == "edgeR") {
	
		
		suppressPackageStartupMessages(library(edgeR)) 
		run_edgeR(exprSet,group_list,opt$contrast)

	}

}

if (opt$style == "signals") {

	if (opt$method == "t.test") {
		suppressPackageStartupMessages(library(pi0)) 
		run_t.test(exprSet,group_list,opt$contrast)
		}

	if (opt$method == "limma") {
		suppressPackageStartupMessages(library(limma)) 
		run_limma(exprSet,group_list,opt$contrast)
		
		
	}

}
