#!/usr/bin/env python
#coding:utf-8
import argparse,os,sys
from tempfile import TemporaryFile
from rpy2.robjects import r
from compiler.ast import flatten
from collections import OrderedDict

parser = argparse.ArgumentParser(description="This Script is used to get deg used DESeq2")
parser.add_argument("-c","--count_tab",type=str,help="The read-count table",required = True)
parser.add_argument("-sg","--sample_group",type=str,help="The all samples and groups rep. like groupA=sample1,sample2,... groupB=sample3,sample4,...",required = True,nargs="+")
parser.add_argument("-o","--output_dir",type=str,help="The output directory",required = True)
parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB ..., ", required = True,nargs="+")
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def write_deg_rscripts(countable,samplesgroups,vs_list,outdir):
    table_col = os.popen("head -1 %s"%countable).read().strip().split()[1:]
    group_dict = OrderedDict()
    for rep in samplesgroups:
        k = rep.split("=")[0]
        v = rep.split("=")[-1].split(",")
        group_dict[k] = v

    samples = set(flatten(group_dict.values()))
    assert len(samples) == len(table_col), "sample names(%s) are not match the countable first line(%s)"%(list(samples),table_col)
    s = str(table_col).strip("[]")
    fr = TemporaryFile()
    fr.write('suppressMessages(library("DESeq2"))\nsuppressMessages(library("edgeR"))\n')
    fr.write('mydata <- read.table("%s",sep="\\t",row.names=1,header=T,as.is=T)\n'%countable)
    fr.write('colnames(mydata) <- c(%s)\n'%str(table_col).strip("[]"))
    fr.write('rnaseqMatrix <- round(as.matrix(mydata))\n')
    fr.write('condition <- factor(c(%s))\n'%s)
    fr.write('coldata <- data.frame(row.names = colnames(rnaseqMatrix), condition)\n')
    fr.write('dds <- DESeqDataSetFromMatrix(rnaseqMatrix,colData = coldata,design=~condition)\n')
    fr.write('dds <- estimateSizeFactors(dds)\nnorcounts <- counts(dds, normalized=T)\n')
    fr.write('write.table(cbind(data.frame(gene_id=rownames(norcounts)),norcounts),file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,os.path.splitext(os.path.basename(countable))[0] + ".normalized" + os.path.splitext(os.path.basename(countable))[-1]))


    for vs in vs_list:
        if "," in vs:
            if vs.endswith("paired"):
                fr.write("\n# diff gene of paired group %s\n"%vs.split(",")[0])
            else:
                print "'paired' string must after the 'vs' info"
                sys.exit(1)
        else:
            fr.write("\n# diff gene of group %s\n"%vs)
        case = vs.split(",")[0].split("-")[0]
        control = vs.split(",")[0].split("-")[-1]
        casestr = str(group_dict[case]).strip("[]")
        constr = str(group_dict[control]).strip("[]")
        g = []
        for i in table_col:
            if i in group_dict[case]:
                g.append(case)
            elif i in group_dict[control]:
                g.append(control)
            else:
                #g.append([ig for ig in group_dict if i in group_dict[ig]][0])
                g.append(i)
        g = str(g).strip("[]")
        if "," in vs:
            vs = vs.split(",")[0]
            if len(group_dict[case]) != len(group_dict[control]):
                print "number of sample must be paired if paired analysis"
                sys.exit(1)
            fr.write('suppressMessages(library("limma"))\n')
            fr.write('limmaMatrix <- mydata[,c(%s)]\n'%str(group_dict[case]+group_dict[control]).strip("[]"))
            fr.write('paired <- factor(rep(seq(1,%d),2))\n'%len(group_dict[case]))
            fr.write('condition <- factor(c(rep(2,%d),rep(1,%d)))\n'%(len(group_dict[case]),len(group_dict[control])))
            fr.write('design <- model.matrix(~paired+condition)\n')
            fr.write('v1 <- voom(limmaMatrix,design,plot=TRUE,normalize="quantile")\n')
            fr.write('fit <- lmFit(v1,design,method="ls")\n')
            fr.write('fit2 <- eBayes(fit)\n')
            fr.write('Output <- topTable(fit2, coef=ncol(design), adjust.method="BH", n=Inf)\n')
            fr.write('res <- data.frame(id=rownames(Output),baseMean=rowMeans(norcounts[rownames(Output),]),baseMean_%s=rowMeans(norcounts[rownames(Output),c(%s)]),baseMean_%s=rowMeans(norcounts[rownames(Output),c(%s)]),log2FoldChange=Output$logFC,pval=Output$P.Value,padj=Output$adj.P.Val,row.names=rownames(Output))\n'%(case,casestr,control,constr))
            fr.write('res <- res[order(res$pval),]\n')
            fr.write('write.table(res,file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,vs,vs+".txt"))
            continue
        fr.write('condition <- factor(c(%s))\n'%g)
        if len(group_dict[case]) == len(group_dict[control]) == 1:
            fr.write('suppressMessages(library("DESeq"))\n')
            fr.write('dds <- newCountDataSet(rnaseqMatrix, condition)\n')
            fr.write('dds <- estimateSizeFactors(dds)\n')
            #fr.write('norcounts <- counts(dds, normalized=T)\n')
            fr.write('keep <- rowSums(rnaseqMatrix)>0\n')      ### count(dds)>0 返回bool值，将数据框中的所有元素全部替换为TRUE或FALSE, 可以代表1和0，然后用rowSums进行求和，有两个以上True的则保留。该条件至少在1个样本中其count值都大于0的基因，保留。
            fr.write('dds <- dds[keep,]\n')
            fr.write('dds <- estimateSizeFactors(dds)\n')
            fr.write('dds <- estimateDispersions(dds, method="blind", sharingMode="fit-only", fitType="local")\n')
            fr.write('res <- nbinomTest(dds, "%s", "%s")\n'%(control,case))
            fr.write('res <- res[order(res$pval),-5]\n')
            fr.write('res[3:4] <- res[4:3]\n')
            fr.write('colnames(res)[3:4] <- c("baseMean_%s","baseMean_%s")\n'%(case,control))
            #fr.write('res <- na.omit(res)\n')
            fr.write('res$baseMean <- rowMeans(norcounts[res$id,])\n')
            fr.write('res <- subset(res,baseMean_%s > 0 | baseMean_%s > 0)\n'%(case,control))
            fr.write('write.table(res,file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,vs,vs+".txt"))

            #fr.write('egdr_d <- DGEList(counts=rnaseqMatrix, group=condition)\n')
            #fr.write('egdr_d <- calcNormFactors(egdr_d)\n')
            #fr.write('sizeFactors(dds) <- egdr_d$samples$norm.factors\n')
            #fr.write('write.table(cbind(data.frame(gene_id=rownames(norcounts)),norcounts),file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,os.path.splitext(os.path.basename(countable))[0] + ".normalized" + os.path.splitext(os.path.basename(countable))[-1]))
        else:
            fr.write('coldata <- data.frame(row.names = colnames(rnaseqMatrix), condition)\n')
            fr.write('dds <- DESeqDataSetFromMatrix(rnaseqMatrix,colData = coldata,design=~condition)\n')
            fr.write('dds <- estimateSizeFactors(dds)\nnorcounts <- counts(dds, normalized=T)\n')
            fr.write('keep <- rowSums(counts(dds,normalized=T)>0)>=2\n')      ### count(dds)>0 返回bool值，将数据框中的所有元素全部替换为TRUE或FALSE, 可以代表1和0，然后用rowSums进行求和，有两个以上True的则保留。该条件至少在2个样本中其count值都大于0的基因，保留。
            fr.write('dds <- dds[keep,]\n')
            fr.write('dds2 <- DESeq(dds)\n')

            fr.write('diff_gene_deseq2 <- results(dds2,contrast = c("condition","%s","%s"))\n'%(case,control))
            #fr.write('diff_gene_deseq2 <- na.omit(diff_gene_deseq2)\n')
            fr.write('diff_gene_deseq2 <- diff_gene_deseq2[order(diff_gene_deseq2$pvalue),]\n')
            fr.write('diff_info <- data.frame(id=rownames(diff_gene_deseq2),baseMean=diff_gene_deseq2$baseMean,baseMean_%s=rowMeans(norcounts[rownames(diff_gene_deseq2),c(%s)]),baseMean_%s=rowMeans(norcounts[rownames(diff_gene_deseq2),c(%s)]),log2FoldChange=diff_gene_deseq2$log2FoldChange,pval=diff_gene_deseq2$pvalue,padj=diff_gene_deseq2$padj,row.names=rownames(diff_gene_deseq2))\n'%(case,casestr,control,constr))
            fr.write('diff_info <- subset(diff_info,baseMean_%s > 0 | baseMean_%s > 0)\n'%(case,control))
            fr.write('write.table(diff_info,file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,vs,vs+".txt"))           

            #fr.write('write.table(cbind(data.frame(gene_id=rownames(norcounts)),norcounts),file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,os.path.splitext(os.path.basename(countable))[0] + ".normalized" + os.path.splitext(os.path.basename(countable))[-1]))
    fr.seek(0)
    a = fr.readlines()
    fr.close();return a

def main():
    countab = os.path.abspath(args.count_tab)
    if not os.path.exists(countab):
        print "%s file not exists!"%countab
        sys.exit(0)
    outdir = os.path.abspath(args.output_dir)
    vs_list = args.vs
    mkdir(outdir,[i.split(",")[0] for i in vs_list])
    rlist = write_deg_rscripts(countab,args.sample_group,vs_list,outdir)
    #fo = open(os.path.join(outdir,"diff_gene.r"),"w");fo.writelines(rlist);fo.close()
    for line in rlist:
        r(line.strip())
    
if __name__ == "__main__":
    main()
    
    
