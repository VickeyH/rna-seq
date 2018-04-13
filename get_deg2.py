#!/usr/bin/env python
#coding:utf-8
#### cd /lustre/work/yongdeng/Project/160761_test/deseq2 && python ~/get_deg.py -c gene_count_matrix.txt -fc 2 -s 2641,2642,2643,2720IB,2821,3136IB,3328NT,3404IB,3551 -g IR,IR,NT,IB,IR,IB,NT,IB,NT -o deg -vs IB-VS-NT
import argparse,os,sys

parser = argparse.ArgumentParser(description="This Script is used to get deg used DESeq2")
parser.add_argument("-c","--count_tab",type=str,help="The read-count table",required = True)
parser.add_argument("-s","--sample",type=str,help="The all samples name in read-count table,the order must be same as group name",required = True)
parser.add_argument("-g","--group",type=str,help="The all group name in read-count table,the order must be same as sample name",required = True)
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

def write_deg_rscripts(rscripts,countable,samples,groups,vs_list,outdir):
    group_dict = {}
    table_col = os.popen("head -1 %s"%countable).read().strip().split()[1:]
    assert len(samples) == len(table_col), "sample names(%s) are not match the countable first line(%s)"%(samples,table_col)
    d = dict(zip(samples,groups))
    s = str(table_col).strip("[]")
    g = str([d[i] for i in table_col]).strip("[]")
    for n,m in enumerate(groups): 
        group_dict.setdefault(m,[]).append(samples[n])    
    fr = open(rscripts,"w")
    fr.write('library("DESeq2")\nlibrary("edgeR")\n')
    fr.write('mydata <- read.table("%s",sep="\\t",row.names=1,header=T,as.is=T)\n'%countable)
    fr.write('colnames(mydata) <- c(%s)\n'%str(table_col).strip("[]"))
    fr.write('condition <- factor(c(%s))\n'%g)
    fr.write('rnaseqMatrix <- round(as.matrix(mydata))\n')

    for vs in vs_list:
        fr.write("\n# diff gene of group %s\n"%vs)
        case = vs.split("-")[0]
        control = vs.split("-")[-1]
        casestr = str(group_dict[case]).strip("[]")
        constr = str(group_dict[control]).strip("[]")
        if len(group_dict[case]) == len(group_dict[control]) == 1:
            fr.write('library("DESeq")\n')
            fr.write('dds <- newCountDataSet(rnaseqMatrix, condition)\n')
            fr.write('dds <- estimateSizeFactors(dds)\n')
            fr.write('norcounts <- counts(dds, normalized=T)\n')
            fr.write('keep <- rowSums(rnaseqMatrix)>=2\n')      ### count(dds)>0 返回bool值，将数据框中的所有元素全部替换为TRUE或FALSE, 可以代表1和0，然后用rowSums进行求和，有两个以上True的则保留。该条件至少在1个样本中其count值都大于0的基因，保留。
            fr.write('dds <- dds[keep,]\n')
            fr.write('dds <- estimateSizeFactors(dds)\n')
            fr.write('dds <- estimateDispersions(dds, method="blind", sharingMode="fit-only", fitType="local")\n')
            fr.write('res <- nbinomTest(dds, "%s", "%s")\n'%(control,case))
            fr.write('res <- res[order(res$pval),-5]\n')
            fr.write('res[3:4] <- res[4:3]\n')
            fr.write('colnames(res)[3:4] <- c("baseMean_%s","baseMean_%s")\n'%(case,control))
            #fr.write('res <- na.omit(res)\n')
            fr.write('write.table(res,file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,vs,vs+".txt"))

            #fr.write('egdr_d <- DGEList(counts=rnaseqMatrix, group=condition)\n')
            #fr.write('egdr_d <- calcNormFactors(egdr_d)\n')
            #fr.write('sizeFactors(dds) <- egdr_d$samples$norm.factors\n')
            fr.write('write.table(cbind(data.frame(gene_id=rownames(norcounts)),norcounts),file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,os.path.splitext(os.path.basename(countable))[0] + ".normalized" + os.path.splitext(os.path.basename(countable))[-1]))
        else:
            fr.write('coldata <- data.frame(row.names = colnames(rnaseqMatrix), condition)\n')
            fr.write('dds <- DESeqDataSetFromMatrix(rnaseqMatrix,colData = coldata,design=~condition)\n')
            fr.write('dds <- estimateSizeFactors(dds)\nnorcounts <- counts(dds, normalized=T)\n')
            fr.write('keep <- rowSums(counts(dds)>0)>=2\n')      ### count(dds)>0 返回bool值，将数据框中的所有元素全部替换为TRUE或FALSE, 可以代表1和0，然后用rowSums进行求和，有两个以上True的则保留。该条件至少在2个样本中其count值都大于0的基因，保留。
            fr.write('dds <- dds[keep,]\n')
            fr.write('dds2 <- DESeq(dds)\n')

            fr.write('diff_gene_deseq2 <- results(dds2,contrast = c("condition","%s","%s"))\n'%(case,control))
            #fr.write('diff_gene_deseq2 <- na.omit(diff_gene_deseq2)\n')
            fr.write('diff_gene_deseq2 <- diff_gene_deseq2[order(diff_gene_deseq2$pvalue),]\n')
            fr.write('diff_info <- data.frame(id=rownames(diff_gene_deseq2),baseMean=diff_gene_deseq2$baseMean,baseMean_%s=rowMeans(norcounts[rownames(diff_gene_deseq2),c(%s)]),baseMean_%s=rowMeans(norcounts[rownames(diff_gene_deseq2),c(%s)]),log2FoldChange=diff_gene_deseq2$log2FoldChange,pval=diff_gene_deseq2$pvalue,padj=diff_gene_deseq2$padj,row.names=rownames(diff_gene_deseq2))\n'%(case,casestr,control,constr))
            fr.write('write.table(diff_info,file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,vs,vs+".txt"))           

            fr.write('write.table(cbind(data.frame(gene_id=rownames(norcounts)),norcounts),file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,os.path.splitext(os.path.basename(countable))[0] + ".normalized" + os.path.splitext(os.path.basename(countable))[-1]))
    fr.close()

def main():
    countab = os.path.abspath(args.count_tab)
    if not os.path.exists(countab):
        print "%s file not exists!"%countab
        sys.exit(0)
    outdir = os.path.abspath(args.output_dir)
    samplename = [i.strip() for i in args.sample.strip().split(",")]
    groupname = [i.strip() for i in args.group.strip().split(",")]
    vs_list = args.vs
    mkdir(outdir,vs_list)
    write_deg_rscripts(os.path.join(outdir,"diff_gene_DESeq2.r"),countab,samplename,groupname,vs_list,outdir)
    
if __name__ == "__main__":
    main()
    
    
