#!/usr/bin/env python
#coding:utf-8
#### cd /lustre/work/yongdeng/Project/160761_test/deseq2 && python ~/get_deg.py -c gene_count_matrix.txt -fc 2 -s 2641,2642,2643,2720IB,2821,3136IB,3328NT,3404IB,3551 -g IR,IR,NT,IB,IR,IB,NT,IB,NT -o deg -vs IB-VS-NT
import argparse,os,sys

parser = argparse.ArgumentParser(description="This Script is used to get deg used DESeq2")
parser.add_argument("-c","--count_tab",type=str,help="The read-count table",required = True)
parser.add_argument("-fc","--fold_change",type=float,help="The foldchange Threshold for diff_gene, default: 2,[log(fc,2)]",default = 2)
parser.add_argument("-s","--sample",type=str,help="The sample name,the order must be same as group name",required = True)
parser.add_argument("-g","--group",type=str,help="The group name,the order must be same as sample name",required = True)
parser.add_argument("-o","--output_dir",type=str,help="The output directory",required = True)
parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB,CasegroupC-VS-ControlgroupD..., '-' symbol must not in group name", required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def write_deg_rscripts(rscripts,countable,samples,groups,vs_list,fc,outdir):
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
    fr.write('info <- data.frame(row.names=c(%s),name=c(%s),group=c(%s))\n'%(s,s,g))
    fr.write('dds <- DESeqDataSetFromMatrix(as.matrix(mydata),colData = info,design=~group)\n')
    fr.write('keep <- rowSums(cpm(counts(dds))>0)>0\n')      ### count(dds)>2 返回bool值，将数据框中的所有元素全部替换为TRUE或FALSE, 可以代表1和0，然后用rowSums进行求和，有两个以上True的则保留。该条件至少在2个样本中其count值都大于0的基因，保留。
    fr.write('dds <- dds[keep,]\ndds2 <- DESeq(dds)\n')
    fr.write('dds <- estimateSizeFactors(dds)\nnorcounts <- counts(dds, normalized=T)\n')
    fr.write('rld <- rlogTransformation(dds,blind=FALSE)\n')
    fr.write('write.table(cbind(data.frame(gene_id=rownames(rld)),assay(rld)),file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,"gene.rlog2.counts.table"))
    fr.write('write.table(cbind(data.frame(gene_id=rownames(norcounts)),norcounts),file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,"gene.rpm.count.table"))
    for vs in vs_list:
        fr.write("\n# diff gene of group %s\n"%vs)
        case = vs.split("-")[0]
        control = vs.split("-")[-1]
        casestr = str(group_dict[case]).strip("[]")
        constr = str(group_dict[control]).strip("[]")
        fr.write('diff_gene_deseq2 <- results(dds2,contrast = c("group","%s","%s"))\n'%(case,control))
        fr.write('diff_gene_deseq2 <- diff_gene_deseq2[order(diff_gene_deseq2$padj),]\n')
        #fr.write('diff_gene_deseq2 <- na.omit(diff_gene_deseq2)\n')
        #fr.write('diff_gene_deseq2 <- subset(diff_gene_deseq2,padj <= 0.05 & (log2FoldChange >= log(%f,2) | log2FoldChange <= -log(%f,2)))\n'%(fc,fc))
        #fr.write('up_down <- c();for(n in diff_gene_deseq2$log2FoldChange){if (n >= log(%f,2)) up_down <- c(up_down,"up") else if (n <= -log(%f,2)) up_down <- c(up_down,"down") else up_down <- c(up_down,"")}\n' %(fc,fc))
        fr.write('up_down <- c()\nfor (n in seq(1,nrow(diff_gene_deseq2))){\n\tif (is.na(diff_gene_deseq2$padj[n])){\n\t\tup_down <- c(up_down,"no")\n\t\tnext\n\t}\n\tif (diff_gene_deseq2$padj[n] <= 0.05 & diff_gene_deseq2$log2FoldChange[n] >= log(%f,2))\n\t\tup_down <- c(up_down,"up")\n\telse if (diff_gene_deseq2$padj[n] <= 0.05 & diff_gene_deseq2$log2FoldChange[n] <= -log(%f,2))\n\t\tup_down <- c(up_down,"down")\n\telse\n\t\tup_down <- c(up_down,"no")\n}\n'%(fc,fc))
        fr.write('diff_info <- data.frame(gene_id=rownames(diff_gene_deseq2),baseMean=diff_gene_deseq2$baseMean,%s=rowMeans(norcounts[rownames(diff_gene_deseq2),c(%s)]),%s=rowMeans(norcounts[rownames(diff_gene_deseq2),c(%s)]),Log2FC=diff_gene_deseq2$log2FoldChange,P_value=diff_gene_deseq2$pvalue,Q_value=diff_gene_deseq2$padj,up_down=up_down,row.names=rownames(diff_gene_deseq2))\n'%(case,casestr,control,constr))
        fr.write('write.table(diff_info,file="%s",quote=F,sep = "\\t",row.names=F)\n'%os.path.join(outdir,vs,vs+".txt"))           
    fr.close()

def main():
    countab = os.path.abspath(args.count_tab)
    if not os.path.exists(countab):
        print "%s file not exists!"%countab
        sys.exit(0)
    outdir = os.path.abspath(args.output_dir)
    fc = args.fold_change
    samplename = [i.strip() for i in args.sample.strip().split(",")]
    groupname = [i.strip() for i in args.group.strip().split(",")]
    vs_list = [i.strip() for i in args.vs.strip().split(",")]
    mkdir(outdir,vs_list)
    write_deg_rscripts(os.path.join(outdir,"diff_gene_DESeq2.r"),countab,samplename,groupname,vs_list,fc,outdir)
    
if __name__ == "__main__":
    main()

    


    
    
    
