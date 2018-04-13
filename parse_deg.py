#!/usr/bin/env python
#coding:utf-8
  
import os,sys,argparse,re
from math import log

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to parse deg file from DEG analysis. The deg file must be 'ID baseMean baseMean_group1   baseMean_group2 log2FoldChange pval padj' colnums")
    parser.add_argument("-i","--input",type=str,help="The input deg result file from DEG analysis",required = True,nargs="+")
    parser.add_argument("-c","--count",type=str,help="The input count matrix file, only for total genes count",required = True)
    parser.add_argument("-p","--pvalue",type=float,help="The pvalue threshold, default:0.05",default=0.05)
    parser.add_argument("-q","--qvalue",type=float,help="The padj value threshold, default:1",default=1)
    parser.add_argument("-fc","--fold_change",type=str,help="The fold change threshold, default:2",default=2)
    #parser.add_argument("-o","--output_dir",type=str,help="The output directory")
    parser.add_argument("--padj",action="store_false",help="Used pval or padj for valcano plot, default: p value. if set, padj value will be used",default=True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
        

def dofile(filename,fc=2,pval=0.05,padj=1,pval_valcano=True):
    outdir = os.path.dirname(filename)
    key = os.path.basename(os.path.dirname(os.path.abspath(filename)))
    #with open(filename) as fi, open(os.path.join(outdir,key+".plot_valcano.file"),"w") as valcano, open(outdir+"/deg.up.txt","w") as up, open(outdir+"/deg.down.txt","w") as down, open(outdir+"/deg.up_down.list","w") as updown, open(os.path.join(outdir,key+".plot_scatter.file"),"w")as scatter:
    with open(filename) as fi, open(os.path.join(outdir,key+".plot_valcano.file"),"w") as valcano, open(outdir+"/deg.up.txt","w") as up, open(outdir+"/deg.down.txt","w") as down, open(outdir+"/deg.up_down.txt","w") as updown:
        if pval_valcano:
            valcano.write("ID\tlog2(fold change)\t-log10(p value)\t%s\n"%key)
        else:
            valcano.write("ID\tlog2(fold change)\t-log10(padj value)\t%s\n"%key)
        header = fi.next()
        up.write(header.strip()+"\tUP_DOWN\n")
        down.write(header.strip()+"\tUP_DOWN\n")
        updown.write(header.strip().split()[0] + "\t" + "Pvalue" + "\tUP_DOWN\n")
        #scatter.write("ID\t%s fpkm value\t%s fpkm valuet\t%s\n"%(key.split("-")[0],key.split("-")[-1],key))
        up_num,down_num = 0,0
        expg1 = set()
        expg2 = set()
        totaldiff = 0
        for line in fi:
            if not line.strip():continue
            totaldiff += 1
            n = line.strip().split("\t")
            if n[-2] == "NA" or n[-1]== "NA":
                #valcano.write(n[0] + "\t" + str(n[-3]) + "\t" + str(-log(float(n[-2]),10)) + "\t\n")
                continue
            n[1:] = map(float,n[1:])
            group1,group2,log2fc,pvalue,qvalue=n[-5:]
            if float(group1) > 0:expg1.add(n[0])
            if float(group2) > 0:expg2.add(n[0])
            if abs(log2fc) >= log(fc,2) and pvalue <= pval and qvalue<= padj:
                if log2fc > 0:
                    up_num += 1
                    up.write("\t".join(map(str,n)) + "\tup\n")
                    updown.write(n[0] + "\t" + str(pvalue) + "\tup\n")
                    #scatter.write(n[0] + "\t" + str(group1) + "\t" + str(group2) + "\tup\n")
                    if pval_valcano:
                        valcano.write(n[0] + "\t" + str(log2fc) + "\t" + str(-log(pvalue,10)) + "\tup\n")
                    else:
                        valcano.write(n[0] + "\t" + str(log2fc) + "\t" + str(-log(qvalue,10)) + "\tup\n")
                else:
                    down_num += 1
                    down.write("\t".join(map(str,n)) + "\tdown\n")
                    updown.write(n[0] + "\t" + str(pvalue) + "\tdown\n")
                    #scatter.write(n[0] + "\t" + str(group1) + "\t" + str(group2) + "\tup\n")
                    if pval_valcano:
                        valcano.write(n[0] + "\t" + str(log2fc) + "\t" + str(-log(pvalue,10)) + "\tdown\n")
                    else:
                        valcano.write(n[0] + "\t" + str(log2fc) + "\t" + str(-log(qvalue,10)) + "\tdown\n")
            else:
                if pval_valcano:
                    valcano.write(n[0] + "\t" + str(log2fc) + "\t" + str(-log(pvalue,10)) + "\t\n")
                else:
                    valcano.write(n[0] + "\t" + str(log2fc) + "\t" + str(-log(qvalue,10)) + "\t\n")
                #scatter.write(n[0] + "\t" + str(group1) + "\t" + str(group2) + "\t\n")
    with open(outdir+"/deg.num.txt","w") as deg_num:
        deg_num.write(key+"\tUP\t%d\n"%up_num)
        deg_num.write(key+"\tDOWN\t%d\n"%down_num)
    with open(os.path.join(outdir,key+".diff.stat"),"w") as ds:
        ds.write("Discription\tNumber\tRatio(%)\n")
        Total_Genes = len(totalgene)
        Expressed_Genes = len(expg1.union(expg2))
        Expressed_In_g1 = len(expg1)
        Expressed_In_g2 = len(expg2)
        ExpressedBoth = len(expg1.intersection(expg2))
        TotalDiffGenes = totaldiff
        updiff = up_num
        downdiff = down_num
        ds.write('Total Genes\t%d\t100.00\n'%Total_Genes)
        ds.write('Expressed Genes\t%d\t%.2f\n'%(Expressed_Genes,float(Expressed_Genes)/Total_Genes*100))
        ds.write('Expressed In %s\t%d\t%.2f\n'%(key.split("-")[0],Expressed_In_g1,float(Expressed_In_g1)/Total_Genes*100))
        ds.write('Expressed In %s\t%d\t%.2f\n'%(key.split("-")[-1],Expressed_In_g2,float(Expressed_In_g2)/Total_Genes*100))
        ds.write('Expressed In Both\t%d\t%.2f\n'%(ExpressedBoth,float(ExpressedBoth)/Total_Genes*100))
        ds.write('Expressed Only In %s\t%d\t%.2f\n'%(key.split("-")[0],Expressed_In_g1-ExpressedBoth,float(Expressed_In_g1-ExpressedBoth)/Total_Genes*100))
        ds.write('Expressed Only In %s\t%d\t%.2f\n'%(key.split("-")[-1],Expressed_In_g2-ExpressedBoth,float(Expressed_In_g2-ExpressedBoth)/Total_Genes*100))
        ds.write('Total Diff Expressed Genes\t%d\t%.2f\n'%(TotalDiffGenes,float(TotalDiffGenes)/Total_Genes*100))
        ds.write('Up Diff Expressed Genes\t%d\t%.2f\n'%(updiff,float(updiff)/Total_Genes*100))
        ds.write('Down Diff Expressed Genes\t%d\t%.2f\n'%(downdiff,float(downdiff)/Total_Genes*100))

def get_total_gene(count):
    geneset = set()
    with open(count) as fi:
        fi.next()
        for line in fi:
            geneset.add(line.split("\t")[0])
    return geneset
                
        
def main():
    #if not os.path.exists(args.output_dir): os.makedirs(args.output_dir)
    if not args.padj and args.qvalue == 1:
        print "Error: used padj value as valcano plot, the padj value threshold must be set, 0.05 or 0.1 as recommendation" 
        sys.exit(1)
    #dofile(args.input,args.output_dir,fc=args.fold_change,pval=args.pvalue,padj=args.qvalue,pval_valcano=args.padj)
    for f in args.input:
        dofile(f,fc=args.fold_change,pval=args.pvalue,padj=args.qvalue,pval_valcano=args.padj)
    degnumfile = " ".join([os.path.dirname(f)+"/deg.num.txt" for f in args.input])
    os.system("cat %s > "%degnumfile + os.path.dirname(os.path.dirname(f)) + "/deg.num.txt")
    expgene = os.popen
        
if __name__ == "__main__":
    args = parseArg()
    totalgene = get_total_gene(args.count)
    main()
    
