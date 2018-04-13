#!/usr/bin/env python
#coding:utf-8

import os,argparse,re
#from pathlib2 import Path
from os.path import join,abspath

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to do sample correlation analysis")
    parser.add_argument("-fpkm","--fpkm",type=str,help="The imput fpkm file",required=True)
    parser.add_argument("-o","--output_dir",type=str,help="The sample correlation analysis output directory",required = True)
    parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def write_sample_cor_makeflow(makeflow,fpkm,outdir):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=filter_fpkm\n")
    fo.write(join(outdir,"gene.corfpkm.txt ") + join(outdir,"gene.anyfpkm.txt ") + join(outdir,"gene.all0.txt : ") + fpkm + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/filter.count.py %s %s &> %s\n\n"%(fpkm,outdir,join(outdir,"filter.count.log")))
    fo.write("CATEGORY=Corheatmap\n")
    fo.write(join(outdir,"genes.fpkm.cor.result ") + join(outdir,"genes.fpkm.cor.png ") + join(outdir,"genes.fpkm.cor.pdf : ") + join(outdir,"gene.corfpkm.txt\n"))
    fo.write("\tRscript /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/Corheatmap.r %s %s %s %s &> %s\n\n"%(join(outdir,"gene.corfpkm.txt"),join(outdir,"genes.fpkm.cor.result"),join(outdir,"genes.fpkm.cor.png"),join(outdir,"genes.fpkm.cor.pdf"),join(outdir,"Corheatmap.log")))
    fo.close()
    
def main():
    args = parseArg()
    fpkm = abspath(args.fpkm)
    out = abspath(args.output_dir)
    if not os.path.isdir(out):
        os.makedirs(out)
    write_sample_cor_makeflow(args.makeflow,fpkm,out)
    
if __name__ == "__main__":
    main()
