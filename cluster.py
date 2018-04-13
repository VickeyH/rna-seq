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
    parser = argparse.ArgumentParser(description="This Script is used to do DEG cluster analysis")
    parser.add_argument("-fpkm","--fpkm",type=str,help="The imput fpkm file",required=True)
    parser.add_argument("-sg","--sample_group",type=str,help="The all samples and groups rep. like groupA=sample1,sample2,... groupB=sample3,sample4,...",required = True,nargs="+")
    parser.add_argument("-o","--output_dir",type=str,help="The sample correlation analysis output directory",required = True)
    parser.add_argument("-deg","--deg",type=str,help="The deg info, 'like group1-vs-group2,deg.up_down.txt'", required=True,nargs="+")
    parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def write_makeflow(makeflow,fpkm,sg,outdir,deg):
    fo = open(makeflow,"a+")
    vslist = [gene.split(",")[0] for gene in deg]
    fo.write("CATEGORY=deg_cluster\n")
    fo.write(join(outdir,"DEG_cluster.input.txt ") + " ".join([join(outdir,vs,vs.split("-")[0] + "_vs_" + vs.split("-")[-1] + ".list") for vs in vslist])  + " " + " ".join([join(outdir,vs,"DEG_cluster.input.txt") for vs in vslist])+ " : " + fpkm + " " + " ".join([gene.split(",")[-1] for gene in deg])+ "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/deg_cluster.py -fpkm %s -sg %s -deg %s -o %s &> %s\n\n"%(fpkm," ".join(sg)," ".join(deg),outdir,join(outdir,"deg_cluster.log")))
    fo.write("CATEGORY=deg_cluster_ng\n")
    fo.write(" ".join([join(outdir,"DEG_cluster.input.ng." + i) for i in ["CPM.log2.centered.dat","CPM.log2.centered.genes_vs_samples_heatmap.pdf","CPM.log2.centered.genes_vs_samples_heatmap.txt","CPM.log2.centered.heatmap.result","CPM.log2.dat","CPM.log2.result","CPM.log2.sample_cor_matrix.pdf","R"]]) + " : " + join(outdir,"DEG_cluster.input.txt\n"))
    fo.write("\tperl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/DifferentialExpression/PtR.pl --CPM --matrix %s --log2 --heatmap --min_colSums 10 --min_rowSums 10 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --output %s &> %s\n\n"%(join(outdir,"DEG_cluster.input.txt"),join(outdir,"DEG_cluster.input.ng"),join(outdir,"DEG_cluster.ng.log")))
    for gene in deg:
        vs = gene.split(",")[0]
        g1 = vs.split("-")[0]
        g2 = vs.split("-")[-1]
        degfile = gene.split(",")[1]        
        fo.write("CATEGORY=deg_cluster_ng_%s\n"%vs)
        fo.write(" ".join([join(outdir,vs,"DEG_cluster.input.ng." + i) for i in ["CPM.log2.centered.dat","CPM.log2.centered.genes_vs_samples_heatmap.pdf","CPM.log2.centered.genes_vs_samples_heatmap.txt","CPM.log2.centered.heatmap.result","CPM.log2.dat","CPM.log2.result","CPM.log2.sample_cor_matrix.pdf","R"]]) + " : " + join(outdir,vs,"DEG_cluster.input.txt\n"))
        fo.write("\tperl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/DifferentialExpression/PtR.pl --CPM --matrix %s --log2 --heatmap --min_colSums 10 --min_rowSums 10 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --output %s &> %s\n\n"%(join(outdir,vs,"DEG_cluster.input.txt"),join(outdir,vs,"DEG_cluster.input.ng"),join(outdir,vs,"DEG_cluster.ng.log")))
        fo.write("CATEGORY=deg_cluster_%s\n"%vs)
        fo.write(" ".join([join(outdir,vs,"DEG_cluster.input." + i) for i in ["CPM.log2.centered.dat","CPM.log2.centered.genes_vs_samples_heatmap.pdf","CPM.log2.centered.genes_vs_samples_heatmap.txt","CPM.log2.centered.heatmap.result","CPM.log2.dat","CPM.log2.result","CPM.log2.sample_cor_matrix.pdf","R"]]) + " : " + join(outdir,vs,"%s_vs_%s.list "%(g1,g2))   + join(outdir,vs,"DEG_cluster.input.txt\n"))
        fo.write("\tperl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/DifferentialExpression/PtR.pl --CPM --samples %s --matrix %s --log2 --heatmap --min_colSums 10 --min_rowSums 10 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --output %s &> %s\n\n"%(join(outdir,vs,"%s_vs_%s.list"%(g1,g2)),join(outdir,vs,"DEG_cluster.input.txt"),join(outdir,vs,"DEG_cluster.input"),join(outdir,vs,"DEG_cluster.log")))
    fo.close()
    
def main():
    args = parseArg()
    fpkm = abspath(args.fpkm)
    out = abspath(args.output_dir)
    if not os.path.isdir(out):
        os.makedirs(out)
    deg = [i.split(",")[0] + "," + os.path.abspath(i.split(",")[-1]) for i in args.deg]
    write_makeflow(args.makeflow,fpkm,args.sample_group,out,deg)
    
if __name__ == "__main__":
    main()
