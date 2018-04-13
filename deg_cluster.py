#!/usr/bin/env python
#coding:utf-8
  
import os,argparse
from subprocess import call
from compiler.ast import flatten

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to get deg cluster result, all 0 fpkm row will be removed")
    parser.add_argument("-fpkm","--fpkm",type=str,help="The input fpkm table",required = True)
    parser.add_argument("-sg","--sample_group",type=str,help="The all samples and groups rep. like groupA=sample1,sample2,... groupB=sample3,sample4,...",required = True,nargs="+")
    parser.add_argument("-deg","--deg",type=str,help="The deg info, 'like group1-vs-group2,deg.file'", required=True,nargs="+")
    parser.add_argument("-o","--output",type=str,help="The output dir",required=True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def parse_group(samplegroup):
    d = {}
    for sg in samplegroup:
        k = sg.split("=")[0]
        v = sg.split("=")[1].split(",")
        d[k] = v
    return d
    
def parse_fpkm(fpkm):
    d = {}
    with open(fpkm) as fi:
        header = fi.next().strip().split("\t")[1:]
        for line in fi:
            line = line.strip().split("\t")
            gene = line[0]
            if any(map(float,line[1:])):
                d[gene] = dict(zip(header,line[1:]))
    return d
        
    
def main():
    args = parseArg()
    sgdict = parse_group(args.sample_group)
    samples = sorted(set(flatten(sgdict.values())))
    fpkm = parse_fpkm(args.fpkm)
    out = args.output
    alldeglist = []
    for deg in args.deg:   
        vs = deg.split(",")[0]
        if not os.path.isdir(os.path.join(out,vs)):os.makedirs(os.path.join(out,vs))
        g1 = vs.split("-")[0]
        g2 = vs.split("-")[-1]
        degfile = deg.split(",")[1]
        with open(os.path.join(out,vs,"%s_vs_%s.list"%(g1,g2)),"w") as fo:
            fo.write("".join([g1+"\t"+i+"\n" for i in sgdict[g1]] + [g2+"\t"+i+"\n" for i in sgdict[g2]]))
        with open(degfile) as fdeg,open(os.path.join(out,vs,"DEG_cluster.input.txt"),"w") as fo:
            header = fdeg.next()
            fo.write("ID\t%s"%"\t".join(sgdict[g1] + sgdict[g2]) + "\n")
            for line in fdeg:
                geneid = line.split("\t")[0]
                alldeglist.append(geneid)
                fo.write(geneid+"\t" + "\t".join([fpkm[geneid][s] for s in sgdict[g1] + sgdict[g2]]) + "\n")
    alldeglist = sorted(set(alldeglist),key = alldeglist.index)
    with open(os.path.join(out,"DEG_cluster.input.txt"),"w") as fo:
        fo.write("ID\t"+"\t".join(samples) + "\n")
        for gene in alldeglist:
            fo.write(gene +"\t" + "\t".join([fpkm[gene][s] for s in samples]) + "\n")
        #print('perl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/DifferentialExpression/PtR.pl --CPM --matrix %s --log2 --heatmap --min_colSums 10 --min_rowSums 10 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --output %s'%(os.path.join(out,vs,"DEG_cluster.input.txt"),os.path.join(out,vs,"DEG_cluster.input.ng")))
        #call('perl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/DifferentialExpression/PtR.pl --CPM --matrix %s --log2 --heatmap --min_colSums 10 --min_rowSums 10 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --output %s'%(os.path.join(out,vs,"DEG_cluster.input.txt"),os.path.join(out,vs,"DEG_cluster.input.ng")),shell=True)
        #print('perl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/DifferentialExpression/PtR.pl --CPM --samples %s --matrix %s --log2 --heatmap --min_colSums 10 --min_rowSums 10 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --output %s'%(os.path.join(out,vs,"%s_vs_%s.list"%(g1,g2)),os.path.join(out,vs,"DEG_cluster.input.txt"),os.path.join(out,vs,"DEG_cluster.input")))
        #call('perl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/DifferentialExpression/PtR.pl --CPM --samples %s --matrix %s --log2 --heatmap --min_colSums 10 --min_rowSums 10 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --output %s'%(os.path.join(out,vs,"%s_vs_%s.list"%(g1,g2)),os.path.join(out,vs,"DEG_cluster.input.txt"),os.path.join(out,vs,"DEG_cluster.input")),shell=True)
    
if __name__ == "__main__":
    main()
        
