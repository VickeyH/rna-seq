#!/usr/bin/env python
#coding:utf-8

import os,argparse,re
#from pathlib2 import Path
from os.path import join,abspath
from collections import OrderedDict
from compiler.ast import flatten
from collections import defaultdict

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to do DEG analysis")
    parser.add_argument("-fpkm","--fpkm",type=str,help="The read count file",required=True)
    parser.add_argument("-sg","--sample_group",type=str,help="The all samples and groups rep. like groupA=sample1,sample2,... groupB=sample3,sample4,...",required = True,nargs="+")
    parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB ..., ", required = True,nargs="+")
    parser.add_argument("-back","--back",type=str,help="the back anno file",required = True)
    parser.add_argument("-degdir","--degdir",type=str,help="the deg analysis dir",required = True)
    parser.add_argument("-anno","--anno",type=str,help="the track_ref anno file",required = True)
    parser.add_argument("-o","--outputfile",type=str,help="The output file",required = True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def get_enrichment(enrich):
    d = {}
    with open(enrich) as fi:
        header = fi.next()
        for line in fi:
            k = line.split("\t")[0]
            v = line.split("\t")[1:]
            d[k] = v
    return d
  
def get_fpkm(fpkmfile):
    d = OrderedDict()
    with open(fpkmfile) as fi:
        header = fi.next().strip().split("\t")[1:]
        for line in fi:
            k = line.split("\t")[0]
            v = dict(zip(header,line.strip().split("\t")[1:]))
            d[k] = v
    return d
    
def gd(sg):
    d = {}
    for i in sg:
        k = i.split("=")[0]
        v = i.split("=")[1].split(",")
        d[k]=v
    return d
 
def get_anno(ano):
    d = {}
    with open(ano) as fi:
        header = fi.next()
        for line in fi:
            k = line.split("\t")[0]
            v = line.strip().split("\t")[1:]
            d[k] = v
    return d
    
def deginfo(degdir,vslist):
    info = defaultdict(dict)
    for vs in vslist:        
        with open(join(degdir,vs,"deg.up_down.txt")) as updown:
            h = updown.next()
            for line in updown: 
                k = line.split("\t")[0]
                v = line.strip().split("\t")[-1]
                info[k].setdefault(vs,{})["mark"] = v
    for vs in vslist:
        with open(join(degdir,vs,vs+".txt")) as vstxt:
            h = vstxt.next()
            for line in vstxt:
                line = line.strip().split("\t")
                k = line[0]
                info[k].setdefault(vs,{})["p"] = line[-2]
                info[k].setdefault(vs,{})["q"] = line[-1]
                info[k].setdefault(vs,{})["fc"] = line[-3]
                if not info.has_key(k):
                    info[k].setdefault(vs,{})["mark"] = "no"
    return info
   
def main():
    args = parseArg()
    fpkm = get_fpkm(args.fpkm)
    enrich = get_enrichment(args.back)
    h = os.popen("head -1 %s"%args.back).read().split("\t")[1:]
    anno = get_anno(args.anno)
    groupdict = gd(args.sample_group)
    deg3info = deginfo(args.degdir,args.vs)
    samples = list(set(flatten(groupdict.values())))
    degd = abspath(args.degdir)
    degene = os.popen("cut -f1 %s |sort -u|grep -v '^id$'"%" ".join([join(degd,vs,vs+".txt") for vs in args.vs])).readlines()
    fo = open(args.outputfile,"w")
    fo.write("ID\t" + "\t".join(samples) + "\t" + "\t".join([vs+i for vs in args.vs for i in [".log2FoldChange",".pvalue",".padj_value",".Mark"]]) + "\t" + "\t".join(["Position","Ref_gene","Ref_gene_Pos","Entrez","Uniprot","Symbol","Description"]) + "\t" )
    fo.write("\t".join(h))
    for gene in degene:
        gene=gene.strip()
        fo.write(gene+"\t" + "\t".join([fpkm[gene][s] for s in samples]) + "\t")
        for vs in args.vs:
            if vs in deg3info[gene]:
                vs_str = "\t".join([deg3info[gene][vs].get(x,"no") for x in ["fc","p","q","mark"]])
                fo.write(vs_str + "\t")
            else:
                fo.write("\t".join(["-"]*4) + "\t")
        fo.write("\t".join(anno.get(gene,["-"]*7)) + "\t")
        entriz = anno.get(gene,["-"]*7)[3]
        fo.write("\t".join(enrich.get(entriz,["-"]*len(h))) + "\n")
            
            
if __name__ == "__main__":
    main()
        
    
