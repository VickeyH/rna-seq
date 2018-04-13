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
    parser = argparse.ArgumentParser(description="This Script is used to get diff info for propro analysis")
    parser.add_argument("-vsdir","--vsdir",type=str,help="The deg vs dir", required=True,nargs="+")
    return parser.parse_args()
    
def main():
    args = parseArg()
    for vs in args.vsdir:
        v = os.path.basename(os.path.normpath(vs))
        with open(join(vs,"deg.up.anno.txt")) as fi, open(join(vs,v+".diff_info"),"w") as fo:
            header = fi.next().split("\t")
            pi = header.index("pvalue")
            ri = header.index("Ref_gene")
            for line in fi:
                ref = line.split("\t")[ri]
                pvalue = line.split("\t")[pi]
                if ref == "-":
                    fo.write("-\t%s\tup\n"%pvalue)
                else:
                    fo.writelines([i +"\t"+pvalue + "\tup\n" for i in ref.split("|")])
        with open(join(vs,"deg.down.anno.txt")) as fi, open(join(vs,v+".diff_info"),"a+") as fo:
            header = fi.next()
            for line in fi:
                ref = line.split("\t")[ri]
                pvalue = line.split("\t")[pi]
                if ref == "-":
                    fo.write("-\t%s\tdown\n"%pvalue)
                else:
                    fo.writelines([i +"\t"+pvalue + "\tdown\n" for i in ref.split("|")])

if __name__ == "__main__":
    main()
