#!/usr/bin/env python
#coding:utf-8
import os,argparse,re,sys,datetime
from os.path import join,abspath
from socket import gethostname

def mkdir(dir,samples=None):
    if samples:
        for n in samples:
            d = os.path.join(dir,n)
            if not os.path.isdir(d):
                os.makedirs(d)
            else:
                continue
    else:
        if not os.path.isdir(dir):
            os.makedirs(dir)
    
     
def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to do gokegg analysis")
    parser.add_argument("-deg","--degdir",type=str,help="The deg analysis dir",required = True)
    parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB ..., ", required = True,nargs="+")
    parser.add_argument("-dt","--database_type",type=str,default="/K/G",help="enrichment database_type, default: /K/G")
    parser.add_argument("-it","--idtype",type=str,default="id:ncbigene",choices =["id:ncbigi", "id:uniprot", "id:ensembl","id:ncbigene"],help='gene type of back list and deg list, can be chosen from ["id:ncbigi", "id:uniprot", "id:ensembl","id:ncbigene"],default: id:ncbigene',metavar = "idtype")
    parser.add_argument("-abbr","--abbr",type=str,help="species abbr",required = True)
    parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
    parser.add_argument("-o","--outputdir",type=str,help="The project out dir, analysis directory will be created under this dir",required = True)
    return parser.parse_args()
    
def write_enrichment(makeflow,abbr,idtype,dbtype,outdir,vslist,degdir):
    fo = open(makeflow,"a+")
    analysis = os.path.dirname(degdir)
    fo.write("CATEGORY=enrichment\n")
    fo.write(join(outdir,"enrichment.all.log : ") + join(analysis,"stringtie","back.list ") + " ".join([join(degdir,vs,"deg.up_down.list") for vs in vslist]) + "\n")
    cmdgokegg = "perl /lustre/work/ranlv/program/kobas3.0.3/enrichment_kobas3.0.pl -specise %s -database_type %s -num 20  -deghead -degId_col 1 -upDown_col 2 -pvalue_col 2 -backlist %s -backId_col 1 -idtype %s -od %s -DEG_overall -dir_filename %s,deg.up_down.list"%(abbr,dbtype,join(analysis,"stringtie","back.list"),idtype,outdir,degdir)
    if gethostname().startswith("hpc"):
        cmdgokegg += " -task_assigned_node hpc"
    fo.write("\t" + cmdgokegg + " &> " + join(outdir,"enrichment.all.log") +"\n\n")
    fo.close()
    
def main():
    args=parseArg()
    mkdir(abspath(args.outputdir))
    write_enrichment(args.makeflow,args.abbr,args.idtype,args.database_type,args.outputdir,args.vs,args.degdir)
    
if __name__ == "__main__":
    main()
    
