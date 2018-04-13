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
    parser.add_argument("-abbr","--abbr",type=str,help="species abbr")
    parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
    parser.add_argument("-o","--outputdir",type=str,help="The project out dir, analysis directory will be created under this dir",required = True)
    parser.add_argument("--nokobas",help="for species which do not in kobas database, used blast for enrichment. if set, species_type, fasta, gtf must be given",action='store_true',default=False)
    parser.add_argument("-st","--species_type",type=str,help="if '--nokobas', species_type must be set, chose from [kobas_animal,kobas_plants,kobas_bacteria,kobas_fungi,kobas_protists] and [uniprot_fungi,uniprot_plants,uniprot_animal,uniprot_microorganism], separated by comma symbol and no space")
    parser.add_argument("-gtf","--gtf",type=str,help="if '--nokobas', total gtf file must be given.")
    parser.add_argument("-fa","--fasta",type=str,help="if '--nokobas', the ref fasta file must be given.")
    return parser.parse_args()
    
def write_enrichment(makeflow,abbr,idtype,dbtype,outdir,vslist,degdir):
    fo = open(makeflow,"a+")
    analysis = os.path.dirname(degdir)
    fo.write("CATEGORY=enrichment\n")
    fo.write(join(outdir,"enrichment.all.log : ") + join(analysis,"stringtie","back.list ") + " ".join([join(degdir,vs,"deg.up_down.list") for vs in vslist]) + "\n")
    cmdgokegg = "perl /lustre/work/ranlv/program/kobas3.0.3/enrichment_kobas3.0.pl -specise %s -database_type %s -num 20  -deghead -degId_col 1 -upDown_col 2 -pvalue_col 3 -backlist %s -backId_col 1 -idtype %s -od %s -DEG_overall -dir_filename %s,deg.up_down.list"%(abbr,dbtype,join(analysis,"stringtie","back.list"),idtype,outdir,degdir)
    if gethostname().startswith("hpc"):
        cmdgokegg += " -task_assigned_node hpc"
    fo.write("\t" + cmdgokegg + " &> " + join(outdir,"enrichment.all.log") +"\n\n")
    fo.close()

def write_nokobas_enrichment(makeflow,gtf,fa,outdir,vslist,degdir,st):
    fo=open(makeflow,"a+")
    analysis = os.path.dirname(degdir)
    fo.write("CATEGORY=enrichment_blast_backfa\n")
    fo.write(join(outdir,"back.fa : ") + gtf + " " + fa + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/python_scripts/extract_longest_trans_seq.py %s %s %s &> %s\n\n"%(fa,gtf,join(outdir,"back.fa"),join(outdir,"extract_backfa.log")))
    fo.write("CATEGORY=enrichment_blast\n")
    fo.write(join(outdir,"enrichment.all.log : ") + " ".join([join(degdir,vs,"deg.up_down.list") for vs in vslist]) + "\n")
    cmdgokegg = 'perl /lustre/work/ranlv/program/kobas3.0.3/fa_enrichment_KG_v3.0.pl -specise_type %s -backfa %s -backfa_type nul -num 20 -deghead -degId_col 1 -upDown_col 2 -pvalue_col 3 -od %s -DEG_overall -dir_filename %s,deg.up_down.list'%(st,join(outdir,"back.fa"),outdir,degdir)
    if gethostname().startswith("hpc"):cmdgokegg += " -task_assigned_node hpc"
    fo.write("\t" + cmdgokegg + " &> " + join(outdir,"enrichment_blast.all.log") +"\n\n")
    fo.close()
    
def main():
    args=parseArg()
    mkdir(abspath(args.outputdir))
    if args.nokobas:
        if args.database_type.rstrip("/") != "/K/G":
            print "Error:if nokabas specise, only /K/G database_type allowed"
            sys.exit(1)
        if not args.fasta or not args.gtf or not args.species_type:
            print "Error:if nokabas specise, fa and gtf and species_type must be given"
            sys.exit(1)
        write_nokobas_enrichment(args.makeflow,abspath(args.gtf),abspath(args.fasta),abspath(args.outputdir),args.vs,abspath(args.degdir),args.species_type)
    else:
        write_enrichment(args.makeflow,args.abbr,args.idtype,args.database_type,args.outputdir,args.vs,args.degdir)
    
if __name__ == "__main__":
    main()
    
