#!/usr/bin/env python
#coding:utf-8
import os,argparse,re,sys,datetime
#from pathlib2 import Path
from os.path import join,abspath
from glob import glob1
from math import ceil
from commands import getstatusoutput
from subprocess import Popen

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue
     
def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to split rawdata")
    parser.add_argument("-raw","--rawdata",type=str,help="The input raw fastq seq files",required=True,nargs="+")
    parser.add_argument("-len","--length",type=int,help="The avg length of input rawdata seq file, default: 150",default=150)
    parser.add_argument("-G","--g_base",type=float,help="The total base(G) count you need, default: 6",default=6.0)
    parser.add_argument("-p","--parallel",type=int,help="the cpu number used",required = False, default = 8)
    #parser.add_argument("-o","--outfile",type=str,help="The out file",required=True,nargs="+")
    return parser.parse_args()
    
def main():
    args=parseArg()
    if args.parallel <= 1:
        print "no more than 2 cpu allowed"
        sys.exit(1)
    if len(args.rawdata) == 2:
        o1 = args.rawdata[0] .replace(".raw.gz","")
        o2 = args.rawdata[1] .replace(".raw.gz","")
        persize = ceil(args.g_base*1024*1024*1024/2.0/args.length) + 2
        (status,output) = getstatusoutput("python /lustre/work/yongdeng/software/python_scripts/extract_random_seq.py -fq %s -size %d -o %s"%(" ".join(args.rawdata),persize," ".join([o1,o2])))
        if status != 0:
            print output
            sys.exit(1)
        else:
            if os.path.isfile(o1) and os.path.isfile(o2):
                p1 = Popen("/mnt/icfs/work/yongdeng/rnacocktail-0.2.2/r.txt.6138.qsub/soft_rnacocktail/bin/pigz-2.3.1/pigz -p %d %s && rm -fr %s"%(args.parallel/2,o1,args.rawdata[0]),shell=True)
                p2 = Popen("/mnt/icfs/work/yongdeng/rnacocktail-0.2.2/r.txt.6138.qsub/soft_rnacocktail/bin/pigz-2.3.1/pigz -p %d %s && rm -fr %s"%(args.parallel/2,o2,args.rawdata[1]),shell=True)
                p1.wait();p2.wait()
            else:
                os.system('mv %s %s'%(args.rawdata[0],args.rawdata[0].replace(".raw","")))
                os.system('mv %s %s'%(args.rawdata[1],args.rawdata[1].replace(".raw","")))
    else:
        print "please use paired seq file"
        sys.exit(1)
        
if __name__ == "__main__":
    main()
             
