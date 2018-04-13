#!/usr/bin/env python
#coding:utf-8
import os,argparse,re,sys
from os.path import join,abspath

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
            
def get_gtf(samplegtf):
    d = {}
    for i in samplegtf:
        k = i.split("=")[0]
        v = i.split("=")[-1]
        d[k]=v
    return d
            

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to do DEG analysis")
    parser.add_argument("-gtf","--gtf",type=str,help="The gtf file, corresponding to sample name, like sample1=/path/to/gtf1 sample2=/path/to/gtf2",required=True,nargs="+")
    parser.add_argument("-hdrs","--hdrs",type=str,help="The fa.hdrs file",required = True)
    parser.add_argument("-o","--output_dir",type=str,help="The DEG output directory",required = True)
    parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')    
    return parser.parse_args()
    
def write_asp_makeflow(makeflow,samplegtf,outputdir,hdrs):
    fo = open(makeflow,"a+")
    for sn in samplegtf:
        gtf = samplegtf[sn]
        fo.write("CATEGORY=asprofile_%s\n"%sn)
        fo.write(join(outputdir,sn,"ASprofile.as : ") + samplegtf[sn] + " " + hdrs +"\n")
        fo.write("\t/lustre/software/target/asprofile-1.0.4/extract-as %s %s > %s 2> %s\n\n"%(gtf,hdrs,join(outputdir,sn,"ASprofile.as"),join(outputdir,sn,"extract-as.err")))
        fo.write("CATEGORY=asprofile_%s\n"%sn)
        fo.write(join(outputdir,sn,"ASprofile.as.summary ") + join(outputdir,sn,"ASprofile.as.nr : ") + join(outputdir,sn,"ASprofile.as\n"))
        fo.write("\tperl /lustre/software/target/asprofile-1.0.4/summarize_as.pl %s %s -p %s &> %s\n\n"%(gtf,join(outputdir,sn,"ASprofile.as"),join(outputdir,sn,"ASprofile"),join(outputdir,sn,"summarize_as.log")))
        fo.write("CATEGORY=asprofile_%s\n"%sn)
        fo.write(join(outputdir,sn,"ASprofile.as.types_stat ") + join(outputdir,sn,"ASprofile.as.info : ")+ join(outputdir,sn,"ASprofile.as\n"))
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/flow/count_as.pl %s &> %s && python /lustre/work/yongdeng/software/protokaryon/flow/addheader.py %s\n\n"%(join(outputdir,sn,"ASprofile.as"),join(outputdir,sn,"count_as.log"),join(outputdir,sn,"ASprofile.as.info")))
    fo.write("CATEGORY=asprofile_state\n")
    fo.write(join(outputdir,"asprofile.state.txt : ") + " ".join([join(outputdir,sn,"ASprofile.as.types_stat") for sn in samplegtf]) + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/stat_asprofile.py -dir %s -o %s &> %s\n\n"%(outputdir,join(outputdir,"asprofile.state.txt"),join(outputdir,"stat_asprofile.log")))
    fo.close()
    
def main():
    args = parseArg()
    gtf = [i.split("=")[0] + "=" + abspath(i.split("=")[-1]) for i in args.gtf]
    out = abspath(args.output_dir)
    samplegtf = get_gtf(gtf)
    mkdir(out,samplegtf.keys())
    write_asp_makeflow(args.makeflow,samplegtf,out,abspath(args.hdrs))
    
if __name__ == "__main__":
    main()
        
