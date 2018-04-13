#!/usr/bin/env python
#coding:utf-8
## cd /lustre/work/yongdeng/Project && python ~/cat_rawdata.py -i /mnt/icfs/rawdata_sequence/rna_seq_ref/160761-01/ -s shiregu,shireye -mf 160761_test/cat_rawdata.mf -o 160761_test/rawdata/

import argparse,os,re,sys
import RNA
from itertools import combinations

parser = argparse.ArgumentParser(description="This Script is used to combine all raw fq file to one fq file per sample. only used for pair-end sequences")
parser.add_argument("-i","--input_dir",type=str,help="The input rawdata directory, please make sure there is no comma symbol in fastq filename under the directory",required = True,nargs = "*",metavar="dir")
parser.add_argument("-o","--output_dir",type=str,help="The output rawdata directory",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-s","--sample",type=str,help="all final samples name which will be uesd for downstream analysis",required = True,nargs='*')
parser.add_argument("-a","--addition",type=str,help="the supplyment sample name file, first colnum is the sample name that will mathch the rawdata filename, second colnum is the final sample name which will be uesd for downstream analysis")
parser.add_argument("-len","--length",type=int,help="The avg length of input rawdata seq file, default: 150",default=150)
parser.add_argument("-G","--g_base",type=float,help="The total base(G) count you need, default: 0 means do not split rawdata",default=0)
parser.add_argument("-p","--parallel",type=int,help="the cpu number used, default:8",required = False, default = 8)
#parser.add_argument("-g","--data_G",type=int,help="the all sequences base-pair data(GB) include pair1 and pair2, 1Gb=1024Mb=1024*1024kb, only used for reads number of raw-date larger than need. default: 10",default=10,metavar="data_G")
#parser.add_argument("-l","--len",type=int,help="the read length. each reads should have be the same lenght",default=150)
#parser.add_argument("-e","--extract",action="store_true",help="weather to extract sequence from raw reads, default:no",default=False)
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def listdir(pathlist):
    rel_dir = []
    for path in pathlist:
        for a,b,c in os.walk(path):
            for i in c:
                if i.endswith("fastq.gz") or i.endswith("fq.gz"):
                    rel_dir.append(os.path.join(a,i))
    return rel_dir

def map_sample_and_raw(sample_list,pathlist,sup):
    name_dict = {}
    if sup:
        with open(sup) as n:
          for line in n:
            if not line.strip():
                continue
            k = line.strip().split()[1]
            v = line.strip().split()[0]
            if k != v: name_dict[k] = v
    new_sample_list = [name_dict.get(i,i) for i in sample_list]
    drpsnset =  set()
    for x in [list(i) for i in list(combinations(new_sample_list,2)) if i[0] in i[1] or i[1] in i[0]]:
       drpsnset.add(x[0])
       drpsnset.add(x[1])
    drpsnset = sorted(list(drpsnset),key=len,reverse=True)
    raw_d = listdir(pathlist)
    raw_n = map(lambda x:re.sub("\W+","_",x),map(os.path.basename,raw_d))
    tmp_d  = dict(zip(raw_n,raw_d))
    sample_raw={}
    drp_raw_n = set()
    for sn in new_sample_list:
      if sn not in drpsnset:
        for p in raw_n:
            if re.search('_?%s_'%sn,p):
                sample_raw.setdefault(sn,[]).append(tmp_d[p])
        if not sample_raw.has_key(sn):
            print "no %s in rawdata directory,please check" %sn
      else:
        for p in raw_n:
            if re.search('_?%s_'%sn,p):
                drp_raw_n.add(p)
    if len(drp_raw_n):
        has = []       
        drp_raw_n = list(drp_raw_n)
        for sn in drpsnset:
            for p in drp_raw_n:
                if re.search('_?%s_'%sn,p) and p not in has:
                    sample_raw.setdefault(sn,[]).append(tmp_d[p])
                    has.append(p)
            
    for k in sample_raw:
        kv =  []
        for n in sample_raw[k]:
            bn = re.sub("\W",lambda m:"\\"+m.group(0),os.path.basename(n))
            bn = bn.replace("\\.",".")
            kv.append(os.path.join(os.path.dirname(n),bn))
        kv.sort()
        sample_raw[k] = kv
    new_sample_raw = {}
    for i,s in enumerate(sample_list):
        new_sample_raw[s] = sample_raw[new_sample_list[i]]
    return new_sample_raw

def write_cat_makeflow(sample_dict,out_dir,makeflow,cpu,avglen,gdata):
    cutadapt = RNA.get_bin_abspath("cutadapt")
    if os.path.exists(makeflow):
        print "%s file exists, please use a new name or remove this file"%makeflow
        sys.exit(1)
    fo =open(makeflow,"a+")
    for s in sample_dict:
        if gdata == 0:
            for r in "12":
                fo.write("CATEGORY=%s_cat_raw\n"%s)
                fo.write(os.path.abspath(os.path.join(out_dir,s)) + '/R%s.fastq.gz : %s\n' %(r," ".join(sample_dict[s][int(r)-1::2])))
                fo.write("\tcat %s >> %s\n\n" %(" ".join(sample_dict[s][int(r)-1::2]),os.path.abspath(os.path.join(out_dir,s)) + '/R%s.fastq.gz'%r))
            continue
        for r in "12":
            fo.write("CATEGORY=%s_cat_raw\n"%s)
            fo.write(os.path.abspath(os.path.join(out_dir,s)) + '/R%s.fastq.raw.gz : %s\n' %(r," ".join(sample_dict[s][int(r)-1::2])))
            fo.write("\tcat %s >> %s\n\n" %(" ".join(sample_dict[s][int(r)-1::2]),os.path.abspath(os.path.join(out_dir,s)) + '/R%s.fastq.raw.gz'%r))
        fo.write("CATEGORY=%s_split_raw\n"%s)
        fo.write(os.path.abspath(os.path.join(out_dir,s,"R1.fastq.gz ")) + os.path.abspath(os.path.join(out_dir,s,"R2.fastq.gz : "))+ os.path.abspath(os.path.join(out_dir,s,"R1.fastq.raw.gz ")) + os.path.abspath(os.path.join(out_dir,s,"R2.fastq.raw.gz\n")))
        fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/split_raw.py -G %f -len %d -p %d -raw %s &> %s\n\n"%(gdata,avglen,cpu,os.path.abspath(os.path.join(out_dir,s,"R1.fastq.raw.gz ")) + os.path.abspath(os.path.join(out_dir,s,"R2.fastq.raw.gz")),os.path.abspath(os.path.join(out_dir,s,"split.raw.log"))))

    for s in sample_dict:
        fo.write("CATEGORY=%s_rm_adapter\n"%s)
        fo.write(os.path.join(out_dir,s,"R1.fq.gz ") + os.path.join(out_dir,s,"R2.fq.gz : ") + os.path.join(out_dir,s,"R1.fastq.gz ") + os.path.join(out_dir,s,"R2.fastq.gz\n"))
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        fo.write("\tpython %s -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -m 50 -f fastq --info-file %s -o %s -p %s %s %s &> %s && rm -fr %s %s\n\n" %(cutadapt,os.path.join(out_dir,s,"rmAdapter.info"),os.path.join(out_dir,s,"R1.fq.gz"),os.path.join(out_dir,s,"R2.fq.gz"),os.path.join(out_dir,s,"R1.fastq.gz"),os.path.join(out_dir,s,"R2.fastq.gz"),os.path.join(out_dir,s,"rmAdapter.log"),os.path.join(out_dir,s,"R1.fastq.gz"),os.path.join(out_dir,s,"R2.fastq.gz")))        
    fo.close()

def main():
    assert len(args.sample) == len(set(args.sample)), "duplication sample name"
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.makeflow))): os.makedirs(os.path.dirname(os.path.abspath(args.makeflow)))
    pathlist = map(os.path.abspath,args.input_dir)
    sample_raw = map_sample_and_raw(args.sample,pathlist,args.addition)
    mkdir(args.output_dir,sample_raw.keys())
    outdir = os.path.abspath(args.output_dir)
    write_cat_makeflow(sample_raw,outdir,args.makeflow,args.parallel,args.length,args.g_base)  
 
if __name__=="__main__":
    main()
    
