#!/usr/bin/env python
#coding:utf-8

import sys,os,argparse,re,datetime
from multiprocessing import Pool,cpu_count
from commands import getstatusoutput

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description= "This script is used to summarize the sequence length infomation or your file fastly!")
parser.add_argument("-fa","--fasta",type=str,help="the fasta input files, '.gz' or '.bz2' file can be allowed",nargs="+")
parser.add_argument("-fq","--fastq",type=str,help="the fastq input files, '.gz' or '.bz2' file can be allowed",nargs="+")
#parser.add_argument("-p","--parallel",type=int,help="the cpu number used, default: 8",required = False, default = 8)
args = parser.parse_args()

def sums(f):
    start = datetime.datetime.now()
    fs = os.path.abspath(f)
    cmd = "/lustre/shared_softwares/vsearch/bin/vsearch --fastx_subsample " if os.path.exists("/lustre/shared_softwares/vsearch/bin/vsearch") else "/lustre/work/yongdeng/software/vsearch/vsearch/bin/vsearch --fastx_subsample "
    cmd += fs + " "
    cmd += "--sample_size %d "%sys.maxint
    if f.endswith(".gz"):cmd += "--gzip_decompress "
    if f.endswith(".bz2"):cmd += "--bzip2_decompress "
    cmd += "--fastaout /dev/null"
    s,tf1 = getstatusoutput(cmd)
    if re.search('\n.+?nt in.+?\n',tf1):
        a = re.findall("\d+",re.search('\n(.+?nt in.+?)\n',tf1).group(0))
        print "File path: " + os.path.abspath(f) + "\n" + "sequence number: "+a[1] + "\n" + "base number: " + a[0] + " bp" + "\n" + 'min reads length: ' + a[2] + " bp\n" + "max reads length: " + a[3] + " bp\n" + "avg reads length: " + a[4] + " bp\n" + "开始时间: %s,\t耗时: %f Seconds\n"%(start.strftime("%x %X"),(datetime.datetime.now()-start).total_seconds())
    else:
        print "Can't get sequence info for %s file!\n"%os.path.abspath(f)
        return

def main():
    if (not args.fasta) and (not args.fastq):
        print "Error: -fa or -fq must be defined!"
        sys.exit(1)
    inputfiles = args.fasta if args.fasta else args.fastq
    a = all([os.path.exists(f) for f in inputfiles])
    if not a:
        print "Error:some input file not exists"
        sys.exit(1)
    print
    pool = Pool(min(len(inputfiles),cpu_count()))
    rl =pool.map(sums,inputfiles)
    pool.close();pool.join()

if __name__ == "__main__":
    main()
            
    
