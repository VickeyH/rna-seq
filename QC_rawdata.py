#!/usr/bin/env python
#coding:utf-8
### 质控步骤
## cd /lustre/work/yongdeng/Project && python ~/QC_rawdata.py -i 160761_test/rawdata/ -o 160761_test/QC -mf 160761_test/QC_rawdata.mf
import argparse,os,sys
import RNA

parser = argparse.ArgumentParser(description="This Script is used to QC rawdata of all samples.")
parser.add_argument("-i","--input_dir",type=str,help="The input rawdata directory",required = True)
parser.add_argument("-m","--methods",type=str,help="The quality control methods, default: IlluQC",required = False, default = "IlluQC",choices=["IlluQC","trimmomatic"])
parser.add_argument("-p","--parallel",type=int,help="the cpu number used, default:8",required = False, default = 8)
parser.add_argument("-o","--output_dir",type=str,help="The output cleandata directory",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def write_trimmomatic_qc_makeflow(samples,makeflow,input_d,output_d,thread):
    cleansoft = RNA.get_soft_abspath("trimmomatic")
    fo =open(makeflow,"a+")
    for sn in samples:
        r1 = os.path.join(input_d,sn,"R1.fq.gz")
        r2 = os.path.join(input_d,sn,"R2.fq.gz")
        od = os.path.join(output_d,sn)
        logfile = os.path.join(od,"%s.trim.log"%sn)
        fo.write("CATEGORY=QC_%s\n"%sn)
        fo.write(os.path.join(od,"R1.clean.fq.gz ") + os.path.join(od,"R1.clean.unpaired.fq.gz ") + os.path.join(od,"R2.clean.fq.gz ") + os.path.join(od,"R2.clean.unpaired.fq.gz ") + logfile + " : " + r1 + " " + r2 + "\n")
        fo.write("\tjava -jar %s PE -threads %d -phred33 -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:/lustre/work/yongdeng/software/protokaryon/soft/QC/Trimmomatic/adapters/adapter:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50 " %(cleansoft,thread,logfile,r1,r2,os.path.join(od,"R1.clean.fq.gz"),os.path.join(od,"R1.clean.unpaired.fq.gz"),os.path.join(od,"R2.clean.fq.gz"),os.path.join(od,"R2.clean.unpaired.fq.gz")) + "> " + os.path.join(od,'trim.out') + " 2> " + os.path.join(od,'trim.err') + "\n\n")       
    fo.close()
    
def write_IlluQC_qc_makeflow(samples,makeflow,input_d,output_d,thread):
    cleansoft = RNA.get_bin_abspath("IlluQC_PRLL_CapitalBio_Adapter.pl")
    fo =open(makeflow,"a+")
    for sn in samples:
        r1 = os.path.join(input_d,sn,"R1.fq.gz")
        r2 = os.path.join(input_d,sn,"R2.fq.gz")
        od = os.path.join(output_d,sn)
        fo.write("CATEGORY=QC_%s\n"%sn)
        fo.write(" ".join([os.path.join(od,i) for i in ["R1.clean.fq.gz","R2.clean.fq.gz","summary.png","Filter_stat","R1_avgQual.png","R2_avgQual.png"]]) + " : " + r1 + " " + r2 + "\n")
        fo.write("\tperl %s -pe %s %s /lustre/work/yongdeng/software/protokaryon/soft/QC/adapter 5 -c %d -z g -o %s > %s 2> %s\n\n"%(cleansoft,r1,r2,thread,od,os.path.join(od,"IlluQC_N.out"),os.path.join(od,"IlluQC_N.err")))
        #fo.write("\tperl %s -pe %s %s /lustre/work/yongdeng/software/protokaryon/soft/QC/adapter 5 -c %d -z g -o %s > %s 2> %s"%(cleansoft,r1,r2,thread,od,os.path.join(od,"IlluQC_N.out"),os.path.join(od,"IlluQC_N.err")))
        #fo.write(' && python /lustre/work/yongdeng/software/protokaryon/flow/checkpaired.py %s %s >> %s\n\n'%(os.path.join(od,"R1.clean.fq.gz"),os.path.join(od,"R2.clean.fq.gz"),os.path.join(od,"IlluQC_N.err")))
    fo.close()
        
def main():
    input_dir = os.path.abspath(args.input_dir)
    thread = args.parallel
    output_dir = os.path.abspath(args.output_dir)
    samples = sorted(map(os.path.basename,filter(os.path.isdir,[os.path.join(input_dir,x) for x in os.listdir(input_dir)])))
    mkdir(output_dir,samples)
    no_mean = write_trimmomatic_qc_makeflow(samples,args.makeflow,input_dir,output_dir,thread) if args.methods == "trimmomatic" else write_IlluQC_qc_makeflow(samples,args.makeflow,input_dir,output_dir,thread)
 
if __name__=="__main__":
    main()
