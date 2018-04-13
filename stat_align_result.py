#!/usr/bin/env python
#coding:utf-8
### python stat_align_result.py -a align_dir -s sample1,sample2,sample3.... -o stat_file
import argparse,os,re,sys,commands,glob

parser = argparse.ArgumentParser(description="This Script is used to stat align result of all sample's align file")
parser.add_argument("-a","--align_dir",type=str,help="The align output directory",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-o","--output",type=str,help="The output state file",required = True)
parser.add_argument("-s","--sample",type=str,help="all samples name, the name must be separate by comma without any space. eg: sample1,sample2,sample3....",required = True)

args = parser.parse_args()

def stat_align(align_file):
    total,unique,mult=0,0,0
    with open(align_file) as ai:
        for line in ai:
            if "reads; of these" in line:
                total += int(line.strip().split()[0])*2
            elif "aligned concordantly exactly 1 time" in line or "aligned discordantly 1 time" in line:
                unique += int(line.strip().split()[0])*2
            elif "aligned exactly 1 time" in line:
                unique += int(line.strip().split()[0])
            elif "aligned concordantly >1 times" in line:
                mult += int(line.strip().split()[0])*2
            elif "aligned >1 times" in line:
                mult += int(line.strip().split()[0])
            else:
                continue
    return total,unique,mult
    
    
def main():
    samples = sorted(args.sample.split(","))
    align_dir = os.path.abspath(args.align_dir)
    align_err_file = os.path.basename(glob.glob(os.path.join(align_dir,samples[0],"*.err"))[0])
    fo = open(args.output,"w")
    fo.write("Samples\tTotal reads\tMapped reads\tUniquely mapped reads\tMultiple mapped reads\n")
    for sn in samples:
        align_file = os.path.join(args.align_dir,sn,align_err_file)
        total,unique,mult = stat_align(align_file)
        m = float(unique+mult)/total * 100
        u = float(unique)/total * 100
        mu = float(mult)/total * 100
        fo.write("%s\t%d(100.00%%)\t%d(%.2f%%)\t%d(%.2f%%)\t%d(%.2f%%)\n"%(sn,total,unique+mult,m,unique,u,mult,mu))
    fo.close()
 
if __name__=="__main__":
    main()       
