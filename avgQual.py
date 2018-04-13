#!/usr/bin/env python
#coding:utf-8
## 统计序列中每个位置上的平均碱基质量
import gzip,sys,argparse
from Bio import SeqIO
from collections import OrderedDict

parser = argparse.ArgumentParser(description="This Script is used to get the average quality along the reads position in fastq files.")
parser.add_argument("-i","--input",type=str,help="The input fastq file",required = True)
parser.add_argument("-o","--output",type=str,help="The output quality matrix file",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
args = parser.parse_args()

def average(seq): 
	return float(sum(seq)) / len(seq)

def avgqual(handle): 
    seqs = 	SeqIO.parse(handle,"fastq")
    qc=OrderedDict()
    seq_num,seq_base,q20,q30 = 0,0,0,0
    for i in seqs:
        seq_num += 1
        phred_quality = i.letter_annotations['phred_quality']
        seq_base += len(phred_quality)
        q20 += len([q for q in phred_quality if q >= 20])
        q30 += len([q for q in phred_quality if q >= 30])
        for index,j in enumerate(phred_quality):
            index += 1
            qc.setdefault(index,[]).append(j)
    return qc,seq_num,seq_base,q20,q30
    
def main():
    if args.input.endswith(".gz"):
        handle = gzip.open(args.input)
    else:
        handle = open(args.input)
    m,seq_num,seq_base,q20,q30 = avgqual(handle)
    with open(args.output,"w") as fo:
        fo.write("reads number\t%d\nbases number\t%d\nQ20\t%d\nQ30\t%d\n"%(seq_num,seq_base,q20,q30))
        for x,y in m.iteritems():
            fo.write("%d\t%.2f\n"%(x,average(y)))
  
if __name__=="__main__":
    main()
