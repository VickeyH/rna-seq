#!/usr/bin/env python
#coding:utf-8
import sys,re,argparse,os
import pandas as pd
import numpy as np
from collections import defaultdict
from itertools import combinations
from operator import itemgetter
parser = argparse.ArgumentParser(description="This Script is used to get the gene expression metrix or transcripts expression metrix by the gene abundance txt file or transcripts.txt file.")
parser.add_argument("-c","--count",type=str,help="the read count matrix with header of all sample",required = True, metavar="count_matrix")
parser.add_argument("-r","--ref_file",type=str,help="The ref gtf file or feature:length file for count feature length",required = True,metavar="gtf|length_file")
parser.add_argument("-v","--value",type=str,help="The expression quality value type. 'TPM','FPKM','FPKM_UQ' can be chosen. default: 'TPM FPKM FPKM_UQ'",default = "TPM FPKM FPKM_UQ",nargs="+")
parser.add_argument("-t","--type",type=str,help="The type of count table, default: 'gene'",default = "gene",choices = ["gene","transcript"],metavar = "gene|transcript")
parser.add_argument("-o","--out",type=str,help="The output directory. if not exists, will creat, default: the dirname of count-table",metavar="outdir")
args = parser.parse_args()

def overlap(alist,blist):  ### 必须确保 alist和blist中第一个元素小于等于第二个元素
    if alist[0]<=blist[0]<=alist[1]+1:return [alist[0],max(blist[1],alist[1])]
    if blist[0]<=alist[0]<=blist[1]+1:return [blist[0],max(blist[1],alist[1])]
    else:return

def remove_overlap(alist):
    if len(alist)==1: return alist
    alist.sort(key=lambda x:x[0])
    b,c,s,r= [],[],[],[]
    for i in list(combinations(alist,2)):b.append(overlap(i[0],i[1]))
    for i in b: c.extend(i) if isinstance(i[1],list) else c.append(i)
    c.sort(key=lambda x:x[0])
    for k in set(map(itemgetter(0),c)):s.append(sorted([i for i in c if i[0]==k],key=lambda x:x[-1])[-1])
    for v in set(map(itemgetter(1),s)):r.append([i for i in s if i[1]==v][0])
    length = len(r)
    return r

def remove_overlap2(alist):
    if len(alist)==1: return alist
    alist.sort(key=lambda x:x[0])
    b,c= [],[]
    for i in list(combinations(alist,2)):b.append(overlap(i[0],i[1]))
    for i in b: c.extend(i) if isinstance(i[1],list) else c.append(i)
    d={k:v for k,v in sorted(c,key=lambda x:(x[0],x[1]))}
    dl = sorted([[d[i],i] for i in d],key=lambda x:(x[0],x[1]),reverse=True)
    e = {k:v for k,v in dl}
    return [[e[k],k] for k in e]

def remove_overlap3(alist):
    if len(alist)==1: return alist
    newlist = list(combinations(alist,2))
    for i in newlist:
        if overlap(i[0],i[1]):
            if i[0] in alist:alist.remove(i[0])
            if i[1] in alist:alist.remove(i[1])
            alist.append(overlap(i[0],i[1]))
    return alist

def overlap2(atuple,btuple):  ## 必须保证每个元组中第一个元素小于等于第二个元素
    #atuple = (atuple[0],atuple[1]) if atuple[0] <= atuple[1] else (atuple[1],atuple[0])
    #btuple = (btuple[0],btuple[1]) if btuple[0] <= btuple[1] else (btuple[1],btuple[0])
    if atuple[0]<=btuple[0]<=atuple[1]+1:return (atuple[0],max(btuple[1],atuple[1]))
    elif btuple[0]<=atuple[0]<=btuple[1]+1:return (btuple[0],max(btuple[1],atuple[1]))
    else:return 
    
def remove_overlap4(aset):
    if len(aset)==1: return aset
    newlist = list(combinations(aset,2))
    for i in newlist:
        if overlap2(i[0],i[1]):
            if i[0] in aset:aset.remove(i[0])
            if i[1] in aset:aset.remove(i[1])
            aset.add(overlap2(i[0],i[1]))
    return aset

def get_gene_length(gtf):
    lengthGeneDict = {}
    with open(gtf) as gf:
        for line in gf:
            if line.startswith("#"):
                continue
            if '\tgene\t' in line:
                gene = re.search('gene_id "(.+?)";',line).group(1)
                length = int(line.split("\t")[4]) - int(line.split("\t")[3]) + 1
                lengthGeneDict[gene] = length
    return lengthGeneDict
    
def get_gene_exon_length(gtfile):      ### 计算基因长度时，只计算外显子长度
    geneLenDict = defaultdict(set)
    #geneLenDict = defaultdict(int)           ### 不考虑外显子的overlap,直接进行求和计算，不处理overlap
    with open(gtfile) as gtf:
        for line in gtf:
            if line.startswith("#") or not line.split("\t")[2]=="exon":
                continue
            line_list = line.strip().split("\t")
            geneid = re.search('gene_id "(.+?)";',line_list[-1]).group(1)
            geneLenDict[geneid].add((min(int(line_list[3]),int(line_list[4])),max(int(line_list[3]),int(line_list[4]))))
            #geneLenDict[geneid].sort(key=lambda x:x[0])   ## 若使用remove_overlap4函数，调用overlap2函数，则无需排序
            #geneLenDict[geneid] += abs(int(line_list[3])-int(line_list[4])) + 1  ### 不考虑外显子的overlap,直接进行求和计算，不处理overlap
    for k in geneLenDict:
        if len(geneLenDict[k])==1:
            s = list(geneLenDict[k])[0][0]
            e = list(geneLenDict[k])[0][1]
            geneLenDict[k]=abs(e-s)+1
            continue
        else:
            count = [0,]
            while True:
                if len(remove_overlap4(geneLenDict[k])) == len(remove_overlap4(remove_overlap4(geneLenDict[k]))) == count[-1]:
                #if len(remove_overlap4(geneLenDict[k])) == len(remove_overlap4(remove_overlap4(geneLenDict[k]))):
                    geneLenDict[k] = remove_overlap4(geneLenDict[k])
                    break
                geneLenDict[k] = remove_overlap4(geneLenDict[k])
                count.append(len(geneLenDict[k]))
            geneLenDict[k] = sum([abs(i[0]-i[1]) for i in geneLenDict[k]]) + len(geneLenDict[k])
    return geneLenDict      

def get_gene_exon_length2(gtfile):
    geneLenDict = defaultdict(list)
    with open(gtfile) as gtf:
        for line in gtf:
            if line.startswith("#") or not line.split("\t")[2]=="exon":continue
            line_list = line.strip().split("\t")
            geneid = re.search('gene_id "(.+?)";',line_list[-1]).group(1)
            geneLenDict[geneid].append([min(int(line_list[3]),int(line_list[4])),max(int(line_list[3]),int(line_list[4]))])

    for k in geneLenDict:
        a = geneLenDict[k]
        a.sort()
        b = [a[0]]
        for i in range(1,len(a)):
            if a[i][0] <= b[-1][1] + 1:
                b[-1][1] = max(a[i][1],b[-1][1])
            else:
                b.append(a[i])
        geneLenDict[k] = sum([abs(i[0]-i[1]) for i in b]) + len(b)
    return geneLenDict

def get_isoform_exon_length(gtfile):   ### 计算转录本长度时，只计算外显子长度 对于有overlap的外显子进行合并处理
    transLenDict = defaultdict(set)
    #transLenDict = defaultdict(int)  ### 不考虑外显子的overlap,直接进行求和计算，不处理overlap
    with open(gtfile) as gtf:
        for line in gtf:
            if line.startswith("#") or not line.split("\t")[2]=="exon":
                continue
            line_list = line.strip().split("\t")
            transcript_id = re.search('transcript_id "(.+?)";',line_list[-1]).group(1)
            transLenDict[transcript_id].add((min(int(line_list[3]),int(line_list[4])),max(int(line_list[3]),int(line_list[4]))))
            #transLenDict[transcript_id].sort(key=lambda x:x[0])
            #transLenDict[transcript_id] += abs(int(line_list[3])-int(line_list[4])) + 1   ### 不考虑外显子的overlap,直接进行求和计算，不处理overlap
    for k in transLenDict:
        if len(transLenDict[k])==1:
            s = list(transLenDict[k])[0][0]
            e = list(transLenDict[k])[0][1]
            transLenDict[k]=abs(e-s)+1
            continue
        else:
            count = [0,]
            while True:
                if len(remove_overlap4(transLenDict[k])) == len(remove_overlap4(remove_overlap4(transLenDict[k])))  == count[-1]:
                    transLenDict[k] = remove_overlap4(transLenDict[k])
                    break
                transLenDict[k] = remove_overlap4(transLenDict[k])
                count.append(len(transLenDict[k]))
            transLenDict[k] = sum([abs(i[0]-i[1]) for i in transLenDict[k]]) + len(transLenDict[k])
    return transLenDict 

def get_isoform_exon_length2(gtfile):
    transLenDict = defaultdict(list)
    with open(gtfile) as gtf:
        for line in gtf:
            if line.startswith("#") or not line.split("\t")[2]=="exon":continue
            line_list = line.strip().split("\t")
            transcript_id = re.search('transcript_id "(.+?)";',line_list[-1]).group(1)
            transLenDict[transcript_id].append((min(int(line_list[3]),int(line_list[4])),max(int(line_list[3]),int(line_list[4]))))

    for k in transLenDict:
        a = transLenDict[k]
        a.sort()
        b = [a[0]]
        for i in range(1,len(a)):
            if a[i][0] <= b[-1][1] + 1:
                b[-1][1] = max(a[i][1],b[-1][1])
            else:
                b.append(a[i])
        transLenDict[k] = sum([abs(i[0]-i[1]) for i in b]) + len(b)
    return transLenDict





def get_isoform_length(gtf):
    lengthTranDict = {}
    with open(gtf) as gf:
        for line in gf:
            if line.startswith("#"):
                continue
            if re.search('transcript_id "(.+?)"',line):
                transcript = re.search('transcript_id "(.+?)"',line).group(1)
                length = int(line.split("\t")[4]) - int(line.split("\t")[3]) + 1
                lengthTranDict.setdefault(transcript,length)
    return lengthTranDict
    
def get_len_dict(lengthfile):
    lengthDict = {}
    with open(lengthfile) as gf:
        for line in gf:
            if line.startswith("#"):
                continue
            else:
                lengthDict[line.split()[0]] = line.split()[-1]
    return lengthDict 
    
def main():
    rc = pd.read_csv(args.count,sep="\t",header = 0,index_col=0)
    outputdir = args.out if args.out else os.path.dirname(args.count)
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    if args.type == "gene":
        geneLen = get_gene_exon_length2(args.ref_file)
        if not len(geneLen):
            geneLen = get_len_dict(args.ref_file)
            if not len(geneLen):
                print "Error: can not get gene length dict from %s file" %args.ref_file
                sys.exit(1)
        length = []
        for i in rc.index:
            length.append(geneLen.get(i,np.nan))
        rc["length"] = length
        rc = rc.dropna()
        rc = rc.astype(float)
        if "TPM" in args.value:
            out = os.path.join(outputdir,"gene_TPM_matrix.txt")
            new_count = rc.apply(lambda x:x[:-1]/x[-1]*1000,axis=1)  ### reads per kilobase (RPK)
            new_count.fillna(value=0,inplace=True)
            new_count = new_count.apply(lambda x:x/np.sum(x)*1000000,axis=0)
            #new_count["length"] = rc["length"].astype(int)
            new_count.to_csv(out,sep="\t")
        if "FPKM" in args.value:
            out = os.path.join(outputdir,"gene_FPKM_matrix.txt")
            new_count1 = rc.iloc[:,:-1].apply(lambda x:x/np.sum(x),axis=0)             
            new_count1["length"] = rc["length"]
            new_count1 = new_count1.apply(lambda x:10**9*x[:-1]/x[-1],axis=1)
            new_count1.fillna(value=0,inplace=True)
            #new_count1["length"] = rc["length"].astype(int)
            new_count1.to_csv(out,sep="\t")
        if "FPKM_UQ" in args.value:
            out = os.path.join(outputdir,"gene_FPKM_UQ_matrix.txt")
            new_count2 = rc.iloc[:,:-1].apply(lambda x:x/np.percentile(x,75),axis=0)
            new_count2["length"] = rc["length"]
            new_count2 = new_count2.apply(lambda x:10**9*x[:-1]/x[-1],axis=1)
            new_count2.fillna(value=0,inplace=True)
            #new_count2["length"] = rc["length"].astype(int)
            new_count2.to_csv(out,sep="\t")
        
    if args.type == "transcript":
        transLen = get_isoform_exon_length2(args.ref_file)
        if not len(transLen):
            transLen = get_len_dict(args.ref_file)
            if not len(transLen):
                print "Error: can not get transcript length dict from %s file" %args.ref_file
                sys.exit(1)
        length = []
        for i in rc.index:
            length.append(transLen.get(i,np.nan))
        rc["length"] = length
        rc = rc.dropna()
        rc = rc.astype(float)
        if "TPM" in args.value:
            out = os.path.join(outputdir,"transcript_TPM_matrix.txt")
            new_count = rc.apply(lambda x:x[:-1]/x[-1]*1000,axis=1)  ### counts per kilobase
            new_count.fillna(value=0,inplace=True)
            new_count = new_count.apply(lambda x:x/np.sum(x)*1000000,axis=0)
            #new_count["length"] = rc["length"].astype(int)
            new_count.to_csv(out,sep="\t")            
        if "FPKM" in args.value:
            out = os.path.join(outputdir,"transcript_FPKM_matrix.txt")
            new_count1 = rc.iloc[:,:-1].apply(lambda x:x/np.sum(x),axis=0) 
            new_count1["length"] = rc["length"]
            new_count1 = new_count1.apply(lambda x:10**9*x[:-1]/x[-1],axis=1)
            new_count1.fillna(value=0,inplace=True)
            #new_count1["length"] = rc["length"].astype(int)
            new_count1.to_csv(out,sep="\t")
        if "FPKM_UQ" in args.value:
            out = os.path.join(outputdir,"transcript_FPKM_UQ_matrix.txt")
            new_count2 = rc.iloc[:,:-1].apply(lambda x:x/np.percentile(x,75),axis=0)
            new_count2["length"] = rc["length"]
            new_count2 = new_count2.apply(lambda x:10**9*x[:-1]/x[-1],axis=1)
            new_count2.fillna(value=0,inplace=True)
            #new_count2["length"] = rc["length"].astype(int)
            new_count2.to_csv(out,sep="\t")

if __name__ == "__main__":
    main()
  
    
